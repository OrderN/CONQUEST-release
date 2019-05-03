      subroutine wrtoud(iun,emin,iconn,ityp,coo,q,iopt)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxcon=10)
      parameter (numres=50000)
      common /athlp/ iatoms, mxnat
      integer*2 ityp
      character*3 chtmp,ambstr
      character*2 gffstr
      common /ffstr/ ambstr(mxamb), gffstr(mxgff)
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /residu/  qtot,ihsres,
     &                 ires(numres),ibeg(numres),iend(numres)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /prot/ ireswr
      integer taskid
      common /mpih/ nproc,taskid
      character*3 hetz
      common /resstr/ hetz(numres)
      dimension icnn(mxcon)
      dimension coo(3,*),q(*),iconn(mxcon+1,*),ityp(*),iopt(*)

      if (nproc.gt.1.and.taskid.ne.0) return

c write tinker .xyz file

      if (box) then
         write(iun,'(a,3(f9.3,1x))') 
     &      '[AMBFOR] box ',abc(1),abc(2),abc(3)
      else
         if (cell) then
            call getcel(abc(1),abc(2),abc(3),alpha,beta,gamma)
            write(iun,'(a,6(f9.3,1x))') 
     &      '[AMBFOR] cell ',abc(1),abc(2),abc(3),alpha,beta,gamma
         else
            write(iun,*) '[AMBFOR]'
         endif
      endif

      write(iun,'(i5,a,f15.3,a)') iatoms,' emin ',emin,
     &   '  AMBFOR/AMBMD generated .xyz (amber/gaff param.)'

      do i=1,iatoms

            ibnds = 0
            do j=1,iconn(1,i)
               if (iconn(j+1,i).gt.0) then
                  ibnds = ibnds + 1
                  icnn(ibnds) = iconn(j+1,i)
               endif
            end do

            itypi = ityp(i)

            if (itypi.gt.0) then

c Amber
                if (itypi.le.mxamb) then
                  chtmp = ambstr(itypi)
                  if (itypi.gt.mxamb) itypi = itypi + (2000-mxamb)
                endif

                write(iun,'(i6,2x,a3,1x,3(f10.5),1x,i4,8(1x,i6))')
     &            iopt(i),chtmp,(coo(j,i),j=1,3),
     &            itypi,(icnn(j),j=1,ibnds)

            elseif (ityp(i).lt.0) then

c GAFF
                if (iabs(itypi).le.mxgff) then
                  chtmp = gffstr(iabs(itypi))//' '
                endif

                write(iun,
     &            '(i6,2x,a3,1x,3(f10.5),1x,i4,1x,f6.3,8(1x,i6))')
     &            iopt(i),chtmp,(coo(j,i),j=1,3),
     &            itypi,q(i),(icnn(j),j=1,ibnds)
            endif


      end do

      if (iun.eq.iun3.or.(iun.eq.iun4.and.ireswr.eq.0)) then

         write(iun,'(a)') '[RESIDUES]'

         do i=1,ihsres
            if (hetz(i).ne."   ") then
               write(iun,'(i6,1x,i6,1x,a3)') ires(i),ibeg(i),hetz(i)
            else
               write(iun,'(i6,1x,i6)') ires(i),ibeg(i)
            endif
         end do

         if (iun.eq.iun4) ireswr = 1

      endif

      return
      end

      subroutine wrtbid(iun,emin,coo)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /athlp/ iatoms, mxnat
      logical box,cell,fast,addbox
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      dimension coo(3,*)

c write tinker coordinates binary form

      write(iun) iatoms,emin

      do i=1,iatoms
           write(iun) (coo(j,i),j=1,3)
      end do

      if (cell) then
           call getcel(abc(1),abc(2),abc(3),alpha,beta,gamma)
           write(iun) abc(1),abc(2),abc(3),alpha,beta,gamma
      endif

      return
      end

      subroutine opfiles
      implicit real (a-h,o-z)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character fniun*132
      common /fnam/ fniun, lenf
      common /opts/ idebug,imon,iarc,ilog,iout,osingl,dolbfgs,oqscal
      character fniunt*132
      logical opfil, osingl,dolbfgs,oqscal

      fniunt = fniun(1:lenf)//'.xyz'
      lenft = lenf + 4

      if (opfil(48,fniunt,lenft,1,1,0)) then
         iun2 = 48
         if (iout.eq.1) then
            fniunt = fniun(1:lenf)//'_opt.xyz'
            lenft = lenf + 7
            if (opfil(50,fniunt,lenft,1,0,0)) then
               iun3 = 50
            else
               print*,'Cant open output file !'//fniunt(1:lenft)
               stop
            endif
         endif

         if (iarc.eq.1) then
            fniunt = fniun(1:lenf)//'.arc'
            lenft = lenf + 4
            if (opfil(60,fniunt,lenft,1,0,0)) then
               iun4 = 60
            else
               print*,'Cant find/open arc file !'//fniunt(1:lenft)
               stop
            endif
         endif

         if (ilog.eq.1) then
            fniunt = fniun(1:lenf)//'.log'
            lenft = lenf + 4
            if (opfil(70,fniunt,lenft,1,0,0)) then
               iun5 = 70
            else
               print*,'Cant find/open log file !'//fniunt(1:lenft)
               stop
            endif
         else
            iun5 = 6
         endif

      else
         print*,'Cant find/open file !'//fniunt(1:lenft)
         stop
      endif

      return
      end

      logical function opfil(iun,filenm,lenfn,iform,iold,isil)
      implicit real (a-h,o-z)
      character*(*) filenm
      character*7 stat

      opfil = .true.

      if (iold.eq.1) then
          stat = 'old'
      else
          stat = 'unknown'
      endif

      if (lenfn.eq.0.or.filenm.eq.' ') then
          print*, 'Invalid Filename !'
          opfil = .false.
      else
          close(iun)
          if (iform.eq.1) then
             open(unit=iun,form='formatted',file=filenm,
     &               status=stat,err=100)
          else
             open(unit=iun,form='unformatted',file=filenm,
     &               status=stat,err=100)
          endif
      endif

      return

100   if (isil.eq.0) print*,'Error Opening File !'

      opfil = .false.
      return
      end

      subroutine wrmod(ncyc,emin)
      implicit real (a-h,o-z), integer (i-n)
      character fniun*132
      common /fnam/ fniun, lenf
      common /opts/ idebug,imon,iarc,ilog,iout,osingl,dolbfgs,oqscal
      integer taskid
      common /mpih/ nproc,taskid
      character*20 ambfil
      character*3 tstr
      character*4 ttstr
      character*5 tttstr
      logical opfil, osingl,dolbfgs,oqscal

      if (nproc.gt.1.and.taskid.ne.0) return

      jlen = lenf
      if (ncyc.lt.1000) then
         call zerstr(ncyc,tstr,3)
         ambfil = fniun(1:lenf)//'.'//tstr
         jlen = jlen + 1 + 3
      elseif (ncyc.lt.10000) then
         write(ttstr,'(i4)') ncyc
         ambfil = fniun(1:lenf)//'.'//ttstr
         jlen = jlen + 1 + 4
      elseif (ncyc.lt.100000) then
         write(tttstr,'(i5)') ncyc
         ambfil = fniun(1:lenf)//'.'//tttstr
         jlen = jlen + 1 + 5
      endif


      if (imon.eq.2.and.ncyc.ne.1) then
         if (opfil(51,ambfil,jlen,0,0,0)) then
             call wrtbin(51,emin)
             close(51)
         endif
      else
         if (opfil(51,ambfil,jlen,1,0,0)) then
             call wrtout(51,emin)
             close(51)
         endif
      endif


      return
      end

      subroutine wrtesc
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /passe/ emin
      
      call wrtout(iun3,emin)

      return
      end

      subroutine wmdesc
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      
      close(iun4)

      return
      end

