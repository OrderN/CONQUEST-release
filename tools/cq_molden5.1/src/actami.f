      subroutine actamd(inum,ikleur,iopt,idosrf,
     &                  iaton,iatclr,iresid,isurf,
     &            icalf,ncalf,ianf,islu,nchain,iamino,ihet,
     &            iclhet,reson,iams,isal)
      implicit double precision (a-h,o-z)
      parameter (mxheta=150)
      common /athlp/ iatoms, mxnat
      common /surf/ natorg,nosncd
      integer reson
      parameter (mxres=42)
      parameter (mxchai=50)
      common /clfhlp/ isndcl(4),iamicl(mxres),ichcol(mxchai)
      integer srfmap,srfloc
      common /hlpsrf/ npts1,npts2,npts3,srfmap,srfloc,ifogl,itsrf
      parameter (mxaln=20)
      common /palign/  nalign,istst(mxaln),istch(mxaln),istres(mxaln),
     &                 istcol(mxaln),istptr(mxaln)
      logical bbone
      dimension iaton(*),iatclr(*),iresid(*),isurf(*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),ihet(*),
     &          iclhet(*),reson(*),iams(*),isal(*)

      if (inum.le.0.and.inum.ge.-3) return

      if (inum.lt.-3) then
         iht = abs(inum) + 1
      endif

      iactkl = ikleur
      if (ikleur.eq.0.and.inum.gt.0) then
          if (nalign.gt.0) then
             do i=1,nalign
                if (inum.le.istres(i)) then
                    iactkl = istcol(i) + 1
                    goto 10
                endif
             end do
10           continue
          else
             if (iamino(inum).ne.0) iactkl = iamicl(iamino(inum))
          endif
      endif
c      if (ikleur.eq.0.and.inum.lt.-3) iactkl = 1
      if (ikleur.eq.0.and.inum.lt.-3) then
          if (iht.le.mxheta) then
             iactkl = iclhet(iht)
          endif
      endif

      if (inum.gt.0) then
          if (iopt.eq.0) then
              reson(inum) = 0
          else
              if (ifogl.eq.1.and.idosrf.eq.1) iams(inum) = 1
              reson(inum) = 1
          endif
      elseif (inum.lt.-3) then
          if (iht.le.mxheta) then
             ihet(iht) = iopt
             if (iopt.ne.0) iclhet(iht) = iactkl
          endif
      endif

      do i=1,iatoms
          if (iresid(i).eq.inum) then
             if (iopt.eq.1) then
                    iaton(i) = 1
                    iatclr(i) = iactkl
                    if (idosrf.eq.1) isurf(i) = 1
             else
                if (iamino(abs(inum)).gt.23) then
                    n = 6
                else
                    n = 3
                endif
                bbone = .false.
                do j=1,n
                    if (i.eq.icalf(j,abs(inum))) bbone = .true.
                end do
                if (bbone.and.ihet(isal(abs(inum))+1).eq.0) then
                     do k=1,nchain
                         do j=ianf(k),islu(k)
                            do l=1,n
                               if (i.eq.icalf(l,j)) 
     &                         iatclr(icalf(l,j)) = ichcol(k)
                            end do
                         end do
                     end do
                else
                     iaton(i) = 0
                endif
             endif
          endif
      end do


      return
      end

      subroutine actexd(que,lque,ikleur,idosrf,ncalf,iamino)
      implicit double precision (a-h,o-z)
      dimension iamino(*)
      logical oque
      integer que(20)

      do i=1,ncalf
         oque = .true.
         do j=1,lque
            iq = que(j)
            ia = iamino(i+j-1)
            if (iq.ne.21.and.iq.ne.ia) then
                if(.not.(((iq.eq.9.or.iq.eq.10).and.ia.eq.21) .or.
     &            ((iq.eq.13.or.iq.eq.14).and.ia.eq.22) )) then
                   oque = .false.
                endif
            endif
         end do
         if (oque) then
            do j=1,lque
               call actami(i+j-1,ikleur,1,idosrf)
            end do
         endif
      end do

      return
      end

      subroutine aacod(vrad,ires,str,nstr,istsrf,iresid,isurf,
     &                 ncalf,ianf,islu,nchain,iamino,ihet,reson,
     &                 isal,irsnr,achain);
c Process commands from the residue commands window
      implicit double precision (a-h,o-z)
      parameter (mxheta=150)
      parameter (mxres=42)
      common /athlp/ iatoms, mxnat
      integer reson
      common /hlpsrf/ npts1,npts2,npts3,srfmap,srfloc,ifogl,itsrf
      character*1 achain
      common /aachlp/ scndst(3),onelet(21)
      common /amino/  aminos(mxres)
      common /helpar/ helrod, ihtype
      common /surf/ natorg,noscnd

      integer que(20)
      logical doque,oclr,oami,oorg,ohasc,oamis(mxres)
      integer srfmap,srfloc
      character keywrd*160
      character*3 aminos,scndst
      character*20 amique
      character*1 onelet,ach
      character*(*) str
      dimension iresid(*),isurf(*)
      dimension ianf(*),islu(*),iamino(*),ihet(*),reson(*),
     &          isal(*),irsnr(*),achain(*)


      data scndst/'HEL','BET','BCK'/

      data onelet/'G','A','S','C','T','I','V','M','D','N','L','K',
     &            'E','Q','P','R','H','F','Y','W','X'/


      toang  = 0.52917706d0
      keywrd = str(1:nstr)
      call tocap(keywrd,nstr)

      if (index(keywrd,'SURON').ne.0) then
          istsrf = 1
          if (index(keywrd,'MAPPED').ne.0) srfmap = 3
          if (index(keywrd,'LOCAL').ne.0) srfloc = 1
          if (index(keywrd,'GLOBAL').ne.0) srfloc = 0
          call inferr('Switching on Surface Drawing',0)
          return
      endif
      if (index(keywrd,'SUROFF').ne.0) then
          istsrf = 0
          srfmap = 0
          srfloc = 1
          call inferr('Switching off Surface Drawing',0)
          return
      endif

      do i=1,mxres
         oamis(i) = .false.
      end do

      doque = .false.
      ibr1 = index(keywrd,'(')
      ibr2 = index(keywrd,')')
      if (ibr1.ne.0.and.ibr2.ne.0.and.ibr2.gt.ibr1) then
          nque = ibr2-ibr1-1
          if (nque.gt.20) then
             call inferr('(Too long) !',0)
          else
             doque = .true.
             amique = keywrd(ibr1+1:ibr2-1)
             do i=ibr1,ibr2
                keywrd(i:i) = ' '
             end do
             lque = 0
             do i=1,nque
                kque = 0
                do j=1,21
                   if (amique(i:i).eq.onelet(j)) then
                      lque = lque + 1
                      que(lque) = j
                      kque = 1
                   endif
                end do
                if (amique(i:i).ne.' '.and.kque.eq.0) 
     &          call inferr('Unrecognised Symbol(s) ignored !',0)
             end do
          endif
      endif

      oami = .false.
      oclr = (index(keywrd,'NOT').ne.0)
      if (oclr) oami = .true.
      if (index(keywrd,'ALL').ne.0.or.
     &    index(keywrd,'CLEAR').ne.0) then
         oami = .true.
         i1 = 1
         i2 = mxres
         if (index(keywrd,'AMINO').ne.0) i2 = 23
         if (index(keywrd,'NUCL').ne.0) i1 = 24
         do i=i1,i2
            oamis(i) = .true.
         end do
      endif
      if (index(keywrd,'CLEAR').ne.0) oclr = .true.

      if (istsrf.eq.1) then
          idosrf = 1
      else
          idosrf = 0
      endif
      if (index(keywrd,'SURF').ne.0) idosrf = 1

      if (idosrf.eq.1.and.oclr) then
          call clrsrf
          return
      endif

      if (idosrf.eq.1.and.oami.and.ifogl.ne.1) then
          call allsrf(1,0,1,4)
          return
      endif

c ires may alternatively contain color

      ikleur = ires
      if (ikleur.lt.1.or.ikleur.gt.15) ikleur = 0
      ii = index(keywrd,'COL')
      if (ii.ne.0) then
          keywrd(ii:ii+2) = '   '
      endif

      if (index(keywrd,'NEIGH').ne.0) then
         rdist = 3.0d0
         ires1 = index(keywrd,'=') 
         if (ires1.ne.0) 
     &     rdist = reada(keywrd,ires1+1,len(keywrd))*1
         rdist = rdist / toang
         call proxim(iresid(ires),rdist,0,idosrf)
      endif

      if (index(keywrd,'ROD').ne.0) then
         ihtype = 1
         vrad = .0d0
      endif
      if (index(keywrd,'BAL').ne.0) then
         ihtype = 0
         vrad = .9d0 / toang
      endif

      if (index(keywrd,'-').eq.0) then
         do ii=1,3
            if (index(keywrd,scndst(ii)).ne.0) then
               inclbb = 0
               if (index(keywrd,'INCL').ne.0) inclbb = 1
               if (oclr) then
                  ihet(ii) = 0
                  call acthel(0,ii-1,ikleur,inclbb)
               else
                  ihet(ii) = 1
                  call acthel(1,ii-1,ikleur,inclbb)
               endif
            endif
         end do
      endif

      if (index(keywrd,'CHA').ne.0.or.
     &    index(keywrd,'POS').ne.0.or.
     &    index(keywrd,'POL').ne.0) then
         oami = .true.
         oamis(12) = .true.
         oamis(16) = .true.
         oamis(17) = .true.
      endif
      if (index(keywrd,'CHA').ne.0.or.index(keywrd,'NEG').ne.0.or.
     &    index(keywrd,'POL').ne.0) then
         oami = .true.
         oamis( 9) = .true.
         oamis(13) = .true.
      endif
      if (index(keywrd,'POL').ne.0) then
         oami = .true.
         oamis(10) = .true.
         oamis(14) = .true.
         oamis(3) = .true.
         oamis(5) = .true.
      endif
      if (index(keywrd,'ARO').ne.0) then
         oami = .true.
         oamis(18) = .true.
         oamis(19) = .true.
         oamis(20) = .true.
      endif
      if (index(keywrd,'ALI').ne.0) then
         oami = .true.
         oamis(1) = .true.
         oamis(2) = .true.
         oamis(6) = .true.
         oamis(7) = .true.
         oamis(11) = .true.
         oamis(15) = .true.
      endif

      do i =1,mxres
         if (index(keywrd,aminos(i)).ne.0) then
            oami = .true.
            oamis(i) = .true.
         endif
      end do

      if (oami) then
         do i=1,mxres
            if (oamis(i)) then
               do j=1,islu(nchain)
                  if (iamino(j).eq.i) then
                     if (oclr) then
                        call actami(j,ikleur,0,idosrf)
                     else
                        call actami(j,ikleur,1,idosrf)
                     endif
                  endif
               end do
            endif
         end do
      else
         iscnd = -1
         if (index(keywrd,'HEL').ne.0) then
             iscnd = 0
         elseif (index(keywrd,'BET').ne.0) then
             iscnd = 1
         elseif (index(keywrd,'RND').ne.0) then
             iscnd = 3
         endif
         ires1 = reada(keywrd,1,len(keywrd))*1
         ires2 = index(keywrd,'-') 
         istar = index(keywrd,'*')
         ichp = index(keywrd,':')
         if (ires2.ne.0) then
            ires2 = reada(keywrd,ires2+1,len(keywrd))*1
         else
            ires2 = ires1
         endif
         oorg = ohasc(ncalf,achain)
         if (ichp.ne.0) then
            ach = keywrd(ichp+1:ichp+1)
            iach = ichar(ach)
            if (iach.ge.65.and.iach.le.90) then
               ifnd = 0
               do j=1,ncalf
                  if (achain(j).eq.ach) ifnd = 1
               end do
               if (ifnd.eq.0) then
                  call inferr('Chain '//ach//' NOT found!',0)
                  ich = 1
                  ires1 = 0
                  oorg = .false.
               endif
            elseif (iach.eq.32) then
               ach = achain(1)
               ich = 1
            else
               ich = reada(keywrd,ichp+1,len(keywrd))*1
               if (ich.lt.1.or.ich.gt.nchain) then
                   call inferr('Invalid Chain number',0)
               else
                   ach = achain(ich)
               endif
            endif
         else
            ach = ' '
            ich = 1
         endif
         if (istar.ne.0.and.ichp.ne.0) then
            if (istar.lt.ichp) then
              if (.not.oorg) then
                ires1 = irsnr(ianf(ich))
                ires2 = irsnr(islu(ich))
              else
                ires1 = 0
                ires2 = 0
                do j=1,ncalf
                   if (achain(j).eq.ach) then
                      if (ires1.eq.0) then
                         ires1 = j
                      endif
                      ires2 = j
                   endif
                end do
              endif
            endif
         endif
         if (ires1.ne.0) then
            if (iscnd.eq.-1) then
               if (index(keywrd,'PRINT').ne.0) then
                  do i=ires1,ires2
                     do j=1,ncalf
                        if (irsnr(j).eq.i)
     &                    print*,i,' ',isal(j)
                     end do
                  end do
               else
                  do i=ires1,ires2
                     if (oorg.and.ach.ne.' ') then
                        do j=1,ncalf
                           if (irsnr(j).eq.i.and.achain(j).eq.ach)
     &                        call actami(j,ikleur,1,idosrf)
                        end do
                     else
                        do j=1,ncalf
                           if (irsnr(j).eq.i)
     &                        call actami(j,ikleur,1,idosrf)
                        end do
                     endif
                  end do
               endif
            else
               iatoms = noscnd
               do i=1,iatoms
                  if (isurf(i).eq.1) idosrf = 1
               end do
               natorg = 0
               do i=ires1,ires2
                  isal(ianf(ich)+i-1) = iscnd
               end do
               call actcal(1)
               do i=1,ncalf
                  reson(i) = 0
               end do
               do i=4,mxheta
                  ihet(i) = 0
               end do
               call ribbon(0,0,1,0,0,0,0)
               call ribbon(1,0,1,0,0,0,0)
               call ribbon(2,0,1,0,0,0,0)
               call ribbon(3,0,1,0,0,0,0)
               do i=1,4
                  ihet(i) = 1
               end do
            endif
         endif
      endif

      if (doque) call actexp(que,lque,ikleur,idosrf)

      if (idosrf.eq.1.and..not.oclr) then
          if (ifogl.eq.1) then
             call initsrf
          else
             call connlp(1.0d0,0,4)
          endif
      endif

      return
      end

      logical function ohasc(ncalf,achain)
      implicit double precision (a-h,o-z)
      character*1 achain
      dimension achain(*)

      ohasc = .true.
      ihasc = 0
      do i=1,ncalf
         if (achain(i).eq.' ') ihasc = ihasc + 1
      end do
      if (ihasc.eq.ncalf) ohasc = .false.

      return
      end

