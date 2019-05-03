      subroutine rdmdf(idebug,istat,coo,ianz,iatclr,iconn,qat,ityp,
     &                 nat,norg,icent,inorm,ncon,nspg,kz,ichx,
     &                 nopr,ir,it,
     &                 xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxmsf=235)
      parameter (numatm=2000)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /athlp/iatoms, mxnat
      integer*2 ir,it
      character*137 line,str
      common /curlin/ line
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /msft/ ltyp(mxmsf)
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      integer*2 ityp
      common /types/ iff
      character*80 col
      character*200 header
      character*100 fname
      character*6 atnam
      character*10 version,label
      character*4 type
      character*4 c4dum
      character*6 c6dum
      integer*2 i2dum,ii2dum,i22dum
      real rdum,rrdum
      dimension rrdum(numat1),
     &          ii2dum(numat1),i22dum(8,numat1)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*),qat(*),
     &          ityp(*),ir(3,3,192),it(3,192)
      data ltyp /4*1,99,4,5,1,1,20*6,10*7,20*8,10*15,10*16,3,11,12,
     &           19,20,25,26,30,37,55,14,13,9,17,35,53,29,23,24,27,
     &           28,33,34,38,39,40,41,42,44,45,46,47,48,50,51,56,74,
     &           76,78,79,80,82,83,57,58,59,89,90,92,52,84,85,43,
     &           21,22,31,32,49,72,73,77,81,87,88,60,61,62,63,64,65,
     &           66,67,68,69,70,71,91,93,94,95,96,97,98,5*99,
     &           2,10,18,36,54,86,99,99,75,14,2*99,7*7,3*99,10*6,17,
     &           4*6,6*7,3*8,6*99,3*7,7*99,3*6,15,7,99/

      if (idebug.eq.1) print*,'routine rdmsf'

      istat = 1
      toang = 0.52917706d0
      ncon = 0

      read(iun2,err=90,end=90) ns,nr,iatoms,version,header

      if (nr.eq.0.or.iatoms.eq.0) goto 90
c
c  read titles
c

      col(1:3) = '   '
      do while (col(1:3).ne.'END')
         read(iun2,err=90,end=90) col
      end do
c
c read segment info
c
      read(iun2,err=90,end=90) (c4dum,i=1,ns)
      read(iun2,err=90,end=90) (i2dum,i=1,ns)
      read(iun2,err=90,end=90) (i2dum,i=1,ns)

c
c read residue information
c
      read(iun2,err=90,end=90) (c6dum,i=1,nr)
      read(iun2,err=90,end=90) (c4dum,i=1,nr)
      read(iun2,err=90,end=90) (i2dum,i=1,nr)
      read(iun2,err=90,end=90) (i2dum,i=1,nr)
      read(iun2,err=90,end=90) (i2dum,i=1,nr)

      do j=1,4
         read(iun2,err=90,end=90) (rdum,i=1,nr)
      end do
c
c read atom info
c
      do j=1,3
         read(iun2,err=90,end=90) (rrdum(i),i=1,iatoms)
         do i=1,iatoms
            coo(j,i) = dble(rrdum(i))/toang
         end do
      end do

      do i=1,iatoms
         iatclr(i) = 1
         iconn(1,i) = 0
      end do

      if (idebug.eq.1) then
         do i=1,iatoms
            print*,i,' ',(coo(j,i),j=1,3)
         end do
      endif
 
      read(iun2,err=90,end=90) (atnam,i=1,iatoms)

      read(iun2,err=90,end=90) (i2dum,i=1,iatoms)

      read(iun2,err=90,end=90) (ii2dum(i),i=1,iatoms)

      if (idebug.eq.1) then
         print*,'atom types:'
         do i=1,iatoms
            print*,i,' ',ii2dum(i)
         end do
      endif

      do i=1,iatoms
         if (ii2dum(i).gt.mxmsf) then
            ianz(i) = 99
            ityp(i) = mxmsf
         else
            ianz(i) = ltyp(int(ii2dum(i)))
            ityp(i) = ii2dum(i)
         endif
      end do

      iff = 6

      read(iun2,err=90,end=90) (rrdum(i),i=1,iatoms)
      do i=1,iatoms
          qat(i) = dble(rrdum(i))
      end do
      ihasq = 1

c
c read the extra data
c

10    read(iun2,err=90,end=100) label,type,nitems,fname

      if (type.eq.'REAL'.or.type.eq.'ASTR') then

        read(iun2,end=90,err=90) (rdum,n=1,nitems)

      elseif (type.eq.'REA3') then

        do i=1,3
            read(iun2,end=90,err=90) (rrdum(j),j=1,nitems)
            do j=1,nitems
c               coo(i,j) = dble(rrdum(j))
            end do
        end do

      elseif (type.eq.'INT4') then

        read(iun2,end=90,err=90) (i4dum,n=1,nitems)

      elseif (type.eq.'INT2') then

        read(iun2,end=90,err=90) (i2dum,n=1,nitems)

      elseif (type.eq.'BOND') then

        if (label.eq.'CONNECT   ') then

          read(iun2,end=90,err=90) 
     &      (ii2dum(i),(i22dum(j,i),j=1,ii2dum(i)),i=1,nitems)
          do i=1,nitems
             iconn(1,i) = int(ii2dum(i))
             if (iconn(1,i).gt.mxcon) iconn(1,i) = mxcon
             do j=1,iconn(1,i)
                iconn(1+j,i) = int(i22dum(j,i))
             end do
          end do
          ncon = 1

          call redcon

        elseif (label.eq.'ORDER     ') then

          read(iun2,end=90,err=90)
     &      (i2dum,(i2dum,j=1,2),i=1,nitems)

        endif

      elseif (type.eq.'SYMM') then

        kz = 0
        idone = 0
        do while(idone.eq.0)
           read(iun2,end=90,err=90) col(1:80)
           if (col(1:4).eq.'CELL'.and.iatoms.lt.1000) then
              line(1:80) = col(1:80)
              ktype = nxtwrd(str,nstr,itype,rtype)
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.3) a = rtype
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.3) b = rtype
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.3) c = rtype
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.3) alpha = rtype
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.3) beta = rtype
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.3) gamma = rtype
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.2) then
                 if (itype.ne.1) goto 100
              else 
                 goto 100
              endif
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.2) nspg = itype

              call prcell(nspg,a,b,c,alpha,beta,gamma)
              call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
              call cpmol(nat,norg,a,b,c,alpha,beta,gamma,
     &                   coo,ianz,iatclr,iconn)
              call cprot(nspg,nopr,icent,ir,it,.false.)
              istat = 2
              ichx = 1
           elseif (col(1:3).eq.'END') then
              idone = 1
           endif

        end do

      endif

      goto 10
c
90    call inferr('Error reading MSF',0)
      iatoms = 0
      istat = 0
      iff = 0
      return

100   continue
      iftyp = 6
      call cooxyz(ianz,iatoms)
      if (iatoms.gt.100) call haszm(.true.)
      return
      end
