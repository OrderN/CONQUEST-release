      subroutine convzmzz(cc,ianc,natoms,igo,ico,ido,
     &              bl,alph,bet,ibl,ialph,ibet,imap,ianz,iz,
     &              c,cz,alpha,beta,ian)

c this is really convzmat
c
c This routines reads a zmatrix from a gamess/gaussian output
c and converts it to cartesian coordinates (stoc)
c

      implicit double precision (a-h,o-z)
      logical ottest,oerror,zreadg,zreado
c
      common /zmfrst/ ihaszm, nz, mxzat
      dimension ian(*),c(3,*),cz(3,*),alpha(*),beta(*)
      dimension cc(3,*),ianc(*)
      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),ianz(*),iz(4,*)
c
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

cd     write(iun3,*)'enter subroutine convzmat'
      do i=1,3
         imap(i) = i
         do j=1,4
          iz(j,i) = 0
         end do
      end do

      maxnz = mxzat
      if (igo.eq.1) then
         if (.not.zreadg(nz,ianz,iz,bl,alph,bet)) then
            call haszm(.false.)
cd           write(iun3,*)'leave subroutine convzmat'
            return
         endif
      elseif (igo.eq.2) then
         if (.not.zreado(nz,ianz,iz,bl,alph,bet)) then
            call haszm(.false.)
            return
         endif
      else
         call zread(nz,ianz,iz,bl,alph,bet)
      endif

      if (ido.eq.1) then
         do i=1,3
            do l=1,3
               cz(l,i) = cc(l,imap(i))
            end do
         end do
      else
         do i=1,3
            do l=1,3
               cz(l,i) = 0.0d0
            end do
         end do
      endif

      ottest=.true.
      call stoc(maxnz,nz,0,0,0,ianz,iz,bl,alph,bet,ottest,natoms,
     &     ian,c,cz,imap,alpha,beta,oerror,.true.,.false.)
      if (.not.oerror) then
         ihaszm = 1
         do i=1,nz
            ibl(i) = 1
            ialph(i) = 1
            ibet(i) = 1
         end do
         if (ico.eq.1) then
            do i=1,natoms
               do j=1,3
                  cc(j,i) = c(j,i)
               end do
               ianc(i) = ian(i)
            end do
         endif
      endif

cd     write(iun3,*)'leave subroutine convzmat'
      return
      end
