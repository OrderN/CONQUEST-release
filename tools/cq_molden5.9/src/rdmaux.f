      subroutine rdmadd(idebug,istatio,istats,
     &                 vectrs,vectrb,focc,focb,eiga,eigb,
     &                 averag,p,halfs,psi,
     &                 nocc,nocb,ncols,ncolb,
     &                 formax,forrms,dismax,disrms,epoints,isav,q)

c THIS IS REALLY rdmaux
 
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (MAXITER=1000)
      common /athlp/ iatoms, mxnat
      common /moldat/ natoms,norbs,nelecs,nat(numatm) 
      common /coord / xyz(3,numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      common /convrg/ convg1(MAXITER),convg2(MAXITER),cnvmax,cnvmin,
     &jstrt1,jend1,jstrt2,jend2,icvav1,icvav2
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      character*137 line,str
      character*2 elt
      character*2 elemnt,tocapf
      common /elem/elemnt(mxel)
      common /mopac/  nfirst(numatm),nlast(numatm),
     &                npq(numatm),pqn(numatm),emus(numatm),
     &                emup(numatm),emud(numatm),consts(numatm),
     &                constp(numatm),constd(numatm),npqref(54)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /curlin/ line
      integer getlin
      logical gnreal
      real eiga,eigb
      dimension vectrs(*),vectrb(*),focc(*),focb(*),eiga(*),eigb(*),
     &          averag(*),p(*),halfs(*),psi(*)
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*),q(*)

      istatio = 1
      toang = 0.52917706d0
      isgau = 0
      istats = 0
      iatfnd = 1
      nelecs = 0
      relecs = 0.0d0
      norbs = 0
      irtype = 0
      ngeoms = 0

      if (idebug.eq.1) call inferr('subroutine rdmaux',0)

      call rewfil

      call search(line,'ATOM_EL',istat)
      if (istat.eq.1) then
         natoms = intlin(line,istat)
         if (istat.eq.0) goto 100

         n40 = natoms / 40
         if (natoms.gt.n40*40) n40 = n40 + 1

         itel = 0
         do j=1,n40

            if (getlin(0).ne.1) goto 100

            do i=1,40

                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.1) then
                   itel = itel + 1
                   elt = str(1:2)
                   if (nstr.eq.1) elt = ' '//str(1:1)
                   do k=1,mxel
                     if (tocapf(elt).eq.tocapf(elemnt(k))) nat(itel) = k
                   end do
                elseif (ktype.eq.0) then
                   goto 5
                else
                   print*,line
                   call inferr('string expected',1)
                   goto 100
                endif

            end do

         end do

      else
         call inferr('No ATOM_EL card found',1)
         goto 100
      end if

5     continue

      call search(line,'AO_ATOMINDEX',istat)
      if (istat.eq.1) then
         norbs = intlin(line,istat)
         if (istat.eq.0) goto 100
         if (norbs.gt.mxorb) then
             nsiz = norbs/256 
             lsiz = norbs - nsiz*256
             if (lsiz.gt.0) nsiz = nsiz + 1
             nsiz = nsiz*256
             call allorb(nsiz,0)
             call rewfil
             istats = -1
             return
         endif

         itel = 0
         iold = 0

         do while (.true.)

50          if (getlin(0).ne.1) goto 100
            if (index(line,'=').ne.0) goto 15

            do while (.true.)

                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.2) then
                   itel = itel + 1
                   if (itype.ne.iold) then
                       nfirst(itype) = itel
                       if (itype.gt.1) nlast(itype-1) = itel - 1
                   endif
                   iold = itype
                elseif (ktype.eq.0) then
                   goto 50
                else
                   print*,line
                   call inferr('integer expected',1)
                   goto 100
                endif

            end do

         end do

      else
         call inferr('No AO_ATOMINDEX card found',1)
         goto 100
      end if

15    nlast(itype) = itel

      call search(line,'AO_ZETA',istat)
      if (istat.eq.1) then
         n10 = norbs / 10
         if (norbs.gt.n10*10) n10 = n10 + 1

         itel = 0
         do j=1,n10

            if (getlin(0).ne.1) goto 100

            do i=1,10
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.3) then
                   itel = itel + 1
                   do k=1,natoms
                      n1 = nfirst(k)
                      n2 = nlast(k)
                      if (itel.eq.n1  .and.itel.le.n2) emus(k) = rtype
                      if (itel.eq.n1+1.and.itel.le.n2) emup(k) = rtype
                      if (itel.eq.n1+4.and.itel.le.n2) emud(k) = rtype
                   end do
                elseif (ktype.eq.0) then
                   goto 25
                else
                   print*,line
                   call inferr('real expected',1)
                   goto 100
                endif
            end do

         end do

      else
         call inferr('No AO_ZETA card found',1)
         goto 100
      end if

25    continue


      if (idebug.eq.1) then
         print*,'AO_ZETA: atom  zetas  zetap zetad'
         do i=1,natoms
            print*,i,emus(i),emup(i),emud(i)
         end do
      endif

      call search(line,'NUM_ELECTRONS=',istat)
      if (istat.eq.1) then
         i1 = index(line,'=')
         nelecs = dint(reada(line,i1+1,len(line)))
      else
         call inferr('No NUM_ELECTRONS card found',1)
         goto 100
      end if

      call search(line,'Final SCF results',istat)

      cf = 5.2917715d-11*1.601917d-19/3.33564d-30

      call search(line,'DIP_VEC',istat)
      if (istat.eq.1) then
         if (idebug.eq.1) call inferr('DIP_VEC card found',0)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (gnreal(dipo,3,.false.)) then
             ihsdp = 1
             do i=1,3
                dipo(i) = dipo(i)/cf
             end do
         endif

      else
         call inferr('No DIP_VEC card found',1)
         goto 100
      end if

      call search(line,'ATOM_X',istat)
      if (istat.eq.1) then
         if (idebug.eq.1) call inferr('ATOM_X card found',0)
      else
         call inferr('No ATOMS_X card found',1)
         iatfnd = 0
         goto 100
      end if
 
      call gtapnt(natoms,xyz,istat)

      if (idebug.eq.1) then
         print*,'ATOM_X'
         do i=1,natoms
             print*,(xyz(j,i),j=1,3)
         end do
      endif

      call search(line,'ATOM_CHARGES',istat)
      if (istat.eq.1) then
         if (idebug.eq.1) call inferr('ATOM_CHARGES card found',0)

         n10 = natoms / 10
         if (natoms.gt.n10*10) n10 = n10 + 1

         itel = 0
         iold = 0

         do j=1,n10

            if (getlin(0).ne.1) goto 100

            do i=1,10

                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.3) then
                   itel = itel + 1
                   q(itel) = rtype
                elseif (ktype.eq.0) then
                   goto 30
                else
                   print*,line
                   call inferr('ATOM_CHARGES: real expected',1)
                   goto 100
                endif

            end do

         end do

      else
         call inferr('No ATOM_CHARGES card found',1)
         goto 100
      end if

30    ihasq = 1

      linear = (norbs*(norbs+1))/2

      call search(line,'OVERLAP_MATRIX',istat)
      if (istat.eq.1) then
         if (getlin(0).ne.1) goto 100
         n10 = linear / 10
         if (linear.gt.n10*10) n10 = n10 + 1

         jtel = 0

         do j=1,n10

            if (getlin(0).ne.1) goto 100

            do i=1,10
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.3) then

                   jtel = jtel + 1
                   vectrs(jtel) = rtype

                elseif (ktype.eq.0) then
                   goto 35
                else
                   print*,line
                   call inferr('OVERLAP_MATRIX: real expected',1)
                   goto 100
                endif
            end do

         end do

      else
         call inferr('No OVERLAP_MATRIX card found',1)
         goto 100
      end if

35    continue

      call rsp(vectrs,norbs,norbs,psi,p)

      do i=1,norbs
         psi(i) = 1.d0/dsqrt(dabs(psi(i)))
      end do

      do i=1,norbs
         do j=1,i

            sum = 0.d0

            do k=1,norbs
               sum = sum + p(i+(k-1)*norbs)*psi(k)
     &                    *p(j+(k-1)*norbs)
            end do

            halfs((i-1)*norbs+j) = sum
            halfs((j-1)*norbs+i) = sum

         end do
      end do

      call search(line,'EIGENVECTORS[',istat)
      if (istat.eq.1) then
         if (index(line,'ALPHA').ne.0) iuhf = 1
         n1 = index(line,'EIGENVECTORS[')
         n2 = index(line,']')
         if (n1.gt.0.and.n2.gt.0) then
            n10 = dint(reada(line,n1+13,n2-1))
            n10 = n10 / 10
         else
            n10 = norbs*norbs / 10
            if (norbs*norbs.gt.n10*10) n10 = n10 + 1
         endif

         itel = 1
         jtel = 0

         do j=1,n10

            if (getlin(0).ne.1) goto 100

            do i=1,10
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.3) then

                   jtel = jtel + 1
                   if (jtel.gt.norbs) then
                      jtel = 1
                      itel = itel + 1
                   endif
                   p((itel-1)*mxorb+jtel) = rtype

                elseif (ktype.eq.0) then
                   goto 45
                else
                   print*,line
                   call inferr('EIGENVECTORS: real expected',1)
                   goto 100
                endif
            end do

         end do

      else
         call inferr('No EIGENVECTORS card found',1)
         goto 100
      end if

45    continue

      call search(line,'EIGENVALUES[',istat)
      if (istat.eq.1) then
         n1 = index(line,'EIGENVALUES[')
         n2 = index(line,']')
         if (n1.gt.0.and.n2.gt.0) then
            n10 = dint(reada(line,n1+12,n2-1))
            n10 = n10 / 10
         else
            n10 = norbs / 10
            if (norbs.gt.n10*10) n10 = n10 + 1
         endif

         itel = 0

         do j=1,n10

            if (getlin(0).ne.1) goto 100

            do i=1,10
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.3) then

                   itel = itel + 1
                   eiga(itel) = rtype

                elseif (ktype.eq.0) then
                   goto 55
                else
                   print*,line
                   call inferr('EIGENVALUES: real expected',1)
                   goto 100
                endif
            end do

         end do

      else
         call inferr('No EIGENVALUES card found',1)
         goto 100
      end if

55    continue

      do i=1,norbs
         focc(i) = 0.0d0
      end do
     
      call search(line,'_OCCUPANCIES[',istat)
      if (istat.eq.1) then
       
         if (index(line,'ALPHA').ne.0) then
            iuhf = 1
            nt = 40
         else
            nt = 10
         endif

         n1 = index(line,'_OCCUPANCIES[')
         n2 = index(line,']')
         if (n1.gt.0.and.n2.gt.0) then
            nc = dint(reada(line,n1+13,n2-1))
            nc = nc / 10
         else
            nc = norbs / nt
            if (norbs.gt.nc*nt) nc = nc + 1
         endif

         itel = 0

         do j=1,nc

            if (getlin(0).ne.1) goto 100

            do i=1,nt
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.3) then

                   itel = itel + 1
                   focc(itel) = rtype
                   if (itel.eq.norbs) goto 65

                elseif (ktype.eq.2) then

                   itel = itel + 1
                   focc(itel) = dble(itype)
                   if (itel.eq.norbs) goto 65

                elseif (ktype.eq.0) then
                   goto 65
                else
                   print*,line
                   call inferr('OCCUPANCIES: real expected',1)
                   goto 100
                endif
            end do

         end do

      else
         call inferr('No OCCUPANCIES card found',1)
         goto 100
      end if

65    nocc = 0
      do i=1,norbs
         if (focc(i).gt.0.0d0) nocc = nocc + 1
      end do
      ncols = norbs


      if (iuhf.eq.1) then
         call averab(averag,p,focc)
      else
         call averg(averag,p)
      endif

      call mulpxs(p,halfs,psi,vectrs)
    
      if (iuhf.eq.1) then
         call rewfil

         call search(line,'BETA_EIGENVECTORS[',istat)
         if (istat.eq.1) then
            n10 = norbs*norbs / 10
            if (norbs*norbs.gt.n10*10) n10 = n10 + 1
   
            itel = 1
            jtel = 0
   
            do j=1,n10
   
               if (getlin(0).ne.1) goto 100
   
               do i=1,10
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.3) then
   
                      jtel = jtel + 1
                      if (jtel.gt.norbs) then
                         jtel = 1
                         itel = itel + 1
                      endif
                      p((itel-1)*mxorb+jtel) = rtype
   
                   elseif (ktype.eq.0) then
                      goto 75
                   else
                      print*,line
                      call inferr('EIGENVECTORS: real expected',1)
                      goto 100
                   endif
               end do
   
            end do
   
         else
            call inferr('No EIGENVECTORS card found',1)
            goto 100
         end if
   
75       continue
   
         call search(line,'BETA_EIGENVALUES[',istat)
         if (istat.eq.1) then
            n10 = norbs / 10
            if (norbs.gt.n10*10) n10 = n10 + 1
   
            itel = 0
   
            do j=1,n10
   
               if (getlin(0).ne.1) goto 100
   
               do i=1,10
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.3) then
   
                      itel = itel + 1
                      eigb(itel) = rtype
   
                   elseif (ktype.eq.0) then
                      goto 85
                   else
                      print*,line
                      call inferr('EIGENVALUES: real expected',1)
                      goto 100
                   endif
               end do
   
            end do
   
         else
            call inferr('No EIGENVALUES card found',1)
            goto 100
         end if
   
85       continue
   
         do i=1,norbs
            focb(i) = 0.0d0
         end do
        
         call search(line,'BETA_MOLECULAR_ORBITAL_OCCUPANCIES[',istat)
         if (istat.eq.1) then
            if (index(line,'BETA').ne.0) then
               nt = 40
            else
               nt = 10
            endif

            nc = norbs / nt
            if (norbs.gt.nc*nt) nc = nc + 1
   
            itel = 0
   
            do j=1,nc
   
               if (getlin(0).ne.1) goto 100
   
               do i=1,nt
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.3) then
   
                      itel = itel + 1
                      focb(itel) = rtype
                      if (itel.eq.norbs) goto 95
   
                    elseif (ktype.eq.2) then

                      itel = itel + 1
                      focb(itel) = dble(itype)
                      if (itel.eq.norbs) goto 95

                   elseif (ktype.eq.0) then
                      goto 95
                   else
                      print*,line
                      call inferr('OCCUPANCIES: real expected',1)
                      goto 100
                   endif
               end do
   
            end do
   
         else
            call inferr('No OCCUPANCIES card found',1)
            goto 100
         end if
   
95       nocb = 0
         do i=1,norbs
            if (focb(i).gt.0.0d0) nocb = nocb + 1
         end do
         ncolb = norbs
   
         do i=1,norbs
            vectrb(i) = averag(i)
         end do

         call averab(averag,p,focb)

         do i=1,norbs
            averag(i) = averag(i) + vectrb(i)
         end do

         call mulpxs(p,halfs,psi,vectrb)

      endif

      call setcst

      call rewfil

      nepnts = 0
      do while (.true.)
         call search(line,'HEAT_OF_FORM_UPDATED',istat)
         if (istat.eq.1) then
            ieav = 1
            nepnts = nepnts + 1
            isav(nepnts) = 1
            epoints(nepnts) = rallin(line,istat)
            call search(line,'GRADIENT_UPDATED',istat)
            if (istat.eq.1) then
               ifmxav = 1
               formax(nepnts) = rallin(line,istat)
            endif
         else
            goto 110
         endif
      end do

110   if (nepnts.gt.0) then
         ngeoms = nepnts
         igcvav = 1
      endif

      if (idebug.eq.1) then
         print*,'HEAT OF FORMATION'
         do i=1,nepnts
            print*,epoints(i)
         end do
      endif

      if (idebug.eq.1) call inferr('AUX: MO data read',0)

      istats = 1
      return

100   continue

c
c      call rewmf
c
c      call srchmf(line,'[FREQ]',istat)
c
c      if (istat.eq.0) then
c
c         call inferr('No FREQ card found',0)
c
c      else
c         irtype = 4
c         if (idebug.eq.1) call inferr('FREQ card found',0)
c         if (iatfnd.ne.1) then
c             call getfra(istat)
c             call getint(istat)
c         endif
c      endif
c

      return
      end


      integer function intlin(line,istat)
      implicit double precision (a-h,o-z)
      character*137 line

      istat = 0
      intlin = 0
      i1 = index(line,'[')
      i2 = index(line,']')
      if (i1.eq.0.or.i2.eq.0) return

      intlin = reada(line,i1+1,i2-1)
      istat = 1

      return
      end

      double precision function rallin(line,istat)
      implicit double precision (a-h,o-z)
      character*137 line

      istat = 0
      rallin = 0
      i1 = index(line,'=')
      if (i1.eq.0.or.i2.eq.0) return

      rallin = reada(line,i1+1,len(line))
      istat = 1

      return
      end

      subroutine rsp(a,n,nvect,root,vect)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(*), root(n), vect(n,n)
      dimension rwork(8*20000), iwork(20000)

      call evvrsp(-1,n,nvect,n*n, n,a,rwork,iwork,root,vect,0,i)

      return
      end

      subroutine evvrsp (msgfl,n,nvect,lena,nv,a,b,ind,root,vect,iorder
     &                 ,ierr)

c       finds   (all) eigenvalues    and    (some or all) eigenvectors
c       of a real symmetric packed matrix.
c
c       the method as presented in this routine consists of four steps:
c       first, the input matrix is reduced to tridiagonal form by the
c       householder technique (orthogonal similarity transformations).
c       second, the roots are located using the rational ql method.
c       third, the vectors of the tridiagonal form are evaluated by the
c       inverse iteration technique.  vectors for degenerate or near-
c       degenerate roots are forced to be orthogonal.
c       fourth, the tridiagonal vectors are rotated to vectors of the
c       original array.
c
c    on entry -
c       msgfl  - integer (logical unit no.)
c                file where error messages will be printed.
c                if msgfl is 0, error messages will be printed on lu 6.
c                if msgfl is negative, no error messages printed.
c       n      - integer
c                order of matrix a.
c       nvect  - integer
c                number of vectors desired.  0 .le. nvect .le. n.
c       lena   - integer
c                dimension of  a  in calling routine.  must not be less
c                than (n*n+n)/2.
c       nv     - integer
c                row dimension of vect in calling routine.   n .le. nv.
c       a      - working precision real (lena)
c                input matrix, rows of the lower triangle packed into
c                linear array of dimension n*(n+1)/2.  the packed order
c                is a(1,1), a(2,1), a(2,2), a(3,1), a(3,2), ...
c       b      - working precision real (n,8)
c                scratch array, 8*n elements
c       ind    - integer (n)
c                scratch array of length n.
c       iorder - integer
c                root ordering flag.
c                = 0, roots will be put in ascending order.
c                = 2, roots will be put in descending order.
c
c    on exit -
c       a      - destoryed.  now holds reflection operators.
c       root   - working precision real (n)
c                all eigenvalues in ascending or descending order.
c                  if iorder = 0, root(1) .le. ... .le. root(n)
c                  if iorder = 2, root(1) .ge. ... .ge. root(n)
c       vect   - working precision real (nv,nvect)
c                eigenvectors for root(1), ..., root(nvect).
c       ierr   - integer
c                = 0 if no error detected,
c                = k if iteration for k-th eigenvalue failed,
c                = -k if iteration for k-th eigenvector failed.
c                (failures should be very rare.  contact c. moler.)
c
c

      implicit double precision (a-h,o-z), integer (i-n)
      logical goparr,dskwrk,maswrk

      common /par/ me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      dimension a(lena),b(n,8),root(n),vect(nv,*),ind(n)

  900 format(26h0*** evvrsp parameters ***/
     +       14h ***      n = ,i8,4h ***/
     +       14h ***  nvect = ,i8,4h ***/
     +       14h ***   lena = ,i8,4h ***/
     +       14h ***     nv = ,i8,4h ***/
     +       14h *** iorder = ,i8,4h ***/
     +       14h ***   ierr = ,i8,4h ***)
  901 format(37h value of lena is less than (n*n+n)/2)
  902 format(39h eqlrat has failed to converge for root,i5)
  903 format(18h nv is less than n)
  904 format(41h einvit has failed to converge for vector,i5)
  905 format(51h value of iorder must be 0 (smallest root first) or
     *      ,23h 2 (largest root first))
  906 format(' value of n is less than or equal zero')

      lmsgfl = msgfl

      if (msgfl.eq.0) lmsgfl=6

      ierr = n - 1
      if (n.le.0) goto 800
      ierr = n + 1

      if ((n*n+n)/2.gt.lena) goto 810

c        reduce real symmetric matrix a to tridiagonal form

      call etred3(n,lena,a,b(1,1),b(1,2),b(1,3))

c        find all eigenvalues of tridiagonal matrix

      call eqlrat(n,b(1,1),b(1,2),b(1,3),root,ind,ierr,b(1,4))

      if (ierr.ne.0) goto 820

c         check the desired order of the eigenvalues

      b(1,3) = float(iorder)

      if (iorder .eq. 0) go to 300

         if (iorder .ne. 2) go to 850

c         order roots in descending order (largest first)...
c         turn root and ind arrays end for end

         do i = 1, n/2
            j = n+1-i
            t = root(i)
            root(i) = root(j)
            root(j) = t
            l = ind(i)
            ind(i) = ind(j)
            ind(j) = l
         end do

c           find i and j marking the start and end of a sequence
c           of degenerate roots

         i=0
         do while(.true.)
            i = i+1
            if (i.gt.n) goto 300
            do j=i,n
               if (root(j).ne.root(i)) goto 240
            end do
            j = n + 1
  240       continue
            j = j - 1
            if (j.ne.i) then

c                    turn around ind between i and j

               jsv = j
               klim = (j-i+1)/2
               do k=1,klim
                  l = ind(j)
                  ind(j) = ind(i)
                  ind(i) = l
                  i = i+1
                  j = j-1
               end do
               i = jsv
            endif

         end do
c
  300 continue
c
      if (nvect.le.0) return
      if (nv.lt.n) goto 830

c        find eigenvectors of tri-diagonal matrix via inverse iteration

      ierr = lmsgfl
      call einvit(nv,n,b(1,1),b(1,2),b(1,3),nvect,root,ind,
     &            vect,ierr,b(1,4),b(1,5),b(1,6),b(1,7),b(1,8))
      if (ierr.ne.0) goto 840

c        find eigenvectors of symmetric matrix via back transformation

  400 continue
      call etrbk3(nv,n,lena,a,nvect,vect)
      return

c        error message section

  800 if (lmsgfl .lt. 0) return
      if (maswrk) write(lmsgfl,906)
      go to 890

  810 if (lmsgfl .lt. 0) return
      if (maswrk) write(lmsgfl,901)
      go to 890

  820 if (lmsgfl .lt. 0) return
      if (maswrk) write(lmsgfl,902) ierr
      go to 890

  830 if (lmsgfl .lt. 0) return
      if (maswrk) write(lmsgfl,903)
      go to 890

  840 continue
      if ((lmsgfl .gt. 0).and.(maswrk)) write(lmsgfl,904) -ierr
      go to 400

  850 ierr=-1
      if (lmsgfl .lt. 0) return
      if (maswrk) write(lmsgfl,905)
      go to 890

  890 continue
      if (maswrk) write(lmsgfl,900) n,nvect,lena,nv,iorder,ierr
      return
      end

      subroutine freda(l,d,a,e)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(*),d(l),e(l)

      jk = 1

c     .......... form reduced a ..........

      do j = 1, l
         f = d(j)
         g = e(j)
 
         do k = 1, j
            a(jk) = a(jk) - f * e(k) - g * d(k)
            jk = jk + 1
         end do
 
      end do

      return
      end

      subroutine einvit(nm,n,d,e,e2,m,w,ind,z,ierr,rv1,rv2,rv3,rv4,rv6)

c       this routine finds those eigenvectors of a tridiagonal
c       symmetric matrix corresponding to specified eigenvalues.
c
c    method -
c       inverse iteration.
c
c    on entry -
c       nm     - integer
c                must be set to the row dimension of two-dimensional
c                array parameters as declared in the calling routine
c                dimension statement.
c       n      - integer
c       d      - w.p. real (n)
c                contains the diagonal elements of the input matrix.
c       e      - w.p. real (n)
c                contains the subdiagonal elements of the input matrix
c                in its last n-1 positions.  e(1) is arbitrary.
c       e2     - w.p. real (n)
c                contains the squares of corresponding elements of e,
c                with zeros corresponding to negligible elements of e.
c                e(i) is considered negligible if it is not larger than
c                the product of the relative machine precision and the
c                sum of the magnitudes of d(i) and d(i-1).  e2(1) must
c                contain 0.0 if the eigenvalues are in ascending order,
c                or 2.0 if the eigenvalues are in descending order.
c                if tqlrat, bisect, tridib, or imtqlv
c                has been used to find the eigenvalues, their
c                output e2 array is exactly what is expected here.
c       m      - integer
c                the number of specified eigenvectors.
c       w      - w.p. real (m)
c                contains the m eigenvalues in ascending
c                or descending order.
c       ind    - integer (m)
c                contains in first m positions the submatrix indices
c                associated with the corresponding eigenvalues in w --
c                1 for eigenvalues belonging to the first submatrix
c                from the top, 2 for those belonging to the second
c                submatrix, etc.
c       ierr   - integer (logical unit number)
c                logical unit for error messages
c
c    on exit -
c       all input arrays are unaltered.
c       z      - w.p. real (nm,m)
c                contains the associated set of orthonormal
c                eigenvectors. any vector which which fails to converge
c                is left as is (but normalized) when iterating stopped.
c       ierr   - integer
c                set to
c                zero    for normal return,
c                -r      if the eigenvector corresponding to the r-th
c                        eigenvalue fails to converge in 5 iterations.
c                        (only last failure to converge is reported)
c
c       rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
c
c       rv1    - w.p. real (n)
c                diagonal elements of u from lu decomposition
c       rv2    - w.p. real (n)
c                super(1)-diagonal elements of u from lu decomposition
c       rv3    - w.p. real (n)
c                super(2)-diagonal elements of u from lu decomposition
c       rv4    - w.p. real (n)
c                elements defining l in lu decomposition
c       rv6    - w.p. real (n)
c                approximate eigenvector
c

      implicit double precision (a-h,o-z), integer (i-n)

      logical convgd,goparr,dskwrk,maswrk
      integer group,p,q,r,s,submat,tag
      double precision norm

      common /par/ me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk

      dimension d(n),e(n+1),e2(n),w(m),z(nm,m),ind(m),
     &          rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)

  001 format(' eigenvector routine einvit did not converge for vector'
     *      ,i5,'.  norm =',1p,e10.2,' performance index =',e10.2/
     *      ' (an error halt will occur if the pi is greater than 100)')

c-----------------------------------------------------------------------

      zero   = 0.0d0
      one    = 1.0d0
      grptol = 0.001d0
      epscal = 0.5d0
      hundrd = 100.0d0
      ten    = 10.0d0

      luemsg = ierr
      ierr   = 0
      x0     = zero
      uk     = zero
      norm   = zero
      eps2   = zero
      eps3   = zero
      eps4   = zero
      group  = 0
      tag    = 0
      order  = one - e2(1)
      q      = 0

      do submat=1,n

         p = q + 1

c        .......... establish and process next submatrix ..........

         do q=p,n-1
            if (e2(q+1).eq.zero) goto 140
         end do

         q = n
  140    continue


c        .......... find vectors by inverse iteration ..........

         tag = tag + 1
         anorm = zero
         s = 0

         do r=1,m

            if (ind(r).eq.tag) then

               its = 1
               x1 = w(r)

               if (s.eq.0) then

c        .......... check for isolated root ..........

                  xu = one
                  if (p.eq.q) then
                     rv6(p) = one
                     convgd = .true.
                     goto 860
                  endif

                  norm = abs(d(p))

                  do i=p+1,q
                     norm = max( norm, abs(d(i)) + abs(e(i)) )
                  end do

c        .......... eps2 is the criterion for grouping,
c                   eps3 replaces zero pivots and equal
c                   roots are modified by eps3,
c                   eps4 is taken very small to avoid overflow .........

                  eps2 = grptol * norm
                  eps3 = epscal * epslon(norm)
                  uk = q - p + 1
                  eps4 = uk * eps3
                  uk = eps4 / sqrt(uk)
                  s = p
                  group = 0

               else

c        .......... look for close or coincident roots ..........

                  if (abs(x1-x0).ge.eps2) then

c                    roots are seperate

                     group = 0

                  else

c                    roots are close

                     group = group + 1

                     if (order*(x1 - x0).le.eps3) 
     &                   x1 = x0 + order * eps3

                  endif

c        .......... elimination with interchanges and
c                   initialization of vector ..........

               endif

               u = d(p) - x1
               v = e(p+1)
               rv6(p) = uk

               do i=p+1,q

                  rv6(i) = uk

                  if (abs(e(i)).gt.abs(u)) then

c                    exchange rows before elimination

c                     *** warning -- a divide check may occur here if
c                      e2 array has not been specified correctly .......

                     xu = u / e(i)
                     rv4(i) = xu
                     rv1(i-1) = e(i)
                     rv2(i-1) = d(i) - x1
                     rv3(i-1) = e(i+1)
                     u = v - xu * rv2(i-1)
                     v = -xu * rv3(i-1)

                  else

c                       straight elimination

                     xu = e(i) / u
                     rv4(i) = xu
                     rv1(i-1) = u
                     rv2(i-1) = v
                     rv3(i-1) = zero
                     u = d(i) - x1 - xu * v
                     v = e(i+1)

                  end if

               end do

               if (abs(u).le.eps3) u = eps3

               rv1(q) = u
               rv2(q) = zero
               rv3(q) = zero

c                 do inverse iterations

               convgd = .false.

               do its=1,5

                  if (its.ne.1) then
 
c                    .......... forward substitution ..........
 
                     if (norm.eq.zero) then

                        rv6(s) = eps4
                        s = s + 1

                        if (s.gt.q) s = p

                     else

                        xu = eps4 / norm
                        call dscal(q-p+1,xu,rv6(p),1)

                     endif

c                     ... elimination operations on next vector

                     do i=p+1,q

                        u = rv6(i)

c                         if rv1(i-1) .eq. e(i), a row interchange
c                         was performed earlier in the
c                         triangularization process ..........

                        if (rv1(i-1).eq.e(i)) then

                           u = rv6(i-1)
                           rv6(i-1) = rv6(i)

                        else

                           u = rv6(i)

                        endif

                        rv6(i) = u - rv4(i) * rv6(i-1)

                     end do

                  endif

c           .......... back substitution

                  rv6(q) = rv6(q) / rv1(q)
                  v = u
                  u = rv6(q)
                  norm = abs(u)

                  do i=q-1,p,-1

                     rv6(i) = (rv6(i) - u*rv2(i) - v*rv3(i)) / rv1(i)
                     v = u
                     u = rv6(i)
                     norm = norm + abs(u)

                  end do

                  if (group.ne.0) then

c                 ....... orthogonalize with respect to previous
c                         members of group ..........

                     j = r

                     do jj=1,group

                        do while(ind(j).ne.tag)
                           j = j - 1
                        end do

                        call daxpy(q-p+1,-ddot(q-p+1,rv6(p),1,z(p,j),1),
     &                             z(p,j),1,rv6(p),1)
                     end do

                     norm = dasum(q-p+1,rv6(p),1)

                  endif

                  if (convgd) goto 840
                  if (norm.ge.one) convgd = .true.

c                 end its loop
               end do

c        .......... normalize so that sum of squares is
c                   1 and expand to full order ..........

  840          continue

               xu = one / dnrm2(q-p+1,rv6(p),1)

  860          continue

               do i=1,p-1
                  z(i,r) = zero
               end do

               do i=p,q
                  z(i,r) = rv6(i) * xu
               end do

               do i=q+1,n
                  z(i,r) = zero
               end do

               if (.not.convgd) then
   
                  rho = estpi1(q-p+1,x1,d(p),e(p),z(p,r),anorm)

                  if (rho.ge.ten .and. luemsg.gt.0 .and. maswrk)
     &                write(luemsg,001) r,norm,rho

c               *** set error -- non-converged eigenvector ..........

                  if (rho.gt.hundrd) ierr = -r

               endif

               x0 = x1

c              endif  (ind(r).eq.tag)
            endif

c           end do r=1,m
         end do

         if (q.eq.n) return

c        end do submat=1,n
      end do

      return
      end

      subroutine elau(hinv,l,d,a,e)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(*),d(l),e(l)


      zero = 0.0d+00
      half = 0.5d+00

      jl = l
      e(1) = a(1) * d(1)
      jk = 2

      do j=2,jl

         f = d(j)
         g = zero
         jm1 = j - 1

         do k=1,jm1

            g = g + a(jk) * d(k)
            e(k) = e(k) + a(jk) * f
            jk = jk + 1

         end do

         e(j) = g + a(jk) * f
         jk = jk + 1

      end do

c     .......... form p ..........

      f = zero

      do j=1,l

         e(j) = e(j) * hinv
         f = f + e(j) * d(j)

      end do

c     .......... form q ..........

      hh = f * half * hinv

      do j=1,l
         e(j) = e(j) - hh * d(j)
      end do

      return
      end

c
c
c    authors -
c       this routine was taken from eispack edition 3 dated 4/6/83
c       this version is by s. t. elbert, ames laboratory-usdoe nov 1986
c
c    purpose -
c       estimate unit roundoff in quantities of size x.
c
c    on entry -
c       x      - working precision real
c                values to find epslon for
c
c    on exit -
c       epslon - working precision real
c                smallest positive value such that x+epslon .ne. zero
c
c    qualifications -
c       this routine should perform properly on all systems
c       satisfying the following two assumptions,
c          1.  the base used in representing floating point
c              numbers is not a power of three.
c          2.  the quantity  a  in statement 10 is represented to
c              the accuracy used in floating point variables
c              that are stored in memory.
c       the statement number 10 and the go to 10 are intended to
c       force optimizing compilers to generate code satisfying
c       assumption 2.
c       under these assumptions, it should be true that,
c              a  is not exactly equal to four-thirds,
c              b  has a zero for its last bit or digit,
c              c  is not exactly equal to one,
c              eps  measures the separation of 1.0 from
c                   the next larger floating point number.
c       the developers of eispack would appreciate being informed
c       about any systems where these assumptions do not hold.
c
c    differences from eispack 3 -
c       use is made of parameter statements and intrinsic functions
c       --no executeable code changes--
c
c    external routines - none
c    intrinsic functions - abs
c
c    note -
c       questions and comments concerning eispack should be directed to
c       b. s. garbow, applied math. division, argonne national lab.
c
c
      double precision function epslon(x)
      implicit double precision (a-h,o-z), integer (i-n)

      a = 4.0d0/3.0d0
      eps = 0.0d0

      do while(eps.eq.0.0d0)
         b = a - 1.0d0
         c = b + b + b
         eps = abs(c - 1.0d0)
      end do

      epslon = eps*abs(x)

      return
      end

      subroutine eqlrat(n,diag,e,e2in,d,ind,ierr,e2)

c       this routine finds the eigenvalues of a symmetric
c       tridiagonal matrix
c
c    on entry -
c       n      - integer
c                the order of the matrix.
c       d      - w.p. real (n)
c                contains the diagonal elements of the input matrix.
c       e2     - w.p. real (n)
c                contains the squares of the subdiagonal elements of
c                the input matrix in its last n-1 positions.
c                e2(1) is arbitrary.
c
c     on exit -
c       d      - w.p. real (n)
c                contains the eigenvalues in ascending order.  if an
c                error exit is made, the eigenvalues are correct and
c                ordered for indices 1,2,...ierr-1, but may not be
c                the smallest eigenvalues.
c       e2     - w.p. real (n)
c                destroyed.
c       ierr   - integer
c                set to
c                zero       for normal return,
c                j          if the j-th eigenvalue has not been
c                           determined after 30 iterations.
c

      implicit double precision (a-h,o-z), integer (i-n)
      dimension d(n),e(n),e2(n),diag(n),e2in(n),ind(n)

      scale  = 1.0d0/64.0d0
      ierr   = 0
      d(1)   = diag(1)
      ind(1) = 1
      k      = 0
      itag   = 0

      if (n.eq.1) return

      do i=2,n

         d(i)    = diag(i)
         e2(i-1) = e2in(i)

      end do

      f = 0.0d0
      t = 0.0d0
      b = epslon(1.0d0)
      c = b*b
      b = b*scale
      e2(n) = 0.0d0

      do l=1,n

         h = abs(d(l)) + abs(e(l))

         if (t.lt.h) then

            t = h
            b = epslon(t)
            c = b*b
            b = b*scale

         endif

c     .......... look for small squared sub-diagonal element ..........

         m = l 

         do while(e2(m).gt.c)
            m = m + 1
         end do

c     .......... e2(n) is always zero, so there is an exit
c                from the loop ..........

         if (m.gt.k) then

            if (m.ne.n) e2in(m+1) = 0.0d0
            k = m
            itag = itag + 1

         endif

         if (m.eq.l) goto 100

c           iterate

         do j=1,30

c              .......... form shift ..........

            l1 = l + 1
            s = sqrt(e2(l))
            g = d(l)
            p = (d(l1) - g) / (2.0d0*s)
            r = sqrt(p*p+1.0d0)
            d(l) = s / (p + sign(r,p))
            h = g - d(l)

            do i=l1,n
               d(i) = d(i) - h
            end do

            f = f + h

c              .......... rational ql transformation ..........

            g = d(m) + b
            h = g
            s = 0.0d0

            do i=m-1,l,-1

               p = g * h
               r = p + e2(i)
               e2(i+1) = s*r
               s = e2(i) / r
               d(i+1) = h + s*(h + d(i))
               g = d(i) - e2(i) / g   + b
               h = g*p / r

            end do

            e2(l) = s*g
            d(l)  = h

c              .......... guard against underflow in convergence test

            if (h.eq.0.0d0 .or. abs(e2(l)).le.abs(c/h) ) goto 100
            e2(l) = h*e2(l)
            if (e2(l).eq.0.0d0) goto 100

         end do

c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........

         ierr = l
         return


c     converged

100      p = d(l) + f

c           .......... order eigenvalues ..........

         i = 1

         if (l.ne.1) then

            if (p.ge.d(1)) then

               i = l - 1

c           .......... loop to find ordered position

               do while(p.lt.d(i))
                  i = i - 1
               end do

               i = i + 1

            endif

            if (i.ne.l) then

               do ii=l,i+1,-1
                  d(ii) = d(ii-1)
                  ind(ii) = ind(ii-1)
               end do

            endif


         endif

         d(i) = p
         ind(i) = itag

      end do

      return
      end

      subroutine etrbk3(nm,n,nv,a,m,z)

c       this routine forms the eigenvectors of a real symmetric
c       matrix by back transforming those of the corresponding
c       symmetric tridiagonal matrix determined by  etred3.
c
c    method -
c       the calculation is carried out by forming the matrix product
c          q*z
c       where  q  is a product of the orthogonal symmetric matrices
c                q = prod(i)[1 - u(i)*.transpose.u(i)*h(i)]
c       u  is the augmented sub-diagonal rows of  a  and
c       z  is the set of eigenvectors of the tridiagonal
c       matrix  f  which was formed from the original symmetric
c       matrix  c  by the similarity transformation
c                f = q(transpose) c q

      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(nv),z(nm,m)

      if (m.eq.0 .or. n.le.2) return

      ii = 3

      do i=3,n

         iz = ii + 1
         ii = ii + i
         h  = a(ii)

         if (h.ne.0.0d0) then

            im1 = i - 1

            do j=1,m
               s = -( ddot(im1,a(iz),1,z(1,j),1) * h) * h
               call daxpy(im1,s,a(iz),1,z(1,j),1)
            end do

         endif

      end do

      return
      end

      subroutine etred3(n,nv,a,d,e,e2)

c       this routine reduces a real symmetric (packed) matrix, stored
c       as a one-dimensional array, to a symmetric tridiagonal matrix
c       using orthogonal similarity transformations, preserving the
c       information about the transformations in  a.
c
c    method -
c       the tridiagonal reduction is performed in the following way.
c       starting with j=n, the elements in the j-th row to the
c       left of the diagonal are first scaled, to avoid possible
c       underflow in the transformation that might result in severe
c       departure from orthogonality.  the sum of squares  sigma  of
c       these scaled elements is next formed.  then, a vector  u  and
c       a scalar
c                      h = u(transpose) * u / 2
c       define a reflection operator
c                      p = i - u * u(transpose) / h
c       which is orthogonal and symmetric and for which the
c       similiarity transformation  pap  eliminates the elements in
c       the j-th row of  a  to the left of the subdiagonal and the
c       symmetrical elements in the j-th column.
c
c       the non-zero components of  u  are the elements of the j-th
c       row to the left of the diagonal with the last of them
c       augmented by the square root of  sigma  prefixed by the sign
c       of the subdiagonal element.  by storing the transformed sub-
c       diagonal element in  e(j)  and not overwriting the row
c       elements eliminated in the transformation, full information
c       about  p  is save for later use in  etrbk3.
c
c       the transformation sets  e2(j)  equal to  sigma  and  e(j)
c       equal to the square root of  sigma  prefixed by the sign
c       of the replaced subdiagonal element.
c
c       the above steps are repeated on further rows of the
c       transformed  a  in reverse order until  a  is reduced to tri-
c       diagonal form, that is, repeated for  j = n-1,n-2,...,3.
c
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(nv+2),d(n+1),e(n+1),e2(n+1)

      if (n.le.2) goto 100

      iz0 = (n*n+n)/2
      aiimax = abs(a(iz0))

      do i=n,3,-1

         l = i - 1
         iia = iz0
         iz0 = iz0 - i
         aiimax = max(aiimax, abs(a(iia)))
         scale = dasum(l,a(iz0+1),1)

         if (scale.eq.abs(a(iia-1)) .or. aiimax+scale.eq.aiimax) then

c           this row is already in tri-diagonal form

            d(i) = a(iia)

            if (aiimax+d(i).eq.aiimax) d(i) = 0.0d0

            e(i) = a(iia-1)

            if (aiimax+e(i).eq.aiimax) e(i) = 0.0d0

            e2(i) = e(i)*e(i)
            a(iia) = 0.0d0

         else

            scalei = 1.0d0 / scale
            call dscal(l,scalei,a(iz0+1),1)
            hroot = dnrm2(l,a(iz0+1),1)

            f = a(iz0+l)
            g = -sign(hroot,f)
            e(i) = scale * g
            e2(i) = e(i)*e(i)
            h = hroot*hroot - f * g
            a(iz0+l) = f - g
            d(i) = a(iia)
            a(iia) = 1.0d0 / sqrt(h)

c           .......... form p then q in e(1:l) ..........

            call elau(1.0d0/h,l,a(iz0+1),a,e)

c           .......... form reduced a ..........

            call freda(l,a(iz0+1),a,e)

         endif

      end do

100   e(1) = 1.0d0
      e2(1)= 1.0d0
      d(1) = a(1)
      e(2) = a(2)
      e2(2)= a(2)*a(2)
      d(2) = a(3)

      return
      end

      function estpi1(n,eval,d,e,x,anorm)

c       evaluate symmetric tridiagonal matrix performance index
c       for 1 eigenvector

c    method -
c       this routine forms the 1-norm of the residual matrix a*x-x*eval
c       where  a  is a symmetric tridiagonal matrix stored
c       in the diagonal (d) and sub-diagonal (e) vectors, eval is the
c       eigenvalue of an eigenvector of  a,  namely  x.
c       this norm is scaled by machine accuracy for the problem size.
c       all norms appearing in the comments below are 1-norms.
c
c    on entry -
c       n      - integer
c                the order of the matrix  a.
c       eval   - w.p. real
c                the eigenvalue corresponding to vector  x.
c       d      - w.p. real (n)
c                the diagonal vector of  a.
c       e      - w.p. real (n)
c                the sub-diagonal vector of  a.
c       x      - w.p. real (n)
c                an eigenvector of  a.
c       anorm  - w.p. real
c                the norm of  a  if it has been previously computed.
c
c    on exit -
c       anorm  - w.p. real
c                the norm of  a, computed if initially zero.
c       estpi1 - w.p. real
c          !!a*x-x*eval!! / (epslon(10*n)*!!a!!*!!x!!);
c          where epslon(x) is the smallest number such that
c             x + epslon(x) .ne. x
c
c          estpi1 .lt. 1 == satisfactory performance
c                 .ge. 1 and .le. 100 == marginal performance
c                 .gt. 100 == poor performance
c          (see lect. notes in comp. sci. vol.6 pp 124-125)
c
      implicit double precision (a-h,o-z), integer (i-n)
      dimension d(n), e(n), x(n)

      estpi1 = 0.0d0

      if (n.le.1) return

      size = 10 * n

      if (anorm.eq.0.0d0) then

c        compute norm of  a

         anorm = max( abs(d(1)) + abs(e(2))
     &               ,abs(d(n)) + abs(e(n)))
         do i=2,n-1
            anorm = max( anorm, abs(e(i))+abs(d(i))+abs(e(i+1)))
         end do

         if (anorm.eq.0.0d0) anorm = 1.0d0

      endif

c     compute norms of residual and eigenvector

      xnorm = abs(x(1)) + abs(x(n))
      rnorm = abs( (d(1)-eval)*x(1) + e(2)*x(2))
     &       +abs( (d(n)-eval)*x(n) + e(n)*x(n-1))

      do i=2,n-1
         xnorm = xnorm + abs(x(i))
         rnorm = rnorm + abs(e(i)*x(i-1) + (d(i)-eval)*x(i)
     &                                     + e(i+1)*x(i+1))
      end do

      estpi1 = rnorm / (epslon(size)*anorm*xnorm)

      return
      end

      subroutine prtxyz
      parameter (numatm=2000)
      implicit double precision (a-h,o-z), integer (i-n)
      common /moldat/ natoms,norbs,nelecs,nat(numatm) 
      common /coord / xyz(3,numatm)
         print*,'test'
         do i=1,natoms
             print*,(xyz(j,i),j=1,3)
         end do

      return
      end

      subroutine gtapnt(natoms,xyz,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      character*137 line,str
      integer getlin
      common /curlin/ line
      dimension xyz(3,*)

      istat = 0
      do i=1,natoms

          if (getlin(0).ne.1) goto 100

          do j=1,3
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                xyz(j,i) = rtype
             else
                print*,line
                call inferr('ATOM_X: real expected',1)
                goto 100
             endif
          end do

      end do

      istat = 1

100   return
      end

      subroutine mult(c,s,vecs,n)
      implicit double precision (a-h,o-z)
      dimension c(n,*), s(n,*), vecs(n,*)

c     operation:-
c                                vecs=back-transformed eigenvectors
c     vecs  =  c*s               c   =un-back-transformed vectors
c                                s   =1/sqrt(overlap matrix)

      do i=1,n
         do j=1,n

            sum=0.d0

            do k=1,n
               sum = sum + c(k,i)*s(j,k)
            end do

            vecs(j,i) = sum

         end do
      end do

      return
      end

      subroutine averab(averag,p,focc)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms,norbs,nelecs,nat(numatm)
      common /mopac/  nfirst(numatm),nlast(numatm),
     &                npq(numatm),pqn(numatm),emus(numatm),
     &                emup(numatm),emud(numatm),consts(numatm),
     &                constp(numatm),constd(numatm),npqref(54)
      dimension averag(*),p(*),focc(*)

      do i=1,norbs
        x = 0.d0
        do j=1,norbs
           x = x + p((j-1)*mxorb+i)**2*focc(j)*0.5d0
        end do
        averag(i) = x*2
      end do

c
c loop to calculate spherical-average atomic orbital occupancy
c
      do i=1,natoms

           il = nfirst(i)
           iu = nlast(i)
           ir = iu - il + 1

           goto (120,130,130,130,140,140,140,140,140),ir
140            x = 0.d0

               do j=1,5
                   ji = j + il + 3
                   x = x + averag(ji)
               end do

               x = x*0.2d0

               do j=1,5
                   ji = j + il + 3
                   averag(ji) = x
               end do

130            x = 0.d0

               do j=1,3
                   ji = j + il
                   x = x + averag(ji)
               end do

               x = x*0.333333d0

               do j=1,3
                   ji = j + il
                   averag(ji) = x
               end do

120       continue
      end do

      return
      end

