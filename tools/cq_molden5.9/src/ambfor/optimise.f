      subroutine optimisd(emin,gtol,nsd,par,par1,g,g1,zr,y,yt,pt,s,
     &                    iopt,ityp)
      implicit real (a-h,o-z), integer (i-n)

c     BFGS method requires stores of Inverse Hessian, 
c                 only for small systems, for big systems:
c     conjugate gradient energy minimization

      parameter (maxbfgs=1000)
      parameter (mxunch=25)
      common /athlp/   iatoms, mxnat
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /opts/ idebug,imon,iarc,ilog,iout,osingl,dolbfgs,oqscal
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw
      logical reset, osingl, dolbfgs,oqscal
      integer*2 ityp
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /cell/ xa,ya,yb,za,zb,zc
      common /convg/ icyc

      dimension par(*),g(*),par1(*),zr(*),g1(*),
     &          y(*),pt(*),yt(*),s(*),iopt(*),ityp(*)
      dimension hessin(maxbfgs*maxbfgs),idx(maxbfgs),
     &          dg(maxbfgs),hdg(maxbfgs),elast(mxunch)
      
c iwrite: write structure file every iwrite iterations
c
c nsd   : maximum number of iterations
c
c deprls: When the gradient-norm is satisfied, test if energy does not change
c         more than DEPRLS upon making a new pair list
c
c derest: When the energy changes more than DEREST upon a making a new
c         pair list, the algoritm is restarted with the gradient vector.
c
c STEP0:  initial step-size in  line-minimization (bv. 1.0E-5)
c         Choose such that NITI (number of iterations in the 
c         line-minimization in the first cycle becomes 1 or 2


      iwrite = 1
      imeth  = 0
      derest = 0.2e0      
      deprls = 0.05e0
      step0  = 0.05e0
      n      = iatoms*3

      ilast = 0
      do i=1,mxunch
         elast(i) = 0.0e0
      end do

      write(iun5,*) 'Gradient tolerance: ',gtol
      write(iun5,*) 'Maximum iterations: ',nsd

      nflx = 0
      do i=1,iatoms
         if (iopt(i).eq.1) then
             if (nflx*3.le.maxbfgs) then
                nflx = nflx + 1
                idx((nflx-1)*3+1) = (i-1)*3 + 1
                idx((nflx-1)*3+2) = (i-1)*3 + 2
                idx((nflx-1)*3+3) = (i-1)*3 + 3
             endif
         endif
      end do
      nflx3 = nflx*3

      if (cell) then
         par(n+1) = xa
         par(n+2) = ya
         par(n+3) = yb
         par(n+4) = za
         par(n+5) = zb
         par(n+6) = zc

         nflx = nflx + 1

         idx((nflx-1)*3+1) = n + 1
         idx((nflx-1)*3+2) = n + 2
         idx((nflx-1)*3+3) = n + 3
         idx((nflx-1)*3+4) = n + 4
         idx((nflx-1)*3+5) = n + 5
         idx((nflx-1)*3+6) = n + 6

         nflx = nflx + 1
         nflx3 = nflx*3

         n = n + 6

      endif

      if (nflx3.gt.maxbfgs.or.cell) then
          imeth = 1
          write(iun5,'(a,i6,a)') 'Number of atoms times 3: ',nflx3,
     &                           ' > MAXBFGS using PB     '
          write(iun5,*) 'Powell-Beale conjugate gradients           '
      else 
          write(iun5,*) 'BFGS Variable Metric                      '

      endif

      nrstrt = nsd
      if (imeth.eq.1) nrstrt = n

      d0   = step0
      ncyc = 0
      icyc = 0
      npow = 0

c     Restart search direction to the steepest descent direction

      call enegrd(eold,par,g)
      if (osingl) then
          write(iun5,'(a,f12.3)') 'Estat=',eold
          return
      endif

      call resdir(istat,imeth,n,nflx3,d0,gtol,eold,g,
     &                  npl,step,d1g1,zr,hessin)
      if (istat.eq.0) then
         emin = eold
         return
      endif
      
c     BEGIN LOOP

      do while(.true.)

         tya = 0.0
         do i=1,n
            tya = tya - zr(i)*g(i)
         end do

         if (tya.ge.0.0e0) then

c           No energy gain in the search direction         

c            write(iun5,*) ' ***** ya positive: ',tya
            reset = .true.
            goto 100

         endif

         call linmin(istat,n,npl,step,d1g1,d0,eold,emin,
     &               zr,par,par1,g,g1,gtol)
         call restr(par1)
         if (istat.eq.0) then
            call enegrd(eold,par,g)
            reset = .true.
            goto 100
         endif

c  We have found a lower energy.
      
         ncyc = ncyc + 1

c  Write out the minimized structure

         if (mod(ncyc,iwrite) .eq. 0) then

c     arc file
             if (iarc.eq.1) call wrtout(iun4,emin)

c     intermediate file

             icyc = icyc + 1
             if (imon.ne.0) call wrmon(icyc,emin)

             if (usesw) then
                if (mod(icyc,10).eq.0) call bldlst
             endif

         end if

c     Calculate RMS gradient, update parameters, calc G1G for new restart
c     Calculate parameter shift s() (for MQN method)

         g1g1 = 0.0e0
         do i=1,n
            g1g1   = g1g1 + g1(i)*g1(i)
            s(i)   = par1(i) - par(i) 
            par(i) = par1(i)
         end do

         rmsgrd = sqrt(g1g1/float(n))

c     Check whether we are done:

c        1. Too many cycles performed:

         if (ncyc.gt.nsd)  then
            write(iun5,*) 'Exceeded maximum interations ',nsd
            return
         endif
      
         if (ilast.lt.mxunch) then
            elast(ilast+1) = emin
         else
            do i=2,mxunch
               elast(i-1) = elast(i)
            end do
            elast(mxunch) = emin
         endif

         ilast = ilast + 1

         eavg = 0.0e0
         emn = 1000000.0e0
         emx = -1000000.0e0
         do i=1,min(ilast,mxunch)
            if (elast(i).lt.emn) emn = elast(i)
            if (elast(i).gt.emx) emx = elast(i)
         end do
         eavg = abs((emx - emn) / 2.0e0)
         if (ilast.gt.mxunch.and.eavg.lt.0.002e0) then
               write(iun5,*) 
     &          'Energy change in last ',mxunch,' steps < 0.002'
               return
         endif

c        print energy

         call prtene(icyc,emin,rmsgrd)

c        2. Gradient vector sufficiently small:

         if (rmsgrd.lt.gtol) then

c           GTOL: stop earlier if RMSgrd gradient < GTOL
c           Check whether the pair list was still good:

            call enegrd(emin2,par,g1)
            if (abs(emin2-emin).lt.deprls) then
               emin = emin2
               write(iun5,*) 'Optimisation Converged'
               return
            endif


            if (abs(emin2-emin).gt.derest) then 

c              Restart search direction if the energy changed too much

               eold = emin2

               do k=1,n
                  g(k) = g1(k)
               end do

               reset = .true.
               goto 100

            else 

c              Continue normally with next search direction

               emin = emin2

            endif

         endif

         if (imeth.eq.0) then
   
            do i=1,nflx3
               dg(i)  = -( g1(idx(i)) - g(idx(i)) )            
               hdg(i) = 0.0e0
            end do
   
         else if (imeth.eq.1) then
   
            yg = 0.0e0         
            py = 0.0e0         
   
            do i=1,n
               y(i) = g(i) - g1(i)
               yg   = yg   - y(i)*g1(i)
               py   = py   + y(i)*zr(i)
            end do
   
         else if (imeth.eq.2) then
c MQN
   
            yy = 0.0e0
            sy = 0.0e0
            yg = 0.0e0
            sg = 0.0e0
   
            do i=1,n
               s(i) = par1(i) - par(i)
               y(i) = g(i)    - g1(i)
               yy   = yy      + y(i)*y(i)
               sy   = sy      + s(i)*y(i)
               yg   = yg      - y(i)*g1(i)
               sg   = sg      - s(i)*g1(i)
            end do
   
         endif
      
c For all methods:
c For restart criterium 
c Update derivatives;

         g1g  = 0.0e0
         g1g1 = 0.0e0

         do i=1,n
            g1g1 = g1g1 + g1(i)*g1(i)
            g1g  = g1g  + g1(i)*g(i)
            g(i) = g1(i)
         end do

         dnew   = g1g1
         rmsgrd = sqrt(g1g1/float(n))
         gnorm  = sqrt(g1g1)

c If RMS gradient is extremely high, perform another Search Dir step

         if (rmsgrd.gt.500.0e0) then
            eold  = emin
            reset = .true.
            goto 100
         endif

c Update inverse hessian by BFGS algoritm

         if (imeth.eq.0) call updhes(nflx3,hessin,hdg,dg,zr,idx)

         npl    = npl + 1
         npow   = npow + 1

         reset  = .false.
         
c Set new Search Direction in zr():

         if (imeth.eq.0) then
   
            do i=1,iatoms
               if (iopt(i).eq.0) then
                 i3 = (i-1)*3
                 do j=1,3
                    zr(i3+j) = 0.0e0
                 end do
               endif
            end do

            do i=1,nflx3

               zr(idx(i)) = 0.0e0

               do j=1,nflx3
                  ij    = (i-1)*nflx3 + j
                  zr(idx(i)) = zr(idx(i)) + hessin(ij)*g(idx(j))
               end do

            end do

            if (npl.eq.nrstrt) reset = .true.
   
         else if (imeth.eq.1) then
   
            gamma = yg/py

            if (abs(g1g).ge. 0.2e0*g1g1 .or. npl.eq.1 .or.
     &           npow.eq.nrstrt) then

               if (npl.ne.1) write(iun5,*) 'POWELL restart'

               gamma1 = 0.0e0
               ptyt   = 0.0e0
               npow   = 0

               do i=1,n
                  yt(i) = y(i)
                  pt(i) = zr(i)
                  ptyt  = ptyt + yt(i)*pt(i)
               end do

            else 

               ytg = 0.0e0

               do i=1,n
                  ytg = ytg - yt(i)*g(i)
               end do

               gamma1 = ytg/ptyt


            endif

            do i=1,n
               zr(i) = g(i) + gamma*zr(i) + gamma1*pt(i)
            end do
   
         else if (imeth.eq.2) then
c MQN
   
            do i=1,n
               zr(i) = g(i) + (y(i)*sg + s(i)*yg)/sy
     &              - (1.0e0 + yy/sy)*(s(i)*sg/sy)
            end do

            if (npl.eq.nrstrt)          reset = .true.
            if (abs(g1g).ge. 0.2*g1g1) reset = .true.
   
         endif
         
         d1g1   = 0.0e0
         zrnorm = 0.0e0

         do i=1,n
            d1g1   = d1g1   - g(i) *zr(i)
            zrnorm = zrnorm + zr(i)*zr(i)
         end do
   
         eold = emin
         step = d0/sqrt(zrnorm)
         if (step.lt.0.00001e0) step = 0.00001e0
         if (step.ne.step) write(iun5,*) 'warning step underflow'

100      if (reset) then

c           Restart search direction to the steepest descent direction

            d0 = step0
            call resdir(istat,imeth,n,nflx3,d0,gtol,eold,g,
     &                        npl,step,d1g1,zr,hessin)
            if (istat.eq.0) then
               emin = eold
               return
            endif
         endif

      end do

c     END LOOP

      return
      end

      subroutine updhes(n,hessin,hdg,dg,zr,idx)
      implicit real (a-h,o-z), integer (i-n)
      dimension hessin(*),hdg(*),dg(*),zr(*),idx(*)

      do j=1,n
         do i=1,n         
            ij = (i-1)*n + j
            hdg(i) = hdg(i) + hessin(ij)*dg(j)
         end do
      end do

      fac = 0.0e0
      fae = 0.0e0

      do i=1,n
         fac = fac + dg(i)*zr(idx(i))
         fae = fae + dg(i)*hdg(i)
      end do

      fac = 1.0e0/fac
      fad = 1.0e0/fae

      do i=1,n
         dg(i) = fac*zr(idx(i)) - fad*hdg(i)
      end do

      do j=1,n
         do i=1,n
            ij = (i-1)*n + j
            hessin(ij) = hessin(ij) + fac*zr(idx(i)) * zr(idx(j))
     &                              - fad*hdg(i) *hdg(j)
     &                              + fae*dg(i)  *dg(j)
         end do
      end do

      return
      end

      subroutine resdir(istat,imeth,n,nflx3,d0,tol,eold,g,
     &                  npl,step,d1g1,zr,hessin)
      implicit real (a-h,o-z), integer (i-n)
c     Restart search direction to the steepest descent direction
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      dimension hessin(*),g(*),zr(*)

      istat  = 1
      npl    = 0

      dold   = 0.0e0
      do i=1,n
         zr(i) = g(i)
         dold  = dold + zr(i)*zr(i)
      end do
      
c Initialize inverse hessian for BFGS variable metric      
      
      if (imeth.eq.0) then

         do i=1,nflx3
            do j=1,nflx3
               ij = (i-1)*nflx3 + j
               if (i.eq.j) then
                  hessin(ij) = 1.0e0
               else
                  hessin(ij) = 0.0e0
               endif
            end do
         end do

      endif

      d1g1  = -dold
      gnorm = sqrt(dold)
      step  = d0/gnorm
      if (step.ne.step) write(iun5,*) 'warning step underflow'

      if (eold.ge.50000.0e0) step = step*100.0e0

      rmsgrd = sqrt(dold/float(n))

      if (rmsgrd.lt.tol) then
         write(iun5,*) ' ***** Gradient norm satisfied at restart CG'
         istat = 0
      endif

      return
      end

      subroutine linmin(istat,n,npl,step,d1g1,
     &                  d0,eold,emin,zr,par,parn,g,gn,gtol)
      implicit real (a-h,o-z), integer (i-n)
      logical cont, updene
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      dimension zr(*),par(*),parn(*),gn(*),g(*)

c     Start LINE MINIMIZATION.
c     Calculate energies (EA, EB) and derivatives (YA, YB) at points
c     A and B along the search direction (ZR), starting with A = 0 and B = STEP.
c     If necessary, adjust A and B to find an optimum search range.
c     Find a minimum X in the search range between A and B by determining
c     the minimum of the third-degree polynomial fixed by EA, EB, YA and YB.

c     FDIV:  Divide step by FDIV, every cycle ( 1.05)
c     FMULT: Multiply step by FMULT if NITI is not zero.( 2.0)
c     eps1,eps2: values for Shanno's convergence criteria

      eps1 = 0.9e0
      eps2 = 0.0001e0
      fdiv = 1.05e0
c      fmult= 1.2e0
      fmult= 2.0e0
      istat = 1
      updene = .true.

      do while (.true.)

         nstijg = 0

         ya = 0.0
         do i=1,n
            ya = ya - zr(i)*g(i)
         end do

         if (ya.ge.0.0e0) then

c           No energy gain in the search direction         

            write(iun5,*) ' ***** ya positive: ',ya
            istat = 0
            return

         endif

         a  = 0.0e0
         ea = eold
         b  = step
         niti = 0

100      cont = .true.
         do while (cont)

            if (updene) then
               do i=1,n
                  parn(i) = par(i) + b*zr(i)
               end do

               call enegrd(eb,parn,gn)

               yb = 0.0
               do i=1,n
                  yb = yb - zr(i)*gn(i)
               end do
            endif

            updene = .true.
            cont = .false.

            if (yb.lt.0.0e0 .and. eb.lt.ea) then

c           There must be a lower minimum further away than B, so shift and
c           double the search range.

               w  = b - a
               a  = b
               b  = a + 2.0e0*w
c               b = 2.0e0*b
               ea = eb
               ya = yb
               niti = niti + 1
               cont = .true.

            endif

            if (eb.gt.ea .and. abs(eb/ea).gt.1000.0e0) then      

c           Very asymmetric energy function, decrease search range

               bt = b
               b = (b+a)/2.0e0
               cont = .true.
            endif
      
         end do

c     Now we can only get here if eb > ea
c     Now there must be a minimum X in the search range between A and B. 
c     Find it and calculate the energy there.

         if (niti.ne.0) then
            step = fmult*step
            if (step.ne.step) write(iun5,*) 'step underflow'
            d0   = fmult*d0
         endif

         z = 3.0e0*(ea-eb)/(b-a) + ya + yb
         w = z*z - ya*yb
   
         if (w.le.0.0e0) then
c            write(iun5,*) 
c     &        ' **** LINMIN: w < 0!!! This should never happen'
            x = 0.0e0
         else
            w = sqrt(w)
            x = b - (yb+w-z)/(yb-ya+2.0e0*w)*(b-a)
         endif
      
c     apply the next guess for the scaling (x) of the search direction 
c     to create new parameters (parn)

         do i=1,n
            parn(i) = par(i) + x*zr(i)
         end do

c     calculate the energy associated with the scaling (x)

         call enegrd(emin,parn,gn)

         d1gs = 0.0e0

         do i=1,n
            d1gs = d1gs - zr(i)*gn(i)
         end do

c     Is this the best scaling of the search direction ?
c     we have three energies now ea, eb and emin

         if (niti.lt.30.and.(emin.gt.eold+d1g1*x*eps2) 
     &     .or.(abs(d1gs).gt.eps1*abs(d1g1)) .or.
     &   (((emin.gt.ea) .or. (emin.gt.eb)).and.((x.gt.1.01*a).and.
     &     (x.lt.0.99e0*b)))) then

c           NOT converged yet

            nstijg = nstijg + 1
   
            if ((x*d1gs).ge.0) then
               b  = x
               eb = emin
               yb = d1gs
            else
               a  = x
               ea = emin
               ya = d1gs
            endif

            updene = .true.

            if (nstijg.lt.5) then

c              Try a next interation within the chosen range

               updene = .false.
               goto 100

            else

               if (npl.eq.0) then

c                 we had already reset search direction, so now
c                 decrease step and try again

                  step = step/10.0e0
                  if (step.ne.step) write(iun5,*) 'step underflow'
                  d0   = d0/10.0e0

                  if (step.lt.1.0e-05) then
                     write(iun5,*) ' STEP in LINMIN too small'
		     d1gs = 0.0e0
                     do i=1,n
                        d1gs = d1gs + gn(i)*gn(i)
                     end do
                     d1gs = sqrt(d1gs/float(n))
                     gtol  = d1gs*1.1e0
                     istat = 1
                     return
                  endif

               else

c                 Restart search direction to the steepest descent direction

                  istat = 0
                  return

               endif
            endif

         else

c Line minimization converged

c            d0 = d0/fdiv      
            return

         endif

      end do
 
      return
      end

      subroutine prtgrd(par,grd)
      implicit real (a-h,o-z), integer (i-n)
      common /athlp/   iatoms, mxnat
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      dimension par(3,*),grd(3,*)
      
      write(iun5,*) 'grd'
      do i=1,iatoms
         write(iun5,'(i5,3f8.3,a,3f8.3)') 
     &       i,(par(j,i),j=1,3),' = ',(grd(j,i),j=1,3)
      end do

      return
      end

      subroutine prtarr(zr,n,str)
      implicit real (a-h,o-z), integer (i-n)
      common /athlp/   iatoms, mxnat
      character*(*) str
      dimension zr(*)

      print*,str

      n3 = n/3
      do i=1,n3
         print*,i,' ',(zr((i-1)*3+j),j=1,3)
      end do

      return
      end

