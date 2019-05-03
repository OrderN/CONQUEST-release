      subroutine lbfgd(emin,ggtol,nsd,x,g,diag,w)
      implicit real (a-h,o-z), integer (i-n)
      common /lb3/gtol,stpmin,stpmax
      common /athlp/   iatoms, mxnat
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      logical osingl,dolbfgs,oqscal
      common /opts/ idebug,imon,iarc,ilog,iout,osingl,dolbfgs,oqscal
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /cell/ xa,ya,yb,za,zb,zc
      common /convg/ icyc
      dimension x(*),g(*),diag(*),w(*)

      write(iun5,*) 'L-BFGS method'
      write(iun5,*) ' '

      gtol = 0.9e0
      stpmin = 1.0e-20
      stpmax = 1.0e+20
      iwrite = 1

      n = iatoms*3

      if (cell) then
         x(n+1) = xa
         x(n+2) = ya
         x(n+3) = yb
         x(n+4) = za
         x(n+5) = zb
         x(n+6) = zc

         n = n + 6
      endif

      m = 5


      eps    = ggtol

      xtol   = 1.0e-16
      write(iun5,*) 'Gradient tolerance: ',eps
      write(iun5,*) 'Maximum iterations: ',nsd
      write(iun5,*) ' '

      ncyc  = 0
      icyc  = 0

      iflag  = 0

      do while (.true.)

         call enegrd(emin,x,g)
         call rmsg(n,g,rmsgrd)

         if (osingl) then
             write(iun5,'(a,f12.3)') 'Estat=',emin
             return
         endif

         do i=1,n
            g(i) = -g(i)
         end do

         call dlbfgs(n,m,x,emin,g,diag,eps,xtol,w,iflag)
         if (iflag.eq.0) return
         if (iflag.lt.0) then
                return
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

         if (ncyc.gt.nsd) then
             write(iun5,*) 'Exceeded maximum interations ',nsd
             return
         endif

      end do

      return
      end

      subroutine dlbfgs(n,m,x,f,g,diag,eps,xtol,w,iflag)
      implicit real (a-h,o-z), integer (i-n)
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      dimension x(*),g(*),diag(*),w(*)

c        limited memory bfgs method for large scale optimization
c                          jorge nocedal

      common /lb3/gtol,stpmin,stpmax
      common /rdwr/ iun1,iun2,iun3,iun4,iun5

      integer bound,cp,point
      logical finish
      save

c     initialize

      if (iflag.eq.1)  goto 172

      iter = 0
      
      if (n.le.0.or.m.le.0) goto 196

      if (gtol.le.1.e-04) then
        write(iun5,245)
        gtol = 9.e-01
      endif

      nfun   = 1
      point  = 0
      finish = .false.

      do i=1,n
         diag(i) = 1.0e0
      end do

      ispt = n+2*m
      iypt = ispt+n*m     

      do i=1,n
         w(ispt+i) = -g(i)*diag(i)
      end do

      gnorm = sqrt(ddot(n,g,g))
      gpnorm = sqrt(ddot(n,g,g)/dble(n))
      stp1  = 1.0e0/gnorm

c     parameters for line search routine
     
      ftol   = 1.0e-4
      maxfev = 20

      call prtene(nfun,f,gpnorm)

c    --------------------
c     main iteration loop
c    --------------------

100      continue

         iter  = iter+1
         info  = 0
         bound = iter-1

         if (iter.ne.1) then
            if (iter.gt.m) bound = m

            ys = ddot(n,w(iypt+npt+1),w(ispt+npt+1))
            yy = ddot(n,w(iypt+npt+1),w(iypt+npt+1))

            do i=1,n
               diag(i) = ys/yy
            end do

            cp = point
            if (point.eq.0) cp = m
            w(n+cp) = 1.0e0/ys

            do i=1,n
               w(i) = -g(i)
            end do

            cp = point
   
            do i=1,bound
               cp = cp-1
               if (cp.eq. -1) cp = m-1
               sq = ddot(n,w(ispt+cp*n+1),w)
               inmc = n + m + cp + 1
               iycn = iypt + cp*n
               w(inmc) = w(n+cp+1)*sq
               call daxpy(n,-w(inmc),w(iycn+1),w)
            end do

            do i=1,n
               w(i) = diag(i)*w(i)
            end do

            do i=1,bound

               yr   = ddot(n,w(iypt+cp*n+1),w)
               beta = w(n+cp+1)*yr
               inmc = n + m + cp + 1
               beta = w(inmc)-beta
               iscn = ispt + cp*n

               call daxpy(n,beta,w(iscn+1),w)

               cp   = cp + 1
               if (cp.eq.m) cp = 0

            end do

            do i=1,n
               w(ispt+point*n+i) = w(i)
            end do

         endif

         nfev = 0
         stp  = 1.0e0
         if (iter.eq.1) stp = stp1

         do i=1,n
            w(i) = g(i)
         end do

 172     continue

         call mcsrch(n,x,f,g,w(ispt+point*n+1),stp,ftol,
     &               xtol,maxfev,info,nfev,diag)
         if (info .eq. -1) then
            iflag = 1
            return
         endif

         if (info.ne.1) goto 190

         nfun = nfun + nfev

c     compute the new step and gradient change 

         npt = point*n
         do i=1,n
            w(ispt+npt+i) = stp*w(ispt+npt+i)
            w(iypt+npt+i) = g(i)-w(i)
         end do

         point = point + 1
         if (point.eq.m) point = 0

c     termination test

         gnorm = sqrt(ddot(n,g,g))
         gpnorm = sqrt(ddot(n,g,g)/dble(n))
         xnorm = sqrt(ddot(n,x,x))
         xnorm = max(1.0e0,xnorm)

         if (gpnorm .le. eps) finish = .true.

         call prtene(nfun,f,gpnorm)

         if (finish) then
            write(iun5,*) 'Optimisation Converged'
            iflag = 0
            return
         endif

      goto 100

c     end of main iteration loop. error exits.

 190  iflag=-1
      write(iun5,200) info
      return

 195  iflag=-2
      write(iun5,235) i
      return

 196  iflag= -3
      write(iun5,240)

c     formats

 200  format(/' iflag= -1 ',/' line search failed. see',
     &        ' documentation of routine mcsrch',/' error return',
     &        ' of line search: info= ',i2,/
     &        ' possible causes: function or gradient are incorrect',/,
     &        ' or incorrect tolerances')
 235  format(/' iflag= -2',/' the',i5,'-th diagonal element of the',/,
     &       ' inverse hessian approximation is not positive')
 240  format(/' iflag= -3',/' improper input parameters (n or m',
     &       ' are not positive)')
 245  format(/'  gtol is less than or equal to 1.e-04',
     &       / ' it has been reset to 9.e-01')

      return
      end

      subroutine daxpy(n,da,dx,dy)
      implicit real (a-h,o-z), integer (i-n)
      dimension dx(*),dy(*)

      if (n.le.0.or.da.eq.0.0e0) return

      m = mod(n,4)

      if (m.ne.0) then
         do i= 1,m
           dy(i) = dy(i) + da*dx(i)
         end do

         if (n.lt.4) return

      endif

      mp1 = m + 1

      do i= mp1,n,4
         dy(i)   = dy(i)   + da*dx(i)
         dy(i+1) = dy(i+1) + da*dx(i+1)
         dy(i+2) = dy(i+2) + da*dx(i+2)
         dy(i+3) = dy(i+3) + da*dx(i+3)
      end do

      return
      end

      real function ddot(n,dx,dy)

c     forms the dot product of two vectors.

      implicit real (a-h,o-z), integer (i-n)
      dimension dx(*),dy(*)

      ddot = 0.0e0
      dtemp = 0.0e0

      if (n.le.0) return


      m = mod(n,5)
      if (m.ne.0) then

          do i=1,m
            dtemp = dtemp + dx(i)*dy(i)
          end do

          if (n.lt.5) goto 60
      endif

      mp1 = m + 1

      do i = mp1,n,5
         dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) +
     &   dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
      end do


60    ddot = dtemp

      return
      end


      subroutine mcsrch(n,x,f,g,s,stp,ftol,xtol,maxfev,info,nfev,wa)

c     line search routine

      implicit real (a-h,o-z), integer (i-n)
      common /lb3/gtol,stpmin,stpmax
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      save
      logical brackt,stage1

      dimension x(*),g(*),s(*),wa(*)

      if (info.ne.-1) then

         infoc = 1


         if (n .le. 0 .or. stp .le. 0.0e0 .or. ftol .lt. 0.0e0 .or.
     &    gtol .lt. 0.0e0 .or. xtol .lt. 0.0e0 .or. stpmin .lt. 0.0e0
     &    .or. stpmax .lt. stpmin .or. maxfev .le. 0) return

         dginit = 0.0e0

         do j = 1, n
            dginit = dginit + g(j)*s(j)
         end do

         if (dginit .ge. 0.0e0) then
            write(iun5,'(a)')
     &       '  the search direction is not a descent direction'
            return
         endif


         brackt = .false.
         stage1 = .true.
         nfev = 0
         finit = f
         dgtest = ftol*dginit
         width = stpmax - stpmin
         width1 = width/0.5e0

         do j = 1, n
            wa(j) = x(j)
         end do

         stx = 0.0e0
         fx  = finit
         dgx = dginit
         sty = 0.0e0
         fy  = finit
         dgy = dginit
   
      endif

c     start of iteration.

      do while (.true.)

         if (info.ne.-1) then
            if (brackt) then
               stmin = min(stx,sty)
               stmax = max(stx,sty)
            else
               stmin = stx
               stmax = stp + 4.0e0*(stp - stx)
            end if

            stp = max(stp,stpmin)
            stp = min(stp,stpmax)


            if ((brackt .and. (stp.le.stmin .or. stp.ge.stmax))
     &      .or. nfev.ge.maxfev-1 .or. infoc.eq.0
     &      .or. (brackt .and. stmax-stmin.le.xtol*stmax)) stp = stx

            do j = 1, n
               x(j) = wa(j) + stp*s(j)
            end do

            info = -1
            return

         endif

         info = 0
         nfev = nfev + 1
         dg = 0.0e0

         do j = 1, n
            dg = dg + g(j)*s(j)
         end do

         ftest1 = finit + stp*dgtest

c        test for convergence.

         if ((brackt .and. (stp.le.stmin .or. stp.ge.stmax))
     &      .or. infoc.eq.0) info = 6

         if (stp.eq.stpmax .and.
     &       f.le.ftest1 .and. dg.le.dgtest) info = 5

         if (stp.eq.stpmin .and.
     &       (f.gt.ftest1 .or. dg.ge.dgtest)) info = 4

         if (nfev.ge.maxfev) info = 3

         if (brackt .and. stmax-stmin.le.xtol*stmax) info = 2

         if (f.le.ftest1 .and. abs(dg).le.gtol*(-dginit)) info = 1

c        check for termination.

         if (info.ne.0) return


         if (stage1 .and. f.le.ftest1 .and.
     &       dg.ge.min(ftol,gtol)*dginit) stage1 = .false.

         if (stage1 .and. f.le.fx .and. f.gt.ftest1) then

            fm   = f - stp*dgtest
            fxm  = fx - stx*dgtest
            fym  = fy - sty*dgtest
            dgm  = dg - dgtest
            dgxm = dgx - dgtest
            dgym = dgy - dgtest

            call mcstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,
     &                 brackt,stmin,stmax,infoc)

            fx = fxm + stx*dgtest
            fy = fym + sty*dgtest
            dgx = dgxm + dgtest
            dgy = dgym + dgtest

         else

            call mcstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg,
     &                 brackt,stmin,stmax,infoc)
         end if

         if (brackt) then

            if (abs(sty-stx).ge.0.66e0*width1)
     &          stp = stx + 0.5e0*(sty - stx)
            width1 = width
            width  = abs(sty-stx)

         end if

      end do

      return
      end

      subroutine mcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
     &                  stpmin,stpmax,info)
      implicit real (a-h,o-z), integer (i-n)
      logical brackt,bound

c     the purpose of mcstep is to compute a safeguarded step for
c     a linesearch and to update an interval of uncertainty for
c     a minimizer of the function.
 

      info = 0

c     check the input parameters for errors.

      if ((brackt .and. (stp .le. min(stx,sty) .or.
     &    stp .ge. max(stx,sty))) .or.
     &    dx*(stp-stx) .ge. 0.0e0 .or. stpmax .lt. stpmin) return

c     determine if the derivatives have opposite sign.

      sgnd = dp*(dx/abs(dx))

c     first case. a higher function value.
c     the minimum is bracketed. if the cubic step is closer
c     to stx than the quadratic step, the cubic step is taken,
c     else the average of the cubic and quadratic steps is taken.

      if (fp .gt. fx) then

              info = 1
              bound = .true.
              theta = 3*(fx - fp)/(stp - stx) + dx + dp
              s = max(abs(theta),abs(dx),abs(dp))
              gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
              if (stp .lt. stx) gamma = -gamma
              p = (gamma - dx) + theta
              q = ((gamma - dx) + gamma) + dp
              r = p/q
              stpc = stx + r*(stp - stx)
              stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)

              if (abs(stpc-stx) .lt. abs(stpq-stx)) then
                 stpf = stpc
              else
                stpf = stpc + (stpq - stpc)/2
              end if

              brackt = .true.

c     second case. a lower function value and derivatives of
c     opposite sign. the minimum is bracketed. if the cubic
c     step is closer to stx than the quadratic (secant) step,
c     the cubic step is taken, else the quadratic step is taken.

      else if (sgnd .lt. 0.0) then

              info = 2
              bound = .false.
              theta = 3*(fx - fp)/(stp - stx) + dx + dp
              s = max(abs(theta),abs(dx),abs(dp))
              gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
              if (stp .gt. stx) gamma = -gamma
              p = (gamma - dp) + theta
              q = ((gamma - dp) + gamma) + dx
              r = p/q
              stpc = stp + r*(stx - stp)
              stpq = stp + (dp/(dp-dx))*(stx - stp)

              if (abs(stpc-stp) .gt. abs(stpq-stp)) then
                 stpf = stpc
              else
                 stpf = stpq
              end if

              brackt = .true.

c     third case. a lower function value, derivatives of the
c     same sign, and the magnitude of the derivative decreases.
c     the cubic step is only used if the cubic tends to infinity
c     in the direction of the step or if the minimum of the cubic
c     is beyond stp. otherwise the cubic step is defined to be
c     either stpmin or stpmax. the quadratic (secant) step is also
c     computed and if the minimum is bracketed then the the step
c     closest to stx is taken, else the step farthest away is taken.

      else if (abs(dp) .lt. abs(dx)) then

              info = 3
              bound = .true.
              theta = 3*(fx - fp)/(stp - stx) + dx + dp
              s = max(abs(theta),abs(dx),abs(dp))

c        the case gamma = 0 only arises if the cubic does not tend
c        to infinity in the direction of the step.

              gamma = s*sqrt(max(0.0e0,(theta/s)**2 - (dx/s)*(dp/s)))

              if (stp .gt. stx) gamma = -gamma

              p = (gamma - dp) + theta
              q = (gamma + (dx - dp)) + gamma
              r = p/q

              if (r .lt. 0.0e0 .and. gamma .ne. 0.0e0) then
                 stpc = stp + r*(stx - stp)
              else if (stp .gt. stx) then
                 stpc = stpmax
              else
                 stpc = stpmin
              endif

              stpq = stp + (dp/(dp-dx))*(stx - stp)

              if (brackt) then

                 if (abs(stp-stpc) .lt. abs(stp-stpq)) then
                    stpf = stpc
                 else
                    stpf = stpq
                 end if

              else

                 if (abs(stp-stpc) .gt. abs(stp-stpq)) then
                    stpf = stpc
                 else
                    stpf = stpq
                 end if

              end if

c     fourth case. a lower function value, derivatives of the
c     same sign, and the magnitude of the derivative does
c     not decrease. if the minimum is not bracketed, the step
c     is either stpmin or stpmax, else the cubic step is taken.

      else

              info = 4
              bound = .false.

              if (brackt) then
                 theta = 3*(fp - fy)/(sty - stp) + dy + dp
                 s = max(abs(theta),abs(dy),abs(dp))
                 gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
                 if (stp .gt. sty) gamma = -gamma
                 p = (gamma - dp) + theta
                 q = ((gamma - dp) + gamma) + dy
                 r = p/q
                 stpc = stp + r*(sty - stp)
                 stpf = stpc
              else if (stp .gt. stx) then
                 stpf = stpmax
              else
                 stpf = stpmin
              end if

      endif

c     update the interval of uncertainty. this update does not
c     depend on the new step or the case analysis above.

      if (fp .gt. fx) then

         sty = stp
         fy = fp
         dy = dp

      else

         if (sgnd.lt.0.0e0) then
            sty = stx
            fy = fx
            dy = dx
         end if

         stx = stp
         fx = fp
         dx = dp

      endif

c     compute the new step and safeguard it.

      stpf = min(stpmax,stpf)
      stpf = max(stpmin,stpf)
      stp = stpf

      if (brackt .and. bound) then

         stt = stx+0.66*(sty-stx)

         if (sty .gt. stx) then
            stp = min(stt,stp)
         else
            stp = max(stt,stp)
         endif

      endif

      return
      end

      subroutine rmsg(n,g1,rmsgrd)
      implicit real (a-h,o-z), integer (i-n)
      dimension g1(*)

      g1g1 = 0.0e0
      do i=1,n
         g1g1   = g1g1 + g1(i)*g1(i)
      end do

      rmsgrd = sqrt(g1g1/dble(n))

      return
      end

