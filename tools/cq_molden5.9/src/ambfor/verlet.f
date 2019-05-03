      subroutine velind(coo,g,v,a,mass)
      implicit real (a-h,o-z), integer (i-n)
      real mass
      common /athlp/   iatoms, mxnat
      common /mdopt/ dt,temp,tau,nstep,ndump,nfree,idumv
      common /mdconv/ convert, gasconst, boltzmann, pi
      integer taskid
      common /mpih/ nproc,taskid
cmpi      include 'mpif.h'
cmpi      integer status(mpi_status_size)
      logical opfil
      dimension vec(3),coo(3,*),v(3,*),a(3,*),g(3,*),mass(*)

      call enegrd(ene,coo,g)

cmpi      call mpi_comm_rank(mpi_comm_world, taskid, ierr)


      do i=1,iatoms

            beta = sqrt(mass(i) / (2.0e0*boltzmann*temp))
            rho = random()
            xs  = erfinv(rho) / beta
            rho = random()
            ys  = erfinv(rho) / beta
            rho = random()
            zs  = erfinv(rho) / beta
            speed = sqrt(xs**2 + ys**2 + zs**2)

            call runitv(vec)

            do j=1,3
                v(j,i) = speed * vec(j)
            end do

      end do

cmpi      if (nproc.gt.1.and.taskid.eq.0) then
         if (idumv.eq.1) then
            if (opfil(53,'velos',5,1,0,0)) then
                do i=1,iatoms
                   write(53,'(3(f13.7,1x))') (v(j,i),j=1,3)
                end do
                close(53)
            endif
         else if (idumv.eq.-1) then
            if (opfil(53,'velos',5,1,1,0)) then
                do i=1,iatoms
                   read(53,'(3(f13.7,1x))') (v(j,i),j=1,3)
                end do
                close(53)
            endif
         endif
cmpi      endif

      do i=1,iatoms
         do j=1,3
             a(j,i) = convert * g(j,i) / mass(i)
         end do
      end do

cmpi      if (taskid.eq.0) then
cmpi         do it=1, nproc - 1
cmpi            call mpi_send(v, 3*iatoms, mpi_real, it, 15,
cmpi     &                    mpi_comm_world, ierr)
cmpi         end do
cmpi      else
cmpi         call mpi_recv(v, 3*iatoms, mpi_real, 0, 15,
cmpi     &                    mpi_comm_world, status, ierr)
cmpi      endif

cmpi      if (taskid.eq.0) then
cmpi         do it=1, nproc - 1
cmpi            call mpi_send(a, 3*iatoms, mpi_real, it, 16,
cmpi     &                    mpi_comm_world, ierr)
cmpi         end do
cmpi      else
cmpi         call mpi_recv(a, 3*iatoms, mpi_real, 0, 16,
cmpi     &                    mpi_comm_world, status, ierr)
cmpi      endif


      return
      end

      subroutine verled(istep,coo,g,v,a,mass,ityp)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxion=2000)
      real mass
      integer*2 ityp
      logical osingl,dolbfgs,oqscal
      common /athlp/  iatoms, mxnat
      common /mdopt/ dt,temp,tau,nstep,ndump,nfree,idumv
      common /opts/ idebug,imon,iarc,ilog,iout,osingl,dolbfgs,oqscal
      common /mdconv/ convert, gasconst, boltzmann, pi
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /h2oer/  numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      common /cyctmp/ icyc
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension coo(3,*),v(3,*),g(3,*),a(3,*),mass(*),ityp(*)

      dt_2  = dt / 2.0e0
      dt2_2 = dt * dt_2
 
c     store the current atom positions, then find new atom
c     positions and half-step velocities via Verlet recursion
 
      do i=1,iatoms

         do j=1,3
            coo(j,i) = coo(j,i) + v(j,i)*dt + a(j,i)*dt2_2
         end do

         do j=1,3
            v(j,i) = v(j,i) + a(j,i)*dt_2
         end do

      end do

      if (box.and.(mod(istep,2).eq.0.or.istep.eq.nstep)) 
     &        call appbnd(coo,ityp)

      if (box.and.fast) then
          if (mod(istep,6).eq.0) call watlst(niwat)
      endif

      call enegrd(ene,coo,g)
 
c     find the full-step velocities using the Verlet recursion
 
      do i=1,iatoms
          do j=1,3
c             a(j,i) = -convert * g(j,i) / mass(i)
             a(j,i) = convert * g(j,i) / mass(i)
             v(j,i) = v(j,i) + a(j,i)*dt_2
          end do
      end do
 
c     calculate the kinetic energy
 
      call kinetic(emv2,v,mass)
 
c     from kinetic energy calculate instantaneous temperature
c     make temperature corrections via berendsen thermostat
 
      call thermst(emv2,v)

      etot = emv2 + ene

      if (mod(istep,ndump).eq.0) then
         call wrtout(iun4,etot)
         icyc = icyc + 1
         if (imon.ne.0) call wrmon(icyc,etot)
      endif

      return
      end

      subroutine kinetic(emv2,v,mass)
      implicit real (a-h,o-z), integer (i-n)
      real mass
      common /mdconv/ convert, gasconst, boltzmann, pi
      common /athlp/  iatoms, mxnat
      dimension ekin(3,3),v(3,*),mass(*)

      emv2 = 0.0e0

      do i = 1, 3
         do j = 1, 3
            ekin(j,i) = 0.0e0
         end do
      end do

      do i=1,iatoms

          term = 0.5e0 * mass(i) / convert

          do j=1,3

             do k=1,3
                value = term * v(j,i) * v(k,i)
                ekin(k,j) = ekin(k,j) + value
             end do

          end do

      end do

      emv2 = ekin(1,1) + ekin(2,2) + ekin(3,3)

      return
      end


      subroutine thermst(emv2,v)
      implicit real (a-h,o-z), integer (i-n)
      common /mdopt/ dt,temp,tau,nstep,ndump,nfree,idumv
      common /mdconv/ convert, gasconst, boltzmann, pi
      common /athlp/  iatoms, mxnat
      dimension v(3,*)

c     scale velocities to satisfy Berendsen thermostat
c     instantaneous temperature from the kinetic energy

      tmpi = 2.0e0 * emv2 / (dble(nfree) * gasconst)

c     scale velocities

      if (tmpi.ne.0.0e0) then
         scale = sqrt(1.0e0 + (dt/tau)*(temp/tmpi-1.0e0))
      else
         scale = 1.0e0
      endif

      do i=1,iatoms

         do j=1,3
            v(j,i) = scale * v(j,i)
         end do

      end do

      return
      end

      subroutine assmad(ityp,mass)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      real mass
      integer*2 ityp
      dimension mass(*),ityp(*)
      common /masses/  ambmas(mxamb),gffmas(mxgff)
      common /athlp/   iatoms, mxnat

      do i=1,iatoms
         it = abs(ityp(i))
         if (ityp(i).gt.0) then
            mass(i) = ambmas(it)
         endif
         if (ityp(i).lt.0) then
            mass(i) = gffmas(it)
         endif
      end do

      return
      end

