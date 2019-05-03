      subroutine enegdd(energy,coo,forces,fintr,fwat,fx,fy,fz,
     &                  potscl,iopt,iresid,
     &                  nac,iac,nad,iad,
     &                  nbnd,ibnd,bl,bk,
     &                  nang,iang,ango,ak,
     &                  nt,it,trs1,trs2,trs3,trs4,
     &                  nti,iti,trsi1,trsi2,
     &                  q,iconn,ityp,iwtpr,nlst,lst)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (mxad=9*mxcon)
      parameter (mxion=2000)
      common /athlp/  iatoms, mxnat
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /h2oer/  numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw
      logical osingl,dolbfgs,oqscal
      common /opts/ idebug,imon,iarc,ilog,iout,osingl,dolbfgs,oqscal
      common /passe/ emin
      logical doint
      common /intr/ ecold,evold,ewold,iprog,doint
      common /wesc/ iwesc
      common /mpih/ nproc,idt
      integer*2 ityp
      real evs,ecs,fx,fy,fz
      dimension forces(3,*),potscl(*),iopt(*),iresid(*)
      dimension fintr(3,*),fwat(3,*),fx(*),fy(*),fz(*)
      dimension nac(*),nad(*),iac(mxac,*),iad(mxad,*)
      dimension ibnd(2,*), bl(*), bk(*)
      dimension iang(3,*),ango(*),ak(*)
      dimension it(4,*),trs1(4,*),trs2(4,*),trs3(4,*),trs4(4,*)
      dimension iti(4,*),trsi1(4,*),trsi2(4,*)
      dimension coo(3,*),q(*),iconn(mxcon+1,*),ityp(*),iwtpr(*)
      dimension nlst(*),lst(*)
      dimension cellder(6)

      eb  = 0.0e0
      ea  = 0.0e0
      eit = 0.0e0
      et  = 0.0e0
      ev  = 0.0e0
      ec  = 0.0e0
      ew  = 0.0e0
      euin = 0.0e0

      do i=1,iatoms
         do j=1,3
            forces(j,i) = 0.0e0
            if (doint) then
               fwat(j,i)  = 0.0e0
               fintr(j,i) = 0.0e0
            endif
         end do
         
      end do

      if (cell) then
          call var2cl(coo,iatoms)
          
c          call uinner(euin,cellder,coo,forces,q,iatoms)
      endif

      call bond (eb,nbnd,ibnd,bl,bk,coo,forces)

      call angle(ea,nang,iang,ango,ak,coo,forces)
      call tors (eit,1,
     &                nt,it,trs1,trs2,trs3,trs4,
     &                nti,iti,trsi1,trsi2,
     &                coo,forces)
      call tors (et,0,
     &                nt,it,trs1,trs2,trs3,trs4,
     &                nti,iti,trsi1,trsi2,
     &                coo,forces)

      if (cell) then
          call qvc(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
     &             q,forces,cellder,iopt,ityp)
      else
          if (usecut) then
             if (usesw) then
                 call qenvdw(ec,ev,nac,iac,nad,iad,iresid,
     &                  coo,iconn,q,forces,iopt,nlst,lst,ityp)
             else
                 call qvnoe(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
     &                  q,forces,iopt,ityp,iwtpr)
             endif
          else
          
          
                 if (doint) then
                    call qvdwrd(ec,ev,fintr,
     &                           nac,iac,nad,iad,iresid,coo,iconn,
     &                           q,iopt,ityp)
                    ecold = ec
                    evold = ev
                 else
                    ec = ecold
                    ev = evold
                 endif

                 do i=1,iatoms
                    do j=1,3
                       forces(j,i) = forces(j,i) + fintr(j,i)
                    end do
                 end do

          endif

          if (fast.and.box) then
              if (doint) then
                 call fstwat(ew,coo,fwat,q)
                 ewold = ew
              else
                 ew = ewold
              endif
              do i=1,iatoms
                 do j=1,3
                    forces(j,i) = forces(j,i) + fwat(j,i)
                 end do
              end do
          endif

      endif

      energy = eb + ea + eit + et + ev + ec + euin
      if (fast.and.box) energy = energy + ew

      emin = energy

      if (idebug.eq.1) then
        print*,'eb=',eb,'ea=',ea,'eit=',eit,'et=',et,'ev=',ev,'ec=',ec
        if (fast.and.box) print*,'ew=',ew
        print*,'energy=',energy
      endif

      do i=1,iatoms
         if (iopt(i).eq.1) then
            do j=1,3
               forces(j,i) = -forces(j,i)
            end do
         else
            do j=1,3
               forces(j,i) = 0.0e0
            end do
         endif
      end do

      if (cell) then
         call grd2var(cellder,forces,iatoms)
         call getcel(a,b,c,alpha,beta,gamma)
c         call prcldr(cellder,coo)
      endif

      if (iwesc.eq.1) then
         call wexit
      endif

      return
      end
