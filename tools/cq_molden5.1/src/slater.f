      subroutine slater(x,y,z,psi)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common/coord/ xyz(3,numatm)
      common /moldat/ natoms,norbs,nelecs,nat(numatm)
      common /mopac/ nfirst(numatm),nlast(numatm),
     &               npq(numatm),pqn(numatm),emus(numatm),
     &               emup(numatm),emud(numatm),consts(numatm),
     &               constp(numatm),constd(numatm),npqref(54)
      common /orbhlp/ mxorb,iuhf,ispd

      dimension psi(*)
c********************** slater type orbitals *********************
      do i=1,mxorb
          psi(i) = 0.0d0
      end do
      ico = 0

      do l=1,natoms

          xk = x - xyz(1,l)
          yk = y - xyz(2,l)
          zk = z - xyz(3,l)

          r2 = xk*xk + yk*yk + zk*zk
          r2 = r2 + 1.d-10
          r1 = dsqrt(r2)

          cons = consts(l)*r1**npq(l)*dexp(-emus(l)*r1)
          conp = constp(l)*r1**(npq(l)-1)*dexp(-emup(l)*r1)
          cond = constd(l)*r1**(pqn(l)-2)*dexp(-emud(l)*r1)

c          print*,'consts ',consts(l),' npq ',npq(l),' emus ',emus(l)
c          print*,'constp ',constp(l),' npq ',npq(l),' emup ',emup(l)
c          print*,'constd ',constd(l),' pqn ',pqn(l),' emud ',emud(l)

          i = nlast(l) - nfirst(l) + 1

          goto  (3,4,4,4,7,7,7,7,7),  i
 
 
    7     psi(ico+9)=               cond*xk*yk
          psi(ico+8)=               cond*yk*zk
          psi(ico+7)=0.2886751346d0*cond*(2.d0*zk*zk-yk*yk-xk*xk)
          psi(ico+6)=               cond*xk*zk
          psi(ico+5)=0.50d0*cond*(xk*xk-yk*yk)
   4      psi(ico+4)=conp*zk
          psi(ico+3)=conp*yk
          psi(ico+2)=conp*xk
   3      psi(ico+1)=cons
 
 
          ico=ico+i
      end do

      return
      end
