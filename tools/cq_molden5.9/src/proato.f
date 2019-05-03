      subroutine proato
c     make projection of atoms onto the plot plane
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /eulx/   ca,cb,sa,sb,cc,sc
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      common /atopla/ xsym(numatm),ysym(numatm),zsym(numatm),
     &                isym(numatm)

      do l=1,natoms
         isym(l) = 0
         xk      = xyz(1,l) - px
         yk      = xyz(2,l) - py
         zk      = xyz(3,l) - pz
c inproduct of vector atom-planecenter with outvector of plane
         xinp    = xk*cx + yk*cy + zk*cz
         xinp    = xinp * xinp
         if (xinp.lt.1.d-10) then
           isym(l) = 1
         endif
         xsym(l) =  - ((ca*xk + sa*yk)*cb - zk*sb)
         ysym(l) =  -  (ca*yk - sa*xk)
         zsym(l) = ((ca*xk + sa*yk)*sb + zk*cb)
         x1 = (cc*xsym(l) + sc*ysym(l))
         y1 = (-sc*xsym(l) + cc*ysym(l))
         xsym(l) = x1
         ysym(l) = y1  
      end do

      return
      end
