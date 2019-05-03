      subroutine tessa(axis,cosi,sinu,dist,cori,crot,natoms,tol)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      dimension axis(3),dist(numat1),cori(3,numat1),
     &          crot(3,numat1)
      dimension vect1(3),vect2(3),vect3(3)

      if (vlen(axis).lt.tol) stop 'tessa: rotation axis = 0 vector'

      iatoms = natoms
      if (iatoms.gt.numat1) iatoms = numat1

      do k=1,iatoms

         if (dist(k).gt.tol) then

            do j=1,3
               vect1(j) = axis(j)
            end do

            call impsc(vect1,cori(1,k),cott)

            if (dabs(dabs(cott)-1.0d0).gt.tol) then

               call crprod(vect1,cori(1,k),vect2)
               call crprod(vect2,vect1,vect3)
               call impsc(vect3,cori(1,k),sitt)
               call vsc1(vect1,dist(k)*cott,tol)
               call vsc1(vect2,dist(k)*sitt,tol)
               call vsc1(vect3,dist(k)*sitt,tol)
                        
               do j=1,3
                  crot(j,k) = vect1(j) + sinu*vect2(j) + cosi*vect3(j)
               end do

            else

               do j=1,3
                  crot(j,k) = cori(j,k)
               end do

            endif
         else

            do j=1,3
               crot(j,k) = cori(j,k)
            end do

         endif

      end do

      return
      end
