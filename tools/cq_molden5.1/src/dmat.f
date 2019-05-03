
      subroutine dmat(p,vectrs,vectrb,focc,focb)

c new style atomic density matrix goes with densmat 

      implicit double precision (a-h,p-z),integer (i-n)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      common /coord / xyz(3,numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      dimension p(*),vectrs(*),vectrb(*),focc(*),focb(*)
c
c   construct a density matrix from orbital occupancies
c
      do i=1,norbs
          do j=i,norbs
              sum = 0.d0
              do k=1,norbs
                  sum = sum + focc(k)*vectrs((k-1)*mxorb+i)
     &                               *vectrs((k-1)*mxorb+j)
                  if (iuhf.eq.1) then
                    if (ispd.eq.0) then
                       sum = sum + focb(k)*vectrb((k-1)*mxorb+i)
     &                                    *vectrb((k-1)*mxorb+j)
                    else
                       sum = sum - focb(k)*vectrb((k-1)*mxorb+i)
     &                                    *vectrb((k-1)*mxorb+j)
                    endif
                  endif
              end do
              p((i-1)*mxorb+j) = sum
              p((j-1)*mxorb+i) = sum
          end do
      end do

      return
      end

c OLD style atomic density matrix creator, goes with densmto

      subroutine dmato(v,pop,norbs,p)
      implicit double precision (a-h,p-w),integer   (i-n)
      common /orbhlp/ mxorb,iuhf,ispd
      dimension v(*),p(*),pop(*)
c
c   construct a density matrix from orbital occupancies
c
      do i=1,norbs
         do j=i,norbs
            sum = 0.d0
            do k=1,norbs
               sum = sum + pop(k)*v((k-1)*mxorb+i)*v((k-1)*mxorb+j)
            end do
            p((i-1)*norbs+j) = sum
            p((j-1)*norbs+i) = sum
         end do
      end do

      return
      end
