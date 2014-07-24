!
! GPF3D The basic 3D complex FFT.
!

!
! Arguments
!
! COMPLEX C(ID,NN(2),NN(3))
! ID first dimension of data array.  
! NN(3) Dimensions of 3D complex data.
! IS Forward (+1)/Reverse (-1)
!
      subroutine GPF3D(c,trigs,size,id,nn,is)
      parameter (NMAX=256)
      complex*16 c(id,id,id)
      integer nn(3)
      integer id,is,size, off
      real*8 trigs(size,3)

      real*8, allocatable, dimension(:) :: a,b

      integer inc, jump

      allocate(a(id*id*id),b(id*id*id))
      n=0
      do i=1,id
         do j=1,id
            do k=1,id
               n=n+1
               a(n) = real(c(k,j,i))
               b(n) = aimag(c(k,j,i))
            end do
         end do
      end do

      lot=nn(2)*nn(3)
      inc=1!2
      jump=id!2*id
      call GPFA(a,b,trigs(1,1),inc,jump,nn(1),lot,-is)
      inc=id!2*id
      jump=1!2
      lot=nn(1)
      do i=0,nn(3)-1
         off = id*nn(2)*i
         call GPFA(a(off+1),b(off+1),trigs(1,2),inc,jump,nn(2),lot,-is)
      end do
      lot=id*nn(2)
      inc=id*nn(2)!2*id*nn(2)
      jump=1
      call GPFA(a,b,trigs(1,3),inc,jump,nn(3),lot,-is)
      n=0
      do i=1,id
         do j=1,id
            do k=1,id
               n=n+1
               c(k,j,i) = dcmplx(a(n),b(n))   ! MIZUHO-IR
            end do
         end do
      end do
      deallocate(a,b)
      return
      end
