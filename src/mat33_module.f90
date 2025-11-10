module mat33
 use numbers
 implicit none
 contains
  subroutine mat33_prod(aa,bb,cc)
  implicit none
  real(double), dimension(3,3),intent(in) :: aa,bb
  real(double), dimension(3,3),intent(out) :: cc
  integer :: ii, jj, kk

  cc = zero
   do ii=1,3
    do jj=1,3
     do kk = 1,3
       cc(ii,jj)=cc(ii,jj)+aa(ii,kk)*bb(kk,jj)
     enddo
    enddo
   enddo
  end subroutine mat33_prod

!inver3
 subroutine inver3(a,ainv)
  implicit none
  integer :: i,j,k
  real(double)::a(3,3),ainv(3,3),c,d
  real(double)::cmax

  d=deter3(a)
  if(abs(d) < very_small) then
    write(*,9001) abs(d), very_small
    9001 format(' ***** ERROR: (inver3) abs(det(a))',d9.2,'< ',e8.3)
    stop ' ***** ERROR: (inver3) abs(det(a))<very_small'
  endif
  ainv(1,1)=(a(2,2)*a(3,3)-a(3,2)*a(2,3))/d
  ainv(2,1)=(a(2,3)*a(3,1)-a(3,3)*a(2,1))/d
  ainv(3,1)=(a(2,1)*a(3,2)-a(3,1)*a(2,2))/d
  ainv(1,2)=(a(3,2)*a(1,3)-a(1,2)*a(3,3))/d
  ainv(2,2)=(a(3,3)*a(1,1)-a(1,3)*a(3,1))/d
  ainv(3,2)=(a(3,1)*a(1,2)-a(1,1)*a(3,2))/d
  ainv(1,3)=(a(1,2)*a(2,3)-a(2,2)*a(1,3))/d
  ainv(2,3)=(a(1,3)*a(2,1)-a(2,3)*a(1,1))/d
  ainv(3,3)=(a(1,1)*a(2,2)-a(2,1)*a(1,2))/d
!check
  cmax=zero
  do i=1,3
   do j=1,3
        c=zero
    do k=1,3
        c=c+ainv(i,k)*a(k,j)
    enddo
    if(i.eq.j) c=c-one
    cmax=max(cmax,abs(c))
   enddo
  enddo
  if(cmax > very_small) then
   write(*,9002) cmax, very_small
   9002 format(' ***** ERROR: (inver3) cmax=',d10.3,'>',e8.3)
   stop ' ***** ERROR: (inver3) cmax>very_small'
  endif
  return
 end subroutine inver3
! -- deter3 ----------
! calculate determinant of rank 3
 real(double) function deter3(a)
  implicit none
  real(double),intent(in):: a(3,3)
  deter3=a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3)) &
        +a(2,1)*(a(3,2)*a(1,3)-a(1,2)*a(3,3)) &
        +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
  return
 end function deter3
 
! write_mat33
 subroutine print_mat33(a)
 implicit none
 real(double), intent(in) :: a(3,3)
 integer :: ii
 do ii=1,3
  write(*,101) a(ii,1:3)
  101 format(10x,3f15.10)
 enddo
 return
 end subroutine print_mat33
end module mat33
