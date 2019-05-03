      subroutine rotc(gamma)
      implicit double precision (a-h,o-z)
      common/rotxyz/amat(3,3),bmat(3,3),cmat(3,3),rmat(3,3),todeg

c--   beta expected in degrees
      radc = gamma / todeg
      cmat(2,2) =  dcos(radc)
      cmat(2,3) =  dsin(radc)
      cmat(3,2) = -dsin(radc)
      cmat(3,3) =  dcos(radc)

      return
      end
