      subroutine rota(alfa)
      implicit double precision (a-h,o-z)
      common/rotxyz/amat(3,3),bmat(3,3),cmat(3,3),rmat(3,3),todeg

c--   alfa expected in degrees
      rada = alfa / todeg
      amat(1,1) =  dcos(rada)
      amat(1,2) =  dsin(rada)
      amat(2,1) = -dsin(rada)
      amat(2,2) =  dcos(rada)

      return
      end
