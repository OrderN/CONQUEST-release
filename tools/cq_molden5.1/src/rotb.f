      subroutine rotb(beta)
      implicit double precision (a-h,o-z)
      common/rotxyz/amat(3,3),bmat(3,3),cmat(3,3),rmat(3,3),todeg

c--   beta expected in degrees
      radb = beta / todeg
      bmat(1,1) =  dcos(radb)
      bmat(1,3) =  dsin(radb)
      bmat(3,1) = -dsin(radb)
      bmat(3,3) =  dcos(radb)

      return
      end
