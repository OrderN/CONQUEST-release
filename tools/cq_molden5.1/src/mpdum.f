      subroutine mpolefit(idip,nesp,esp,connl,dx,dy,dz,iz,dmachg,
     &                    keywrd,ichadd)
      implicit double precision (a-h,o-z)
      logical dmachg
      dimension esp(*), connl(3,*)
      character*(*) keywrd

      print*,'!!!!!!!!!!! This executable has NO MPFIT !!!!!!'
      call inferr('This executable has NO MPFIT',0)
      idum = 1

      return
      end
