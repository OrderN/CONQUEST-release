      subroutine parstp
      implicit double precision (a-h,o-z)
      logical ostep,fine,oscal
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      dimension setval(20)

      setval(1)  = 1.d-4
      setval(2)  = 2.5d-4
      setval(3)  = 5.d-4
      setval(4)  = 1.d-3
      setval(5)  = 2.5d-3
      setval(6)  = 5.d-3
      setval(7)  = 1.d-2
      setval(8)  = 2.5d-2
      setval(9)  = 5.d-2
      setval(10) = 1.d-1
      setval(11) = 2.5d-1
      setval(12) = 5.d-1
      setval(13) = 1.d 0
      setval(14) = 2.5d 0
      setval(15) = 5.d 0
      setval(16) = 1.d 1
      setval(17) = 2.5d 1
      setval(18) = 5.d 1
      setval(19) = 1.d 2
      setval(20) = 2.5d 2

      do 215 i=1,20
  215        if(pmax*cut.lt.setval(i)) goto 216
      if(i.eq.21) i=20
  216 step=setval(i)*0.05d0

      if (ostep) step = stepf

      if(fine) then
          dellh=0.249999d0
      else
          dellh=1.d0
      endif
      step=step*dellh

      if (pmax.eq.0.0d0) then
         scale =1.0d0
      else
         scale = -0.5d0/pmax
      endif
      if (oscal) scale = scalf

      return
      end
