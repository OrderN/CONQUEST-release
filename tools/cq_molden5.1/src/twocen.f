      subroutine twocen(t2,fac1,strt,t,fac2,ijtyp)
      implicit double precision (a-h,o-z)
      dimension t(45),t2(9)
      dimension h(9)

      h(1) = strt
      if (ijtyp.le.1) goto 100
      h(2) = fac2*fac1*h(1)
      if (ijtyp.eq.2) goto 90
      h(3) = fac2*(fac1*h(2) - h(1))
      if (ijtyp.eq.3) goto 80
      h(4) = fac2*(fac1*h(3) - 2.d0*h(2))
      if (ijtyp.eq.4) goto 70
      h(5) = fac2*(fac1*h(4) - 3.d0*h(3))
      if (ijtyp.eq.5) goto 60
      h(6) = fac2*(fac1*h(5) - 4.d0*h(4))
      if (ijtyp.eq.6) goto 50
      h(7) = fac2*(fac1*h(6) - 5.d0*h(5))
      if (ijtyp.eq.7) goto 40
      h(8) = fac2*(fac1*h(7) - 6.d0*h(6))
      if (ijtyp.eq.8) goto 30
      h(9) = fac2*(fac1*h(8) - 7.d0*h(7))

      t2(9) = h(1)*t(37)+h(3)*t(39)+h(5)*t(41)+h(7)*t(43)+h(9)*t(45)
30    t2(8) = h(2)*t(30) + h(4)*t(32) + h(6)*t(34) + h(8)*t(36)
40    t2(7) = h(1)*t(22) + h(3)*t(24) + h(5)*t(26) + h(7)*t(28)
50    t2(6) = h(2)*t(17) + h(4)*t(19) + h(6)*t(21)
60    t2(5) = h(1)*t(11) + h(3)*t(13) + h(5)*t(15)
70    t2(4) = h(2)*t(8)  + h(4)*t(10)
80    t2(3) = h(1)*t(4)  + h(3)*t(6)
90    t2(2) = h(2)*t(3)
100   t2(1) = h(1)

      return
      end
