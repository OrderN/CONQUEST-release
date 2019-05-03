      subroutine calct(t,cset,ijtyp)
      implicit double precision (a-h,o-z)
      dimension t(45)

      if (ijtyp.le.1) return
      t(3) = cset
      if (ijtyp.eq.2) goto 10
      t(4) = t(3)
      t(6) = t(3) * t(3)
      if (ijtyp.eq.3) goto 10
      t(8) = t(3) * (2.d0 * t(3) + t(4))
      t(10) = t(3) * t( 6)
      if (ijtyp.eq.4) goto 10
      t(11) = 3.d0 * t(3) * t( 4)
      t(13) = t(3) * (3.d0 * t( 6) + t( 8))
      t(15) = t(3) * t(10)
      if (ijtyp.eq.5) goto 10
      t(17) = t(3) * (4.d0 * t( 8) + t(11))
      t(19) = t(3) * (4.d0 * t(10) + t(13))
      t(21) = t(3) * t(15)
      if (ijtyp.eq.6) goto 10
      t(22) = 5.d0 * t(3) * t(11)
      t(24) = t(3) * (5.d0 * t(13) + t(17))
      t(26) = t(3) * (5.d0 * t(15) + t(19))
      t(28) = t(3) * t(21)
      if (ijtyp.eq.7) goto 10
      t(30) = t(3) * (6.d0 * t(17) + t(22))
      t(32) = t(3) * (6.d0 * t(19) + t(24))
      t(34) = t(3) * (6.d0 * t(21) + t(26))
      t(36) = t(3) * t(28)
      if (ijtyp.eq.8) goto 10
      t(37) = 7.d0 * t(3) * t(22)
      t(39) = t(3) * (7.d0 * t(24) + t(30))
      t(41) = t(3) * (7.d0 * t(26) + t(32))
      t(43) = t(3) * (7.d0 * t(28) + t(34))
      t(45) = t(3) * t(36)

  10  continue

      return
      end
