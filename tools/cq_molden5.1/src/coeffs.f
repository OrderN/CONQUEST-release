      subroutine coeffs(cc,ap,bp,itype1,jtype1)
      implicit double precision (a-h,o-z)
      dimension cc(192)
c
      cc(1) = 1.0d0
      if (jtype1.eq.1) goto 200
      cc(2) = bp
      cc(3) = cc(1)
      if (jtype1.eq.2) goto 200
      cc(4) = bp*cc(2)
      cc(5) = bp + cc(2)
      cc(6) = cc(3)
      if (jtype1.eq.3) goto 200
      cc(7) = bp*cc(4)
      cc(8) = bp*cc(5) + cc(4)
      cc(9) = bp + cc(5)
      cc(10) = cc(6)
c
  200 if (itype1.eq.1) return
      cc(11) = ap
      cc(12) = cc(1)
      if (jtype1.eq.1) goto 300
      cc(13) = ap*cc(2)
      cc(14) = ap + cc(2)
      cc(15) = cc(3)
      if (jtype1.eq.2) goto 300
      cc(16) = ap*cc(4)
      cc(17) = ap*cc(5) + cc(4)
      cc(18) = ap + cc(5)
      cc(19) = cc(6)
      if (jtype1.eq.3) goto 300
      cc(20) = ap*cc(7)
      cc(21) = ap*cc(8) + cc(7)
      cc(22) = ap*cc(9) + cc(8)
      cc(23) = ap + cc(9)
      cc(24) = cc(10)
c
  300 if (itype1.eq.2) return
      cc(25) = ap*cc(11)
      cc(26) = ap + cc(11)
      cc(27) = cc(12)
      if (jtype1.eq.1) goto 400
      cc(28) = ap*cc(13)
      cc(29) = ap*cc(14) + cc(13)
      cc(30) = ap + cc(14)
      cc(31) = 1.0d0
      if (jtype1.eq.2) goto 400
      cc(32) = ap*cc(16)
      cc(33) = ap*cc(17) + cc(16)
      cc(34) = ap*cc(18) + cc(17)
      cc(35) = ap + cc(18)
      cc(36) = 1.0d0
      if (jtype1.eq.3) goto 400
      cc(37) = ap*cc(20)
      cc(38) = ap*cc(21) + cc(20)
      cc(39) = ap*cc(22) + cc(21)
      cc(40) = ap*cc(23) + cc(22)
      cc(41) = ap + cc(23)
      cc(42) = 1.0d0
c
  400 if (itype1.eq.3) return
      cc(43) = ap*cc(25)
      cc(44) = ap*cc(26) + cc(25)
      cc(45) = ap + cc(26)
      cc(46) = 1.0d0
      if (jtype1.eq.1) goto 500
      cc(47) = ap*cc(28)
      cc(48) = ap*cc(29) + cc(28)
      cc(49) = ap*cc(30) + cc(29)
      cc(50) = ap + cc(30)
      cc(51) = 1.0d0
      if (jtype1.eq.2) goto 500
      cc(52) = ap*cc(32)
      cc(53) = ap*cc(33) + cc(32)
      cc(54) = ap*cc(34) + cc(33)
      cc(55) = ap*cc(35) + cc(34)
      cc(56) = ap + cc(35)
      cc(57) = 1.0d0
      if (jtype1.eq.3) goto 500
      cc(58) = ap*cc(37)
      cc(59) = ap*cc(38) + cc(37)
      cc(60) = ap*cc(39) + cc(38)
      cc(61) = ap*cc(40) + cc(39)
      cc(62) = ap*cc(41) + cc(40)
      cc(63) = ap + cc(41)
      cc(64) = 1.0d0
c
  500 if (itype1.eq.4) return
      cc(65) = ap*cc(43)
      cc(66) = ap*cc(44) + cc(43)
      cc(67) = ap*cc(45) + cc(44)
      cc(68) = ap + cc(45)
      cc(69) = 1.0d0
      if (jtype1.eq.1) goto 600
      cc(70) = ap*cc(47)
      cc(71) = ap*cc(48) + cc(47)
      cc(72) = ap*cc(49) + cc(48)
      cc(73) = ap*cc(50) + cc(49)
      cc(74) = ap + cc(50)
      cc(75) = 1.0d0
      if (jtype1.eq.2) goto 600
      cc(76) = ap*cc(52)
      cc(77) = ap*cc(53) + cc(52)
      cc(78) = ap*cc(54) + cc(53)
      cc(79) = ap*cc(55) + cc(54)
      cc(80) = ap*cc(56) + cc(55)
      cc(81) = ap + cc(56)
      cc(82) = 1.0d0
      if (jtype1.eq.3) goto 600
      cc(83) = ap*cc(58)
      cc(84) = ap*cc(59) + cc(58)
      cc(85) = ap*cc(60) + cc(59)
      cc(86) = ap*cc(61) + cc(60)
      cc(87) = ap*cc(62) + cc(61)
      cc(88) = ap*cc(63) + cc(62)
      cc(89) = ap + cc(63)
      cc(90) = 1.0d0
c
  600 if (itype1.eq.5) return
      cc(91) = ap*cc(65)
      cc(92) = ap*cc(66) + cc(65)
      cc(93) = ap*cc(67) + cc(66)
      cc(94) = ap*cc(68) + cc(67)
      cc(95) = ap + cc(68)
      cc(96) = 1.0d0
      if (jtype1.eq.1) goto 700
      cc(97) = ap*cc(70)
      cc(98) = ap*cc(71) + cc(70)
      cc(99) = ap*cc(72) + cc(71)
      cc(100) = ap*cc(73) + cc(72)
      cc(101) = ap*cc(74) + cc(73)
      cc(102) = ap + cc(74)
      cc(103) = 1.0d0
      if (jtype1.eq.2) goto 700
      cc(104) = ap*cc(76)
      cc(105) = ap*cc(77) + cc(76)
      cc(106) = ap*cc(78) + cc(77)
      cc(107) = ap*cc(79) + cc(78)
      cc(108) = ap*cc(80) + cc(79)
      cc(109) = ap*cc(81) + cc(80)
      cc(110) = ap + cc(81)
      cc(111) = 1.0d0
      if (jtype1.eq.3) goto 700
      cc(112) = ap*cc(83)
      cc(113) = ap*cc(84) + cc(83)
      cc(114) = ap*cc(85) + cc(84)
      cc(115) = ap*cc(86) + cc(85)
      cc(116) = ap*cc(87) + cc(86)
      cc(117) = ap*cc(88) + cc(87)
      cc(118) = ap*cc(89) + cc(88)
      cc(119) = ap + cc(89)
      cc(120) = 1.0d0
c
  700 if (itype1.eq.6) return

      cc(121) = ap*cc(91)
      cc(122) = ap*cc(92) + cc(91)
      cc(123) = ap*cc(93) + cc(92)
      cc(124) = ap*cc(94) + cc(93)
      cc(125) = ap*cc(95) + cc(94)
      cc(126) = ap + cc(95)
      cc(127) = 1.0d0
      if (jtype1.eq.1) goto 800

      cc(128) = ap*cc(97)
      cc(129) = ap*cc(98) + cc(97)
      cc(130) = ap*cc(99) + cc(98)
      cc(131) = ap*cc(100) + cc(99)
      cc(132) = ap*cc(101) + cc(100)
      cc(133) = ap*cc(102) + cc(101)
      cc(134) = ap + cc(102)
      cc(135) = 1.0d0
      if (jtype1.eq.2) goto 800

      cc(136) = ap*cc(104)
      cc(137) = ap*cc(105) + cc(104)
      cc(138) = ap*cc(106) + cc(105)
      cc(139) = ap*cc(107) + cc(106)
      cc(140) = ap*cc(108) + cc(107)
      cc(141) = ap*cc(109) + cc(108)
      cc(142) = ap*cc(110) + cc(109)
      cc(143) = ap + cc(110)
      cc(144) = 1.0d0
      if (jtype1.eq.3) goto 800

      cc(145) = ap*cc(112)
      cc(146) = ap*cc(113) + cc(112)
      cc(147) = ap*cc(114) + cc(113)
      cc(148) = ap*cc(115) + cc(114)
      cc(149) = ap*cc(116) + cc(115)
      cc(150) = ap*cc(117) + cc(116)
      cc(151) = ap*cc(118) + cc(117)
      cc(152) = ap*cc(119) + cc(118)
      cc(153) = ap + cc(119)
      cc(154) = 1.0d0
c
  800 if (itype1.eq.7) return
      cc(155) = ap*cc(121)
      cc(156) = ap*cc(122) + cc(121)
      cc(157) = ap*cc(123) + cc(122)
      cc(158) = ap*cc(124) + cc(123)
      cc(159) = ap*cc(125) + cc(124)
      cc(160) = ap*cc(126) + cc(125)
      cc(161) = ap + cc(126)
      cc(162) = 1.0d0
      if (jtype1.eq.1) return

      cc(163) = ap*cc(128)
      cc(164) = ap*cc(129) + cc(128)
      cc(165) = ap*cc(130) + cc(129)
      cc(166) = ap*cc(131) + cc(130)
      cc(167) = ap*cc(132) + cc(131)
      cc(168) = ap*cc(133) + cc(132)
      cc(169) = ap*cc(134) + cc(133)
      cc(170) = ap + cc(134)
      cc(171) = 1.0d0
      if (jtype1.eq.2) return

      cc(172) = ap*cc(136)
      cc(173) = ap*cc(137) + cc(136)
      cc(174) = ap*cc(138) + cc(137)
      cc(175) = ap*cc(139) + cc(138)
      cc(176) = ap*cc(140) + cc(139)
      cc(177) = ap*cc(141) + cc(140)
      cc(178) = ap*cc(142) + cc(141)
      cc(179) = ap*cc(143) + cc(142)
      cc(180) = ap + cc(143)
      cc(181) = 1.0d0
      if (jtype1.eq.3) return

      cc(182) = ap*cc(145)
      cc(183) = ap*cc(146) + cc(145)
      cc(184) = ap*cc(147) + cc(146)
      cc(185) = ap*cc(148) + cc(147)
      cc(186) = ap*cc(149) + cc(148)
      cc(187) = ap*cc(150) + cc(149)
      cc(188) = ap*cc(151) + cc(150)
      cc(189) = ap*cc(152) + cc(151)
      cc(190) = ap*cc(153) + cc(152)
      cc(191) = ap + cc(153)
      cc(192) = 1.0d0

      return
      end
