      double precision function distot()
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      common /athlp/ iatoms, mxnat
      common /gamori/imap(numat1),cstand(3,numat1),istand(numat1),
     &               msucc
      dimension b(3,numat1)

      call rotcor(b)

      distot = 0.0d0

      do i=1,iatoms
          d1 = b(1,i) - cstand(1,i)
          d2 = b(2,i) - cstand(2,i)
          d3 = b(3,i) - cstand(3,i)
          dist = d1*d1 + d2*d2 + d3*d3
          distot = distot + dist
      end do

      return
      end
