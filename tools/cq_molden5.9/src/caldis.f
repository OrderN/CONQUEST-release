      subroutine caldis(calds,crot,fdist,cdist,tol,ianz)
      implicit double precision (a-h,o-z)
      parameter (bignum=1.0d12)
      parameter (numat1=20000)
      common /athlp/ iatoms, mxnat
      common /gamori/imap(numat1),cstand(3,numat1),istand(numat1),
     &               msucc
      dimension crot(3,numat1),fdist(numat1),cdist(numat1)
      dimension ianz(*)

c dist = the shortest distance found between a particular atom of crot
c        whilst screening it against all of the atoms of cstand
c dis  = the total sum over all atoms of dist(iatom)

      calds = 0.0d0

      natoms = iatoms
      if (iatoms.gt.numat1) natoms = numat1

      do k=1,natoms
         dist = bignum
         imap(k) = 0
         do l=1,natoms
            if((ianz(k).eq.istand(l)).and.
     &      ((fdist(k)-cdist(l))**2.lt.tol)) then
                distt = 0.0d0
                do j=1,3
                   distt = distt + 
     &             (crot(j,k)-cstand(j,l))**2
                end do
                if (distt.lt.dist) then
                   dist = distt
                   imap(k) = l
                endif
            endif
         end do
         calds = calds + dsqrt(dist)
      end do

      calds = dsqrt(calds)

      return
      end
