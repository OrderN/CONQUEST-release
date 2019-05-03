      subroutine setbas(oooke)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      logical oooke
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      common /coord / xyz(3,numatm)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common/bond/nbas(numatm),ipol(numatm)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      dimension npr(numatm),ltot1(numatm)

c     initialise

      oooke = .true.

      do i=1,natoms
        ltot1(i)=0
        npr(i)=0
      end do
c
c  set ltot1 per atom sum off maximum angular quantum number per shell +1
c      npr   total number off primitives per atom
c
c  TABEL npr and (ltot1)
c
c  nbas   basis    h-he   li-ne  na-cl
c--------------------------------------
c  0      sto3g     3 (1)   6 (3)  9 (5)
c  1      321g      3 (2)   6 (5)  9 (7)
c  2      431g      4       8     12
c  2      431g*     5       9     13
c  3      631g      4      10     16
c  3      631g*     5      11     17
c-------------------------------------
c
      do i=1,nshell
         ltot1(jan(i)) = ltot1(jan(i)) + shellt(i) + 1
         npr(jan(i))   = npr(jan(i)) + shelln(i)
      end do
c
c     determine type off basis per atom
c             nbas  ipol
c     sto3g    0     0
c     321g     1     0
c     431g     2     0
c     431g*    2     1
c     631g     3     0
c     631g*    3     1
      do 40 i=1,natoms
      goto (50,50,60,60,60,60,60,60,60,60,
     +      70,70,70,70,70,70,70) nat(i)
c.....h-he.................................
50    if (npr(i).eq.3) then
         nbas(i) = 0
         ipol(i) = 0
         if (ltot1(i).eq.2) nbas(i) = 1
         goto 40
      endif
      if (npr(i).eq.4) then
         ipol(i) = 0
         nbas(i) = 2
         goto 40
      endif
      if (npr(i).eq.5) then
         ipol(i) = 1
         nbas(i) = 2
         goto 40
      endif
      goto 200 
c.....li-ne............................
60    if (npr(i).eq.6) then
         nbas(i) = 0
         ipol(i) = 0
         if (ltot1(i).eq.5) nbas(i) = 1
         goto 40
      endif
      if (npr(i).eq.8) then
         ipol(i) = 0
         nbas(i) = 2
         goto 40
      endif
      if (npr(i).eq.9) then
         ipol(i) = 1
         nbas(i) = 2
         goto 40
      endif
      if (npr(i).eq.10) then
         ipol(i) = 0
         nbas(i) = 3
         goto 40
      endif
      if (npr(i).eq.11) then
         ipol(i) = 1
         nbas(i) = 3
         goto 40
      endif
      goto 200
c.....na-cl............................
70    if (npr(i).eq.9) then
         nbas(i) = 0
         ipol(i) = 0
         if (ltot1(i).eq.7) nbas(i) = 1
         goto 40
      endif
      if (npr(i).eq.12) then
         ipol(i) = 0
         nbas(i) = 2
         goto 40
      endif
      if (npr(i).eq.13) then
         ipol(i) = 1
         nbas(i) = 2
         goto 40
      endif
      if (npr(i).eq.16) then
         ipol(i) = 0
         nbas(i) = 3
         goto 40
      endif
      if (npr(i).eq.17) then
         ipol(i) = 1
         nbas(i) = 3
         goto 40
      endif
      goto 200
c
40    continue
      return
200   call inferr('When using BONDS only the following basissets',0)
      call inferr('on atoms from H to Cl are allowed ;',0)
      call inferr('Bonds: only sto3g 3-21g 4/6-31g(*/**)',1)

      oooke = .false.

      end
