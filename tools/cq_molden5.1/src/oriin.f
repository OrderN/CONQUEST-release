      subroutine oriin(i)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /keywrd/ keywrd,keyori
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common/orihlp/ori(numatm),iuser(numatm),oalpha(numatm),
     $              obeta(numatm),norien,ibal
      character*320 keywrd,keyori,string
      character*30  cunit

      norien = 0
      i1 = index(keywrd(i+1:320),'(')
      i2 = index(keywrd(i+1:320),')')
      string(1:i2-i1-1) = keywrd(i+i1+1:i+i2-1)
      l = i2 - i1 - 1
      j = 1
2300  if (string(j:j).eq.' ') then
        if (j.le.l-1) then
           string(j:l-1) = string(j+1:l)
        endif
        l = l - 1 
      else
        j = j + 1
      endif
      if (j.ne.l+1) goto 2300
2400  i = index(string(1:l),',')
      if (i.eq.0) then
         cunit(1:l) = string(1:l)
         i = l + 1
         l = 0
      else
         cunit(1:i-1)  = string(1:i-1)
         string(1:l-i) = string(i+1:l) 
         l = l - i
      endif
c---- cunit now contains item between "," or "( ," or ", )"
      ind = index(cunit(1:i-1),'/')
      jnd = ind + index(cunit(ind+1:i-1),'/')
      if (ind.eq.0) then
         k = int(reada(cunit(1:i-1)//' ',1,len(cunit)))
         if (nat(k).ne.8.and.nat(k).ne.9.and.nat(k).ne.16
     &      .and.nat(k).ne.17) 
     &      call inferr('ORIENT: allowed atoms are O,F,S,CL',1)
         norien        = norien + 1
         ori(norien)   = k
         iuser(norien) = 0
      elseif (jnd.gt.ind) then
         k = int(reada(cunit(1:ind-1)//' ',1,len(cunit)))
         if (nat(k).ne.8.and.nat(k).ne.9.and.nat(k).ne.16
     &      .and.nat(k).ne.17) 
     &      call inferr('ORIENT: allowed atoms are O,F,S,CL',1)
         norien         = norien + 1
         ori(norien)    = k
         iuser(norien)  = 1
         oalpha(norien) = reada(cunit(ind+1:jnd-1)//' ',1,len(cunit))
         obeta(norien)  = reada(cunit(jnd+1:i-1)//' ',1,len(cunit))
      else
         call inferr('error in ORIENT directive',1)
      endif
      if (l.ne.0) goto 2400

      return
      end
