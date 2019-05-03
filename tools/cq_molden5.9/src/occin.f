      subroutine occin(i,gocc,ndim)
      implicit double precision (a-h,o-z)
      character*320 keywrd,keyori,string
      common /keywrd/ keywrd,keyori
      character*30  cunit
      dimension gocc(*)

      ind = i + 4
      i = index(keywrd(ind:320),'=')
      i1 = index(keywrd(i+ind:320),'(')
      i2 = index(keywrd(i+ind:320),')')
      string(1:i2-i1-1) = keywrd(i+ind+i1:i+ind+i2-1)
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
         cunit(1:i-1) = string(1:i-1)
         string(1:l-i) = string(i+1:l) 
         l = l - i
      endif
      ind = index(cunit(1:i-1),'/')
      jnd = index(cunit(1:i-1),'-')
      if (jnd.eq.0) then
         ibeg = int(reada(cunit(1:ind-1)//' ',1,len(cunit)))
         iend = ibeg
      else
         ibeg = int(reada(cunit(1:jnd-1)//' ',1,len(cunit)))
         iend = int(reada(cunit(jnd+1:ind-1)//' ',1,len(cunit)))
      endif

      if (iend.lt.ibeg)
     &   call inferr('OCCU: orbital label1 exceeds label2 ',1)

      if (ibeg.gt.ndim.or.iend.gt.ndim) 
     &   call inferr('OCCU: orbital label > No. of Orbitals',1)

      val = reada(cunit(ind+1:i-1)//' ',1,len(cunit))

      do k=ibeg,iend
         gocc(k) = val
      end do

      if (l.ne.0) goto 2400

      return
      end
