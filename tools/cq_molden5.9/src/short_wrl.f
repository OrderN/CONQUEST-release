c This program removes the redundancy in wrl-files from marching cube alg.
c and thus shortens them by about a factor 4. The program is called by 
c       short_wrl  filename(incl. path)
c filename must have extension wrl.
c A new file is written which end on '+.wrl'
c This program was written by 
c Dr. Andreas Klamt, COSMOlogic GmbH&CoKG, Leverkusen, Germany (17/08/2000)
c It may be freely used without removal of the copyright statement.
c Adapted by G.Schaftenaar to include processing of molden VRML orbitals
c     maxpts is the maximum number of non-redundant vertices 
c     Thus the number of redundant vertices may be about 5*maxpts
      parameter (maxpts=500000)
      dimension x(3,maxpts),id(maxpts),r(3)
      logical isz(maxpts)
      character*150 zeile
c
      call getarg(1,zeile)
      ind=index(zeile,'.wrl')
      if (ind.eq.0) stop 'No wrl-file!!!'
      open (7,file=zeile, status='old')
      open (8,file=zeile(:ind)//'+.wrl', status='unknown')
5     read(7,'(a)') zeile
      ind=index(zeile,'                                            ')-1
      write(8,'(a)') zeile(:ind)
      if(index(zeile,'coord Coordinate { point [').eq.0) go to 5
      l=0
      do 30 i=1,500000
        isz(i)=.false.
        k=l+1
        read(7,'(a)') zeile
        read(zeile,*,err=50) (x(ix,k),ix=1,3)
        do 10 j=1,l
          if((x(1,j)-x(1,k))**2 .gt. 1.d-10) go to 10
          if((x(2,j)-x(2,k))**2 .gt. 1.d-10) go to 10
          if((x(3,j)-x(3,k))**2 .gt. 1.d-10) go to 10
          id(i)=j-1
          go to 30
10      continue
        id(i)=l
        l=l+1
        isz(i)=.true.
        write(8,'(3f12.5,1h,)') (x(ix,l),ix=1,3)
30    continue
50    np=i-1
      ind=index(zeile,'                                            ')-1
      write(8,'(a)') zeile(:ind)
51    read(7,'(a)') zeile
      ind=index(zeile,'                                            ')-1
      write(8,'(a)') zeile(:ind)
      if(index(zeile,'normal Normal { vector [').eq.0) go to 51
      do i=1,np
        read(7,*)  r
        if(isz(i)) write(8,'(3f10.6,1h,)') r
      end do
61    read(7,'(a)') zeile
      ind=index(zeile,'                                            ')-1
      write(8,'(a)') zeile(:ind)
      if(index(zeile,'CIndex [').ne.0) go to 72
      if(index(zeile,'color Color { color [').eq.0) go to 61
      do i=1,np
        read(7,*)  r
        if(isz(i)) write(8,'(3f7.3,1h,)') r
      end do
71    read(7,'(a)') zeile
      ind=index(zeile,'                                            ')-1
      write(8,'(a)') zeile(:ind)
      if(index(zeile,'CIndex [').eq.0) go to 71
72    print *,'Vertices reduced from',np,' to',l
      print *,'by Short_wrl, a program from COSMOlogic ',
     &      '(www.cosmologic.de)'
73    read(7,'(a)') zeile 
      if (index(zeile,']').ne.0) then
        ind=index(zeile,
     &    '                                            ')-1
        write(8,'(a)') zeile(:ind)
        go to 81
      endif
      read(zeile,*)  i1,i2,i3
      write(8,'(3(i7,1h,)3h-1,)') id(i1+1),id(i2+1),id(i3+1)
      goto 73
81    read(7,'(a)',end=91) zeile
      ind=index(zeile,'                                            ')-1
      write(8,'(a)') zeile(:ind)
      if(index(zeile,'CIndex [').eq.0) go to 81
82    read(7,'(a)') zeile 
      if (index(zeile,']').ne.0) then
        ind=index(zeile,
     &    '                                            ')-1
        write(8,'(a)') zeile(:ind)
        go to 83
      endif
      read(zeile,*)  i1,i2,i3
      write(8,'(3(i7,1h,)3h-1,)') id(i1+1),id(i2+1),id(i3+1)
      goto 82
83    read(7,'(a)',end=91) zeile
      ind=index(zeile,'                                            ')-1
      write(8,'(a)') zeile(:ind)
      goto 83
91    end
