      program unmullik

c written by g.schaftenaar
c to convert the unformatted GRAPH file to a formatted
c FGRAPH file , which can be read by a VAX

      implicit double precision (a-h,o-z)

      parameter (maxhev=75 , maxlit=75 )
      parameter (numatm=maxhev+maxlit)
      parameter (maxorb=4*maxhev+maxlit)
      parameter (morb2=maxorb**2)
      dimension h(morb2),c(morb2)
      common
     1       /molkst/ numat,nat(numatm),nfirst(numatm),nmidle(numatm),
     2                nlast(numatm), norbs, nelecs,nalpha,nbeta,
     3                nclose,nopen,ndumy,fract
      common /expont/ zs(numatm),zp(numatm),zd(numatm)
      dimension xyz(3,numatm)

c     First do the reading

      call filopn(13,'GRAPH','UNFORMATTED','UNKNOWN')
      rewind(13)

      read(13)numat,norbs,nelecs,((xyz(i,j),j=1,numat),i=1,3)
      if (numat.gt.numatm.or.norbs.gt.maxorb) stop 
     1    'Limits 150 atoms 375 orbitals'
      read(13)(nlast(i),nfirst(i),i=1,numat)
      read(13)(zs(i),i=1,numat),(zp(i),i=1,numat),
     1         (zd(i),i=1,numat),(nat(i),i=1,numat)
      linear=norbs*norbs
      read(13)(c(i),i=1,linear)
      read(13)(h(i),i=1,linear)

c     Now do the writing

      call filopn(14,'FGRAPH','FORMATTED','UNKNOWN')
      write(14,*)numat,norbs,nelecs,((xyz(i,j),j=1,numat),i=1,3)
      write(14,*)(nlast(i),nfirst(i),i=1,numat)
      write(14,*)(zs(i),i=1,numat),(zp(i),i=1,numat),
     1         (zd(i),i=1,numat),(nat(i),i=1,numat)
      write(14,*)(c(i),i=1,linear)
      write(14,*)(h(i),i=1,linear)

      end
      subroutine filopn(iunit, filnam, fmt, stat)
      character*(*) filnam, fmt, stat
      character*80 envfil
C
      call getenv(filnam,envfil)
      i = 80
   10 if (envfil(i:i) .ne. ' ') goto 20
      i = i - 1
      if (i .gt. 0) goto 10
   20 if (i .eq. 0) then
        open(unit=iunit, file=filnam, form=fmt, status=stat)
      else
        open(unit=iunit, file=envfil(1:i), form=fmt, status=stat)
      endif
      return
      end
