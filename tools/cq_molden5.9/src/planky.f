      subroutine planky(npts1,npts2,npts3,keywrd,defau)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxvalc=10)
      logical defau
      logical keyi,keyr,keyiv,keyrv,keyirv
      character*(*) keywrd
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      logical oscal,ostep,fine
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      common /srfhlp/edge,ctval(mxvalc),pxyz(3),nvalc,nspts,istyp
      common /grdhlp/ mx3d,mx3d2
      dimension d1(3),d2(3),c12(3)

      todeg = 45.0d0 / datan(1.0d0)
      iflag = 0
      imxflg = 0

      if (.not.keyrv(keywrd,'CENTER',px,py,pz)) then
          if (keyi(keywrd,'CENTER',i)) then
             if (i.gt.natoms) then
                call inferr('CENTER: Atom Nr. > Nr. of Atoms !',1)
                return
             endif
             iplat = 0
             iflag = 1
             px = xyz(1,i)
             py = xyz(2,i)
             pz = xyz(3,i)
          endif
      else
          iplat = 0
          iflag = 1
      endif

      if (.not.keyrv(keywrd,'LINE',cx,cy,cz)) then
          if (keyi(keywrd,'LINE',i)) then
             if (i.gt.natoms) then
                call inferr('LINE: Atom Nr. > Nr. of Atoms !',1)
                return
             endif
             iplat = 0
             iflag = 1
             cx = xyz(1,i) - px
             cy = xyz(2,i) - py
             cz = xyz(3,i) - pz
          endif
      else
          iplat = 0
          iflag = 1
      endif

      if (keyiv(keywrd,'PLANE',iat1,iat2,iat3)) then
          if (iflag.eq.1) then
             call inferr('CENTER/LINE and Plane EXCLUSIVE!',1)
             return
          endif
          if((iat1.gt.natoms).or.(iat2.gt.natoms).or.(iat3.gt.natoms))
     &    then
             call inferr('PLANE: Atom nr. > number of atoms !',1)
             return
          endif
          call parpla(iat1,iat2,iat3,istat)
          if (istat.eq.1) print*,'Three atoms lie on a line !'
          iflag = 2
      endif

      if (iflag.eq.0.and.defau) call defpc

      if (keyirv(keywrd,'ROT',iata,iatb,degree)) then
          if (iflag.eq.1.or.iflag.eq.0) then
            if (iflag.eq.1) 
     &         call inferr('Center/Line and Rot are EXCLUSIVE',1)
            if (iflag.eq.0) 
     &         call inferr('ROT used before PLANE !',1)
            return
          endif
          if (((iata.ne.iat1).and.(iata.ne.iat2).and.(iata.ne.iat3))
     &    .or.((iatb.ne.iat1).and.(iatb.ne.iat2).and.(iatb.ne.iat3)))
     &    then
             call inferr('ROT: atom used not used with Plane !',1)
             return
          endif
          radian = degree / todeg
          do i=1,3
            d1(i) = (xyz(i,iatb) - xyz(i,iata))
          end do
          d2(1) = cx
          d2(2) = cy
          d2(3) = cz
          px = xyz(1,iata) + 0.5d0*d1(1)
          py = xyz(2,iata) + 0.5d0*d1(2)
          px = xyz(3,iata) + 0.5d0*d1(3)
          d2norm = vlen(d2)
          do i=1,3
            d2(i) = d2(i)/d2norm
          end do
          call crprod(d1,d2,c12)
          c12nrm = vlen(c12)
          do i=1,3
            c12(i) = c12(i)/c12nrm
          end do
          cx = dcos(radian)*d2(1) - dsin(radian)*c12(1)
          cy = dcos(radian)*d2(2) - dsin(radian)*c12(2)
          cz = dcos(radian)*d2(3) - dsin(radian)*c12(3)
          iplat = 0
      endif

      if (keyr(keywrd,'LIFT',rlift)) call dolift(rlift)

      if (.not.keyr(keywrd,'EDGE',r(1))) then
         if (defau) call defrad(.false.)
      else
         r(2) = r(1)
         r(3) = r(1)
      endif

      if (.not.keyr(keywrd,'EDZ',r(3))) then
         r(3) = r(1)
      endif

      if (.not.keyr(keywrd,'EDY',r(2))) then
         r(2) = r(1)
      endif

      if (.not.keyr(keywrd,'EDX',r(1))) idum = 1

      if (r(1).lt.0.01d0) then
         r(1) = 3.0
         r(2) = r(1)
         r(3) = r(1)
      endif

      if (keyi(keywrd,'NPTSX',npts1)) then
         if (npts1.gt.mx3d) imxflg = 1
      endif
      if (keyi(keywrd,'NPTSY',npts2)) then
         if (npts2.gt.mx3d) imxflg = 1
      endif
      if (keyi(keywrd,'NPTSZ',npts3)) then
         if (npts3.gt.mx3d) imxflg = 1
      endif

      if (imxflg.eq.1) then
         nptsmx = npts1
         if (npts2.gt.nptsmx) nptsmx = npts2
         if (npts3.gt.nptsmx) nptsmx = npts3
         call allgrd(nptsmx)
         if (nptsmx.gt.mx3d) then
             call inferr('Exceeding maximum points !',1)
             if (npts1.gt.mx3d) npts1 = mx3d
             if (npts2.gt.mx3d) npts2 = mx3d
             if (npts3.gt.mx3d) npts3 = mx3d
         endif
      endif

      oneold = one
      if (.not.keyr(keywrd,'PHASE',one)) then
         if (index(keywrd,'PHASE').ne.0) one = -1.d0*one
      else
         if (one.eq.0.0d0) one = -1.d0*oneold
      endif

      if (index(keywrd,'ALIGN').ne.0) then
         px = pxyz(1)
         py = pxyz(2)
         pz = pxyz(3)
         iplat = 0
      endif

      write(iun3,10) px, py, pz, cx, cy, cz, r(1), r(2), r(3)
10    format(9x,'  CENTER OF GRAPH           AXIS OF GRAPH       RADII',
     1/9X,      'PX      PY      PZ        CX      CY      CZ',
     24X,      'EDX   EDY   EDY',
     3/,5X,3F8.3,F10.3,2F8.3,3F6.1,///)

      return
      end

      logical function keyi(keyl,keyw,ival)
      implicit double precision (a-h,o-z)
      character*(*) keyl,keyw
      
      keyi = .false.
      l = linlen(keyw)
      i = index(keyl,keyw)
      if (i.ne.0) then
          if ((i.eq.1.or.keyl(i-1:i-1).eq.' ').and.
     &         keyl(i+l:i+l).eq.'=') then
             ival = reada(keyl,i,len(keyl))
             keyi = .true.
          endif
      endif

      return
      end

      logical function keyr(keyl,keyw,rval)
      implicit double precision (a-h,o-z)
      character*(*) keyl,keyw
      
      keyr = .false.
      i = index(keyl,keyw)
      l = linlen(keyw)
      if (i.ne.0) then
          if ((i.eq.1.or.keyl(i-1:i-1).eq.' ').and.
     &         keyl(i+l:i+l).eq.'=') then
             rval = reada(keyl,i,len(keyl))
             keyr = .true.
          endif
      endif

      return
      end

      logical function keyrv(keyl,keyw,r1,r2,r3)
      implicit double precision (a-h,o-z)
      character*(*) keyl,keyw
      
      keyrv = .false.
      kend = len(keyl)

      i = index(keyl,keyw)
      if (i.eq.0) return
      i = i + len(keyw)
      j = index(keyl(i:kend),'(')
      if (j.eq.0) return
      do k=0,j-2
         if (keyl(i+k:i+k).ne.' '.and.keyl(i+k:i+k).ne.'=') return
      end do
      i = i + j 
      r1 = reada(keyl,i,len(keyl))
      j = index(keyl(i:kend),',')
      if (j.eq.0) return
      i = i + j 
      r2 = reada(keyl,i,len(keyl))
      j = index(keyl(i:kend),',')
      if (j.eq.0) return
      i = i + j 
      r3 = reada(keyl,i,len(keyl))

      keyrv = .true.
      return
      end

      integer function keyr3v(keyl,keyw,r,n)
      implicit double precision (a-h,o-z)
      character*(*) keyl,keyw
      dimension r(*)
      
      keyr3v = 0
      kend = len(keyl)

      i = index(keyl,keyw)
      if (i.eq.0) return
      i = i + len(keyw)
      j = index(keyl(i:kend),'(')
      if (j.eq.0) return
      do k=0,j-2
         if (keyl(i+k:i+k).ne.' '.and.keyl(i+k:i+k).ne.'=') return
      end do

      do l=1,n
         i = i + j 
         r(l) = reada(keyl,i,len(keyl))
         keyr3v = keyr3v + 1
         j = index(keyl(i:kend),',')
         if (j.eq.0) return
      end do

      return
      end

      logical function keyiv(keyl,keyw,i1,i2,i3)
      implicit double precision (a-h,o-z)
      character*(*) keyl,keyw
      
      keyiv = .false.
      kend = len(keyl)

      i = index(keyl,keyw)
      if (i.eq.0) return
      i = i + len(keyw)
      j = index(keyl(i:kend),'(')
      if (j.eq.0) return
      do k=0,j-2
         if (keyl(i+k:i+k).ne.' '.and.keyl(i+k:i+k).ne.'=') return
      end do
      i = i + j 
      i1 = reada(keyl,i,len(keyl))
      j = index(keyl(i:kend),',')
      if (j.eq.0) return
      i = i + j 
      i2 = reada(keyl,i,len(keyl))
      j = index(keyl(i:kend),',')
      if (j.eq.0) return
      i = i + j 
      i3 = reada(keyl,i,len(keyl))

      keyiv = .true.
      return
      end

      logical function keyirv(keyl,keyw,i1,i2,r)
      implicit double precision (a-h,o-z)
      character*(*) keyl,keyw
      
      keyirv = .false.
      kend = len(keyl)

      i = index(keyl,keyw)
      if (i.eq.0) return
      i = i + len(keyw)
      j = index(keyl(i:kend),'(')
      if (j.eq.0) return
      do k=0,j-2
         if (keyl(i+k:i+k).ne.' '.and.keyl(i+k:i+k).ne.'=') return
      end do
      i = i + j 
      i1 = reada(keyl,i,len(keyl))
      j = index(keyl(i:kend),',')
      if (j.eq.0) return
      i = i + j 
      i2 = reada(keyl,i,len(keyl))
      j = index(keyl(i:kend),',')
      if (j.eq.0) return
      i = i + j 
      r = reada(keyl,i,len(keyl))

      keyirv = .true.
      return
      end
