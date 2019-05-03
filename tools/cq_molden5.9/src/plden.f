      subroutine plded(ndim1,ndim2,scale,icells,adjus,idisml,
     &                 dens,ix,iy,rz)
      implicit double precision (a-h,o-z)

c THIS IS REALLY plden

      parameter (mxel=100)
      parameter (numatm=2000)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)

      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      common /eul/   ca,cb,sa,sb,ccab,scab,ssab,scba
      common /atopla/ xsym(numatm),ysym(numatm),zsym(numatm),
     &                isym(numatm)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      real xx, yy
      character str*100
      integer*2 ipl1(8)
      integer*2 ipl2(8)
      integer*2 ixx(4)
      dimension temp1(3),temp2(3),temp3(3),zvect(3)
      dimension dens(*),ix(*),iy(*),rz(*)

      npoly = 4
      zvect(1) = -0.5d0
      zvect(2) = -0.5d0
      zvect(3) = -1.5d0

      npts = ndim1*ndim2
      scald = -1.0d0*scale*ndim1
      scle1  = ndim1*1.2d0 * r(2) / r(1)
      scle2  = ndim2*1.2d0
      nmed1 = ndim1/2
      nmed2 = ndim2/2
      if (ndim1-nmed1*2.eq.1) nmed1 = nmed1 + 1
      if (ndim2-nmed2*2.eq.1) nmed2 = nmed2 + 1

      call gethei(ihigh)

      ij = 0
      do i=1,ndim1
         do j=1,ndim2
           ij = ij + 1
           ix(ij) = (ccab*(i-nmed1)/scle1 + scab*(j-nmed2)/scle2 -
     &              sb*dens(ij)*scald/scle1 + 0.5d0)*ihigh
           iy(ij) = (0.6d0-(-sa*(i-nmed1)/scle1 + ca*(j-nmed2)/scle2))
     &              *ihigh
           rz(ij) = (scba*(i-nmed1)+ssab*(j-nmed2)+cb*dens(ij)*scald)
         end do
      end do

      k = 0
      do i=1,ndim1
         do j =1,ndim2
           k = k + 1
           if (i.ne.ndim1.and.j.ne.ndim2) then
              m = k + ndim2
              ipl1(1) = iy(k)
              ipl1(2) = ix(k)
              ipl1(3) = iy(k+1)
              ipl1(4) = ix(k+1)
              ipl1(5) = iy(m)
              ipl1(6) = ix(m)
              ipl1(7) = iy(k)
              ipl1(8) = ix(k)
   
              ipl2(1) = iy(k+1)
              ipl2(2) = ix(k+1)
              ipl2(3) = iy(m)
              ipl2(4) = ix(m)
              ipl2(5) = iy(m+1)
              ipl2(6) = ix(m+1)
              ipl2(7) = iy(k+1)
              ipl2(8) = ix(k+1)
   
              if (icells.ge.256) then
               temp1(1) = (ix(k) - ix(k+1)) *1.0d0
               temp1(2) = (iy(k) - iy(k+1)) *1.0d0
               temp1(3) = ((rz(k) - rz(k+1))/scle1)*ihigh
               temp2(1) = (ix(m) - ix(k+1)) *1.0d0
               temp2(2) = (iy(m) - iy(k+1)) *1.0d0
               temp2(3) = ((rz(m) - rz(k+1))/scle1)*ihigh
               call crprod(temp1,temp2,temp3)
               call impsc(temp3,zvect,cot1)
               temp1(1) = (ix(k+1) - ix(m)) *1.0d0
               temp1(2) = (iy(k+1) - iy(m)) *1.0d0
               temp1(3) = ((rz(k+1) - rz(m))/scle1)*ihigh
               temp2(1) = (ix(m+1) - ix(m)) *1.0d0
               temp2(2) = (iy(m+1) - iy(m)) *1.0d0
               temp2(3) = ((rz(m+1) - rz(m))/scle1)*ihigh
               call crprod(temp1,temp2,temp3)
               call impsc(temp3,zvect,cot2)
               cot2 = -1.0d0*cot2
               ctt1 = 5.0d0*cot1 +0.5
               ctt2 = 5.0d0*cot2 +0.5
               icol1 = 130 + ctt1
               icol2 = 130 + ctt2
               if (icol1.lt.126) icol1 = 126
               if (icol2.lt.126) icol2 = 126
              else
               icol1 = 12
               icol2 = 12
              endif
   
              call drwpol(ipl1,npoly,icol1,1,0,1)
              call drwpol(ipl2,npoly,icol2,1,0,1)

           endif
         end do
      end do

      if (idisml.eq.0) return

      xx=3.0
      call xwin(xx,yy,10,str,nstr,idum1,idum2)
      call xwin(xx,yy,8,str,nstr,idum1,idum2)

      do i=1,natoms
         xa = -xsym(i)*ndim1 / r(1)
         ya = -ysym(i)*ndim2 / r(2)
         za = -zsym(i)*ndim1 / r(1)
         do j=i+1,natoms
              dmaxsq = (vdwr(nat(i))+vdwr(nat(j)))**2
              dijsq = ((xyz(1,i)-xyz(1,j))*adjus)**2
     +               +((xyz(2,i)-xyz(2,j))*adjus)**2
     +               +((xyz(3,i)-xyz(3,j))*adjus)**2
              if (dijsq.lt.dmaxsq) then
                  if (isym(i).eq.1.and.isym(j).eq.1) then
                     if (icells.eq.2) then
                        xx = 0.0
                     else
                        xx = 1.0
                     endif
                     call unstip
                  else
                     xx = 11.0
                     call ststip
                  endif
                  yy= 0.0
                  call xwin(xx,yy,99,str,nstr,idum1,idum2)
                  xb = -xsym(j)*ndim1 / r(1)
                  yb = -ysym(j)*ndim2 / r(2)
                  zb = -zsym(j)*ndim1 / r(1)
                  ixx(2) = ((ccab*xa/scle1+scab*ya/scle2-sb*za/scle1)
     &                     + 0.5d0)*ihigh
                  ixx(1) = (0.6d0-(-sa*xa/scle1+ca*ya/scle2))*ihigh
                  ixx(4) = ((ccab*xb/scle1+scab*yb/scle2-sb*zb/scle1)
     &                     + 0.5d0)*ihigh
                  ixx(3) = (0.6d0-(-sa*xb/scle1+ca*yb/scle2))*ihigh
                  call drawseg(ixx,1,0)
              endif 
         end do
      end do
      call unstip
      xx=1.0
      call xwin(xx,yy,10,str,nstr,idum1,idum2)
      call xwin(xx,yy,8,str,nstr,idum1,idum2)

      return
      end

      subroutine p3dd(iun,scale,ndimx,ndimz,adjus,dens)
      implicit double precision (a-h,o-z)

c THIS IS REALLY p3dv

      common /plrat/ rat(3)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      dimension dens(*)

      hinvx = dble(1.0d0/ndimx)
      hinvz = dble(1.0d0/ndimz)
      write(iun,'(''#VRML V2.0 utf8'')')
      write(iun,*) 'NavigationInfo { type ',char(34),
     &             'EXAMINE',char(34),' }'
      write(iun,*) 'Viewpoint { position -0.5 0 0.5 description ',
     &              char(34),'MoldenVRML',char(34),' }'
      write(iun,*) '#Background {'
      write(iun,*) '#skyColor [0 0.2 0.7,0 0.5 1,1 1 1 ]'
      write(iun,*) '#skyAngle [ 1.309, 1.571 ]'
      write(iun,*) '#groundColor [0.1 0.1 0,0.4 0.25 0.2,0.6 0.6 0.6]'
      write(iun,*) '#groundAngle [1.309,1.571]}'
      write(iun,*) 'Transform {'
      write(iun,*) 'translation -0.5 0 -0.5'
      write(iun,*) 'rotation 0 0 1 3.1415927'
      write(iun,*) 'children ['
      write(iun,*) 'Transform {'
      write(iun,*) '  scale ',hinvx,' ',hinvx,' ',hinvx
c      write(iun,*) '  translation -1.0 0.0 -1.0'
c      write(iun,*) '  translation -1.0 0.0 ',-1.0 + rat(1) / rat(2)
      write(iun,*) '  translation ',-0.5*r(2)/r(1),' 0.0 -0.5'
      write(iun,*) '  children ['
      write(iun,*) '    Shape {'
      write(iun,*) '       appearance Appearance {'
      write(iun,*) '         material Material {'
      write(iun,*) '           diffuseColor 1.0 0.0 1.0'
      write(iun,*) '         }'
      write(iun,*) '       }'
      write(iun,*) '       geometry ElevationGrid {'
      write(iun,*) '        xDimension ',ndimz
      write(iun,*) '        zDimension ',ndimx
      write(iun,*) '        xSpacing ',(dble(ndimx)/dble(ndimz))*
     &                                 (r(2)/r(1))
      write(iun,*) '        zSpacing 1'
      write(iun,*) '        solid FALSE'
      write(iun,*) '        height ['
      ij = 0
      do i=1,ndimz
          do j=1,ndimx
             ij = ij + 1
             write(iun,*) '                ',
     &               dens(ij)*scale*dble(ndimx)
          end do
      end do
      write(iun,*) '               ]'
      write(iun,*) '       }'
      write(iun,*) '    }'
      write(iun,*) '  ]'
      write(iun,*) '}'

      write(iun,*) 'Transform {'
      write(iun,*) ' rotation 1 0 0 1.5708'
      write(iun,*) ' children ['
      call plvmol(iun,ndimx,adjus)
      write(iun,*) '  ]'
      write(iun,*) '}'
      
      write(iun,*) '  ]'
      write(iun,*) '}'

      return
      end
