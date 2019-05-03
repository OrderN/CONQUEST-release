      subroutine eucmol(vdwr,adjus)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      character tk4014*5
      character str*100, el*2
      character esc, etx, gs, us
      real xx, yy
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /atopla/ xsym(numatm),ysym(numatm),zsym(numatm),
     &                isym(numatm)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /coord / xyz(3,numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plot/   iplot,iplwin,icolps
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      dimension vdwr(*)

      etx    = char(3)
      esc    = char(27)
      gs     = char(29)
      us     = char(31)

c
c WRITE ELEMENT NAMES and connect atoms in plot plane
c

      do i=1,natoms

         xa = 0.5d0 + ysym(i) / r(2)
         ya = 0.5d0 + xsym(i) / r(1)
         za = zsym(i) / r(3)

         if(xa.lt.1.0.and.xa.gt.0.0.and.ya.lt.1.0.and.ya.gt.0.0)then

            call eulerh(ya,xa,za,xaa,yaa)
            xaa = 0.5d0 + xaa
            yaa = 0.5d0 + yaa

            if (iplot.eq.0) then
               write(iun4,'(''.m '',f5.3,f6.3)') xaa,yaa
               write(iun4,'(''.to 5'')')
               write(iun4,'(''.sy 4'')')
               write(iun4,'(''.to 3'')')
               write(iun4,'(''.pt '',a2)')  elemnt(nat(i)) 
            endif

            if (iplot.eq.1) then
               write(iun4,'(''PU'',f5.3,f6.3,'';'')')xaa,yaa
               write(iun4,'(''WG0.005,0,360;'')')
               write(iun4,'(''PU'',f5.3,f6.3,'';LB'',a3,a,'';'')')
     +            xaa,yaa,elemnt(nat(i)),etx
            endif

            if (iplot.eq.3) then
c               write(iun4,*)esc//'/0d'
               write(iun4,*)gs//tk4014(int(xaa*3100),int(yaa*3100))
     +                    //us//elemnt(nat(i))
            endif

            if (iplot.eq.4) then
              write(iun4,'(''   0 setgray'')')
              write(iun4,'(''n '',i4,'' '',i4,'' 10 0 360 arc fill'')')
     +        int(xaa*2000),int(yaa*2000+125)
              write(iun4,'(''   1 setgray'')')
              write(iun4,'(''n '',i4,'' '',i4,'' 17 0 360 arc fill'')')
     +        int(xaa*2000+40),int(yaa*2000+138)
              write(iun4,'(''   0 setgray'')')
              write(iun4,'(i4,'' '',i4,'' m (   '',a,'') show'')')
     +        int(xaa*2000),int(yaa*2000+125),elemnt(nat(i))
            endif

            if (iplot.eq.6) then
               el = elemnt(nat(i))
               xx = 5.0
               yy = 0.0
               call xwin(xx,yy,99,str,nstr,idum1,idum2)
               xx = xaa
               yy = yaa
               call xwin(xx,yy,4,el,2,idum1,idum2)
            endif

         endif

         if (iplot.eq.6) then
           xx=3.0
           call xwin(xx,yy,10,str,nstr,idum1,idum2)
           call xwin(xx,yy,8,str,nstr,idum1,idum2)
           call ststip
         endif

         do j=1,natoms

            if(isym(j).eq.1.and.isym(i).eq.1)then
               if (iplot.eq.6) then
                  xx = 15.0
                  yy = 0.0
                  call xwin(xx,yy,99,str,nstr,idum1,idum2)
                  call unstip
               endif
               if (iplot.eq.4) then
                  write(iun4,'(''0 setgray'')')
               endif
            else
               if (iplot.eq.6) then
                  xx = 8.0
                  yy = 0.0
                  call xwin(xx,yy,99,str,nstr,idum1,idum2)
               endif
               if (iplot.eq.4) then
                  write(iun4,'(''0.5 setgray'')')
               endif
            endif

            dmaxsq = (vdwr(nat(i)) + vdwr(nat(j)))**2
            dijsq = ((xyz(1,i) - xyz(1,j))*adjus)**2
     +            + ((xyz(2,i) - xyz(2,j))*adjus)**2
     +            + ((xyz(3,i) - xyz(3,j))*adjus)**2
            if (dijsq.lt.dmaxsq) then
               xb = 0.5d0 + ysym(j) / r(2)
               yb = 0.5d0 + xsym(j) / r(1)
               zb = zsym(j) / r(3)
               if (xa.lt.1.0.and.xa.gt.0.0.and.ya.lt.1.0
     +         .and.ya.gt.0.0.and.xb.lt.1.0.and.xb.gt.0.0
     +         .and.yb.lt.1.0.and.yb.gt.0.0) then
                  if (iplot.eq.3) then
                     idum = 1
                     call plotgh(idum,0.0d0,0.0d0)
                  endif
                  if (iplot.eq.4) then
                     write(iun4,'(''n'')')
                  endif
                  call euler(ya,xa,za,1)
                  call euler(yb,xb,zb,2)
                  if (iplot.eq.4) then
                     write(iun4,'(''s'')')
                  endif
               endif
            endif
         end do

         if (iplot.eq.6) then
           xx=1.0
           call xwin(xx,yy,10,str,nstr,idum1,idum2)
           call xwin(xx,yy,8,str,nstr,idum1,idum2)
           call unstip
         endif

         if (iplot.eq.4) then
           write(iun4,'(''s'')')
           write(iun4,'(''0 setgray'')')
           write(iun4,'(''n'')')
         endif

      end do

      return
      end
