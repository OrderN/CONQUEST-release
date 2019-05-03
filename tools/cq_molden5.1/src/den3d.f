      subroutine dendd(npts1,npts2,scale,dens,edx,edy,iedlog)
      implicit double precision (a-h,o-z)
      real xx,yy
      character str*100

c THIS IS REALLY den3d

      parameter (isize=2000)
      character keywrd*320, keyori*320
      common /hide/   upx(isize),upy(isize),downx(isize),downy(isize),
     &                nup,ndown
      common /keywrd/ keywrd,keyori
      common /plot/   iplot,iplwin,icolps
      dimension dens(*),edx(*),edy(*),iedlog(*)
c
c   make a 3d picture ------
c
      factx = 1.d0 / (npts1-1)
      facty = 1.d0 / (npts2-1)

      if (iplot.eq.6) then
        xx = 15.0
        yy = 0.0
        call xwin(xx,yy,99,str,nstr,idum1,idum2)
      endif

      if (index(keywrd,'TRAN') .ne. 0 ) goto 100

c     make a 3d picture with hidden lines

      do i=0,npts1-1
         ii = i*npts2
         ij = ii + 1
         xvect = 1.d0 - i*factx
         call eulerh(xvect,1.d0,(dens(ij))*scale,edx(ij),edy(ij))
         do j = 2,npts2
           yvect = 1.0d0 - (j-1)*facty
           ij = i*npts2 + j
           call eulerh(xvect,yvect,(dens(ij))*scale,edx(ij),edy(ij))
         end do
      end do

      do i=1,npts2
         ij = i
         xvect = 1.d0 - (i-1)*facty
         call eulerh(1.d0,xvect,(dens(ij))*scale,edx(ij),edy(ij))
         iedlog(ij) = 0
         do j=1,npts1-1
             yvect = j*factx
             ij = j*npts2 + i
             call eulerh(1.d0-yvect,xvect,(dens(ij))*scale
     +       ,edx(ij),edy(ij))
             iedlog(ij) = 0
         end do
      end do

      nup = 1
      nsq = npts1*npts2
      fxt = edx(npts2)
      fyt = edy(npts2)
      upx(1) = fxt
      downx(1) = fxt
      upy(1) = fyt
      downy(1) = fyt
      iedlog(npts2) = 1
      call plotgr(1,fyt,fxt)

      do j=2,npts1
         ij = j*npts2
         nup = nup + 1
         fxt = edx(ij)
         fyt = edy(ij)
         iedlog(ij) = 1
         upx(nup) = fxt
         downx(nup) = fxt
         upy(nup) = fyt
         downy(nup) = fyt
         call plotgr(2,fyt,fxt)
      end do

      fxt = edx(nsq)
      fyt = edy(nsq)
      call plotgr(1,fyt,fxt)
      do i=2,npts2
         ij = nsq - i + 1
         nup = nup + 1
         fxt = edx(ij)
         fyt = edy(ij)
         iedlog(ij) = 1
         upx(nup) = fxt
         downx(nup) = fxt
         upy(nup) = fyt
         downy(nup) = fyt
         call plotgr(2,fyt,fxt)
      end do

      ndown = nup
      do k=2,npts1
         do l=2,npts2
            ij = (npts1-k)*npts2 + npts2 - l
            x1 = edx(ij+2)
            y1 = edy(ij+2)
            x2 = edx(ij+1)
            y2 = edy(ij+1)
            ij = ij + npts2 + 1
            x3 = edx(ij)
            y3 = edy(ij)
            call hidedr(x1,y1,x2,y2,x3,y3,iedlog(ij))
         end do
      end do

      return

100   continue

c     make a transparent 3d picture 

      do i=0,npts1-1
         ii = i*npts2
         ij = ii + 1
         xvect = 1.d0 - i*factx
         call euler(xvect,1.d0,(dens(ij))*scale,1)
         do j=2,npts2
           yvect = 1.0d0 - (j-1)*facty
           ij = i*npts2 + j
           call euler(xvect,yvect,(dens(ij))*scale,2)
         end do
      end do

      do i=1,npts2
         ij=i
         xvect=1.d0-(i-1)*facty
         call euler(1.d0,xvect,(dens(ij))*scale,1)
         do j=1,npts1-1
            yvect = j*factx
            ij = j*npts2 + i
            call euler(1.d0-yvect,xvect,(dens(ij))*scale,2)
         end do
      end do

      return
      end
