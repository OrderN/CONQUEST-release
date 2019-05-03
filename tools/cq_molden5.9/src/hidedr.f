      subroutine hidedr(x1,y1,x2,y2,x3,y3,iovis)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (isize=2000)
      common /hide/ upx(isize),upy(isize),downx(isize),downy(isize)
     &              ,nup,ndown
      dimension jdel(100),xadd(100),yadd(100)
      
      if (nup.gt.isize.or.ndown.gt.isize) return
      ipen = 1
      iadd = 0
      istart1 = 1
      iend1 = 0
      idel = 0
      istart2 = 1
      iend2 = 1
      do i=1,nup-1
         if (x1.ge.upx(i).and.x1.lt.upx(i+1)) istart1 = i
         if (x2.gt.upx(i).and.x2.le.upx(i+1)) iend1 = i
         if (x2.ge.upx(i).and.x2.lt.upx(i+1)) istart2 = i
         if (x3.gt.upx(i).and.x3.le.upx(i+1)) iend2 = i
         if (iend2.ne.1) goto 15
      end do
15    continue
      do i=istart1,iend1
         call above(upx(i),upy(i),upx(i+1),upy(i+1),x1,y1,oa)
         call above(upx(i),upy(i),upx(i+1),upy(i+1),x2,y2,ob)
         if (i.eq.istart1.and.oa) then
           call plotgh(ipen,x1,y1)
           iadd = iadd + 1
           xadd(iadd) = x1
           yadd(iadd) = y1
         endif
         if ((oa.and..not.ob).or.(ob.and..not.oa)) then
           call cross(upx(i),upy(i),upx(i+1),upy(i+1),x1,y1,x2,y2,os,
     &                xs,ys)
           if (os) then
             call plotgh(ipen,xs,ys)
             iadd = iadd + 1
             xadd(iadd) = xs
             yadd(iadd) = ys
           endif
         endif
         if (i.eq.iend1.and.ob) call plotgh(ipen,x2,y2)
         call above(x1,y1,x2,y2,upx(i),upy(i),ou)
         if (upx(i).ge.x1.and.(.not.ou).and.upy(i).ne.y1) then
           idel = idel + 1
           jdel(idel) = i
         endif
      end do

      do i=istart2,iend2
         call above(upx(i),upy(i),upx(i+1),upy(i+1),x2,y2,oa)
         call above(upx(i),upy(i),upx(i+1),upy(i+1),x3,y3,ob)
         if (i.eq.istart2.and.oa) then
           call plotgh(ipen,x2,y2)
           iadd = iadd + 1
           xadd(iadd) = x2
           yadd(iadd) = y2
         endif
         if ((oa.and..not.ob).or.(ob.and..not.oa)) then
           call cross(upx(i),upy(i),upx(i+1),upy(i+1),x2,y2,x3,y3,os,
     &                xs,ys)
           if (os) then
             call plotgh(ipen,xs,ys)
             iadd = iadd + 1
             xadd(iadd) = xs
             yadd(iadd) = ys
             if (x3.eq.xs.and.y3.eq.ys) iovis = 1
           endif
         endif
         if (i.eq.iend2.and.ob) then
           call plotgh(ipen,x3,y3)
           iadd = iadd + 1
           xadd(iadd) = x3
           yadd(iadd) = y3
           iovis = 1
         endif
         call above(x2,y2,x3,y3,upx(i),upy(i),ou)
         if (upx(i).ge.x2.and.(.not.ou).and.upy(i).ne.y2) then
           idel = idel + 1
           jdel(idel) = i
         endif
         if (i.eq.iend2.and.upx(i+1).eq.x3.and.oa) then
           call above(x2,y2,x3,y3,upx(i+1),upy(i+1),ou)
           if (.not.ou) then
             idel = idel + 1
             jdel(idel) = i + 1
           endif
         endif
      end do

      ieff = iadd - idel
      nup = nup + ieff
      iend = iend2 + 1

      do j=1,idel
         do inge = jdel(j)-j+1,iend-j
            upx(inge) = upx(inge+1)
            upy(inge) = upy(inge+1)
         end do
      end do

      if (ieff.eq.0) goto 70
         if (ieff.gt.0) then
           do inge=nup,iend+ieff+1,-1
             upx(inge) = upx(inge-ieff)
             upy(inge) = upy(inge-ieff)
           end do
         else
           do inge=iend+ieff+1,nup             
             upx(inge) = upx(inge-ieff)
             upy(inge) = upy(inge-ieff)
           end do
         endif
70       continue
      itemp = iend - idel
      if (istart1.eq.1) istart1 = 2
      do 80 j=1,iadd
         do i=istart1-1,itemp+j-1
            if (i.gt.0.and.j.gt.0) then
               if (xadd(j).ge.upx(i).and.xadd(j).lt.upx(i+1)) then
                 do k=itemp+j,i+2,-1
                    upx(k) = upx(k-1)
                    upy(k) = upy(k-1)
                 end do
                 upx(i+1) = xadd(j)
                 upy(i+1) = yadd(j)
                 goto 80
               endif
            endif
         end do
80    continue
      ipen = 1
      iadd = 0
      idel = 0
      iend2 = 0

      do i=1,ndown-1
         if (x1.ge.downx(i).and.x1.lt.downx(i+1)) istart1=i
         if (x2.gt.downx(i).and.x2.le.downx(i+1)) iend1=i
         if (x2.ge.downx(i).and.x2.lt.downx(i+1)) istart2=i
         if (x3.gt.downx(i).and.x3.le.downx(i+1)) iend2=i
         if (iend2.ne.0) goto 115
      end do

115   continue

      if (istart1.eq.1) istart1 = 2

      do i=istart1,iend1
         call under(downx(i),downy(i),downx(i+1),downy(i+1),x1,y1,oa)
         call under(downx(i),downy(i),downx(i+1),downy(i+1),x2,y2,ob)
         if (i.eq.istart1.and.oa) then
           call plotgh(ipen,x1,y1)
           iadd = iadd + 1
           xadd(iadd) = x1
           yadd(iadd) = y1
         endif
         if ((oa.and..not.ob).or.(ob.and..not.oa)) then
           call cross(downx(i),downy(i),downx(i+1),downy(i+1),x1,y1,x2
     &                ,y2,os,xs,ys)
           if (os) then
             call plotgh(ipen,xs,ys)
             iadd = iadd + 1
             xadd(iadd) = xs
             yadd(iadd) = ys
           endif
         endif
         if (i.eq.iend1.and.ob) call plotgh(ipen,x2,y2)
         call under(x1,y1,x2,y2,downx(i),downy(i),ou)
         if (downx(i).ge.x1.and.(.not.ou).and.downy(i).ne.y1) then
           idel = idel + 1
           jdel(idel) = i
         endif
      end do

      do i=istart2,iend2
         call under(downx(i),downy(i),downx(i+1),downy(i+1),x2,y2,oa)
         call under(downx(i),downy(i),downx(i+1),downy(i+1),x3,y3,ob)
         if (i.eq.istart2.and.oa) then
           call plotgh(ipen,x2,y2)
           iadd = iadd + 1
           xadd(iadd) = x2
           yadd(iadd) = y2
         endif
         if ((oa.and..not.ob).or.(ob.and..not.oa)) then
           call cross(downx(i),downy(i),downx(i+1),downy(i+1),x2,y2,x3
     &                ,y3,os,xs,ys)
           if (os) then
             call plotgh(ipen,xs,ys)
             iadd = iadd + 1
             xadd(iadd) = xs
             yadd(iadd) = ys
             if (x3.eq.xs.and.y3.eq.ys) iovis = 1
           endif
         endif
         if (i.eq.iend2.and.ob) then
           call plotgh(ipen,x3,y3)
           iadd = iadd + 1
           xadd(iadd) = x3
           yadd(iadd) = y3
           iovis = 1
         endif
         call under(x2,y2,x3,y3,downx(i),downy(i),ou)
         if (downx(i).ge.x2.and.(.not.ou).and.downy(i).ne.y2) then
           idel = idel + 1
           jdel(idel) = i
         endif
         if (i.eq.iend2.and.downx(i+1).eq.x3.and.oa) then
           call under(x2,y2,x3,y3,downx(i+1),downy(i+1),ou)
           if (.not.ou) then
             idel = idel + 1
             jdel(idel) = i + 1
           endif
         endif
      end do

      ieff = iadd - idel
      ndown = ndown + ieff
      iend = iend2 + 1
      if (ndown.gt.isize) ndown = isize

      do j=1,idel
         do inge=jdel(j)-j+1,iend-j
            downx(inge) = downx(inge+1)
            downy(inge) = downy(inge+1)
         end do
      end do

      if (ieff.eq.0) goto 170
         if (ieff.gt.0) then
           do inge=ndown,iend+ieff+1,-1
             downx(inge) = downx(inge-ieff)
             downy(inge) = downy(inge-ieff)
           end do
         else
           do inge=iend+ieff+1,ndown             
             downx(inge) = downx(inge-ieff)
             downy(inge) = downy(inge-ieff)
           end do
         endif
170   continue

      itemp = iend - idel
      if (itemp+iadd.gt.isize) goto 180

      do j=1,iadd
         do i=istart1-1,itemp+j-1
            if (xadd(j).ge.downx(i).and.xadd(j).lt.downx(i+1)) then
              do k=itemp+j,i+2,-1
                 downx(k) = downx(k-1)
                 downy(k) = downy(k-1)
              end do
              downx(i+1) = xadd(j)
              downy(i+1) = yadd(j)
              goto 180
            endif
          end do
      end do

180   continue
      return
      end
