      subroutine cross(x1,y1,x2,y2,x3,y3,x4,y4,oans,xs,ys)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      dimension d(6)
      toler=1.0d-10
      oans=.false.
      if(x1.eq.x3.and.y1.eq.y3)then
        oans=.true.
        xs=x1
        ys=y1
        return
      endif
      if(x1.eq.x4.and.y1.eq.y4)then
        oans=.true.
        xs=x1
        ys=y1
        return
      endif
      if(x2.eq.x3.and.y2.eq.y3)then
        oans=.true.
        xs=x2
        ys=y2
        return
      endif
      if(x2.eq.x4.and.y2.eq.y4)then
        oans=.true.
        xs=x2
        ys=y2
        return
      endif
      d(1)=x2-x1
      d(2)=x3-x1
      d(3)=x4-x3
      d(4)=y2-y1
      d(5)=y3-y1
      d(6)=y4-y3
      rnoem=d(1)*d(6)-d(3)*d(4)
      if ( rnoem.ne.0.0d0 ) then
        rmu=(d(2)*d(4)-d(1)*d(5))/rnoem
        rlambd=(d(2)*d(6)-d(3)*d(5))/rnoem
        xs=x3+rmu*d(3)
        ys=y3+rmu*d(6)
        if ( (rmu.le.1.0d0+toler.and.rmu.ge.0.0d0).and.
     &       (rlambd.le.1.0d0+toler.and.rlambd.ge.0.0d0)
     &       ) oans=.true.
      endif
      return
      end
