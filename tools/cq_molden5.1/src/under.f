      subroutine under(x1,y1,x2,y2,xq,yq,oans)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      oans=.true.
      if (x1.eq.xq.and.y1.eq.yq)then
         oans=.false.
         return
      endif
      if (x2.eq.xq.and.y2.eq.yq)then
         oans=.false.
         return
      endif
      dum=x2-x1
      if ( dum.ne.0.0d0 ) then
         dum=(y2-y1)*(xq-x1)/dum+(y1-yq)
         if ( dum.le.0.0d0 ) oans=.false.
      else
         if ( yq.ge.y1.or.yq.ge.y2 ) oans=.false.
      endif
      return
      end
