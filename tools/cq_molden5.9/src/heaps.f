      subroutine heaps(nwords,arr,iarr)
      implicit double precision (a-h,o-z)
      dimension arr(*),iarr(*)

      do j=1,nwords
          iarr(j)=j
      end do
      l=nwords/2+1
      ir=nwords
10      continue
            if(l.gt.1) then
               l=l-1
               iarrt=iarr(l)
               arrt=arr(iarrt)
            else
               iarrt=iarr(ir)
               arrt=arr(iarrt)
               iarr(ir)=iarr(1)
               ir=ir-1
               if(ir.eq.1)then
                   iarr(1)=iarrt
                   goto 100
               endif
            endif
            i=l
            j=l+l
20          if(j.le.ir)then
               if(j.lt.ir)then
                  if(arr(iarr(j)).lt.arr(iarr(j+1)))j=j+1
               endif
               if(arrt.lt.arr(iarr(j)))then
                  iarr(i)=iarr(j)
                  i=j
                  j=j+j
               else
                  j=ir+1
               endif
           goto 20
           endif
           iarr(i)=iarrt
        goto 10
100     continue

        return
        end
