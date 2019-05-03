      subroutine snypnt(x1,y1,z1,x2,y2,z2,hght,sx,sy,sz)
      implicit double precision (a-h,o-z)
      dimension temp(3)

      temp(1) = x2 - x1
      temp(2) = y2 - y1
      temp(3) = z2 - z1

      rlambd = (hght - z1) / temp(3) 

      sx = x1 + rlambd*temp(1)
      sy = y1 + rlambd*temp(2)
      sz = z1 + rlambd*temp(3)

      return
      end
