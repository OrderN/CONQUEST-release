      subroutine maxmid(npts1,npts2,scale,dens)
      implicit double precision (a-h,o-z)

c THIS IS REALLY maxmin

      logical euclid, yes
      common /cntval/ cntval, euclid, yes, icnt, iflag
      dimension dens(*)

      fact1 = 1.d0 / (npts1-1)
      fact2 = 1.d0 / (npts2-1)

      do 20 i=2,npts1-1
          do 20 j=2,npts2-1
              l = j + (i-1)*npts2
              k = l - npts2
              m = l + npts2
              if (dens(l).ge.dens(l-1)) goto 20
              if (dens(l).ge.dens(l+1)) goto 20
              if (dens(l).ge.dens(k-1)) goto 20
              if (dens(l).ge.dens(k  )) goto 20
              if (dens(l).ge.dens(k+1)) goto 20
              if (dens(l).ge.dens(m-1)) goto 20
              if (dens(l).ge.dens(m  )) goto 20
              if (dens(l).ge.dens(m+1)) goto 20
c
c found a minimum
c
              x = 1.d0 - (i-1)*fact1
              y = 1.d0 - (j-1)*fact2
              if (euclid) then
                 z = 0.d0
              else
                 z = dens(l)*scale
              endif
              call eulstr(x,y,z,dens(l))
  20  continue

      do 40 i=2,npts1-1
          do 40 j=2,npts2-1
              l = j + (i-1)*npts2
              k = l - npts2
              m = l + npts2
              if(dens(l).le.dens(l-1)) goto 40
              if(dens(l).le.dens(l+1)) goto 40
              if(dens(l).le.dens(k-1)) goto 40
              if(dens(l).le.dens(k  )) goto 40
              if(dens(l).le.dens(k+1)) goto 40
              if(dens(l).le.dens(m-1)) goto 40
              if(dens(l).le.dens(m  )) goto 40
              if(dens(l).le.dens(m+1)) goto 40
c
c found a maximum
c
              x = 1.d0 - (i-1)*fact1
              y = 1.d0 - (j-1)*fact2
              if (euclid) then
                 z = 0.d0
              else
                 z = dens(l)*scale
              endif
              call eulstr(x,y,z,dens(l))
  40  continue

      return
      end
