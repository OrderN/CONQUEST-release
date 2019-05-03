      integer function genaos(count)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      integer shella, shelln, shellt, shellc, shladf, aos
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg

      logical count
 
      genaos = 0
      itel = 1
      do 100 i = 1, nshell
c S
          if (shellt(i).eq.0) then
              aos(i) = itel
              itel = itel + 1
              goto 100
          endif

c sp
          if (shellt(i).ne.1.or.shellc(i).eq.1) goto 330
              aos(i) = itel
              itel = itel + 4
              goto 100

c p
330       if (shellt(i).eq.1) then
               aos(i) = itel - 1
               itel = itel + 3
               goto 100
          endif

c spd
          if (shellt(i).ne.2.or.shellc(i).ne.0) goto 350
              aos(i) = itel
              if (ido5d.eq.1) then
                 itel = itel + 9
              else
                 itel = itel + 10
              endif
              goto 100
c d

350       if (shellt(i).eq.2) then
              aos(i) = itel - 4
              if (ido5d.eq.1) then
                 itel = itel + 5
              else
                 itel = itel + 6
              endif
              goto 100
          endif

c          if (shellt(i).ne.3.or.(shellc(i).ne.0.and.shellc(i).ne.3))
c     &        call inferr('error in genaos',1)
          if (shellt(i).ne.3.or.(shellc(i).eq.2)) goto 370
          aos(i) = itel - 10
          if (ido7f.eq.1) then
             itel = itel + 7
          else
             itel = itel + 10
          endif
          goto 100

370       if (shellt(i).ne.4.or.(shellc(i).eq.2))
     &        call inferr('error in genaos',1)
          aos(i) = itel - 20
          if (ido9g.eq.1) then
             itel = itel + 9
          else
             itel = itel + 15
          endif
100       continue

      aos(nshell+1) = itel
      nbasis = itel - 1
      if (count) then
          genaos = nbasis
      else
          if (nbasis.ne.norbs) 
     &       call inferr('genaos: nbasis.ne.norbs',1)
      endif

      return
      end
