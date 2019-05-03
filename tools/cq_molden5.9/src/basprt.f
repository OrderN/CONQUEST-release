      subroutine basprt(iun,gcompl,docrys)
      implicit double precision (a-h,o-z),integer (i-n)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      character*2 elm
      common /elem/elm(mxel)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      logical gcompl,docrys,dowrt
      character*3  namesh(5,5)
      character*3  name
      dimension iatflg(mxel),ich(7)
      data namesh/'s  ','sp ','spd','f  ','g  ',
     &            '???','p  ','d  ','???','g  ',
     &            's  ','sp ','d  ','f  ','g  ',
     &            's  ','sp ','d  ','f  ','g  ',
     &            's  ','sp ','d  ','f  ','g  '/
      data pi32 /5.56832799683170d0/
      data pt187 /1.875d+00/
      data pt5,pt75,pt656 /0.5d0,0.75d0,6.5625d0/
      data ich /2,8,6,10,18,14,18/
 1000 format(' ****')
 1010 format(a2,' atom ',i3)
 1020 format(1x,a3,1x,i2,' 1.00')
 1030 format(4d18.10)
 1040 format('################## Basis set #################')
 1050 format(' ',/)
c
      if (docrys) then
         do i=1,mxel
            iatflg(i) = 0
         end do
         natsh = 1
         do i=2,nshell
            if (jan(i).ne.jan(i-1)) then
               nt = nat(jan(i-1))
               if (iatflg(nt).eq.0) iatflg(nt) = natsh
               natsh = 0
            endif
            natsh = natsh + 1
         end do
         nt = nat(jan(nshell))
         if (iatflg(nt).eq.0) iatflg(nt) = natsh
      endif

      if (.not.gcompl.and..not.docrys) write(iun,1040)
      do i = 1, nshell
          if ((i.gt.1.and.jan(i).ne.jan(i-1)).and..not.docrys) then
              if (gcompl) then
                 write(iun,'(a)') '  '
              else
                 write(iun,1000)
              endif
          endif
          if (i.eq.1.or.jan(i).ne.jan(i-1)) then
              icht = nat(jan(i))
              if (docrys) then
                 dowrt = .false.
                 if (iatflg(nat(jan(i))).ne.0) then
                    write(iun,'(i3,1x,i3)') nat(jan(i)),
     &                        iatflg(nat(jan(i)))
                    iatflg(nat(jan(i))) = 0
                    dowrt = .true.
                 endif
              else
                 if (gcompl) then
                    write(iun,'(i3,a)') jan(i),' 0'
                 else
                    write(iun,1010)elm(nat(jan(i))),jan(i)
                 endif
                 dowrt = .true.
              endif
          endif
          if (.not.dowrt) goto 100

          iscod = 0
          if(shellt(i).eq.0) iscod = 0
          if(shellt(i).eq.1.and.shellc(i).eq.1) iscod = 2
          if(shellt(i).eq.1.and.shellc(i).ne.1) iscod = 1
          if(shellt(i).eq.2.and.shellc(i).eq.0) iscod = 4
          if(shellt(i).eq.2.and.shellc(i).ne.0) iscod = 3
          if(shellt(i).eq.3) iscod = 5
          if(shellt(i).eq.4) iscod = 6

          if (docrys.and.(iscod.ge.4)) then
              call inferr('CRYSTAL: spd/f/g not supported !',0)
              return
          endif

          ichpr = ich(iscod+1)
          if (icht.lt.ichpr) then
             ichpr = icht
             icht = 0
          else
             icht = icht - ichpr
          endif

          if (docrys) then
             write(iun,'(a,i1,a,i2,a,f3.1,a)') 
     &           '0 ',iscod,' ',shelln(i),' ',dfloat(ichpr),' 1.0'
          else
             name = namesh(shellt(i)+1,shellc(i)+1)
             write(iun,1020) name, shelln(i)
          endif

          mm = shella(i)
          mmdf = shladf(i)
          do j = 1, shelln(i)
              ee = exx(mm)+exx(mm)
              facs = pi32/(ee*dsqrt(ee))
              facp = dsqrt(pt5*facs/ee)
              facd = dsqrt(pt75*facs/(ee*ee))
              facf = dsqrt(pt187*facs/(ee**3))
              facg = dsqrt(pt656*facs/(ee**4))
              facs = dsqrt(facs)
              if (iscod.eq.0) write(iun,1030)exx(mm),c1(mm)*facs
              if (iscod.eq.2)
     &            write(iun,1030) exx(mm), c2(mm)*facp
              if (iscod.eq.1)
     &            write(iun,1030) exx(mm), c1(mm)*facs, c2(mm)*facp
              if (iscod.eq.4)
     &            write(iun,1030) exx(mm), c1(mm)*facs, c2(mm)*facp
     &                          , c3(mmdf)*facd
              if (iscod.eq.3)
     &            write(iun,1030) exx(mm), c3(mmdf)*facd
              if (iscod.eq.5) write(iun,1030)exx(mm),c4(mmdf)*facf
              if (iscod.eq.6) write(iun,1030)exx(mm),c5(mmdf)*facg
              mm = mm + 1
              mmdf = mmdf + 1
          end do
100       continue
      end do

      if (docrys) then
         write(iun,'(a)') '99 0'
      else
         write(iun,1050)
      endif

      return
      end
