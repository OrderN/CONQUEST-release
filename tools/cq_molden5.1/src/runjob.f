      subroutine runjob(imod,iqopt,ihaszm)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb
      character*80 glin1,glin2,gtitl,rungam
      character*15 jname,qname
      common /gauopt/ ito,imo,ibo,itotc,imult,ibatch,ihess,itime,
     &                iwxyz,ichh,ichm,ichl,imh,imm,iml,
     &                glin1,glin2,gtitl,jname,qname,rungam
      common /athlp/  iatoms, mxnat
      character*80 tnknm
      integer tnkbg,tnkit,tnkarc,tnkarf,tnkprg
      common /tnkopt/ rmsgrd,tnkbg,tnkit,tnkarc,tnkarf,tnkprg,tnknm,icst
      common /types/ iff
      integer dolabs,fancy,persp,shade,atcol,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      common /pbc/ abc(3),ibox,icell,igfmap
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*80 taroot
      character*80 pwd
      character*11 ename
      character*80 outfil
      character*137 line
      logical mapxyz,opfil,domap
      
      izmtmp = ihaszm

      jlen = linlen(jname)
      iqlen = linlen(qname)
      itnkln = linlen(tnknm)

      ename = 'PWD'
      call getenv(ename,pwd)

      ig94 = 1
      itnk = 0
      if (imod.eq.1) then
         ename = 'TA_ROOT'
      elseif (imod.eq.2) then
         ename = 'QNT_ROOT'
      elseif (imod.eq.3) then
         ename = 'TNK_ROOT'
      elseif (imod.eq.4) then
         ename = 'GAMESS_ROOT'
      elseif (imod.eq.5) then
         ename = 'g94root'
      elseif (imod.eq.6) then
         ename = 'MOPAC_ROOT'
      elseif (imod.eq.7) then
         ename = 'XTINKER'
      elseif (imod.eq.8) then
         ename = 'AMBFOR'
      endif

      call getenv(ename,taroot)
      if (taroot(1:1).eq.' ') then
          if (imod.eq.1) then
            call inferr('Environment variable TA_ROOT not set !',0)
            call zmterr('Environment variable TA_ROOT not set !',1,0,1)
          elseif (imod.eq.2) then
            call inferr('Environment variable QNT_ROOT not set !',0)
            call zmterr('Environment variable QNT_ROOT not set !',1,0,1)
          elseif (imod.eq.3) then
            call inferr('Environment variable TNK_ROOT not set !',0)
            call zmterr('Environment variable TNK_ROOT not set !',1,0,1)
          elseif (imod.eq.4) then
            call inferr('Environment variable GAMESS_ROOT not set !',0)
            call zmterr('Environment variable GAMESS_ROOT not set !',1,
     &                  0,1)
          elseif (imod.eq.5) then

c ig94   0  => g92
c ig94   1  => g94
c ig94   2  => g98
c ig94   3  => g03
c ig94   4  => g09

            ename = 'g92root'
            call getenv(ename,taroot)
            if (taroot(1:1).eq.' ') then
               ename = 'g98root'
               call getenv(ename,taroot)
            else
               ig94 = 0
            endif
            if (taroot(1:1).eq.' ') then
               ename = 'g03root'
               call getenv(ename,taroot)
            else
               ig94 = 2
            endif
            if (taroot(1:1).eq.' ') then
               ename = 'g09root'
               call getenv(ename,taroot)
            else
               ig94 = 3
               goto 10
            endif
            if (taroot(1:1).eq.' ') then
               call inferr(
     &     'Environment variable g92/g94/g98/g03/g09root not set !',0)
               call zmterr(
     &   'Environment variable g92/g94/g98/g03/g09root not set !',1,0,1)
               goto 10
            else
               ig94 = 4
               goto 10
            endif
          elseif (imod.eq.6) then
            call inferr('Environment variable MOPAC_ROOT not set !',0)
            call zmterr('Environment variable MOPAC_ROOT not set !',1,
     &                   0,1)
          elseif (imod.eq.7) then
            ename = 'TNK_ROOT'
            call getenv(ename,taroot)
            if (taroot(1:1).eq.' ') then
               call inferr(
     &         'Environment variable TNK_ROOT/XTINKER not set !',0)
            else
               itnk = 1
               goto 10
            endif
          elseif (imod.eq.8) then
            taroot = '/usr/local/bin/'
            goto 10
          endif
          return
      endif
10    itlen = linlen(taroot)

      iun = 50
      iuncln = 51

      if (imod.eq.1) then
         ixyz = 1
         call wrpnt('molin.mol2',10,4,0,0,0,0,1,0)
      elseif (imod.eq.2) then
         ixyz = 2
         call wrpnt('molin.msf',9,4,0,0,0,0,1,0)
      elseif (imod.eq.3) then
         ixyz = 3
         call wrpnt(tnknm(1:itnkln)//'.xyz',itnkln+4,4,0,0,0,0,1,0)
      elseif (imod.eq.4) then
         call wrpnt(jname(1:jlen)//'.in',1,19,0,0,0,0,1,0)
      elseif (imod.eq.5) then
         call wrpnt(jname(1:jlen)//'.com',1,20,0,0,0,0,1,0)
      elseif (imod.eq.6) then
         call wrpnt(jname(1:jlen)//'.dat',1,21,0,0,0,0,1,0)
      elseif (imod.eq.7) then
         if (iff.ne.1) then
             iff = 1
             call dotyp(1)
         endif
         ixyz = 5
         call wrpnt(tnknm(1:itnkln)//'.cssr',itnkln+5,4,0,0,0,0,1,0)
         if (iqopt.eq.0) then
            call fdat(16,0,0,0,0,0)
         else
            call fdat(17,0,0,0,0,0)
         endif
         call fdat(ifd,0,0,0,0,0)
      elseif (imod.eq.8) then
         nwramb = iatoms
         ixyz = 11
         call wrpnt(tnknm(1:itnkln)//'.xyz',itnkln+4,4,0,0,0,0,1,0)
      endif

      if (imod.ne.8) then
         if (imod.ge.3) then
            if (imod.eq.3.or.imod.eq.7) then
               open(unit=iun,form='formatted',
     &           file=tnknm(1:itnkln)//'.run',
     &           status='unknown',err=100)
               open(unit=iuncln,form='formatted',file='clnjob',
     &           status='unknown',err=100)
            else
               open(unit=iun,form='formatted',
     &           file=jname(1:jlen)//'.run',
     &           status='unknown',err=100)
            endif
         else
            open(unit=iun,form='formatted',file='runjob',
     &           status='unknown',err=100)
         endif
      endif

      if (imod.eq.1) then
      write(iun,'(''#!/bin/sh'')')
      write(iun,*) 'rm molout.mol2'
      write(iun,*) '. $TA_ROOT/lib/.profile'
      write(iun,*) '$TA_ROOT/bin/sybyl << EOF > molout.log'
      write(iun,*) 'tailor set mol file_format mol2 | |'
      write(iun,*) 'set autosave off'
      write(iun,*) 'mol in M1 molin.mol2'
      write(iun,*) 'SETVAR TAILOR!MAXIMIN2!MINIMIZATION_METHOD POWELL'
      write(iun,*) 'SETVAR TAILOR!MAXIMIN2!TERMINATION_OPTION GRADIENT'
      write(iun,*) 'SETVAR TAILOR!MAXIMIN2!MIN_ENERGY_CHANGE  0.050'
      write(iun,*) 'SETVAR TAILOR!MAXIMIN2!RMS_GRADIENT  0.050'
      write(iun,*) 'SETVAR TAILOR!MAXIMIN2!RMS_DISPLACEMENT  0.001'
      write(iun,*) 'SETVAR TAILOR!MAXIMIN2!MAXIMUM_ITERATIONS 1000'
      write(iun,*) 'SETVAR TAILOR!FORCE_FIELD!PARAMETER_SET Tripos'
      write(iun,*) 'SETVAR TAILOR!FORCE_FIELD!NON_BONDED_CUTOFF 8.0'
      write(iun,*) 'SETVAR TAILOR!FORCE_FIELD!DIELECTRIC_CONSTANT 1.0'
      write(iun,*) 
     &  'SETVAR TAILOR!FORCE_FIELD!DIELECTRIC_FUNCTION distance'
      write(iun,*) 
     &  'SETVAR TAILOR!FORCE_FIELD!REVIEW_HS_AND_LPS DO_NOT_CHANGE'
      write(iun,*) 'SETVAR TAILOR!MAXIMIN2!LIST_TERMS NO'
      write(iun,*) 'SETVAR TAILOR!FORCE_FIELD!ONE_FOUR_SCALING 1.0'
      write(iun,*) 'SETVAR TAILOR!FORCE_FIELD!HBOND_RAD_SCALING 0.7'
      write(iun,*) 'MAXIMIN2 M1 Electrostatics IGNORE_ELECTROSTATICS ',
     &             char(92)
      write(iun,*) 'Interesting M1(*-(0)) Boundary_Conditions ',char(92)
      write(iun,*) 'IGNORE_PBCS List DONE INTERACTIVE'
      write(iun,*) 'mol out M1 molout.mol2'
      write(iun,'(''EOF'')')
      elseif (imod.eq.2) then
      write(iun,'(''#!/bin/csh'')')
      write(iun,*) 'unalias mv'
      write(iun,*) 'unalias rm'
      write(iun,*) 'source /compchem/1/msi/quanta96/.setquanta'
      write(iun,*) 'mkdir qnt$$'
      write(iun,*) 'cd qnt$$'
      write(iun,*) 'cat <<EOF > quanta.scr'
      write(iun,*) 'mole init molin.msf'
      write(iun,*) 'charm init'
      write(iun,*) 'charm writ psf'
      write(iun,*) 'charm set mini nstep 1000'
      write(iun,*) 'charm do mini'
      write(iun,*) 'msf save over'
      write(iun,*) 'end'
      write(iun,'(''EOF'')')
      write(iun,*) 'cp $HYD_LIB/default.cst .cst'
      write(iun,*) 'mv ../molin.msf .'
      write(iun,*) 'chmod u+rwx .cst'
      write(iun,*) 'quanta -n -i quanta.scr -l quanta >& /dev/null'
      write(iun,*) 'mv molin.msf ..'
      write(iun,*) 'mv quanta.log ..'
      write(iun,*) 'cd ..'
      write(iun,*) 'unalias rm'
      write(iun,*) 'rm -rf qnt$$'
      elseif (imod.eq.3.or.imod.eq.7) then
      write(iuncln,'(''#!/bin/csh -f'')')
      write(iuncln,*) 'unalias rm'
      write(iuncln,*) 'if ( -e '//tnknm(1:itnkln)//
     &                ' ) mv '//tnknm(1:itnkln)//' '//
     &                tnknm(1:itnkln)//'.bck'

      if (imod.eq.3) then
         write(iuncln,*) 'rm '//tnknm(1:itnkln)//'.xyz_*'
      else
         write(iuncln,*) 'rm '//tnknm(1:itnkln)//'.cssr_*'
      endif
      write(iuncln,*) 'rm '//tnknm(1:itnkln)//'.[0-9]*'
      write(iuncln,*) 'rm '//tnknm(1:itnkln)//'.arc*'
      close(iuncln)
      write(iun,'(''#!/bin/csh -f'')')
      write(iun,*) 'unalias rm'
      write(iun,*) 'rm clnjob'
      write(iun,*) 'setenv AMBER amber'
      write(iun,*) 'setenv CHARMM charmm'
      write(iun,*) 
     &  'if ( -f $TNK_ROOT/params/amber99.prm ) setenv AMBER amber99'
      write(iun,*) 
     &  'if ( -f $TNK_ROOT/params/charmm27.prm ) setenv CHARMM charmm27'
      if (imod.eq.3) then
        write(iun,*) 'cat <<EOF > '//tnknm(1:itnkln)//'.key'
        if (iff.le.1) then
         write(iun,*) 'parameters $TNK_ROOT/params/mm3'
         write(iun,*) 'torsion       5   2   3   6       '//
     &                '0.000 +1     8.000 -2     0.000 +3'
         write(iun,*) 'torsion       2   3   6  24       '//
     &                '0.000 +1     0.500 -2     0.000 +3'
         write(iun,*) 'angle         1    3    77    '//
     &                '0.850    123.500    123.500      0.000'
         write(iun,*) 'angle         1    3    78    '//
     &                '0.850    123.500    123.500      0.000'
         write(iun,*) 'angle         5    1    75    '//
     &                '0.820    110.000    108.900    108.700'
         write(iun,*) 'opbend        3    77         0.650'
         write(iun,*) 'opbend        3    78         0.650'
        elseif (iff.eq.2) then
         write(iun,*) 'parameters $TNK_ROOT/params/$CHARMM'
        elseif (iff.eq.3) then
         write(iun,*) 'parameters $TNK_ROOT/params/$AMBER'
        elseif (iff.eq.4) then
         write(iun,*) 'parameters $TNK_ROOT/params/amoebapro.prm'
        endif
        if (.not.(tnkbg.eq.1.and.tnkarc.eq.0)) 
     &     write(iun,*) 'save-cycle'
        write(iun,*) 'maxiter ',tnkit
        if (tnkprg.eq.5) write(iun,*) 'saddlepoint'
        write(iun,'(''EOF'')')
      else
        write(iun,*) 'cat < tnk.key > '//tnknm(1:itnkln)//'.key'
        write(iun,*) 'rm tnk.key'
        write(iun,*) 'cat <<EOF >>'//tnknm(1:itnkln)//'.key'
        if (.not.(tnkbg.eq.1.and.tnkarc.eq.0)) 
     &     write(iun,*) 'save-cycle'
        write(iun,*) 'maxiter ',tnkit
        write(iun,*) 'writeout 1'
        if (tnkprg.eq.5) write(iun,*) 'saddlepoint'
        write(iun,'(''EOF'')')
      endif

      if (tnkprg.eq.0) then
         write(iun,*) '$TNK_ROOT/bin/minimize '//tnknm(1:itnkln)
     &                //' <<EOF >& '//tnknm(1:itnkln)//'.log'
      elseif (tnkprg.eq.1) then
         write(iun,*) '$TNK_ROOT/bin/optimize '//tnknm(1:itnkln)
     &                //' <<EOF >& '//tnknm(1:itnkln)//'.log'
      elseif (tnkprg.eq.2.or.tnkprg.eq.5) then
         write(iun,*) '$TNK_ROOT/bin/newton '//tnknm(1:itnkln)
     &                //' <<EOF >& '//tnknm(1:itnkln)//'.log'
         write(iun,*) ' '
         write(iun,*) ' '
      elseif (tnkprg.eq.3) then
         write(iun,*) '$TNK_ROOT/bin/dynamic '//tnknm(1:itnkln)
     &                //' <<EOF >& '//tnknm(1:itnkln)//'.log'
         write(iun,*) '10000'
         write(iun,*) ' '
         write(iun,*) ' '
         write(iun,*) ' '
         write(iun,'(''EOF'')')
      elseif (tnkprg.eq.4) then
         if (itnk.eq.1) then
            write(iun,*) '$TNK_ROOT/bin/xxtalmins c '//tnknm(1:itnkln)
     &                //' <<EOF >& '//tnknm(1:itnkln)//'.log'
         else
            write(iun,*) '$XTINKER/bin/xxtalmins c '//tnknm(1:itnkln)
     &                //' <<EOF >& '//tnknm(1:itnkln)//'.log'
         endif
      endif

      if (tnkprg.lt.3.or.tnkprg.eq.4) then
         write(iun,'(f12.8)') rmsgrd
         write(iun,'(''EOF'')')
      endif

      if (tnkbg.eq.1.and.tnkarc.eq.1) then
         if (imod.eq.3) then
            write(iun,*) '$TNK_ROOT/bin/archive '//tnknm(1:itnkln)
     &                   //' <<EOF >>& '//tnknm(1:itnkln)//'.log'
            write(iun,*) ' '
            write(iun,*) ' '
            write(iun,*) '1 1000 ',tnkarf
            write(iun,'(''EOF'')')
         else
            write(iun,*) 'cat '//tnknm(1:itnkln)//'.[0-9]* > '//
     &                    tnknm(1:itnkln)//'.arc'
         endif
         write(iun,*) 'rm '//tnknm(1:itnkln)//'.[0-9]*'
      endif

      elseif (imod.eq.4) then
         write(iun,'(''#!/bin/csh'')')
         write(iun,'(a)') rungam(1:linlen(rungam))
      elseif (imod.eq.5) then
         write(iun,'(''#!/bin/csh'')')
         if (ig94.eq.1) then
            write(iun,*) 'source $g94root/g94/bsd/g94.login'
            if (ibatch.eq.1) then
               write(iun,*) 'g94 '//jname(1:jlen)
            else
               write(iun,*) 'subg94 '//qname(1:iqlen)//' '
     &                      //jname(1:jlen)
            endif
         elseif (ig94.eq.2) then
            write(iun,*) 'source $g98root/g98/bsd/g98.login'
            if (ibatch.eq.1) then
               write(iun,*) 'g98 '//jname(1:jlen)
            else
               write(iun,*) 'subg98 '//qname(1:iqlen)//' '
     &                      //jname(1:jlen)
            endif
         elseif (ig94.eq.3) then
            write(iun,*) 'source $g03root/g03/bsd/g03.login'
            if (ibatch.eq.1) then
               write(iun,*) 'g03 '//jname(1:jlen)
            else
               write(iun,*) 'subg03 '//qname(1:iqlen)//' '
     &                      //jname(1:jlen)
            endif
         elseif (ig94.eq.4) then
            write(iun,*) 'source $g09root/g09/bsd/g09.login'
            if (ibatch.eq.1) then
               write(iun,*) 'g09 '//jname(1:jlen)
            else
               write(iun,*) 'subg09 '//qname(1:iqlen)//' '
     &                      //jname(1:jlen)
            endif
         else
            write(iun,*) 'source $g92root/g92/bsd/g92.login'
            if (ibatch.eq.1) then
               write(iun,*) 'g92 '//jname(1:jlen)
            else
               write(iun,*) 'subg92 '//qname(1:iqlen)//' '//
     &                       jname(1:jlen)
            endif
         endif
         write(iun,*) ' '
         write(iun,*) ' '
         write(iun,*) ' '
      elseif (imod.eq.6) then
         write(iun,'(''#!/bin/csh'')')
         if (ibatch.eq.0) write(iun,*) 'cd '//pwd(1:linlen(pwd))
         write(iun,*) taroot(1:itlen)//'/bin/mopac.csh '//
     &                jname(1:jlen)//
     &                ' >& /dev/null'
      endif

      if (imod.ne.8) close(iun)

c make jobscript executable

      if (imod.ge.3) then
         if (imod.eq.3.or.imod.eq.7) then
            taroot = 'chmod u+x '//tnknm(1:itnkln)//'.run'
            nstr = itnkln+14
         else
            taroot = 'chmod u+x '//jname(1:jlen)//'.run'
            nstr = jlen+14
         endif
      else
         taroot = 'chmod u+x runjob'
         nstr = 16
      endif

      iextmp = 0
      if (imod.ne.8) call exstr(taroot,nstr,iextmp)

c prepare cleaning up

      if (imod.eq.3.or.imod.eq.7) then
         taroot = 'chmod u+x clnjob'
         nstr = 16
         iextmp = 0
         call exstr(taroot,nstr,iextmp)
         taroot = '.'//char(47)//'clnjob'
         nstr = 8
         iextmp = 0
         call exstr(taroot,nstr,iextmp)
      endif

c are you ready ?

      if (imod.ne.3.and.imod.ne.7.and.imod.ne.8) call curs(1)
      call confrm(0,istat)
      if (istat.ne.1) then
          call curs(0)
          return
      endif

      if (imod.ge.3.and.imod.ne.8) then
         if (imod.eq.6.and.ibatch.eq.0) then
            taroot = 'qsub -eo -q '//qname(1:iqlen)//' '//
     &               jname(1:jlen)//'.run'
            nstr = jlen+iqlen+17
         else
            if (imod.eq.3.or.imod.eq.7) then
               taroot = '.'//char(47)//tnknm(1:itnkln)//'.run'
               nstr = itnkln+6
               if (tnkbg.eq.0) then
                   taroot = tnknm(1:itnkln)
                   nstr = itnkln
               endif
            else
               taroot = '.'//char(47)//jname(1:jlen)//'.run'
               nstr = jlen+6
            endif
         endif
      elseif (imod.eq.8) then
         taroot = 
     &      'ambfor '//tnknm(1:itnkln)//'.xyz'
         nstr = 7 + itnkln + 4
      else
        
         taroot = '.'//char(47)//'runjob'
         nstr = 8
      endif

c execute or batching of jobscript

      iextmp = 0
      if (imod.lt.3.or.(imod.ge.4.and.ibatch.eq.0).or.
     &    ((imod.eq.3.or.imod.eq.7).and.tnkbg.eq.0).or.
     &    (imod.eq.8.and.tnkbg.eq.0)) then

         if ((imod.eq.3.or.imod.eq.7.or.imod.eq.8).and.
     &        tnkbg.eq.0) then

c    tinker in the background updating screen; also fork

            if (imod.eq.8) then
               taroot = tnknm(1:itnkln)
               nstr = itnkln
               iextmp = 3
               call exstr(taroot,nstr,iextmp)
            else
               iextmp = 2
               call exstr(taroot,nstr,iextmp)
            endif
         else

c    system

            iextmp = 0
            call exstr(taroot,nstr,iextmp)
         endif

      else

c    fork

         if (imod.eq.8) then
            taroot = tnknm(1:itnkln)
            nstr = itnkln
            iextmp = 4
         else
            iextmp = 1
         endif
         call exstr(taroot,nstr,iextmp)

      endif

      iexst = iextmp

c clean up tinker cycle files of a non detached job

      if ((imod.eq.3.or.imod.eq.7).and.tnkbg.eq.0) then
         open(unit=iuncln,form='formatted',file='clnjob',
     &     status='unknown',err=100)
         write(iuncln,'(''#!/bin/csh -f'')')
         if (tnkarc.eq.1) then
          if (imod.eq.3) then
            write(iuncln,*) '$TNK_ROOT/bin/archive '//tnknm(1:itnkln)
     &                     //' <<EOF >>& '//tnknm(1:itnkln)//'.log'
            write(iuncln,*) ' '
            write(iuncln,*) ' '
            write(iuncln,*) '1 1000 ',tnkarf
            write(iuncln,'(''EOF'')')
          else
            write(iuncln,*) 'cat '//tnknm(1:itnkln)//'.[0-9]* > '//
     &                    tnknm(1:itnkln)//'.arc'
          endif
         endif
         write(iuncln,*) 'unalias rm'
         write(iuncln,*) 'rm '//tnknm(1:itnkln)//'.[0-9]*'
         write(iuncln,*) 'rm clnjob'
         close(iuncln)
         iextmp = 0
         taroot = 'chmod u+x clnjob'
         nstr = 16
         call exstr(taroot,nstr,iextmp)
         taroot = '.'//char(47)//'clnjob'
         nstr = 8
         call exstr(taroot,nstr,iextmp)
      endif

c open to be mapped file

      ierr = 0
      call tomap(imap,iambfr)
      if (imap.eq.0) ierr = 21 

      iform = 1
      domap = .true.
      if (imod.eq.1) then
         nstr = 11
         outfil = 'molout.mol2'
      elseif (imod.eq.2) then
         iform = 0
         nstr = 9
         outfil = 'molin.msf'
      elseif (imod.eq.3.and.tnkbg.eq.0.and.iexst.ne.-1) then
         nstr = itnkln+6
         outfil = tnknm(1:itnkln)//'.xyz_2'
      elseif (imod.eq.7.and.tnkbg.eq.0.and.iexst.ne.-1) then
         nstr = itnkln+7
         outfil = tnknm(1:itnkln)//'.cssr_2'
c      elseif (imod.eq.8.and.tnkbg.eq.0.and.iexst.ne.-1) then
      elseif (imod.eq.8.and.tnkbg.eq.0) then
         nstr = itnkln+8
         outfil = tnknm(1:itnkln)//'_opt.xyz'
         if (igfmap.ne.1.or.icell.eq.1) izmtmp = 0

c check for ambfor md run, does not need to map _opt file:

         call tomap(imap,iambfr)
         if (imap.eq.0) domap = .false.

      else
         domap = .false.
      endif

      if (domap) then
         if (opfil(iun,outfil,nstr,iform,1,0)) then
             rewind(iun)
         else
             ierr = 1
             goto 50
         endif
      else
         goto 50
      endif

c map the file onto the z-matrix, or just read it in

      if (imod.ne.4.and..not.
     &    ((imod.eq.3.or.imod.eq.7.or.imod.eq.8)
     &      .and.tnkbg.eq.1)) then
         ierr = 1
         call tomap(imap,iambfr)
         if (imap.eq.0) ierr = 21 

         if (iexst.ne.-1) iexst = imod
         if (mapxyz(iun,iexst,iff,izmtmp)) ierr = 0
         if (iexst.eq.-1) then
            ierr = 20
            call dohcon(0)
         else
            if (imod.eq.7) call fdat(ifd,0,0,0,0,0)
         endif
      endif
      if (domap) close(iun)

c check for tinker error conditions (non detached job)

50    if ((imod.eq.3.or.imod.eq.7)
     &   .and.tnkbg.eq.0.and.iexst.ne.-1) then
         if (opfil(iun,tnknm(1:itnkln)//'.log',itnkln+4,1,1,0)) then
             iuntmp = iun2
             iun2 = iun
             call searcht(line,'MECHANIC','Incomplete Convergence',
     &                'Too many Parameters',istat)
             if (istat.ne.0) then
                if (index(line,'MECHANIC').ne.0) then
                   ierr = 2
                elseif (index(line,'Incomplete Convergence').ne.0) then
                   ierr = 3
                elseif (index(line,'Too many Parameters').ne.0) then
                   ierr = 4
                endif
             endif
             close(iun)
             iun2 = iuntmp
         endif
      endif
      if ((imod.eq.3.or.imod.eq.7).and.iexst.eq.-1) ierr = 5

      if ((imod.ge.4.and.imod.lt.7).or.
     &   ((imod.eq.3.or.imod.eq.7.or.imod.eq.8).and.tnkbg.eq.1)) 
     &   ierr = 6
 
      if (imod.ne.3.and.imod.ne.7.and.imod.ne.8) call curs(0)

c popup status message

      call messg(ierr)
      call upajob

      return

100   call inferr('Error opening file !',0)
      return
      end

      subroutine tnkpnt(ipnt,iret)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      character*80 tnknm
      integer tnkbg,tnkit,tnkarc,tnkarf,tnkprg
      common /tnkopt/ rmsgrd,tnkbg,tnkit,tnkarc,tnkarf,tnkprg,tnknm,icst
      common /types/ iff
      integer dolabs,fancy,persp,shade,atcol,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*20 tnkfil
      character*20 lckfil
      character*3 tstr
      character*4 ttstr
      character*5 tttstr
      logical opfil

      iret = 0
      idummy = 0
      
      jlen = linlen(tnknm)
      jlent = jlen

      if (ipnt.eq.-1) then
         tnkfil = tnknm(1:jlen)//'.tmp'
         jlen = jlen + 1 + 3
      elseif (ipnt.lt.1000) then
         call zerstr(ipnt,tstr,3,0)
         tnkfil = tnknm(1:jlen)//'.'//tstr
         jlen = jlen + 1 + 3
      elseif (ipnt.lt.10000) then
         write(ttstr,'(i4)') ipnt
         tnkfil = tnknm(1:jlen)//'.'//ttstr
         jlen = jlen + 1 + 4
      elseif (ipnt.lt.100000) then
         write(ttstr,'(i5)') ipnt
         tnkfil = tnknm(1:jlen)//'.'//tttstr
         jlen = jlen + 1 + 5
      endif

      lckfil = tnkfil(jlent+1:jlen)//".ambforw"
      llen = jlen-jlent+8

      if (ixyz.eq.11) then
c ambfor
         do while (islck(ipnt).eq.1)
         end do

         if (ipnt.ne.1) then
            if (opfil(52,tnkfil,jlen,0,1,0)) then
                call rdbin(52,heat)
                iret = 1
                close(52)
            endif
         else

            if (opfil(52,tnkfil,jlen,1,1,0)) then
                iuntmp = iun2
                iun2 = 52
                call tnkfst(igttnk,0)
                if (igttnk.eq.1) iret = 1
                close(52)
                iun2 = iuntmp
            endif
         endif

      else

c tinker
         if (opfil(52,tnkfil,jlen,1,1,0)) then
             iuntmp = iun2
             iun2 = 52
             if (tnkprg.eq.4) then
                call rdchx(0,3,0,0,0,istat,icrtp)
                if (istat.gt.0) then
                   iret = 1
                   call fdat(ifd,0,0,0,0,0)
                endif
             else
                call gettnk(igttnk,0,idummy,iff,iheat,heat)
                if (igttnk.eq.1) iret = 1
             endif
             close(52)
             iun2 = iuntmp
         endif

      endif

      return
      end

      subroutine rdbid(iun,emin,coo,ianz,iatclr,iconn,ireord,
     &                  nat,norg,
     &                  xa,ya,yb,za,zb,zc,
     &                  a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb
      common /pbc/ abc(3),ibox,icell,igfmap
      real tmp,ar,br,cr,alphar,betar,gammar,emin
      dimension ireord(*),tmp(3)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*)

      toang = 0.52917706d0

c read tinker coordinates binary form

      read(iun) natoms,emin

      do i=1,natoms
           read(iun,err=100,end=100) (tmp(j),j=1,3)

           if (i.gt.nwramb) then
              ire = i
           else
              ire = ireord(i)
           endif

           if (ire.gt.0.and.ire.le.mxnat) then
              do j=1,3
                 coo(j,ire) = tmp(j)/toang
              end do
           endif
      end do

      if (icell.eq.1) then
           read(iun,err=100,end=100) ar,br,cr,alphar,betar,gammar
           a = dble(ar)
           b = dble(br)
           c = dble(cr)
           alpha = dble(alphar)
           beta = dble(betar)
           gamma = dble(gammar)
           call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,0)
           call cpmol2(nat,norg,xa,ya,yb,za,zb,zc,
     &                coo,ianz,iatclr,iconn)
           call updc(coo,xa,ya,yb,za,zb,zc)
      endif

100   continue
      return
      end

      subroutine tnkfsd(igttnk,idebug,coo,ianz,iatclr,iconn,ireord,
     &                  ichx,nat,norg,
     &                  xa,ya,yb,za,zb,zc,
     &                  a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,o-z), integer ( i-n)
      integer getlin
      character*137 str, tstr
      parameter (numatm=2000)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*137 line,tline
      common /curlin/ line
      common /pbc/ abc(3),ibox,icell,igfmap
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb
      logical dall,box,gnreal,owat
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*)
      dimension ireord(*),v1(3),v2(3),v3(3)

      igttnk = 1
      dall = .false.
      owat = .true.
      nwat =  0
      igaff = 0
      toang = 0.52917706d0
      box = .false.
      natoms = 0

      if (idebug.eq.1) print*,'subroutine fsttnk'
   
      if (getlin(0).eq.1) then

         tline = line

         if (icdex(line,'[AMBFOR]').ne.0) then
            igaff = 1
            ic = icdex(line,'box')
            if (ic.ne.0) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
               box = .true.
               do i=1,3
                  v1(i) = 0.0d0
                  v2(i) = 0.0d0
                  v3(i) = 0.0d0
               end do
               if (gnreal(abc,3,.false.)) then
                  ibox = 1
                  v1(1) = abc(1)/toang
                  v2(2) = abc(2)/toang
                  v3(3) = abc(3)/toang
               else
                  box = .false.
               endif
            endif

            ic = icdex(line,'cell')
            if (ic.ne.0) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
               icell = 1

               if (gnreal(abc,3,.false.)) then
                  a = abc(1)
                  b = abc(2)
                  c = abc(3)
                  if (gnreal(abc,3,.false.)) then
                     alpha = abc(1)
                     beta  = abc(2)
                     gamma = abc(3)
                  else
                     icell = 0
                  endif
               else
                  icell = 0
               endif
               if (icell.eq.1) then
                  call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,0)
               endif
            endif

            if (getlin(0).eq.1) idum = 1
         endif

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            natoms = itype
            if (natoms.gt.mxnat) then
               dall = .true.
               goto 100
            endif
         else
            goto 100
         endif
         
      else
         goto 100
      endif

      do iat=1,natoms
         if (getlin(0).eq.1) then

           tline = line

           ktype = nxtwrd(str,nstr,itype,rtype)
           ktype = nxtwrd(str,nstr,itype,rtype)

           if (iat.gt.nwramb) then
              ire = iat
              if (ktype.eq.1.and.owat) then
                 if (nstr.eq.2) then
                    if (str(1:2).eq.'OW'.or.str(1:2).eq.'HW') then
                        nwat = iat
                        owat = .false.
                    endif
                 endif
              endif
           else
              ire = ireord(iat)
              if (ire.eq.0) ire = iat
           endif

           if (ire.gt.0.and.ire.le.mxnat) then
              do i=1,3
                 ktype = nxtwrd(str,nstr,itype,rtype)
                 if (ktype.eq.3) then
                    coo(i,ire) = rtype/toang
                 else
                    goto 100
                 endif
              end do
           endif

         else
           print*,'Number of atoms read ',iat,
     &           ' less than specified in the header of file ',
     &           natoms
           goto 100
         endif

      end do

      if (icell.ne.1) then
         if (natoms.gt.nwramb) then
            iatoms = natoms
            call allon(nwramb,nwat,natoms)
            call doscal
            call docent
         endif
      endif


      if (box) then
         call addtbx(v1,v2,v3)
      endif

      if (icell.eq.1) then
         call cpmol2(nat,norg,xa,ya,yb,za,zb,zc,
     &              coo,ianz,iatclr,iconn)
         call updc(coo,xa,ya,yb,za,zb,zc)
      endif

      return

100   if (dall) then
         rewind iun2
         igttnk = -1
         return
      endif

      if (idebug.eq.1) print*,'ERROR:',tline
      igttnk = 0

      return
      end

      subroutine wrtsng
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb
      character*80 tnknm
      integer tnkbg,tnkit,tnkarc,tnkarf,tnkprg
      common /tnkopt/ rmsgrd,tnkbg,tnkit,tnkarc,tnkarf,tnkprg,tnknm,icst

      itnkln = linlen(tnknm)

      ixyz = 11
      call wrpnt(tnknm(1:itnkln)//'.xyz',itnkln+4,4,0,0,0,0,1,0)

      return
      end
