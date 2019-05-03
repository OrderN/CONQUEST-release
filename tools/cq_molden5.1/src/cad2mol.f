C
c 11 jul 01: header fixed; send to Gijs
c  7 jun 01: determine runtype, scftype, dft, mp2 from output
c  6 jun 01: uhf fixed
c  4 jun 01: prepares MOLDEN Format file from CADPAC output for:
c  { rhf | uhf }
c  { ene | opt | hess }
c  { dft | mp2 }

c
c Mariusz Klobukowski 
c cad2mol is a minor modification of Chris Lovallo's hon2mol
c
C This program will take a CADPAC output file and convert it to
C Molden Format for use with Molden (as filename.mdn). The program
C should work for SP energy, geometry optimizations, and Hessians for
C any single-determinant wavefunction

C For the orbitals to be read correctly, NOPRINT OCCVECTORS must not
c be used

      program cad2mol    
      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      character filenm*70
      parameter(inp=1,iout=2,itmp=3,maxbf=2000)
      common /flags/ flgrhf,flguhf,flgopn,flgdft,flggvb,flgmp2,flgmp4
      common /flags2/ flgmc
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      dimension occ(maxbf)

C Open necessary files (input, output, and temp)

      call open_files(inp,iout,itmp,filenm)

C Print title onto files

      call write_files(iout,itmp,'[Molden Format]')

C Read and setup parameters of run

      call find_flags(inp,occ)
      rewind(inp)

C Find and write coordinates and basis set

      call find_coord(inp,iout,itmp)

C Find and write eigenvectors and eigenvalues

      call find_eigen(inp,iout,occ)

C Next (and last for SP energy) is SCF energy convergence

      rewind(inp)
      call find_conv1(inp,iout)

C End if the run is SP energy, skip to frequencies if run is Hessian

      if (irun .eq. 0) goto 9999
      if (irun .eq. 2) goto 10

C Find and write out geometry optimization data

      call find_gopt(inp,iout)
      if (irun .eq. 3) goto 10
      goto 9999

C Find and write frequencies and normal modes (for Hessian runs)

 10   call find_freq(inp,iout,itmp)

C Close files

 9999 call close_files(inp,iout,itmp)

C That's it! It's over! Bye now!

      end
C-----------------------------------------------------------------------
      subroutine open_files(inp,iout,itmp,filenm)
      character filenm*70

C Open input file (CADPAC output file)

      write(6,*) 'Opening necessary files...'
      call getarg(1,filenm)
      if(filenm.eq.' ') then
        write(6,fmt='('' Enter the input file name (with .out): '',$)')
        read (5,fmt='(a70)') filenm
      endif
      open(unit=inp,file=filenm,access='SEQUENTIAL',
     1     status='OLD',form='FORMATTED')

C Open output and temp files

      do i=1,67
        if (filenm(i:i+3) .eq. '.out') filenm(i:i+3)='.mdn'
      enddo
      open(unit=iout,file=filenm,access='SEQUENTIAL',
     1     status='UNKNOWN',form='FORMATTED')
      open(unit=itmp,access='SEQUENTIAL',
     1     status='UNKNOWN',name='TEMP',form='FORMATTED')
c    1     status='SCRATCH',form='FORMATTED')
      return
      end
C-----------------------------------------------------------------------
      subroutine find_flags(inp,occ)
      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      character line*130
c     character kyword*26
      parameter(maxbf=2000)
      common /flags/ flgrhf,flguhf,flgopn,flgdft,flggvb,flgmp2,flgmp4
      common /flags2/ flgmc
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      dimension occ(*)

C This subroutine searches the CADPAC file for various data pertaining
C to the run (type of wavefunction, number of basis functions and
C occupied orbitals, etc.)

      write(6,*) 'Beginning initialization...'

C Initialize flags

      flgrhf=.false.
      flguhf=.false.
      flgopn=.false.
      flgdft=.false.
      flggvb=.false.
      flgmp2=.false.
      flgmp4=.false.
      flgmc= .false.

C Search for type of run and properties - decode from command line
C Determine type of wavefunction - decode from command line
c cad2mol job.out
c  { rhf | uhf | rohf | gvb | mp2 | mp4 }
c  { ene | opt | hess }
c  { dft }

*     inparg=iargc()
*     write(6,'('' inparg:'',i3)') inparg

      irun = -1
      iwfn = -1
*     do narg=2,inparg
*       call getarg(narg,kyword)
*       write(6,'(''   narg:'',i3,'' kyword:'',a)') narg,kyword
*       call up_case(kyword,26)
*       call last_char(kyword,26,loc)
*       if     (kyword(1:loc).eq.'RHF') then
*         flgrhf=.true.
*         iwfn=0
*       else if(kyword(1:loc).eq.'ROHF') then
*         flguhf=.true.
*         flgopn=.true.
*         iwfn=0
*       else if(kyword(1:loc).eq.'GVB') then
*         flguhf=.true.
*         flgopn=.true.
*         iwfn=0
*       else if(kyword(1:loc).eq.'UHF') then
*         flguhf=.true.
*         flgopn=.true.
*         iwfn=0
*       else if(kyword(1:loc).eq.'MP2') then
*         iwfn=5
*         flgmp2=.true.
*       else if(kyword(1:loc).eq.'MP4') then
*         iwfn=6
*         flgmp4=.true.
*       else if(kyword(1:loc).eq.'CAS') then
*         iwfn=1
*       else if(kyword(1:loc).eq.'DFT') then
*         flgdft=.true.
*       else if(kyword(1:loc).eq.'ENE') then
*         irun=0
*       else if(kyword(1:loc).eq.'OPT') then
*         irun=1
*       else if(kyword(1:loc).eq.'HESS') then
*         irun=2
*       else if(kyword(1:loc).eq.'OPHS') then
*         irun=3
*       else
*         write(6,*) 'wrong parameter on command line'
*         call error
*       endif
*     enddo

      rewind(inp)
      ifound=0
      call search(inp,line,'RUN TYPE = ENERGY',ifound,0)
      if(ifound.eq.1) then
        write(6,*) 'RUN TYPE = ENERGY'
        irun=0
      endif
      rewind(inp)
      ifound=0
      call search(inp,line,'RUN TYPE = OPTIMIZE',ifound,0)
      if(ifound.eq.1) then
        write(6,*) 'RUN TYPE = OPTIMIZE'
        irun=1
      endif
      rewind(inp)
      ifound=0
      call search(inp,line,'RUN TYPE = SECDER',ifound,0)
      if(ifound.eq.1.and.irun.eq.-1) then
        irun=2
        write(6,*) 'RUN TYPE = SECDER -- only'
      endif
      if(ifound.eq.1.and.irun.eq.1) then
        irun=3
        write(6,*) 'RUN TYPE = SECDER -- after OPTIMIZE'
      endif
      rewind(inp)
      ifound=0
      call search(inp,line,'RUN TYPE = FORCE',ifound,0)
      if(ifound.eq.1.and.irun.eq.-1) then
        irun=2
        write(6,*) 'RUN TYPE = FORCE -- only'
      endif
      if(ifound.eq.1.and.irun.eq.1) then
        irun=3
        write(6,*) 'RUN TYPE = FORCE -- after OPTIMIZE'
      endif
      rewind(inp)
      ifound=0
      call search(inp,line,'SCF TYPE = CLOSED',ifound,0)
      if(ifound.eq.1) then
        flgrhf=.true.
        iwfn=0
        write(6,*) 'SCF TYPE = CLOSED'
      endif
      rewind(inp)
      ifound=0
      call search(inp,line,'UHF Calculation',ifound,0)
      if(ifound.eq.1) then
        iwfn=0
        flguhf=.true.
        flgopn=.true.
        write(6,*) 'UHF Calculation'
      endif 
      rewind(inp)
      ifound=0
      call search(inp,line,'CLOSED-SHELL KOHN-SHAM',ifound,0)
      if(ifound.eq.1) then
        flgdft=.true.
        write(6,*) 'CLOSED-SHELL KOHN-SHAM'
      endif
      rewind(inp)
      ifound=0
      call search(inp,line,'Unrestricted KOHN-SHAM',ifound,0)
      if(ifound.eq.1) then
        iwfn=0
        flgdft=.true.
        flguhf=.true.
        flgopn=.true.
        write(6,*) 'Unrestricted KOHN-SHAM'
      endif
      rewind(inp)
      ifound=0
      call search(inp,line,'RHF MP2 Calculation',ifound,0)
      if(ifound.eq.1) then
        iwfn=5
        flgmp2=.true.
        write(6,*) 'RHF MP2 Calculation'
      endif
      rewind(inp)
      ifound=0
      call search(inp,line,'UHF MP2 Calculation',ifound,0)
      if(ifound.eq.1) then
        iwfn=5
        flgmp2=.true.
        flguhf=.true.
        flgopn=.true.
        write(6,*) 'UHF MP2 Calculation'
      endif
      
      write(6,fmt='('' irun   ='',i2)') irun
      write(6,fmt='('' iwfn   ='',i2)') iwfn
      write(6,fmt='('' flgrhf ='',L2)') flgrhf
      write(6,fmt='('' flguhf ='',L2)') flguhf
      write(6,fmt='('' flgopn ='',L2)') flgopn
      write(6,fmt='('' flgdft ='',L2)') flgdft
      write(6,fmt='('' flgmp2 ='',L2)') flgmp2

      if (irun.ne.0. and.
     1    irun.ne.1 .and.
     1    irun.ne.2 .and.
     1    irun.ne.3) then
        write(6,*)  'Only RUNTYPs=0-3 (SP energy, geom opt, and Hessian)
     1 are currently implemented!'
        call error
      endif

      if (iwfn .ne. 0 .and. iwfn .ne. 5 .and. iwfn .ne. 1) then
       write(6,*)   'Only single-determinant and MCSCF wavefunctions are
     1 currently implemented!'
       call error
      endif

      rewind (inp)

C Find number of basis functions

      call search(inp,line,'Number of Basis Functions',ifound,-1)
      read(unit=line,fmt='(46x,i5)') nbf
      write(*,*) 'nbf  =',nbf
      if (nbf .ge. maxbf) then
        write(6,*)   'Error code 5! Number of basis functions is greater
     1 than parameter MAXBF!'
        write(6,*) 'Check parameter MAXBF and run program again.'
        call error
      endif


C Determine type of run and orbital occupancies

      if (iwfn .eq. 1 .or. iwfn .eq. 10 .or. iwfn .eq. 12) flgmc=.true.
      if (flgmc .or. iwfn .eq. 2) then
        call search(inp,line,'NUMBER OF ORBITALS',ifound,-1)
        read(unit=line,fmt='(27x,i5)') norb
        call search(inp,line,'# OF CORE   ORBITALS',ifound,-1)
        read(unit=line,fmt='(26x,i5)') ndoc  
        goto 20
      endif

C The program does not search for MP2/4 lines, it sets the
C flags according to iwfn

 10   if (iwfn .eq. 5) flgmp2=.true.
      if (iwfn .eq. 6) flgmp4=.true.
      if(flguhf) then
        call search(inp,line,'occupied alpha orbitals',ifound,-1)
        read(unit=line,fmt='(46x,i5)') ndoc
      else
        call search(inp,line,
     1       'Number of doubly occupied orbitals',ifound,-1)
        read(unit=line,fmt='(46x,i5)') ndoc
      endif
      write(*,*) 'ndoc =',ndoc
      if (.not. flgopn) goto 20
      if (flguhf) then
        call search(inp,line,'occupied beta orbitals',ifound,-1)
        read(unit=line,fmt='(46x,i5)') nsoc
        write(*,*) 'nsoc=',nsoc
        goto 20
      endif

 20   if (.not. flguhf) then
        do i=1,ndoc
          occ(i)=2.0
        enddo
        if (nsoc .gt. 0) then
          do i=ndoc+1,ndoc+1+nsoc
            occ(i)=1.0
          enddo
          do i=ndoc+1+nsoc,nbf
            occ(i)=0.0
          enddo
        else
          do i=ndoc+1,nbf
            occ(i)=0.0
          enddo
        endif
        goto 30
      else
        do i=1,ndoc
          occ(i)=1.0
        enddo
        do i=ndoc+1,nbf
          occ(i)=0.0
        enddo
      endif

C Determine maximum number of iterations (SCF or CI)  

 30   if (iwfn .eq. 2) then
        call search(inp,line,'MAX. NUMB. OF ITERATIONS',ifound,-1)
        read(unit=line,fmt='(29x,i6)') maxit
      else
        if(flguhf) then
         call search(inp,line,'Maximum Number of Iterations',ifound,-1)
        else
         call search(inp,line,'Maximum number of iterations',ifound,-1)
        endif
        read(unit=line,fmt='(32x,i6)') maxit
        write(*,*) 'maxit=',maxit
      endif
      maxit=maxit+1

      write(6,*) 'Initialization done!'
      return
      end
C-----------------------------------------------------------------------
      subroutine find_coord(inp,iout,itmp)
      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      character line*130,lintmp*130,atmnam*8
      parameter(maxat=100)
      dimension atmnam(maxat)

C This subroutine finds and writes out the coordinates of the atoms
C It searches for CHARGE, and then skips the --- and blank lines
C It then calls find_bas to write out the basis set

      write(6,*) 'Reading coordinates...'
      call write_files(iout,itmp,'[Atoms] AU')
      call search(inp,line,'Atom     Nuclear',ifound,-1)
      call skip(inp,2)
      i=1
 10   call read_line(inp,line)
      read(unit=line,fmt='(4x,a8)') atmnam(i)  
      if (atmnam(i) .eq. '--------') goto 20
      read(unit=line,fmt='(4x,a8,2x,f5.1,3f15.8)') atmnam(i),chrg
     1 ,x,y,z
      ichrg=int(chrg)
      write(*,fmt='(1x,a8,i3,i3,3(2x,f12.8))') atmnam(i),i,ichrg,x,y,z
      write(unit=lintmp,fmt='(1x,a8,i3,i3,3(2x,f12.7))') atmnam(i),i
     1      ,ichrg,x,y,z
      call write_files(iout,itmp,lintmp)
      i=i+1
      goto 10
 20   nat=i-1
      write(*,*) 'nat=',nat
      write(6,*) 'Done writing coordinates!'
      call find_bas(inp,iout,itmp,atmnam)
      return
      end
C-----------------------------------------------------------------------
      subroutine find_bas(inp,iout,itmp,atmnam)
      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      character line*130,lintmp*130,temp8*8,temp4*4,atmnam*8,shltyp*4
      parameter(mxprim=30)
      dimension atmnam(*),coef1(mxprim),coef2(mxprim),exp(mxprim)
 
C This subroutine reads and writes out the basis set, including looping
C over the non-symmetry unique atoms
 
      write(6,*) 'Reading basis set...'
      call skip(inp,6)
      call write_files(iout,itmp,'[GTO]')
      do j=1,nat
*     write(*,*) j,atmnam(j)
c ... reset cgtf count for each new atom
        ikount=0
        write(unit=lintmp,fmt='(i3,a5)') j, ' 0   '
        call write_files(iout,itmp,lintmp)
        call read_line(inp,line)
*       write(*,*) '>'//line(1:76)
        read(unit=line,fmt='(1x,a8)') temp8
*       write(*,*) 'temp8'//temp8
        if (temp8 .ne. atmnam(j)) then
C Either an error has occurred, or this is not a symmetry unique atom.
C Checking for an atom with the same name and copying basis set...
          backspace(inp)
          temp8=atmnam(j)
          do k=1,j-1
            if (temp8 .eq. atmnam(k)) then
              rewind(itmp)
              call search(itmp,line,atmnam(k),ifound,-1)
              read(unit=line,fmt='(10x,i3)') n
              write(unit=temp8,fmt='(i3,a5)') n, ' 0   '
              call search(itmp,line,temp8,ifound,-1)
              write(unit=temp8,fmt='(i3,a5)') n+1, ' 0   '
 10           call read_line(itmp,line)
*             write(*,*) '>'//line(1:76)
              if (line(1:8) .ne. temp8) then
                call write_files(iout,0,line)
                goto 10
              endif
              call ffwrd(itmp)
              goto 40
            endif
          enddo
 
C Nope! There was an error! Stopping...
 
          write(6,*) 'Error code 1 reading basis set!'
          call error
        endif
        
        l=1
 20     call read_line(inp,line)
*       write(*,*) '>'//line(1:76)
        if(line(1:6).eq.' Model') call read_line(inp,line)
        read(unit=line,fmt='(12x,i2,5x,a4,6x,f16.6,f15.6)')
     1     jkount,temp4,exp(l),coef1(l)
*       write(*,*) 'l ikount jkount',l,ikount,jkount,exp(l)
        if(jkount.eq.0) then
*          backspace(inp)
           ikount=jkount
           goto 30
        endif
        if(ikount.gt.0 .and. jkount.ne.ikount) then
          ikount=jkount
          backspace(inp)
          goto 30
        else
          ikount=jkount
        endif
        shltyp=temp4
        if (shltyp(1:1) .eq. 'L') then
          read(unit=line,fmt='(19x,a4,6x,f16.6,2f15.6)')
     1         shltyp,exp(l),coef1(l),coef2(l)
        endif
        l=l+1
        goto 20
 30     nprim=l-1
*       write(*,*) 'nprim exp',nprim,exp(nprim)
 
C Write basis set
 
        if (shltyp(1:1) .ne. 'L') then
          write(unit=lintmp,fmt='(a2,i5,a5)') shltyp(1:1),nprim, ' 1.00'
          call write_files(iout,itmp,lintmp)
          do l=1,nprim
            write(unit=lintmp,fmt='(2(1x,d17.10))') exp(l),coef1(l)
            call write_files(iout,itmp,lintmp)
          enddo
        elseif (shltyp(1:1) .eq. 'L') then
          write(unit=lintmp,fmt='(1x,a2,i4,a5)') 'SP',nprim,
     1          ' 1.00'
          call write_files(iout,itmp,lintmp)
          do l=1,nprim
            write(unit=lintmp,fmt='(3(1x,d17.10))') exp(l),coef1(l)
     1            ,coef2(l)
            call write_files(iout,itmp,lintmp)
          enddo
        else
          write(6,*) 'Error code 3 writing basis set!'
          call error
        endif
        if(jkount.eq.0) goto 40
        l=1
        goto 20
 40     call blank_line(iout)
      enddo
      call blank_line(iout)
      write(6,*) 'Done writing basis set!'
      return
      end
C-----------------------------------------------------------------------
      subroutine find_eigen(inp,iout,occ)
      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      common /flags/ flgrhf,flguhf,flgopn,flgdft,flggvb,flgmp2,flgmp4
      common /flags2/ flgmc
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      character line*130,lintmp*130
      parameter(maxbf=2000)
      parameter (zero=0.0d+00, one=1.0d+00)
      dimension occ(*),orbs(maxbf,maxbf),eval(maxbf)

C This subroutine finds and writes out the molecular orbitals and
C orbital energies, along with their occupancies

      write(6,*) 'Reading orbitals and eigenvalues...'
      call write_files(iout,0,'[MO]')
c ... 1st MO output
      first_SCF=.true.
      if(flguhf) then
        call search(inp,line,'Alpha-spin Orbitals',ifound,-1)
      else
        call search(inp,line,'Molecular Orbitals',ifound,-1)
      endif
c ... 2nd MO output - used if opt+hess in the same run
      if (irun .eq. 1 .or. irun.eq.2 .or. irun.eq.3) then
        if(flguhf) then
          call search(inp,line,'Alpha-spin orbitals',ifound,0)
        else
          call search(inp,line,'Molecular orbitals',ifound,0)
        endif
      endif
      first_SCF=.false.
      if(ifound.eq.0) then
        write(*,*) '... only hessian'
        first_SCF=.true.
        rewind(inp)
        if(flguhf) then
          call search(inp,line,'Alpha-spin Orbitals',ifound,-1)
        else
          call search(inp,line,'Molecular Orbitals',ifound,-1)
        endif
      endif
c ... locate the 1st eigenvalue
      call search(inp,line,'     1    ',ifound,-1)
      backspace(inp)
c ... read the number of MOs printed
      found=.false.
      nmos=0
      do while (.not.found)
        call read_line(inp,line)
        if(line(1:20).ne.'                    ' .and.
     1     line(1:13).ne.' Eigenvectors') then
          nmos=nmos+1
        else
          found=.true.
          backspace(inp)
        endif
      enddo
      write(*,*) 'nmos=',nmos
      write(*,*) 'first_SCF ',first_SCF
      if(.not.flguhf) then
        if(irun.eq.0.or.first_SCF) call skip(inp,2)
      endif
      if(flguhf.and.first_SCF) call skip(inp,1)
      if (flggvb) call skip(inp,1)
      j=1
 10   k=j+4
      call skip(inp,2)
      call read_line(inp,line)
      read(unit=line,fmt='(12x,5(f13.6))') (eval(i), i=j,k)
      write(*,fmt='(5(i4,f13.6))') (i,eval(i), i=j,k)
      call skip(inp,1)
      do l=1,nbf
        call read_line(inp,line)
        read(unit=line,fmt='(12x,5(f13.6))') (orbs(i,l), i=j,k)
      enddo
      if (k .ge. nmos) then
        do i=1,nmos
          write(unit=lintmp,fmt='(a5,2x,f9.4)') ' Ene=',eval(i)
          call write_files(iout,0,lintmp)
          call write_files(iout,0,' Spin= Alpha')
          if(flguhf.and.i.gt.ndoc) occ(i)=zero
          write(unit=lintmp,fmt='(a7,f11.6)') ' Occup=',occ(i)
          call write_files(iout,0,lintmp)
          do j=1,nbf
            write(unit=lintmp,fmt='(i4,f11.6)') j,orbs(i,j)
            call write_files(iout,0,lintmp)
          enddo
        enddo
        if (.not. flguhf) goto 20

C Now read beta orbitals (for UHF runs only)

        write(6,*) 'Done writing alpha orbitals...'
        call find_beta(inp,iout,occ)
        write(6,*) 'Done writing beta  orbitals...'
        goto 20
      endif
      j=j+5
      goto 10
 20   write(6,*) 'Done writing orbitals and eigenvalues!' 
      return
      end
C-----------------------------------------------------------------------
      subroutine find_beta(inp,iout,occ)

C This subroutine finds and write the beta molecular orbitals and
C orbital energies for a UHF wavefunction

      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      common /flags/ flgrhf,flguhf,flgopn,flgdft,flggvb,flgmp2,flgmp4
      common /flags2/ flgmc
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      character line*130,lintmp*130
      parameter(maxbf=2000)
      dimension occ(*),orbs(maxbf,maxbf),eval(maxbf)

c ... locate start of current beta's
      call search(inp,line,'Beta-spin',ifound,-1)
c ... locate 1st beta MO
      call search(inp,line,'     1     ',ifound,-1)
      backspace(inp)
c ... read the number of MOs printed
      found=.false.
      nmos=0
      do while (.not.found)
        call read_line(inp,line)
        if(line(1:20).ne.'                    '.and.
     1     line(1:13).ne.' Eigenvectors') then
          nmos=nmos+1
        else
          found=.true.
          backspace(inp)
        endif
      enddo
      write(*,*) 'nmos=',nmos
      call read_line(inp,line)
      if(line(1:13).ne.' Eigenvectors') backspace(inp)
      j=1
 10   k=j+4
      call skip(inp,2)
      call read_line(inp,line)
      read(unit=line,fmt='(12x,5(f13.6))') (eval(i), i=j,k)
      write(*,fmt='(5(i4,f13.6))') (i,eval(i), i=j,k)
      call skip(inp,1)
      do l=1,nbf
        call read_line(inp,line)
        read(unit=line,fmt='(12x,5(f13.6))') (orbs(i,l), i=j,k)
      enddo
      if (k .ge. nmos) then
        do i=1,nmos
          write(unit=lintmp,fmt='(a5,2x,f9.4)') ' Ene=',eval(i)
          call write_files(iout,0,lintmp)
          call write_files(iout,0,' Spin= Beta')
          if(flguhf.and.i.gt.nsoc) occ(i)=zero
          write(unit=lintmp,fmt='(a7,f11.6)') ' Occup=',occ(i)
          call write_files(iout,0,lintmp)
          do j=1,nbf
            write(unit=lintmp,fmt='(i4,f11.6)') j,orbs(i,j)
            call write_files(iout,0,lintmp)
          enddo
        enddo
        goto 20
      endif
      j=j+5
      goto 10
 20   return
      end
C----------------------------------------------------------------------
      subroutine find_eigmc(inp,iout,occ)
      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      common /flags/ flgrhf,flguhf,flgopn,flgdft,flggvb,flgmp2,flgmp4
      common /flags2/ flgmc
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      character line*130,lintmp*130
      parameter(maxbf=2000)
      dimension occ(*),orbs(maxbf,maxbf),eval(maxbf)

C This subroutine finds and writes out the natural orbitals and their
C occupancies (energies are all set to 0 if the Fock energies are not
C printed)

      do i=1,nbf
        eval(i)=0.0
      enddo
      rewind(inp)
      call search(inp,line,'EIGENVALUES OF -MC-',ifound,0)
      if (ifound .eq. 1) then
        call skip(inp,2)
        j=1
 10     k=j+5
        call read_line(inp,line)
        read(unit=line,fmt='(6(10x,f11.4))') (eval(i), i=j,k)
        j=j+6
        if (k .lt. nbf) goto 10
      endif
      rewind(inp)
      call write_files(iout,0,'[MO]')
      if (irun .ne. 1) goto 30
      n=0
 20   call search(inp,line,'CAN.CORE + ) NAT.VAL',ifound,0)
      if (ifound .eq. 1) then
        n=n+1
        goto 20
      endif
      rewind(inp)
      do i=1,n
        call search(inp,line,'CAN.CORE + ) NAT.VAL',ifound,-1)
      enddo
      goto 40
 30   call search(inp,line,'CAN.CORE + ) NAT.VAL',ifound,-1)
 40   call skip(inp,9)
      call read_line(inp,line)
      j=1
 50   k=j+6
      read(unit=line,fmt='(15x,7(f15.10))') (occ(i), i=j,k)
      call skip(inp,2)
      do l=1,nbf
        call read_line(inp,line)
        read(unit=line,fmt='(15x,7(f15.10))') (orbs(i,l), i=j,k)
      enddo
      if (k .ge. norb) then
        do i=1,norb
          write(unit=lintmp,fmt='(a5,2x,f9.4)') ' Ene=',eval(i)
          call write_files(iout,0,lintmp)
          call write_files(iout,0,' Spin= Alpha')
          write(unit=lintmp,fmt='(a7,f11.6)') ' Occup=',occ(i)
          call write_files(iout,0,lintmp)
          do j=1,nbf
            write(unit=lintmp,fmt='(i4,f11.6)') j,orbs(i,j)
            call write_files(iout,0,lintmp)
          enddo
        enddo
        goto 60
      endif
      j=j+7
      call skip(inp,8)
      call read_line(inp,line)
      goto 50
 60   return
      end
C-----------------------------------------------------------------------
      subroutine find_eigci(inp,iout,occ)
      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      common /flags/ flgrhf,flguhf,flgopn,flgdft,flggvb,flgmp2,flgmp4
      common /flags2/ flgmc
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      character line*130,lintmp*130
      parameter(maxbf=2000)
      dimension occ(*),orbs(maxbf,maxbf),eval(maxbf)

C This subroutine finds and writes out the natural orbitals and their
C occupancies (energies are all set to 0)

      do i=1,nbf
        eval(i)=0.0
      enddo
      rewind(inp)
      call write_files(iout,0,'[MO]')
      if (irun .ne. 1) goto 30
      n=0
 20   call search(inp,line,'NATURAL ORBITALS IN ATOMIC',ifound,0)
      if (ifound .eq. 1) then
        n=n+1
        goto 20
      endif
      rewind(inp)
      do i=1,n
        call search(inp,line,'NATURAL ORBITALS IN ATOMIC',ifound,-1)
      enddo
      goto 40
 30   call search(inp,line,'NATURAL ORBITALS IN ATOMIC',ifound,-1)
 40   call skip(inp,6)
      call read_line(inp,line)
      j=1
 50   k=j+6
      read(unit=line,fmt='(15x,7(f15.10))') (occ(i), i=j,k)
      call skip(inp,2)
      do l=1,nbf
        call read_line(inp,line)
        read(unit=line,fmt='(15x,7(f15.10))') (orbs(i,l), i=j,k)
      enddo
      if (k .ge. norb) then
        do i=1,norb
          write(unit=lintmp,fmt='(a5,2x,f9.4)') ' Ene=',eval(i)
          call write_files(iout,0,lintmp)
          call write_files(iout,0,' Spin= Alpha')
          write(unit=lintmp,fmt='(a7,f11.6)') ' Occup=',occ(i)
          call write_files(iout,0,lintmp)
          do j=1,nbf
            write(unit=lintmp,fmt='(i4,f11.6)') j,orbs(i,j)
            call write_files(iout,0,lintmp)
          enddo
        enddo
        goto 60
      endif
      j=j+7
      call skip(inp,5)
      call read_line(inp,line)
      goto 50
 60   return
      end
C-----------------------------------------------------------------------
      subroutine find_conv1(inp,iout)
      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      common /flags/ flgrhf,flguhf,flgopn,flgdft,flggvb,flgmp2,flgmp4
      common /flags2/ flgmc
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      character line*130,lintmp*130,temp4*4
      parameter(maxits=101)
      dimension etot(maxits)

C This subroutine finds and writes out the SCF energy convergence
C (the first one)

      write(6,*) 'Reading energy convergence...' 
      if(flgdft) then
        write(6,*) '... DFT calculations - no energy convergence'
        return
      endif
      if (maxit .ge. maxits) then
        write(6,*)   'Error code 5! Number of iterations is greater than
     1 parameter MAXITS!'
        write(6,*) 'Check parameter MAXITS and run program again.'
        call error
      endif
      if (flgmc) then
        call search(inp,line,' CYCLE  TOTAL ENERGY',ifound,-1)
      else if (iwfn .eq. 2) then
        call search(inp,line,' ITER.    MAX.DEV.',ifound,-1)
        call skip(inp,1)
      else
*       call search(inp,line,'Cycle      Total energy',ifound,-1)
        call search(inp,line,'Cycle     ',ifound,-1)
        if(.not.flguhf) call skip(inp,1)
      endif
      i=0
      finish=.false.
      do while(.not.finish)
        call read_line(inp,line)
        if (iwfn .eq. 2) then
          read(unit=line,fmt='(a4,23x,f16.9)') temp4,etot(i)
        else
          if(line(1:5).ne.'     '.and.line(1:4).ne.' Den'
     1                           .and.line(1:5).ne.' Last'
     1                           .and.line(1:5).ne.' HOMO') then
            i=i+1
            read(unit=line,fmt='(5x,f18.9)') etot(i)
*           write(*,*) 'i e',i,etot(i)
          endif
        endif
        if (line(1:4) .eq. ' Den') then
          finish=.true.
        endif
      enddo
      ncyc=i
 10   call write_files(iout,0,'[SCFCONV]')
      write(unit=lintmp,fmt='(a22,i5)') 'scf-first    1 THROUGH',ncyc
      call write_files(iout,0,lintmp)
      do i=1,ncyc
       write(unit=lintmp,fmt='(f14.6)') etot(i)
        call write_files(iout,0,lintmp)
      enddo
      write(6,*) 'Finished writing energy convergence!'
      return
      end
C-----------------------------------------------------------------------
      subroutine find_gopt(inp,iout)

C This subroutine finds and writes the data pertaining to the geometry
C optimization process.

      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      common /flags/ flgrhf,flguhf,flgopn,flgdft,flggvb,flgmp2,flgmp4
      common /flags2/ flgmc
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      character line*130,lintmp*130,temp4*4,temp8*8,atmnam*8
      parameter(maxstp=200,bohr2a=5.29177249D-01,maxits=101,maxat=100)
      dimension estep(maxstp),xmxfor(maxstp),rmsfor(maxstp)
      dimension ecyc(maxits),atmnam(maxat),geom(3,maxat,maxstp)
*     dimension grdint(5*maxat)

      write(6,*) 'Reading geometry optimization information...'

c ... DFT calculations do not show SCF convergence 
c ... irrelevant for MP2
      if(flgdft .or. flgmp2) goto 25
      rewind(inp)
      nstep=0
 10   if (flgmc) then
        call search(inp,line,' CYCLE  TOTAL ENERGY',ifound,0)
      else if (iwfn .eq. 0) then
        if(flguhf) then
          call search(inp,line,' UHF Calculation',ifound,0)
        else
          call search(inp,line,' CLOSED-SHELL SCF CALCULATION',ifound,0)
        endif
      else
        call search(inp,line,' Cycle     ',ifound,0)
      endif
      if (ifound .eq. 1) then
        nstep=nstep+1
        goto 10
      endif
      write(*,*) 'nstep=',nstep

c ... position to the last SCF
      rewind(inp)
      do i=1,nstep
        if (flgmc) then
          call search(inp,line,' CYCLE  TOTAL ENERGY',ifound,-1)
        else if (iwfn .eq. 0) then
          if(flguhf) then
          call search(inp,line,' UHF Calculation',ifound,0)
          else
          call search(inp,line,' CLOSED-SHELL SCF CALCULATION',ifound,0)
          endif
          call skip(inp,1)
        else
          call search(inp,line,' CYCLE   TOTAL ENERGY',ifound,-1)
          call skip(inp,1)
        endif
      enddo
c ... located last SCF; find 1st scf cycle
      call search(inp,line,'  1  ',ifound,-1)
      backspace(inp)
      i=0
      finish=.false.
      do while(.not.finish)
        call read_line(inp,line)
        if (iwfn .eq. 2) then
          read(unit=line,fmt='(a4,23x,f16.9)') temp4,ecyc(i)
        else
          if(line(1:5).ne.'     '.and.line(1:4).ne.' Den'
     1                           .and.line(1:5).ne.' Last') then
            i=i+1
            read(unit=line,fmt='(5x,f18.9)') ecyc(i)
            write(*,*) 'i e',i,ecyc(i)
          endif
        endif
        if (line(1:4) .eq. ' Den') then
          finish=.true.
        endif
      enddo
      ncyc=i
      write(*,*) 'ncyc=',ncyc
c ... write scf convergence - last SCF
 20   write(unit=lintmp,fmt='(a21,i5)') 'scf-last    1 THROUGH',ncyc
      call write_files(iout,0,lintmp)
      do i=1,ncyc
        write(unit=lintmp,fmt='(f14.6)') ecyc(i)
        call write_files(iout,0,lintmp)
      enddo

  25  continue

C The program only uses those steps that have a progress line at the end
C of them (one with NSERCH, etc.), so nstep above may not be the nstep
C used later (see line 30)

      rewind(inp)
      nstep=0
 11   continue
        call search(inp,line,' RMS of gradient    ',ifound,0)
        if (ifound .eq. 1) then
          nstep=nstep+1
          goto 11
        endif
      write(*,*) 'n RMS=',nstep

      rewind(inp)
      do i=1,nstep
        call search(inp,line,' RMS of gradient',ifound,0)
        if (ifound .eq. 0)  goto 30
        call backs(inp,line,' Molecular Geometry',ifound,0)
        call skip(inp,2)
        do j=1,nat
          call read_line(inp,line)
          read(unit=line,fmt='(6x,a8,2x,3f20.10)') atmnam(j),
     1         geom(1,j,i),geom(2,j,i), geom(3,j,i)
          write(*,*) j,geom(1,j,i),geom(2,j,i), geom(3,j,i)
          do k=1,3
            geom(k,j,i) = geom(k,j,i)*bohr2a
          enddo
        enddo
        call search(inp,line,' Current energy    ',ifound,-1)
        read(unit=line,fmt='(32x,f19.10)') estep(i)
        write(*,*) 'estep(i)',i,estep(i)
        call skip(inp,2)
        call read_line(inp,line)
        read(unit=line,fmt='(32x,f13.6)') rmsfor(i)
        write(*,*) 'rmsfor(i)',rmsfor(i)
        call read_line(inp,line)
        read(unit=line,fmt='(32x,f13.6)') xmxfor(i)
        write(*,*) 'xmxfor(i)',xmxfor(i)
      enddo
      nstep=i
      write(*,*) 'nstep',nstep

C Now write everything out

 30   nstep=i-1
      call write_files(iout,0,'[GEOCONV]')
      call write_files(iout,0,'energy')
      do i=1,nstep
        write(unit=lintmp,fmt='(1x,f13.6)') estep(i)
        call write_files(iout,0,lintmp)
      enddo
      call write_files(iout,0,'max-force')
      do i=1,nstep
        write(unit=lintmp,fmt='(4x,f10.6)') xmxfor(i)
        call write_files(iout,0,lintmp)
      enddo
      call write_files(iout,0,'rms-force')
      do i=1,nstep
        write(unit=lintmp,fmt='(4x,f10.6)') rmsfor(i)
        call write_files(iout,0,lintmp)
      enddo

C Note that the program prints out only the Cartesian geometries 
C (no matter what coordinates were used in the optimization)

      call write_files(iout,0,'[GEOMETRIES] XYZ')
      write(unit=temp8,fmt='(i5)') nat
      do i=1,nstep
        write(unit=lintmp,fmt='(a9,f13.6)') 'scf done:', estep(i)
        call write_files(iout,0,temp8)
        call write_files(iout,0,lintmp)
        do j=1,nat
          write(unit=lintmp,fmt='(1x,a8,3(f11.6,2x))') atmnam(j),
     1          geom(1,j,i),geom(2,j,i), geom(3,j,i)
          call write_files(iout,0,lintmp)
        enddo
      enddo
      write(6,*) 'Finished writing geometry optimization information!'
      return
      end
C-----------------------------------------------------------------------
      subroutine find_freq(inp,iout,itmp)

C This subroutine finds and writes the frequencies and normal modes of
C vibration.

      implicit double precision(a-e,g-h,o-z)
      implicit logical(f)
      common /flags/ flgrhf,flguhf,flgopn,flgdft,flggvb,flgmp2,flgmp4
      common /flags2/ flgmc
      common /rundat/ irun,iwfn,nbf,ndoc,nsoc,maxit,nat,norb
      character line*130,lintmp*130
      parameter(maxat=100)
      dimension efreq(3*maxat),vibcor(3*maxat,3*maxat)

      write(6,*) 'Reading vibrational frequencies and displacements...'
      call write_files(iout,0,'[FREQ]')
      rewind(inp)
      call search(inp,line,'Transformation matrix from',ifound,-1)
      j=1
 10   k=j+5
      call skip(inp,3)
      call read_line(inp,line)
      read(unit=line,fmt='(17x,6f10.3)') (efreq(i), i=j,k)
      write(*,fmt='(6(i3,f9.3))') (i,efreq(i), i=j,k)
      call skip(inp,2)
      do l=1,3*nat
        call read_line(inp,line)
        read(unit=line,fmt='(17x,6f10.5)') (vibcor(i,l),i=j,k)
*       write(*,fmt='(6f10.5)') (vibcor(i,l),i=j,k)
      enddo
      if (k .ge. 3*nat) then
        do i=1,3*nat
          write(unit=lintmp,fmt='(f10.4)') efreq(i)
          call write_files(iout,0,lintmp)
        enddo
        call write_files(iout,0,'[FR-COORD]')
        rewind(itmp)
        call skip(itmp,2)
        do i=1,nat
          call read_line(itmp,line)
          read(unit=line,fmt='(1x,a8,6x,3(2x,f12.7))') atmnam,x,y,z
          write(unit=lintmp,fmt='(a8,3(1x,f12.6))') atmnam,x,y,z
          call write_files(iout,0,lintmp)
        enddo
        call write_files(iout,0,'[FR-NORM-COORD]')
        do i=1,3*nat
          write(unit=lintmp,fmt='(a9,i5)') 'vibration',i
          call write_files(iout,0,lintmp)
          do j=1,3*nat,3
            write(unit=lintmp,fmt='(3(f12.6,1x))') vibcor(i,j),
     1            vibcor(i,j+1),vibcor(i,j+2)
            call write_files(iout,0,lintmp)
          enddo
        enddo
      goto 20
      endif
      j=j+6
      goto 10
 20   write(6,*)           'Finished writing vibrational frequencies and
     1 displacements!' 
      return
      end
C----------------------------------------------------------------------
      subroutine close_files(inp,iout,itmp)  

C Closes all open files (temp file is deleted)

      write(6,*) 'Closing files...'
      close(unit=inp, err=10, status='KEEP')
      close(unit=iout, err=20, status='KEEP')
      close(unit=itmp, err=30, status='DELETE')
      goto 40
 10   write(6,*) 'An error occurred while closing the input file!'
      call error
 20   write(6,*) 'An error occurred while closing the output file!'
      call error
 30   write(6,*) 'An error occurred while closing the temp file!'
      call error
 40   write(6,*) 'All done! Ending program...' 
      return
      end
C-----------------------------------------------------------------------
      subroutine blank_line(iout)

C Writes a blank line to unit=iout.

      write(iout,fmt='(a)') ' '
      return
      end
C----------------------------------------------------------------------
      subroutine error
      write(6,*) 'An error has occurred while processing the file.'
      write(6,*) 'Aborting run...'
      stop
      end
C-----------------------------------------------------------------------
      subroutine ffwrd(iunit)

C Fast-forward unit=iunit to the end of the file

 10   read(iunit,fmt='(a130)',end=20) line
      goto 10
 20   return
      end
C----------------------------------------------------------------------
      subroutine last_char(line,length,loc)
      character*(*) line

C  Given string LINE of length LENGTH locate and return in LOC
C  the position of the last character.

      loc=length
      do i=length,1,-1
        if(line(i:i).ne.' ') then
          loc=i
          return
        endif
      enddo
      return
      end
C-----------------------------------------------------------------------
      subroutine read_line(iunit,line)
      character line*130
      read(iunit,fmt='(a130)') line
      return
      end
C-----------------------------------------------------------------------
      subroutine search(iunit,line,text,ifound,idbug)
      character*(*) text
      character line*130

C Scan unit IUNIT for the next occurence of TEXT.
C If idbug > 0 write line to unit idbug (probably 6).
C If idbug = -1 then stop program if TEXT not found

      ifound=0
      l=len(text)
 10   read(iunit,fmt='(a130)',end=20) line
      if (idbug .gt. 0) then
        call write_files(idbug,0,line)
      endif
      do i=1,130-l+1
        if (line(i:i+l-1) .eq. text) then
          ifound=1
          return
        endif
      enddo
      goto 10
 20   if (idbug .eq. 1) then
        write(6,*) 'Cannot find ', text(1:l+1)
        call error
      endif
      return
      end
C-----------------------------------------------------------------------
      subroutine skip(iunit,nlines)
      character char*1

C Skip next NLINES on unit IUNIT

      do i=1,nlines
        read(iunit,fmt='(1a1)') char
      enddo
      return
      end
C-----------------------------------------------------------------------
      subroutine write_files(iout1,iout2,text)
      character*(*) text

C Writes TEXT to files iout1 and iout2 (usually output and temp files)
C or just output file (iout2=0)

      loc1=len(text)
      call last_char(text,loc1,loc)
      write(iout1,'(a)') text(1:loc)
      if (iout2 .eq. 0) goto 10
      write(iout2,'(a)') text(1:loc)

 10   return
      end
C----------------------------------------------------------------------

      subroutine up_case(string,lenstr)
c convert all lower case characters to upper case
C
      CHARACTER*(*) STRING
      CHARACTER*26 UCASE,LCASE
C
      DATA UCASE /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA LCASE /'abcdefghijklmnopqrstuvwxyz'/
C
C        CONVERTS LOWER CASE TO UPPER CASE IN THE GIVEN STRING.
C
      DO 100 I=1,LENSTR
         IC = INDEX(LCASE,STRING(I:I))
         IF (IC.GT.0) STRING(I:I) = UCASE(IC:IC)
  100 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      subroutine backs(iunit,line,text,ifound,idbug)
      character*(*) text
      character line*130
 
C Scan unit IUNIT backwards for the occurence of TEXT.
C If idbug > 0 write line to unit idbug (probably 6).
C If idbug = -1 then stop program if TEXT not found
 
      ifound=0
      l=len(text)
 10   read(iunit,fmt='(a130)',end=20) line
      if (idbug .gt. 0) then
        call write_files(idbug,0,line)
      endif
      do i=1,130-l+1
        if (line(i:i+l-1) .eq. text) then
          ifound=1
          return
        endif
      enddo
      backspace(iunit)
      backspace(iunit)
      goto 10
 20   if (idbug .eq. 1) then
        write(6,*) 'Cannot find ', text(1:l+1)
        call error
      endif
      return
      end
