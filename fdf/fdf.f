c----------------------------------------------------------
c
c FDF (Flexible Data Format) routine package
c
c Copyright Alberto Garcia, Jose Soler, 1996, 1997, 1998
c
c------------------
c
c     Notes: 
c
c     This package implements fdf Version 0.6 (see file fdf.Standard)
c
c     User callable routines:
c
c     fdf_init(filein,fileout)
c     fdf_inhibit
c     fdf_shutdown
c     fdf_enabled
c     fdf_{integer,single,double,string,boolean} (label,default)
c     fdf_physical(label,default,unit)
c     fdf_block(label,io_unit)
c     fdf_defined(label)
c
c     Implementation notes:
c
c     This version needs the io module for logical unit management.
c
c-------------------------------------------------------------------
c
      block data fdf_data
      implicit none
      include 'fdf.h'
      data fdf_started, fdf_debug, fdf_debug2, fdf_donothing
     $     /.false.,.false.,.false.,.false./
      end
c
      subroutine fdf_inhibit
      implicit none
      include 'fdf.h'
      fdf_donothing = .true.
      end

      subroutine fdf_init(filein,fileout)
c
c     New initialization for fdf. Simplified user interface using
c     the io package.
c
      implicit none

      character*(*) filein, fileout

      include 'fdf.h'

      integer debug_level
c
      integer fdf_integer
      external fdf_integer
c
      if (fdf_donothing) return
c
c     Prevent the user from opening two head files
c
      if (fdf_started) then 
         write(fdf_err,'(a)') 'FDF: Head file already set...'
         stop 'HEAD'
      endif

      call io_geterr(fdf_err)

      ndepth = 0

      call io_assign(fdf_out)
      open(unit=fdf_out,file=fileout,form='formatted',
     $                               status='unknown')
      rewind(fdf_out)

      call fdf_open(filein)
      write(fdf_out,'(/,a,a,a,i3,/)')
     $          '#FDF: Opened ',filein, ' for input. Unit:',fdf_in

      fdf_started = .true.

      debug_level = fdf_integer('fdf-debug',0)
      call fdf_setdebug(debug_level)

      return
      end
c---------------------------------------------------
      subroutine fdf_shutdown
c
c     Closes the 'head' file
c
      implicit none

      include 'fdf.h'

      if (.not. fdf_started) return

      call fdf_refresh
      call io_close(fdf_in)
      close(fdf_out)
      fdf_started = .false.

      return
      end
c--------------------------------------------------
      logical function fdf_enabled()
      implicit none
      include 'fdf.h'
      fdf_enabled = fdf_started
      return
      end
c---------------------------------------------------
      subroutine fdf_setdebug(level)
c
c     Debugging levels: 
c     level <=0: nothing
c     level  =1: standard
c     level >=2: exhaustive
c
      implicit none

      integer level
      include 'fdf.h'
c
      if (level .le. 0) then

         if (fdf_debug) then
            call io_close(fdf_log)
            fdf_debug = .false.
         endif

      else

         if (.not. fdf_debug) then
            call io_assign(fdf_log)
            open(fdf_log,file='FDF.debug',form='formatted',
     $           status='unknown')
            rewind(fdf_log)
            fdf_debug = .true.
         endif
      endif
      
      fdf_debug2 = (level .ge. 2)

      return
      end
c---------------------------------------------------
c
      subroutine fdf_open(filename)
      implicit none
c
c     Opens a file for fdf processing.
c
      character*(*) filename

      include 'fdf.h'

      integer lun
      logical file_exists

      logical leqi
      external leqi
c
      ndepth = ndepth + 1
      if (ndepth .gt. maxdepth) then
         write(fdf_err,'(a)') 'FDF: Too many nested fdf files...'
         stop 'DEPTH'
      endif

      if (leqi(filename,'stdin')) then
         lun = 5
         if (fdf_debug) write(fdf_log,'(a,i1,a)')
     $        '--->Reading from Standard Input [depth:', ndepth,'] '

      else

         call io_assign(lun)

         inquire(file=filename,exist=file_exists)
         if (file_exists) then
            open(unit=lun,file=filename,status='old',form='formatted')
            rewind(lun)
            if (fdf_debug) write(fdf_log,'(a,i1,a,a50)')
     $           '--->Opened [depth:', ndepth,'] ', filename
         else
            write(fdf_err,'(a,a60)')
     $           'FDF: Cannot open ',filename
         endif
      endif

      fdf_stack(ndepth) = lun
      fdf_in = lun
      
      return
      end
c-----------------------------------
      subroutine fdf_close
      implicit none
c
c     Closes currently opened fdf file, except if it is the original one.
c
      include 'fdf.h'

      if (ndepth .gt. 1) then
         call io_close(fdf_in)
         if (fdf_debug) write(fdf_log,'(a,i1,a)')
     $        '--->Closed [depth:', ndepth,']'

         ndepth = ndepth -1
         fdf_in = fdf_stack(ndepth)
      endif

      return
      end
c-------------------------------------
      subroutine fdf_refresh
c
c     Closes all the open files in the stack (except the first).
c     Failure to do so would imply that the next Label is searched 
c     first in the 'deeper' files. fdf_locate calls fdf_refresh 
c     before doing anything. 
c
      implicit none
      
      include 'fdf.h'

      integer i

      do i = ndepth, 1 , -1
         call fdf_close
      enddo
      
      return
      end
c
      subroutine fdf_parse
c
c     Processes the input line looking for meaningful tokens.
c
      implicit none

      include 'fdf.h'
c
      logical intoken, instring

      integer c
      integer stringdel
c
c     Character statement functions
c
      integer i
      logical isdigit, isupper, islower, isalpha,
     $        isalnum, isextra, istokch
      logical iscomment, isdelstr, isspecial
c
      isdigit(i) = (i .ge. 48) .and. (i .le. 57)
      isupper(i) = (i .ge. 65) .and. (i .le. 90)
      islower(i) = (i .ge. 97) .and. (i .le. 122)
      isalpha(i) = isupper(i) .or. islower(i)
      isalnum(i) = isdigit(i) .or. isalpha(i)

c     Extra characters allowed in tokens:  $ % * + & - . / @ ^ _ ~
      isextra(i) = ((i .ge. 36) .and. (i .le. 38))
     $           .or. (i .eq. 42) .or. (i .eq. 43)
     $           .or. (i .eq. 45) .or. (i .eq. 46)
     $           .or. (i .eq. 47) .or. (i .eq. 64) .or. (i .eq. 94)
     $           .or. (i .eq. 95) .or. (i .eq. 126)

      istokch(i) = isalnum(i) .or. isextra(i)
c
c     Comments are signaled by:  !  #  ; 
      iscomment(i) = (i.eq.33) .or. (i.eq.35) .or. (i.eq.59)
c
c     String delimiters: "  '  `
      isdelstr(i) = (i.eq.34) .or. (i.eq.39) .or. (i.eq.96)
c
c     Special characters which are tokens by themselves: <
      isspecial(i) = (i.eq.60)
c
c========================================================
c
      intoken = .false.
      instring = .false.
      ntokens = 0
      stringdel = 0
      
      do i = 1, len(line)
         c = ichar(line(i:i))

         if (iscomment(c)) then
c possible comment...
            if (instring) then
               last(ntokens) = i
            else
               goto 1000
            endif

         else if (istokch(c)) then
c character allowed in a token...
            if (.not. intoken) then
               intoken = .true.
               ntokens = ntokens+1
               first(ntokens) = i
            endif
            last(ntokens) = i

         else if (isspecial(c)) then
c character that forms a token by itself...
            if (.not. instring) then
               ntokens=ntokens+1
               first(ntokens) = i
               intoken = .false.
            endif
            last(ntokens) = i

         else if (isdelstr(c)) then
c string delimiter... make sure it is the right one before
c closing the string.
c If we are currently in a token, the delimiter is appended to it.

            if (instring) then
               if (c.eq.stringdel) then
                  instring = .false.
                  intoken = .false.
                  stringdel = 0
               else
                  last(ntokens) = i
               endif
            else
               if (intoken) then
                  last(ntokens) = i
               else
                  instring = .true.
                  stringdel = c
                  intoken = .true.
                  ntokens = ntokens+1
                  first(ntokens) = i+1
                  last(ntokens) = i+1
               endif
            endif

         else
c token delimiter...

            if (instring) then
               last(ntokens) = i
            else
               if (intoken) intoken=.false.
            endif
         endif

      enddo

 1000 continue

      if (fdf_debug2) then
         write(fdf_log,*) '            ',  ntokens, ' tokens:'
         do i=1,ntokens
            write(fdf_log,*) '                 ',
     $           '|',line(first(i):last(i)),'|'
         enddo
      endif

      return
      end
c------------------------------------------------------------------
c
      integer function fdf_search(label)
c
c     Performs a case-and-punctuation-insensitive search for 'label'
c     among the tokens in a line.
c
      implicit none

      character*(*) label
      include 'fdf.h'
      integer i
      logical labeleq
      external labeleq

      fdf_search = 0
      do i = 1, ntokens
         if (labeleq(label,line(first(i):last(i)))) then
            fdf_search = i
            return
         endif
      enddo

      return
      end
c-----------------------------------------------------------------
c
c
c====================================================================
c
      character *(*) function fdf_string(label,default)
c
c     Returns a string associated with label label, or default if label
c     is not found in the fdf file.
c
      implicit none

      character*(*) label
      character*(*) default

      include 'fdf.h'

      logical fdf_locate
      external fdf_locate
c
      fdf_string = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,a,5x,a)') label, default,
     $        '# Default value'
         return
      endif
c
c     From the second token up...
c
      fdf_string = line(first(2):last(ntokens))
      write(fdf_out,'(a,5x,a)') label, fdf_string

      return

      end
c
      logical function fdf_boolean(label,default)
c
c     Returns true if label appears by itself or in the form
c     label {Yes,true,.true.,T} (case insensitive).
c
c     Returns false if label appears in the form
c     label {No,false,.false.,F} (case insensitive).
c
c     If label is not found in the fdf file, fdf_boolean returns the 
c     logical variable default.
c
      implicit none

      character*(*) label
      logical default

      character valstr*40

      include 'fdf.h'

      logical fdf_locate, leqi
      external fdf_locate, leqi
c
      fdf_boolean = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,l10,5x,a)')
     $                  label, default, '# Default value'
         return
      endif

c
c     If the label appears by itself, we interpret it as .true.
c
      if (ntokens .eq. 1) then
         fdf_boolean = .true.
         write(fdf_out,'(a,5x,l10,5x,a)') label, fdf_boolean,
     $                                    '# Label by itself'
         return
      endif
c
c     Look for second word
c
      valstr=line(first(2):last(2))
c
      if (leqi(valstr,'yes') .or.
     $    leqi(valstr,'true') .or.
     $    leqi(valstr,'.true.') .or.
     $    leqi(valstr,'t') .or.
     $    leqi(valstr,'y'))       then
         
         fdf_boolean = .true.
         write(fdf_out,'(a,5x,l10)') label, fdf_boolean

      else if (leqi(valstr,'no') .or.
     $         leqi(valstr,'false') .or.
     $         leqi(valstr,'.false.') .or.
     $         leqi(valstr,'f') .or.
     $         leqi(valstr,'n'))       then

         fdf_boolean = .false.
         write(fdf_out,'(a,5x,l10)') label, fdf_boolean

      else

         write(fdf_err,*)
     $        'FDF_BOOLEAN: Unexpected fdf logical value ',
     $        label, ' = ', valstr
         stop 

      endif

      return
      end
c
c----------------------------------------------------------------------
c
      integer function fdf_integer(label,default)
c
c     Returns an integer associated with label, or default if label
c     is not found in the fdf file.
c
      implicit none

      character*(*) label
      integer default

      include 'fdf.h'

      logical fdf_locate
      external fdf_locate
c
      character*10 fmtstr
c
      fdf_integer = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,i10,5x,a)')
     $                   label, default, '# Default value'
         return
      endif

      if (ntokens.eq.1) then
         write(fdf_err,*) 'FDF_INTEGER: No value for ', label
         stop
      endif

      write(fmtstr,9000) last(2)-first(2)+1
 9000 format('(i',i2.2,')')
      read(line(first(2):last(2)),fmt=fmtstr) fdf_integer
      write(fdf_out,'(a,5x,i20)') label, fdf_integer

      return

      end
c
      real function fdf_single(label,default)
c
c     Returns a single precision value associated with label label, 
c     or default if label is not found in the fdf file.
c
      implicit none

      character*(*) label
      real default

      include 'fdf.h'

      logical fdf_locate
      external fdf_locate
c
      character*10 fmtstr
c
      fdf_single = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,g20.10,5x,a)')
     $                   label, default, '# Default value'
         return
      endif

      if (ntokens.eq.1) then
         write(fdf_err,*) 'FDF_SINGLE: No value for ', label
         stop
      endif
      write(fmtstr,9000) last(2)-first(2)+1
 9000 format('(g',i2.2,'.0)')
      read(line(first(2):last(2)),fmt=fmtstr) fdf_single
      write(fdf_out,'(a,5x,g20.10)') label, fdf_single

      return

      end
c
      double precision function fdf_double(label,default)
c
c     Returns a double precision value associated with label label, 
c     or default if label is not found in the fdf file.
c
      implicit none

      character*(*) label
      double precision default

      include 'fdf.h'

      logical fdf_locate
      external fdf_locate
c
      character*10 fmtstr
c
      fdf_double = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,g20.10,5x,a)')
     $                   label, default, '# Default value'
         return
      endif

      if (ntokens.eq.1) then
         write(fdf_err,*) 'FDF_DOUBLE: No value for ', label
         stop
      endif
      write(fmtstr,9000) last(2)-first(2)+1
 9000 format('(g',i2.2,'.0)')
      read(line(first(2):last(2)),fmt=fmtstr) fdf_double
      write(fdf_out,'(a,5x,g20.10)') label, fdf_double

      return

      end
c
c------------------------------------------------------
      double precision function fdf_physical(label,default,defunit)
c
c     Returns a double precision value associated with label label, 
c     or default if label is not found in the fdf file. Converts
c     the units to defunit.
c
      implicit none

      character*(*) label, defunit
      double precision default

      character unitstr*10
      double precision value

      include 'fdf.h'

      logical fdf_locate, leqi
      double precision fdf_convfac
      external fdf_locate, fdf_convfac, leqi
c
      character*10 fmtstr
c
      fdf_physical = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,g20.10,1x,a,5x,a)')
     $                   label, default, defunit, '# Default value'
         return
      endif

      if (ntokens.eq.1) then
         write(fdf_err,*) 'FDF_PHYSICAL: No value for ', label
         stop
      endif
      write(fmtstr,9000) last(2)-first(2)+1
 9000 format('(g',i2.2,'.0)')
      read(line(first(2):last(2)),fmt=fmtstr) value
      fdf_physical = value
c
c     Look for unit
c
      if (ntokens.eq.2) then
         write(fdf_err,*) 'FDF_PHYSICAL: No unit specified for ', label
         stop
      endif

      unitstr=line(first(3):last(3))
      if (.not. leqi(unitstr,defunit))
     $     fdf_physical = value * fdf_convfac(unitstr,defunit)
      write(fdf_out,'(a,5x,g20.10,1x,a10)')
     $               label, fdf_physical, defunit
      write(fdf_out,'(a,a,5x,g20.10,1x,a10)')
     $               '# Above item originally: ', label, value, unitstr

      return

      end
c
c-----------------------------------------------------------------------
c
      logical function fdf_block(label,unit)
c
c     Returns "true" and the unit number of the file from which to read
c     the contents of a block if "label" is associated with a block, and
c     false if not (unit is set to -1 in this case).
c
      implicit none
      character*(*) label
      integer unit
      include 'fdf.h'
      
      character*50 token1, filename
      integer iless

      logical fdf_locate, leqi
      integer fdf_search

      fdf_block = .false.
      unit = -1
      if (.not. fdf_locate(label)) return

      token1 = line(first(1):last(1))
      if (.not. leqi(token1,'%block')) then
         write(fdf_err,*) 'FDF_BLOCK: Not a block:',label
c
c        Return instead of stopping
c
         return
      endif

      iless = fdf_search('<')
      if ((iless .ne. 0) .and. (ntokens .gt. iless)) then
c
c           Read block from file
c
         filename = line(first(iless+1):last(iless+1))
         if (fdf_debug) write(fdf_log,'(2a)')
     $        '*Reading block from file', filename
         call fdf_open(filename)
         if (fdf_search('%dump') .ne. 0) call fdf_dumpfile(label)
         fdf_block = .true.
         unit = fdf_in
         return
      endif
c
c     Standard block in fdf file. Dump contents
c
      call fdf_dumpblock(label)
      fdf_block = .true.
      unit = fdf_in
      
      return
      end
c
c
      logical function fdf_block_old(label,unit)
c
c     Returns "true" and the unit number of the file from which to read
c     the contents of a block if "label" is associated with a block, and
c     false if not (unit is set to -1 in this case).
c
c     Wrapper function to get around scope bug in pgf90

      implicit none
      character*(*) label
      integer unit
      logical fdf_block
      external fdf_block

      fdf_block_old = fdf_block(label,unit)
      end
c-----------------------------------------------------------------------
c
      logical function fdf_locate(label)
c
c     Searches for label in the fdf hierarchy. If it appears and it
c     is not part of a comment, the function returns .true. and leaves
c     the file positioned at the next line. Otherwise, it returns .false.
c
c     It supports two kinds of "include" files:
c
c     %include filename  
c     Indicates an unconditional opening of filename for 
c     further fdf processing.
c
c     Label1 Label2 ... < filename  
c     Indicates that filename should be opened only when 
c     searching for any of the labels indicated.
c     'filename' should be an fdf file.
c
      implicit none
      
      character*(*) label

      character*60 token1, filename
      integer ilabel, iless
c
      include 'fdf.h'

      integer fdf_search
      logical leqi, fdf_getline
      external fdf_search, leqi, fdf_getline
c
      fdf_locate = .false.
      if (fdf_donothing) return
c
      call fdf_refresh
      if (fdf_debug) write(fdf_log,'(/,a,1x,a)')
     $        'Looking for ', label

      rewind(fdf_in)

 10   continue

      if (.not. fdf_getline()) then
         if (ndepth .gt. 1) then
            call fdf_close
            goto 10
         endif
         if (fdf_debug) write(fdf_log,'(a,1x,a)')
     $        '*Did not find ', label
         return
      endif
c
      if (ntokens .eq. 0) goto 10
c
      token1 = line(first(1):last(1))
c
      if (leqi(token1,'%include')) then
c
c        Include file
c
         if (ntokens .eq. 1) then
            write(fdf_err,*) 'FDF: No valid filename after %include'
            stop
         endif
         filename = line(first(2):last(2))
         call fdf_open(filename)
         goto 10
      endif

      ilabel = fdf_search(label)
      if (ilabel .ne. 0) then
c
c        Label found...
c
         if (leqi(token1,'%block')) then
            fdf_locate = .true.
            if (fdf_debug) write(fdf_log,'(a,1x,a)')
     $     '*Found ', label
            return
         endif

         iless = fdf_search('<')
         if ((iless .ne. 0) .and. (ntokens .gt. iless)) then
c
c           Continue search in other file
c
            filename = line(first(iless+1):last(iless+1))
            call fdf_open(filename)
            goto 10
         endif
c
c        If we reach this point we must be dealing with a line
c        of the form 'Label Value'. But we are not interested if
c        the string appears in the "Value" section
c
         if (ilabel .eq. 1) then
            fdf_locate = .true.
            if (fdf_debug) write(fdf_log,'(a,1x,a)') '*Found ', label
            return
         else
            goto 10
         endif

      else

         goto 10

      endif

      end
c
c---------------------------------------
c     
      logical function fdf_defined(label)
c
      implicit none
      
      character*(*) label
      include 'fdf.h'

      logical fdf_locate
      external fdf_locate
      
      fdf_defined = fdf_locate(label)
      if (fdf_defined) write(fdf_out,'(a)') label

      return
      end
c
      logical function fdf_getline()
      implicit none

      include 'fdf.h'
      
      read(fdf_in,end=100,err=100,fmt='(a)') line
      fdf_getline = .true.
      if (fdf_debug2)
     $     write(fdf_log,'(a,a76)') '> ', line
      call fdf_parse
      return

 100  continue
      fdf_getline = .false.
      return
      end

c
c------------------------------------------------------
c
      DOUBLE PRECISION FUNCTION fdf_convfac( FROM, TO )

C Returns conversion factor between a subset of physical units
C Written by J.M.Soler. Dec'96.
c Modified by Alberto Garcia, Jan'97.

      IMPLICIT      NONE
      CHARACTER*(*) FROM, TO
      INTEGER       IU, IFROM, ITO, NU

      PARAMETER ( NU = 54 )
      CHARACTER         DIM(NU)*10, NAME(NU)*10
      DOUBLE PRECISION  UNIT(NU)
c
c
c     We allow case variations in the units. This could be dangerous
c     (meV --> MeV !!) in real life, but not in this restricted 
c     field.
c     

      logical leqi
      external leqi

      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=1,10) /
     .  'mass  ', 'Kg      ', 1.D0,
     .  'mass  ', 'g       ', 1.D-3,
     .  'mass  ', 'amu     ', 1.66054D-27,
     .  'length', 'm       ', 1.D0,
     .  'length', 'nm      ', 1.D-9,
     .  'length', 'Ang     ', 1.D-10,
     .  'length', 'Bohr    ', 0.529177D-10,
     .  'time  ', 's       ', 1.D0,
     .  'time  ', 'fs      ', 1.D-15,
     .  'energy', 'J       ', 1.D0/
      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=11,20) /
     .  'energy', 'erg     ', 1.D-7,
     .  'energy', 'eV      ', 1.60219D-19,
     .  'energy', 'meV     ', 1.60219D-22,
     .  'energy', 'Ry      ', 2.17991D-18,
     .  'energy', 'mRy     ', 2.17991D-21,
     .  'energy', 'Hartree ', 4.35982D-18,
     .  'energy', 'K       ', 1.38066D-23,
     .  'energy', 'kcal/mol', 6.94780D-21,
     .  'force ', 'N       ', 1.D0,
     .  'force ', 'eV/Ang  ', 1.60219D-9/
      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=21,30) /
     .  'force ', 'Ry/Bohr ', 4.11943D-8,
     .  'length  ', 'cm      ', 1.d-2,
     .  'time    ', 'ps      ', 1.d-12,
     .  'time    ', 'ns      ', 1.d-9,
     .  'energy  ', 'mHartree', 4.35982D-21,
     .  'energy  ', 'kJ/mol  ', 1.6606d-21,
     .  'energy  ', 'Hz      ', 6.6262d-34,
     .  'energy  ', 'THz     ', 6.6262d-22,
     .  'energy  ', 'cm-1    ', 1.986d-23,
     .  'energy  ', 'cm^-1   ', 1.986d-23/
      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=31,40) /
     .  'pressure', 'Pa      ', 1.d0,
     .  'pressure', 'MPa     ', 1.d6,
     .  'pressure', 'GPa     ', 1.d9,
     .  'pressure', 'atm     ', 1.01325d5,
     .  'pressure', 'bar     ', 1.d5,
     .  'pressure', 'Mbar    ', 1.d11,
     .  'charge  ', 'C       ', 1.d0,
     .  'charge  ', 'e       ', 1.602177d-19,
     .  'dipole  ', 'C*m     ', 1.d0,
     .  'dipole  ', 'D       ', 3.33564d-30/
      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=41,50) /
     .  'dipole  ', 'debye   ', 3.33564d-30,
     .  'dipole  ', 'e*Bohr  ', 8.47835d-30,
     .  'dipole  ', 'e*Ang   ', 1.602177d-29,
     .  'energy  ', 'cm**-1    ', 1.986d-23,
     .  'pressure', 'Ry/Bohr**3', 1.47108d13,
     .  'pressure', 'eV/Ang**3 ', 1.60219d11,
     .  'MomInert', 'Kg*m**2   ', 1.d0,
     .  'MomInert', 'Ry*fs**2  ', 2.17991d-48,
     .  'Efield  ', 'V/m       ', 1.d0,
     .  'Efield  ', 'V/nm      ', 1.d9 /
      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=51,54) /
     .  'Efield  ', 'V/Ang     ', 1.d10,
     .  'Efield  ', 'V/Bohr    ', 1.8897268d10,
     .  'Efield  ', 'Ry/Bohr/e ', 2.5711273d11,
     .  'Efield  ', 'Har/Bohr/e', 5.1422546d11 /
c
      IFROM = 0
      ITO   = 0
      DO 10 IU = 1,NU
        IF (leqi(NAME(IU),FROM)) IFROM = IU
        IF (leqi(NAME(IU),TO))   ITO   = IU
   10 CONTINUE

      IF (IFROM .EQ. 0) THEN
        WRITE(6,*) 'FDF_CONVFAC: Unknown unit = ', FROM
        STOP
      ENDIF
      IF (ITO .EQ. 0) THEN
        WRITE(6,*) 'FDF_CONVFAC: Unknown unit = ', TO
        STOP
      ENDIF

      IF (leqi(DIM(IFROM),DIM(ITO))) THEN
        FDF_CONVFAC = UNIT(IFROM) / UNIT(ITO)
      ELSE
        WRITE(6,*)
     $    'FDF_CONVFAC: Unit''s physical dimensions don''t match: ',
     .     FROM, ', ', TO
        STOP
      ENDIF
      END
c
c-------------------------------------------------------------------
c
      subroutine fdf_dumpblock(label)
c     
c     Dumps block contents starting at the current line
c     
      implicit none

      character*(*) label

      include 'fdf.h'

      integer i, lblock
      character*60 token1

      logical fdf_getline, leqi
      integer fdf_search

      write(fdf_out,'(/,a79)') line
      lblock = 0
 120  continue
      if (fdf_getline()) then
         lblock = lblock + 1
         write(fdf_out,'(a79)') line
         token1 = line(first(1):last(1))
         if (.not. leqi(token1,'%endblock')) goto 120
      else
         write(fdf_err,'(a,a,a)')
     $        'FDF_LOCATE: Block ', label, ' does not end!'
         stop 'FDF'
      endif
      write(fdf_out,*)
c     
c     Sanity check (optional construct %endblock [ Label [Label] ])
c     
      if ((ntokens .gt. 1) .and.
     $     (fdf_search(label) .eq. 0)) then
         write(fdf_err,'(a,a,a)')
     $        'FDF_LOCATE: Block ', label, ' does not end!'
         stop 'FDF'
      endif
c
c     Backspace the lines read
c
      do i=1,lblock
         backspace(fdf_in)
      enddo

      return
      end
c     
c--------------------------------------------------------------
c     
      subroutine fdf_dumpfile(label)
c     
c     Dumps the contents of a file to fdf_out.
c     The lines are embedded in a %block ... %endblock pair.
c     
      implicit none

      include 'fdf.h'

      character*(*) label
      character form*30
      integer length

      logical fdf_getline
      external fdf_getline, chrlen
c     
c     Build the right format
c     
      call chrlen(label,0,length)
      write(form,'(a,i2.2,a)') '(a,a',length,',10x,a)'

      write(fdf_out,*)
      write(fdf_out,form) '%block ', label,
     $     '# Originally in include file' 
c     
      rewind(fdf_in)
 10   continue
      if (fdf_getline()) then
         write(fdf_out,'(a79)') line
         goto 10
      endif

      write(fdf_out,form) '%endblock ', label,
     $     '# Originally in include file' 
      write(fdf_out,*)
      rewind(fdf_in)
c     
      return
      end
c
c----------------------------------------------------------------
c
      logical function labeleq(s1,s2)
c
c     Compares s1 and s2 without regard for case, or appearance
c     of '_', '.', '-'.
c
      implicit none

      include 'fdf.h'

      character*(*) s1, s2
      character*80 n1, n2
      logical leqi

      call fdf_pack(s1,n1)
      call fdf_pack(s2,n2)
      labeleq=leqi(n1,n2)
      if (fdf_debug) then
         if (labeleq .and. .not. leqi(s1,s2))
     $        write(fdf_log,'(a,/,a,/,a)')
     $        '--------- Considered equivalent:', s1, s2
      endif
      return
      end
c-----------------------------
      subroutine fdf_pack(s,n)
      implicit none
      character*(*) s, n
c
c     Removes occurrences of '_ .-'  from s1
c
      character*1 c
      integer i, j
      logical issep
      issep(i) = (i.eq.95) .or. (i.eq.46) .or. (i.eq.45)

      n = ' '
      j = 0
      do i = 1, len(s)
         c = s(i:i)
         if (.not.issep(ichar(c))) then
            j = j+1
            n(j:j) = c
         endif
      enddo
      return
      end
c-------------
      SUBROUTINE CHRLEN(STRING,NCHAR,LCHAR)
C
C  CHRLEN accepts a STRING of NCHAR characters and returns LCHAR,
C  the length of the string up to the last nonblank, nonnull.
C     
      implicit none

      CHARACTER CHAR*1
      CHARACTER STRING*(*)
      integer nchar,lchar
C
      integer ncopy, i
c
      NCOPY=NCHAR
      IF(NCOPY.LE.0)NCOPY=LEN(STRING)
C
      DO 10 I=1,NCOPY
        LCHAR=NCOPY+1-I
        IF(STRING(LCHAR:LCHAR).NE.' '.AND.
     *     STRING(LCHAR:LCHAR).NE.CHAR(0))RETURN
 10     CONTINUE
      LCHAR=0
      RETURN
      END
c
      SUBROUTINE CHRCAP(STRING,NCHAR)
C
C  CHRCAP accepts a STRING of NCHAR characters and replaces
C  any lowercase letters by uppercase ones.
C
      implicit none

      CHARACTER CHAR*1
      integer nchar, ncopy, i, itemp
      LOGICAL   LGE
      LOGICAL   LLE
      CHARACTER STRING*(*)
C
      NCOPY=NCHAR
      IF(NCOPY.LE.0)NCOPY=LEN(STRING)
      DO 10 I=1,NCOPY
C
        IF(LGE(STRING(I:I),'a').AND.LLE(STRING(I:I),'z'))THEN
          ITEMP=ICHAR(STRING(I:I))+ICHAR('A')-ICHAR('a')
          STRING(I:I)=CHAR(ITEMP)
          ENDIF
10      CONTINUE
      RETURN
      END
c
      LOGICAL FUNCTION LEQI(STRNG1,STRNG2)
C
C  Case-insensitive lexical equal-to comparison
C
      implicit none
c
      CHARACTER*1   S1,S2
      CHARACTER*(*) STRNG1
      CHARACTER*(*) STRNG2
C
      integer len1, len2, lenc, i
      LEN1=LEN(STRNG1)
      LEN2=LEN(STRNG2)
      LENC=MIN(LEN1,LEN2)
C
      LEQI=.FALSE.
      DO 10 I=1,LENC
        S1=STRNG1(I:I)
        S2=STRNG2(I:I)
        CALL CHRCAP(S1,1)
        CALL CHRCAP(S2,1)
        IF(S1.NE.S2)RETURN
10      CONTINUE
C 
      IF(LEN1.GT.LENC.AND.STRNG1(LENC+1:LEN1).NE.' ')RETURN
      IF(LEN2.GT.LENC.AND.STRNG2(LENC+1:LEN2).NE.' ')RETURN
      LEQI=.TRUE.
      RETURN
      END
c


