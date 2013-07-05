! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module input_module
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/input_module
!!  NAME
!!   input_module
!!  PURPOSE
!!   This implements some of the functionality of FDF but operating on an character array
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2008/08/20
!!  MODIFICATION HISTORY
!!   2008/09/03 10:43 dave
!!    Possible bug fix: add datatypes as module wide for typing of fdf_double
!!  SOURCE
!!  
module input_module

  use datatypes
  use GenComms, ONLY: inode, ionode

  implicit none

  save
  ! Input array
  integer :: input_lines
  character(len=132), dimension(:), allocatable :: input_array

  ! Parsing
  integer :: current_line, ntokens
  character(len=132) :: line
  integer, dimension(100) :: first, last

  ! Control
  logical :: fdf_debug = .false.
  logical :: inblock = .false.
  integer :: fdf_log, fdf_out
  integer :: block_start, block_end

  ! Files
  integer, parameter :: lun_min = 10
  integer, parameter :: lun_max = 99
  integer, parameter :: nunits = 90!lun_max-lun_min+1
  logical, dimension(lun_min:lun_max) :: free_lun != (/90*.true./)!nunits*.true./)
  data free_lun /nunits*.true./

  character(len=80), private :: RCSid = "$Id: initial_read.module.f90 64 2008-08-07 07:50:31Z astorralba $"
!!***

contains

! ------------------------------------------------------------------------------
! Subroutine 
! ------------------------------------------------------------------------------

!!****f* input_module/load_input *
!!
!!  NAME 
!!   load_input
!!  USAGE
!!   
!!  PURPOSE
!!   Loads the input file into a character array and broadcasts it
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2008/08/20
!!  MODIFICATION HISTORY
!!   2013/07/01 M.Arita
!!    - Bug fix in closing a file
!!  SOURCE
!!  
  subroutine load_input

    use GenComms, ONLY: inode, ionode, gcopy, cq_abort

    implicit none

    integer :: lun, stat, i, l,j
    character(len=132) :: line
    character(len=10) :: slabel
    logical :: good_line, done
    real(double) :: r

    ! Count lines in input file on ionode
    if(inode==ionode) then
       call io_assign(lun)
       open(unit=lun,file='Conquest_input',iostat=stat,status='old')
       if(stat/=0) call cq_abort("We need Conquest_input to run !")
       input_lines = 0
       done = .false.
       do while(.NOT.done)
          ! Need error/eof added here !
          read(lun,fmt='(a)',iostat=stat) line
          if(stat<0) then
             done = .true.
             exit
          end if
          ! Test for comment or blank lines
          good_line = .false.
          do i=1,len(line)
             if(line(i:i)==' ') cycle ! Remove leading blanks
             if(line(i:i)=='#') then ! Comment line
                exit
             else ! Do we want an al-num test here ?
                good_line = .true.
             end if
          end do
          if(good_line) input_lines = input_lines+1
       end do
       !OLD close(lun)
       call io_close(lun)       !01/07/2013 michi
    end if
    ! Broadcast size of array
    call gcopy(input_lines)
    ! Allocate array
    allocate(input_array(input_lines))
    if(inode==ionode) then
       open(unit=lun,file='Conquest_input',iostat=stat,status='old')
       l = 1
       done = .false.
       do while(.NOT.done)
          read(lun,fmt='(a)',iostat=stat) line
          if(stat<0) then
             done = .true.
             exit
          end if
          ! Test for comment or blank lines
          good_line = .false.
          do i=1,len(line)
             if(line(i:i)==' ') cycle ! Remove leading blanks
             if(line(i:i)=='#') then ! Comment line
                exit
             else ! Do we want an al-num test here ?
                good_line = .true.
             end if
          end do
          if(good_line) then
             if(l>input_lines) call cq_abort("Input reading error !",l,input_lines)
             input_array(l) = line
             l = l+1
          end if
       end do
       call io_close(lun)
    end if
    do i=1,input_lines
       call gcopy(input_array(i),132)
    end do
    current_line = 0
    if(inode==ionode) then
       call io_assign(fdf_out)
       open(unit=fdf_out,file='input.log')
       i =  fdf_integer('fdf-debug',0)
       if(i>0) then
          fdf_debug = .true.
          call io_assign(fdf_log)
          open(unit=fdf_log,file='input_debug.log')
       else
          fdf_debug = .false.
       end if
    else
       fdf_debug = .false.
    end if
    return
  end subroutine load_input
!!***
       
  integer function fdf_integer(label,default)
!
!     Returns an integer associated with label, or default if label
!     is not found in the fdf file.
!
    use GenComms, ONLY: cq_abort

    implicit none

    character(len=*) :: label
    integer :: default

    character(len=10) :: fmtstr

    fdf_integer = default

    if (.not. fdf_locate(label)) then
       if(inode==ionode) write(fdf_out,'(a,5x,i10,5x,a)') label, default, '# Default value'
       return
    endif

    if (ntokens==1) then
       if(inode==ionode) write(fdf_out,*) 'FDF_INTEGER: No value for ', label
       call cq_abort("Input problem: integer not found")
    endif

    write(fmtstr,fmt='("(i",i2.2,")")') last(2)-first(2)+1
    read(line(first(2):last(2)),fmt=fmtstr) fdf_integer
    if(inode==ionode) write(fdf_out,'(a,5x,i20)') label, fdf_integer

    return

  end function fdf_integer

  real(double) function fdf_double(label,default)
!
!     Returns an integer associated with label, or default if label
!     is not found in the fdf file.
!
    use GenComms, ONLY: cq_abort

    implicit none

    character(len=*) :: label
    real(double) :: default

    character(len=10) :: fmtstr

    fdf_double = default

    if (.not. fdf_locate(label)) then
       if(inode==ionode) write(fdf_out,'(a,5x,g20.10,5x,a)') label, default, '# Default value'
       return
    endif

    if (ntokens==1) then
       if(inode==ionode) write(fdf_out,*) 'FDF_DOUBLE: No value for ', label
       call cq_abort("Input problem: double not found")
    endif

    write(fmtstr,fmt='("(g",i2.2,".0)")') last(2)-first(2)+1
    read(line(first(2):last(2)),fmt=fmtstr) fdf_double
    if(inode==ionode) write(fdf_out,'(a,5x,g20.10)') label, fdf_double

    return

  end function fdf_double

  function fdf_string(n,label,default)
!
! Returns a string associated with label label, or default if label
!     is not found in the fdf file.
    implicit none

    integer :: n
    character(len=n) :: fdf_string
    character(len=*) :: label, default
    
    fdf_string = default
    
    if (.not. fdf_locate(label)) then
       if(inode==ionode) write(fdf_out,'(a,5x,a,5x,a)') label, default, '# Default value'
       return
    endif
    fdf_string = line(first(2):last(ntokens))
    if(inode==ionode) write(fdf_out,'(a,5x,a)') label, fdf_string
  end function fdf_string

  logical function fdf_boolean(label,default)

    !     Returns true if label appears by itself or in the form
    !     label {Yes,true,.true.,T} (case insensitive).
    !
    !     Returns false if label appears in the form
    !     label {No,false,.false.,F} (case insensitive).
    !
    !     If label is not found in the fdf file, fdf_boolean returns the 
    !     logical variable default.

    use GenComms, ONLY: cq_abort
    implicit none

    character(len=*) :: label
    logical :: default

    character(len=40) :: valstr

    fdf_boolean = default

    if (.not. fdf_locate(label)) then
       if(inode==ionode) write(fdf_out,'(a,5x,l10,5x,a)') label, default, '# Default value'
       return
    endif
    !     If the label appears by itself, we interpret it as .true.

    if (ntokens .eq. 1) then
       fdf_boolean = .true.
       if(inode==ionode) write(fdf_out,'(a,5x,l10,5x,a)') label, fdf_boolean, '# Label by itself'
       return
    endif

    !     Look for second word

    valstr=line(first(2):last(2))

    if (leqi(valstr,'yes') .or. &
         leqi(valstr,'true') .or. &
         leqi(valstr,'.true.') .or. &
         leqi(valstr,'t') .or. &
         leqi(valstr,'y'))       then

       fdf_boolean = .true.
       if(inode==ionode) write(fdf_out,'(a,5x,l10)') label, fdf_boolean

    else if (leqi(valstr,'no') .or. &
         leqi(valstr,'false') .or. &
         leqi(valstr,'.false.') .or. &
         leqi(valstr,'f') .or. &
         leqi(valstr,'n'))       then

       fdf_boolean = .false.
       if(inode==ionode) write(fdf_out,'(a,5x,l10)') label, fdf_boolean

    else
       call cq_abort("FDF_BOOLEAN: Unexpected fdf logical value "//label//" "//valstr)
    endif
    return
  end function fdf_boolean

  logical function fdf_block(label)

    use GenComms, ONLY: cq_abort
    implicit none

    character(len=*) :: label
    integer :: start, end

    character(len=50) :: token1
    logical :: done

    fdf_block = .false.
    if (.not. fdf_locate(label)) return

    token1 = line(first(1):last(1))
    if (.not. leqi(token1,'%block')) then
       if(inode==ionode) write(fdf_log,*) 'FDF_BLOCK: Not a block:',label
!        Return instead of stopping
       return
    endif
    block_start = current_line+1
    fdf_block = .true.
    done = .false.
    do while(.NOT.done)
       if(fdf_getline()) then
          token1 = line(first(1):last(1))
          if (leqi(token1,'%endblock')) done = .true.
       else
          call cq_abort("Block fails to end: "//label)
       end if
    end do
    block_end = current_line-1
    inblock = .true.
    return
  end function fdf_block

  subroutine fdf_endblock

    implicit none

    inblock = .false.
  end subroutine fdf_endblock

  logical function fdf_defined(label)

    implicit none

    character(len=*) :: label

    fdf_defined = fdf_locate(label)
    if(fdf_defined.AND.inode==ionode) write(fdf_out,'(a)') label

    return
  end function fdf_defined

  logical function fdf_locate(label)
!
!     Searches for label in the fdf hierarchy. If it appears and it
!     is not part of a comment, the function returns .true. and leaves
!     the file positioned at the next line. Otherwise, it returns .false.
!
!     It supports two kinds of "include" files:
!
!     %include filename  
!     Indicates an unconditional opening of filename for 
!     further fdf processing.
!
!     Label1 Label2 ... < filename  
!     Indicates that filename should be opened only when 
!     searching for any of the labels indicated.
!     'filename' should be an fdf file.
!
    implicit none

    character(len=*) :: label

    character(len=60) :: token1, filename
    integer :: ilabel, iless
    logical :: done

    done = .false.
    fdf_locate = .false.
    if(inblock) then
       current_line = block_start - 1
    else
       current_line = 0
    end if
    if (fdf_debug.AND.inode==ionode) write(fdf_log,'(/,a,1x,a)') 'Looking for ', label
    do while(.NOT.done) 
       if (.not. fdf_getline()) then ! get a line
          if (fdf_debug.AND.inode==ionode) write(fdf_log,'(a,1x,a)') '*Did not find ', label
          return
       endif
       if (ntokens .eq. 0) cycle
       token1 = line(first(1):last(1))

       ilabel = fdf_search(label)
       if (ilabel .ne. 0) then
          if (leqi(token1,'%block')) then ! We've found a block
             fdf_locate = .true.
             if (fdf_debug.AND.inode==ionode) write(fdf_log,'(a,1x,a)') '*Found ', label
             return
          endif
          !        If we reach this point we must be dealing with a line
          !        of the form 'Label Value'. But we are not interested if
          !        the string appears in the "Value" section
          if (ilabel .eq. 1) then
             fdf_locate = .true.
             if (fdf_debug.AND.inode==ionode) write(fdf_log,'(a,1x,a)') '*Found ', label
             return
          endif
       endif
    end do
  end function fdf_locate

  logical function fdf_getline()

    implicit none

    integer :: maxlines

    if(inblock) then
       maxlines = block_end
    else
       maxlines = input_lines
    end if
    ! point to or copy next line ?
    current_line = current_line+1
    if(current_line<=maxlines) then
       fdf_getline = .true.
       line = input_array(current_line)
       if (fdf_debug.AND.inode==ionode) write(fdf_log,'(a,a76)') '> ', line
       call fdf_parse
       return
    else
       fdf_getline = .false.
       return
    endif
  end function fdf_getline

  subroutine fdf_parse
    !     Processes the input line looking for meaningful tokens.
    implicit none

    logical intoken, instring

    integer c
    integer stringdel

    !     Character statement functions

    integer i
    logical isdigit, isupper, islower, isalpha, isalnum, isextra, istokch
    logical iscomment, isdelstr, isspecial

    isdigit(i) = (i .ge. 48) .and. (i .le. 57)
    isupper(i) = (i .ge. 65) .and. (i .le. 90)
    islower(i) = (i .ge. 97) .and. (i .le. 122)
    isalpha(i) = isupper(i) .or. islower(i)
    isalnum(i) = isdigit(i) .or. isalpha(i)

    !     Extra characters allowed in tokens:  $ % * + & - . / @ ^ _ ~
    isextra(i) = ((i .ge. 36) .and. (i .le. 38)) &
         .or. (i .eq. 42) .or. (i .eq. 43) &
         .or. (i .eq. 45) .or. (i .eq. 46) &
         .or. (i .eq. 47) .or. (i .eq. 64) .or. (i .eq. 94) &
         .or. (i .eq. 95) .or. (i .eq. 126)

    istokch(i) = isalnum(i) .or. isextra(i)

    !     Comments are signaled by:  !  #  ; 
    iscomment(i) = (i.eq.33) .or. (i.eq.35) .or. (i.eq.59)

    !     String delimiters: "  '  `
    isdelstr(i) = (i.eq.34) .or. (i.eq.39) .or. (i.eq.96)

    !     Special characters which are tokens by themselves: <
    isspecial(i) = (i.eq.60)

    !========================================================

    intoken = .false.
    instring = .false.
    ntokens = 0
    stringdel = 0

    do i = 1, len(line)
       c = ichar(line(i:i))
       if (iscomment(c)) then
          ! possible comment...
          if (instring) then
             last(ntokens) = i
          else
             exit !goto 1000
          endif
       else if (istokch(c)) then
          ! character allowed in a token...
          if (.not. intoken) then
             intoken = .true.
             ntokens = ntokens+1
             first(ntokens) = i
          endif
          last(ntokens) = i
       else if (isspecial(c)) then
          ! character that forms a token by itself...
          if (.not. instring) then
             ntokens=ntokens+1
             first(ntokens) = i
             intoken = .false.
          endif
          last(ntokens) = i
       else if (isdelstr(c)) then
          ! string delimiter... make sure it is the right one before closing the string.
          ! If we are currently in a token, the delimiter is appended to it.
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
          ! token delimiter...
          if (instring) then
             last(ntokens) = i
          else
             if (intoken) intoken=.false.
          endif
       endif
    enddo

    if (fdf_debug) then
       if(inode==ionode) write(fdf_log,*) '            ',  ntokens, ' tokens:'
       do i=1,ntokens
          if(inode==ionode) write(fdf_log,*) '                 ', &
               '|',line(first(i):last(i)),'|'
       enddo
    endif

    return
  end subroutine fdf_parse

  integer function fdf_search(label)

    !     Performs a case-and-punctuation-insensitive search for 'label'
    !     among the tokens in a line.

    implicit none

    character(len=*) :: label

    integer i

    fdf_search = 0
    do i = 1, ntokens
       if (labeleq(label,line(first(i):last(i)))) then
          fdf_search = i
          return
       endif
    enddo

    return
  end function fdf_search

  logical function labeleq(s1,s2)
    !
    !     Compares s1 and s2 without regard for case, or appearance
    !     of '_', '.', '-'.
    !
    implicit none

    character(len=*) :: s1, s2
    character(len=80) :: n1, n2
    !logical :: leqi

    call fdf_pack(s1,n1)
    call fdf_pack(s2,n2)
    labeleq=leqi(n1,n2)
    if (fdf_debug) then
       if (labeleq.and.(.not. leqi(s1,s2)).AND.&
            inode==ionode) write(fdf_log,'(a,/,a,/,a)') '--------- Considered equivalent:', s1, s2
    endif
    return
  end function labeleq

  subroutine fdf_pack(s,n)
    implicit none

    character(len=*) s, n
    !
    !     Removes occurrences of '_ .-'  from s1
    !
    character :: c
    integer :: i, j
    logical :: issep
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
  end subroutine fdf_pack

  logical function leqi(strng1,strng2)
!
!  case-insensitive lexical equal-to comparison
!
    implicit none

    character :: s1,s2
    character(len=*) :: strng1
    character(len=*) :: strng2

    integer :: len1, len2, lenc, i

    len1=len(strng1)
    len2=len(strng2)
    lenc=min(len1,len2)

    leqi=.false.
    do  i=1,lenc
       s1=strng1(i:i)
       s2=strng2(i:i)
       call chrcap(s1,1)
       call chrcap(s2,1)
       if(s1.ne.s2) return
    end do

    if(len1.gt.lenc.and.strng1(lenc+1:len1).ne.' ')return
    if(len2.gt.lenc.and.strng2(lenc+1:len2).ne.' ')return
    leqi=.true.
    return
  end function leqi

  subroutine chrcap(string,nchar)
    !
    !  CHRCAP accepts a STRING of NCHAR characters and replaces
    !  any lowercase letters by uppercase ones.
    !
    implicit none

    character :: char
    integer nchar, ncopy, i, itemp
    character(len=*) ::  string

    ncopy=nchar
    if(ncopy.le.0)ncopy=len(string)
    do i=1,ncopy
       if(lge(string(i:i),'a').and.lle(string(i:i),'z'))then
          itemp=ichar(string(i:i))+ichar('A')-ichar('a')
          string(i:i)=char(itemp)
       endif
    end do
    return
  end subroutine chrcap

  subroutine io_assign(lun)

    use GenComms, ONLY: cq_abort

    implicit none

    integer :: lun, iostat
    logical :: used

    do lun = lun_min, lun_max
       if(free_lun(lun)) then
          inquire(unit=lun,opened=used,iostat=iostat)
          if(iostat/=0) used = .true.
          free_lun(lun) = .false.
          if(.NOT.used) return
       end if
    end do
    call cq_abort("Error in io_assign: no free luns between ",lun_min,lun_max)
  end subroutine io_assign

  subroutine io_close(lun)

    implicit none

    integer :: lun
    
    close(lun)
    if(lun>=lun_min.AND.lun<=lun_max) free_lun(lun) = .true.
    return
  end subroutine io_close

end module input_module
