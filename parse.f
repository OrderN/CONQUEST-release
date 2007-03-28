      module parse
!
!     Copyright Alberto Garcia <wdpgaara@lg.ehu.es> (1999, 2000-)
!
!=====================================================================
!
! This module provides a simple yet powerful way to analyze the information
! in a string (such as an input line). 
! 
! Routine, 'digest' takes as input a string 'line' and returns a pointer
! to a derived type 'parsed_line':
! 
!       type parsed_line
!       private
!           integer                            ::  ntokens
!           character(len=max_length), pointer ::  token_arr(:) 
!           character(len=1), pointer          ::  id(:)
!       end type parsed_line
! 
! which holds a list of tokens and token tags. The parsing (splitting of the
! string into tokens) is done by a helper routine 'line_parse' which
! currently behaves according to the FDF standard. Each token is
! classified by helper routine 'line_morphol' and a token id assigned
! in the following way:
! 
! * Tokens that can be read as real numbers are assigned to class
! 'values' and given a token id 'v'. These are further classified as
! 'integers' (id 'i') or 'reals' (id 'r').
! * All other tokens are tagged as 'names' (id 'n').
! 
! The recommended usage follows the outline:
! 
!     use parse
!     character(len=?) line
!     type(parsed_line), pointer :: p
!     ...
!     p=>digest(line)
!     (extract information from p)
!     call destroy(p)
! 
! Note the pointer assignment and the explicit call to a destroyer
! routine that frees the storage associated to p.
! 
! The information is extracted by module procedures that fall into three
! classes: 
! 
! a) Enquiry functions: 'search' and 'match'
! 
! *  'search' determines whether a token in 'line' matches the given
!    string, optionally returning an index. The search is
!    case-insensitive by default, but this can be changed by supplying
!    an extra procedure argument 'eq_func' with interface:
! 
!       interface
!          function eq_func(s1,s2)
!          logical eq_func
!          character(len=*), intent(in) :: s1,s2
!          end function eq_func
!       end interface
! 
!    Example:  if (search(p,'Mary',ind=j)) ...
!    will return in 'ind' the index of the first token that matches
!    "Mary", or -1 if not found. The function itself will return
!    .true. or .false.
! 
!    This function can take an optional keyword 'after=' (see below).
! 
! *  'match' is probably the most powerful routine in the module. It
!    checks whether the token morphology of 'line' conforms to the
!    sequence of characters specified. For example,
! 
!    if (match(p,'nii')) ...
! 
!    returns .true. if 'line' contains at least three tokens and they are
!    a 'name' and two 'integers'. A 'v' is matched by both an 'integer'
!    and a 'real'.
!    This function can take an optional keyword 'after=' (see below).
! 
! b) Number functions: ntokens, nnames, nreals, nintegers, nvalues
! 
!    These functions return the number of tokens of each kind in 'line':
! 
!    number_of_energies = nreals(p)
! 
!    These functions can take an optional keyword 'after=' (see below).
! 
! c) Extraction functions: tokens, names, reals, integers, values
! 
!    These functions return a piece of data which corresponds to a token
!    of the specified kind with sequence number matching the index
!    provided. For example,
! 
!    nlevels = integers(p,2)
! 
!    assigns to variable 'nlevels' the second integer in 'line'.
!    Execution stops in the assignment cannot be made. The user should
!    call the corresponding 'number' routine to make sure there are
!    enough tokens of the given kind.
! 
!    These functions can take an optional keyword 'after=' (see below).
! 
! 
! By default, the routines in the module perform any indexing from the
! beginning of 'line', in such a way that the first token is assigned the
! index 1. It is possible to specify a given token as 'origin' by using
! the 'after=' optional keyword. For example:
! 
!     if (search(p,'P',index=jp)) then
!        if (match(p,'i',after=jp) npol = integers(p,1,after=jp)
!     endif			 
! 
! first checks whether 'P' is found in 'line'. If so, 'match' is used to
! check whether it is followed by at least an 'integer'. If so, its
! valued is assigned to variable 'npol'.
! 
! If the 'after=' optional keyword is used in routine 'search', the
! returned index is absolute, not relative. For example, to get the
! real number coming right after the first 'Q' which appears to the
! right of the 'P' found above:
! 
!     if (search(p,'Q',after=jp,index=jq)) then
!        if (match(p,'r',after=jq) energy = reals(p,1,after=jq)
!     endif			 
! 
!-----------------------------------------------------------------------
! Future enhancements:
! 
! * Configurable parsing routine.
! * Better exception handling.
! * Optional 'reverse=' keyword.
! * Removal of arbitrary length limits (max_length and maxntokens)
!========================================================================

      implicit none

      private
      public parsed_line, digest, destroy
      public match, search

      public nintegers, integers
      public nreals, reals
      public nvalues, values
      public nnames, names
      public ntokens, tokens

      interface destroy
        module procedure destroy_lp
      end interface

      interface search
        module procedure search_line
      end interface

c========================================================================
c
      integer, save         ::  line_log = 0
      logical, save         ::  line_debug = .false.

      integer, parameter    ::  max_length = 132

      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)

      type parsed_line
      ! use 'token_arr' instead of 'tokens' to work around a 
      ! bug in pgf90
      private
          integer                            ::  ntokens
          character(len=max_length), pointer ::  token_arr(:) 
          character(len=1), pointer          ::  id(:)
      end type parsed_line
          
      CONTAINS
!
      subroutine create(p)
      type(parsed_line), pointer     :: p

      if (associated(p)) call destroy_lp(p)
      allocate(p)
      nullify(p%token_arr,p%id)
      end subroutine create

      subroutine destroy_lp(p)
      type(parsed_line), pointer     :: p

      if (.not. associated(p)) return
      if (associated(p%token_arr)) deallocate(p%token_arr)
      if (associated(p%id)) deallocate(p%id)
      deallocate(p)
      end subroutine destroy_lp

!=====================================================================
!     Checks whether the morphology of the line or part of it
!     matches the 'signature' string str.
!
      function match(p,str,after)  !!! not implemented: ,reverse)
      logical match
      type(parsed_line), pointer :: p
      character(len=*), intent(in)  :: str
      integer, intent(in), optional :: after   ! 
!!!      logical, intent(in), optional :: reverse ! Start from end of line

      integer nids, i, shift
      character(len=1) c, id

      shift = 0
      if (present(after)) shift = after
      if (shift.lt.0) call die("Wrong starting position in match")

      match = .false.
      nids = len_trim(str)
      if (p%ntokens - shift .lt. nids) return
      do i = 1, nids
         c = str(i:i)
         id = p%id(shift+i)
         if (leqi(c,'v') .and. 
     $               .not. (leqi(id,'i') .or. leqi(id,'r'))) return
         if (leqi(c,'i') .and. .not. leqi(id,'i')) return
         if (leqi(c,'r') .and. .not. leqi(id,'r')) return
         if (leqi(c,'n') .and. .not. leqi(id,'n')) return
      enddo
      match = .true.
      end function match
!------------------------------------------------------------------------
!     
!========================================================================
!     Return the number of different kinds of tokens.
!
      function nintegers(p,after)
      type(parsed_line), pointer     :: p
      integer, intent(in), optional  :: after
      integer nintegers

      integer i, starting_pos

      starting_pos = 0
      if (present(after)) then
         if (after.lt.0) call die("Wrong starting pos. in nintegers")
         starting_pos = after
      endif
      nintegers = 0
      do i=starting_pos+1,p%ntokens
         if (leqi(p%id(i),'i'))  nintegers = nintegers + 1
      enddo
      end function nintegers
!.......................................................................
      function nreals(p,after)
      type(parsed_line), pointer     :: p
      integer, intent(in), optional  :: after
      integer nreals

      integer i, starting_pos

      starting_pos = 0
      if (present(after)) then
         if (after.lt.0) call die("Wrong starting pos. in nreals")
         starting_pos = after
      endif
      nreals = 0
      do i=starting_pos+1,p%ntokens
         if (leqi(p%id(i),'r'))  nreals = nreals + 1
      enddo
      end function nreals
!.......................................................................
      function nvalues(p,after)
      type(parsed_line), pointer     :: p
      integer, intent(in), optional  :: after
      integer nvalues

      integer i, starting_pos

      starting_pos = 0
      if (present(after)) then
         if (after.lt.0) call die("Wrong starting pos. in nvalues")
         starting_pos = after
      endif
      nvalues = 0
      do i=starting_pos+1,p%ntokens
         if (leqi(p%id(i),'i').or.leqi(p%id(i),'r'))
     $        nvalues = nvalues + 1
      enddo
      end function nvalues
!.......................................................................
      function nnames(p,after)
      type(parsed_line), pointer     :: p
      integer, intent(in), optional  :: after
      integer nnames

      integer i, starting_pos

      starting_pos = 0
      if (present(after)) then
         if (after.lt.0) call die("Wrong starting pos. in nnames")
         starting_pos = after
      endif
      nnames = 0
      do i=starting_pos+1,p%ntokens
         if (leqi(p%id(i),'n'))  nnames = nnames + 1
      enddo
      end function nnames
!.......................................................................
      function ntokens(p,after)
      type(parsed_line), pointer     :: p
      integer, intent(in), optional  :: after
      integer ntokens

      integer starting_pos

      starting_pos = 0
      if (present(after)) then
         if (after.lt.0) call die("Wrong starting pos. in ntokens")
         starting_pos = after
      endif
      ntokens = p%ntokens - starting_pos
      if (ntokens.lt.0) call die("Wrong starting pos. in ntokens")
      end function ntokens
!-----------------------------------------------------------------------
!
!=======================================================================
!     Return a given token, specifying it by its sequence number
!     within a given kind. It is also possible to make the sequence
!     start after a given token number in the line.
!
      function integers(p,ind,after)
      integer  integers
      type(parsed_line), pointer     :: p
      integer, intent(in)            :: ind
      integer, intent(in), optional  :: after

      integer i, starting_pos, j

      starting_pos = 0
      if (present(after)) then
         if (after.lt.0) call die("Wrong starting pos. in integers")
         starting_pos = after
      endif

      j = 0
      do i=starting_pos+1, p%ntokens
         if (leqi(p%id(i),'i'))  j = j + 1
         if (j.eq.ind) then
            integers  = s2i(p%token_arr(i))
            return
         endif
      enddo
      call die("Not enough integers in line")
      end function integers
!.......................................................................
      function reals(p,ind,after)
      real(dp)  reals
      type(parsed_line), pointer     :: p
      integer, intent(in)            :: ind
      integer, intent(in), optional  :: after

      integer i, starting_pos, j

      starting_pos = 0
      if (present(after)) then
         if (after.lt.0) call die("Wrong starting pos. in reals")
         starting_pos = after
      endif

      j = 0
      do i=starting_pos+1, p%ntokens
         if (leqi(p%id(i),'r'))  j = j + 1
         if (j.eq.ind) then
            reals  = s2r(p%token_arr(i))
            return
         endif
      enddo
      call die("Not enough reals in line")

      end function reals
!.......................................................................
      function values(p,ind,after)
      real(dp) values
      type(parsed_line), pointer     :: p
      integer, intent(in)            :: ind
      integer, intent(in), optional  :: after

      integer i, starting_pos, j

      starting_pos = 0
      if (present(after)) then
         if (after.lt.0) call die("Wrong starting pos. in values")
         starting_pos = after
      endif

      j = 0
      do i=starting_pos+1, p%ntokens
         if (leqi(p%id(i),'i').or.leqi(p%id(i),'r'))  j = j + 1
         if (j.eq.ind) then
            values  = s2r(p%token_arr(i))
            return
         endif
      enddo
      call die("Not enough values in line")
      end function values
!.......................................................................
      function names(p,ind,after)
      character(len=max_length)      :: names
      type(parsed_line), pointer     :: p
      integer, intent(in)            :: ind
      integer, intent(in), optional  :: after

      integer i, starting_pos, j

      starting_pos = 0
      if (present(after)) then
         if (after.lt.0) call die("Wrong starting pos. in names")
         starting_pos = after
      endif

      j = 0
      do i=starting_pos+1, p%ntokens
         if (leqi(p%id(i),'n'))  j = j + 1
         if (j.eq.ind) then
            names  = p%token_arr(i)
            return
         endif
      enddo
      call die("Not enough names in line")
      end function names
!.......................................................................
      function tokens(p,ind,after)
      character(len=max_length) tokens
      type(parsed_line), pointer     :: p
      integer, intent(in)            :: ind
      integer, intent(in), optional  :: after

      integer starting_pos

      starting_pos = 0
      if (present(after)) then
         if ((after.lt.0).or.(after.ge.p%ntokens))
     $        call die("Wrong starting pos. in tokens")
         starting_pos = after
      endif

      tokens = p%token_arr(starting_pos+ind)
      end function tokens
!-----------------------------------------------------------------------

!
!    Main processing function
!

      function digest(str) result(pline)
      character(len=*), intent(in)   :: str
      type(parsed_line), pointer     :: pline

      integer, parameter         ::  maxntokens = 50

      character(len=max_length), save   ::  line
      integer, save              ::  ntokens
      integer, save              ::  first(maxntokens)
      integer, save              ::  last(maxntokens)
      character(len=1), save     ::  token_id(maxntokens)
c

      integer i

      line = str
      call line_parse
      call line_morphol

      nullify(pline)
      call create(pline)

      pline%ntokens = ntokens

      if (ntokens.ne.0) then
         allocate(pline%token_arr(ntokens))
         allocate(pline%id(ntokens))
         do i=1,ntokens
            pline%token_arr(i) = line(first(i):last(i))
            pline%id(i) = token_id(i)
         enddo
      endif


      CONTAINS   !!      line_parse, line_morphol, and helpers

      subroutine line_parse
      logical intoken, instring

      integer c
      integer stringdel

!     Character statement functions

      integer i
      logical isdigit, isupper, islower, isalpha,
     $        isalnum, isextra, istokch
      logical iscomment, isdelstr, isspecial

      isdigit(i) = (i .ge. 48) .and. (i .le. 57)
      isupper(i) = (i .ge. 65) .and. (i .le. 90)
      islower(i) = (i .ge. 97) .and. (i .le. 122)
      isalpha(i) = isupper(i) .or. islower(i)
      isalnum(i) = isdigit(i) .or. isalpha(i)

!     Extra characters allowed in tokens:  $ % * + & - . / @ ^ _ ~

      isextra(i) = ((i .ge. 36) .and. (i .le. 38))
     $             .or. (i .eq. 42) .or. (i .eq. 43)
     $             .or. (i .eq. 45) .or. (i .eq. 46) 
     $          .or. (i .eq. 47) .or. (i .eq. 64) .or. (i .eq. 94)
     $             .or. (i .eq. 95) .or. (i .eq. 126)

      istokch(i) = isalnum(i) .or. isextra(i)

!     Comments are signaled by:  !  #  ; 
      iscomment(i) = (i.eq.33) .or. (i.eq.35) .or. (i.eq.59)

!     String delimiters: "  '  `
      isdelstr(i) = (i.eq.34) .or. (i.eq.39) .or. (i.eq.96)

!     Special characters which are tokens by themselves: <
      isspecial(i) = (i.eq.60)

!!========================================================
!
      ntokens = 0

      intoken = .false.
      instring = .false.
      stringdel = 0
      
      do i = 1, len(line)
         c = ichar(line(i:i))

         if (iscomment(c)) then
! possible comment...
            if (instring) then
               last(ntokens) = i
            else
               exit
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
! string delimiter... make sure it is the right one before
! closing the string.
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

      if (line_debug) then
         write(line_log,*) '            ',  ntokens, ' tokens:'
         do i=1,ntokens
            write(line_log,*) '                 ',
     $           '|',line(first(i):last(i)),'|'
         enddo
      endif

      end subroutine line_parse

c-----------------------------------------------------------------
c
      subroutine line_morphol
c
c     Classifies the tokens according to their morphology
c
      integer i, ierr
      real(dp) real_value
      character(len=max_length) token

      do i = 1, ntokens
         token = line(first(i):last(i))
         if (is_value(token)) then
c
c           This read also serves to double check the token for
c           real meaning (for example, ".d0" should give an error)
c
            read(token,fmt=*,iostat=ierr) real_value
            if (ierr.ne.0) then
               write(0,'(a,i3,1x,a,/,a)')
     $              'Numeric conversion error at token number',i,
     $              'in line', line
               stop 'CONV'
            endif

            if (is_integer(token)) then
               token_id(i) = 'i'
            else
               token_id(i) = 'r'
            endif
         else
            token_id(i) = 'n'
         endif
      enddo
      end subroutine line_morphol
c
      logical function is_integer(string)
      character(len=*) string
c
      integer i, length

      logical is_digit, is_sign
      character c*1
      is_digit(c) = (ichar(c).ge.48 .and. ichar(c).le.57)
      is_sign(c) = (c.eq.'+' .or. c.eq.'-')
c
      length = len_trim(string)
      is_integer = .false.
      if (length.eq.0) return
      
      c = string(1:1)
      if (.not. (is_digit(c) .or. is_sign(c))) return
      do i=2,length
         c = string(i:i)
         if (.not. (is_digit(c))) return
      enddo
      is_integer = .true.
      end function is_integer
!
      logical function is_value(string)
      character(len=*) string
c
      integer i, length, exp_mark

      logical is_digit, is_sign, is_dot, is_expmark
      character c*1
      logical dotsok
      is_digit(c) = (ichar(c).ge.48 .and. ichar(c).le.57)
      is_sign(c) = (c.eq.'+' .or. c.eq.'-')
      is_dot(c) = (c.eq.'.' .and. dotsok)
      is_expmark(c) = (c.eq.'e' .or. c.eq.'E' .or. c.eq.'d'
     $                          .or. c.eq.'D')
c

      length = len_trim(string)

      is_value = .false.
      dotsok = .true.
c
c     Find the starting point of a possible exponent
c
      exp_mark = length+1
      do i = 1, length
         c = string(i:i)
         if (is_expmark(c)) exp_mark = i
      enddo
      if (exp_mark .eq. length) return
c
      c = string(1:1)
      if (.not. (is_digit(c) .or. is_sign(c))) then
         if (is_dot(c)) then
            dotsok = .false.
         else
            return
         endif
      endif
      do i=2,exp_mark-1
         c = string(i:i)
         if (.not. (is_digit(c))) then
            if (is_dot(c)) then
               dotsok = .false.
            else
               return
            endif
         endif
      enddo
c 
c     Is the exponent an integer?
c
      if (exp_mark .lt. length) then
         if (.not. is_integer(string(exp_mark+1:length))) return
      endif
c
c     Here we could do some extra checks to see if the string still makes
c     sense... For example, "." and ".d0" pass the above tests but are not
c     readable as numbers. I believe this should be reported by the
c     conversion routine, to warn the user of a mis-typed number, instead
c     of reporting it as a string and break havoc somewhere else.
c
      is_value = .true.

      end function is_value

      end function digest

!--------------------------------------------------------
c-----------------------------------------------------------------
c
c     Search function
c
      function search_line(p,string,ind,after,eq_func) result(found)
      logical found
      type(parsed_line), pointer :: p
      character(len=*) string
      integer, intent(out) :: ind
      integer, intent(in), optional :: after
      optional :: eq_func

      interface
         function eq_func(s1,s2)
         logical eq_func
         character(len=*), intent(in) :: s1,s2
         end function eq_func
      end interface

      integer i, starting_pos

      starting_pos = 0
      if (present(after)) then
         if (after.lt.0) call die("Wrong starting pos. in search_line")
         starting_pos = after
      endif

      ind = -1
      found = .false.
      if (.not. associated(p)) return
!
!     The default comparison routine is 'leqi' (case-insensitive)
!
      if (present(eq_func)) then
         do i = starting_pos+1, p%ntokens
            if (eq_func(string,p%token_arr(i))) then
               ind = i
               found = .true.
               return
            endif
         enddo
      else
         do i = starting_pos+1, p%ntokens
            if (leqi(string,p%token_arr(i))) then
               ind = i
               found = .true.
               return
            endif
         enddo
      endif

      end function search_line
c
c     Examples of eq_func's...
c
      logical function eq_strict(str1,str2)
      character(len=*) str1, str2
      eq_strict = (str1 .eq. str2)
      end function eq_strict
c
c-----------------------------------------------------------------
c     Control routines
c
      subroutine line_setdebug(level)
      integer level
      line_debug = (level .eq. 1)
      end subroutine line_setdebug
c
      subroutine line_setlog(unit)
      integer unit
      line_log = unit
      end subroutine line_setlog
c---------------------------------------------------------------------
c
      SUBROUTINE CHRCAP(STRING,NCHAR)
C
C  CHRCAP accepts a STRING of NCHAR characters and replaces
C  any lowercase letters by uppercase ones.
C
      CHARACTER STRING*(*)
      CHARACTER CHAR*1
      integer nchar, ncopy, i, itemp
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
      END subroutine chrcap
c
      FUNCTION LEQI(STRNG1,STRNG2)
      logical leqi
      character(len=*) strng1, strng2
C
C  Case-insensitive lexical equal-to comparison
C
      CHARACTER*1   S1,S2
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
      END function leqi

!
!
      function s2i(str)
      integer s2i
      character(len=*), intent(in) :: str

      integer ierr
      read(str,fmt=*,iostat=ierr) s2i
      if (ierr.ne.0) call die("Integer conversion error")
      end function s2i

      function s2r(str)
      real(dp) s2r
      character(len=*), intent(in) :: str

      integer ierr
      read(str,fmt=*,iostat=ierr) s2r
      if (ierr.ne.0) call die("Real conversion error")
      end function s2r
!
!     Private copy of "die" (not really MPI-safe yet)
!
      subroutine die(str)
      character(len=*), intent(in), optional   :: str
      if (present(str)) write(6,'(a)') trim(str)
      write(6,'(a)') 'Stopping Program'
      stop
      end subroutine die

      end module parse



