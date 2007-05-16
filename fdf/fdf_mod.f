! 
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2003
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
      module fdf
!
!  Copyright Alberto Garcia, Jose M. Soler (1996-)
!
!=====================================================================
!
!     This module implements an extended Fortran 90/95 interface
!     to the Flexible Data Format library of A. Garcia and J.M. Soler,
!     originally written in Fortran 77.
!
!     It provides interface blocks for the routines in the f77 library,
!     as well as new routines which take advantage of the new facilities
!     in Fortran 90/95.
!
!
!     NEW FEATURES:
!
!     a) Block pointers. 
!
!     Block content can now be flexibly handled by means of a pointer
!     to a derived type "block". Typical usage:
!
!     use fdf
!     type(block), pointer :: bp
!
!     if (fdf_block('SomeBlock',bp)) then
!         loop: do
!               if (.not. fdf_bline(bp,line)) exit loop
!               (process line, possibly with 'digest')
!         enddo loop
!     endif
!     call destroy(bp)
!
!     The generic name 'fdf_block' stands for both the old and the new
!     syntax (In the old syntax: fdf_block('SomeBlock',unit), the function
!     returns in 'unit' a unit number from which to read the contents of
!     the block.)
!
!     Routine fdf_bline returns in 'line' the next non-blank, 
!     non-comment line from the block, unless there are no more
!     lines, in which case it returns .false. and 'line' is undefined.
!     Optionally, by specifying 'any=.true.', fdf_bline will return
!     the next line, even if blank.
!
!     Routine 'backspace' moves an internal pointer to the previous line
!     returned (unless the optional 'last_physical=.true.' is specified, 
!     in which case the internal pointer will move to the previous physical 
!     line in the block (even if a blank or comment line); this can lead 
!     to unpredictable behavior) 
!
!     Routine 'rewind' moves the internal pointer to the beginning of 
!     the block.
!
!     The call to 'destroy' frees the storage associated to the pointer bp.
!
!     Among the advantages of the new interface to block handling, are:
!
!     * Automatic detection of the end of the block.
!     * Block pointers can be kept around for as long as needed.
!     * Will work even if the underlying mechanism for block management
!       is not based on files.
!
!     b) Generic interface to scalar routines.
!
!     The generic function 'fdf_get' can be used instead of any of the
!     old scalar routines. The old names are also accepted, for backwards
!     compatibility. 
!
!     c) New routine returning a "digested" string.
!
!     Function fdf_parsed_string(label,default) returns a pointer to
!     a 'parsed_line' derived type (see module parse).
!
!------------------------------------------------------------------------
!     Future enhancements:
!
!     * Stand-alone f90/95 version.
!     * Better exception handling.
!     * MPI-awareness.
!     * Array primitives.
!     * Nested FDF namespaces.
!     * Non-recursive implementation of 'destroy'
!========================================================================

      implicit none

      private 

!
!     Routines already in the f77 library
!
      public fdf_init
      public fdf_convfac
      public fdf_integer, fdf_single, fdf_double, fdf_physical
      public fdf_string, fdf_boolean
      public fdf_defined, fdf_enabled, fdf_inhibit
!
!     Generic interface for old and new 'block' routines
!
      public fdf_block
!
!     New routines
!
      public block, destroy
      public print_block, backspace, rewind
      public fdf_bline
!
!     Private kind declarations
!
      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)

c     Declarations for fdf procedures

      interface fdf_block
          module procedure fdf_blockf, fdf_blockp
      end interface
      
      interface destroy
        module procedure destroy_bp
      end interface

      interface backspace
        module procedure backspace_fdf_block
      end interface

      interface rewind
        module procedure rewind_fdf_block
      end interface

      interface
!
!        Functions in fdf.f
!
         function fdf_defined(label)
         logical fdf_defined
         character(len=*), intent(in) :: label
         end function fdf_defined

         function fdf_enabled()
         logical fdf_enabled
         end function fdf_enabled

         function fdf_integer(label,default)
         integer fdf_integer
         character(len=*), intent(in) :: label
         integer, intent(in) ::  default
         end function fdf_integer

         function fdf_single(label,default)
         real fdf_single
         character(len=*), intent(in) :: label
         real, intent(in) ::  default
         end function fdf_single

         function fdf_double(label,default)
         real*8 fdf_double
         character(len=*), intent(in) :: label
         real*8, intent(in) ::  default
         end function fdf_double

         function fdf_physical(label,default,unit)
         real*8 fdf_physical
         character(len=*), intent(in) :: label, unit
         real*8, intent(in) ::  default
         end function fdf_physical

         function fdf_boolean(label,default)
         logical fdf_boolean
         character(len=*), intent(in) :: label
         logical, intent(in) ::  default
         end function fdf_boolean

         function fdf_string(label,default)
         character(len=132) fdf_string
         character(len=*), intent(in) :: label
         character(len=*), intent(in) ::  default
         end function fdf_string


         function fdf_convfac(unit1,unit2)
         real(selected_real_kind(14,100)) fdf_convfac
         character(len=*), intent(in) :: unit1, unit2
         end function fdf_convfac

         subroutine fdf_init(filein,fileout)
         character(len=*), intent(in) :: filein, fileout
         end subroutine fdf_init

         subroutine fdf_inhibit
         end subroutine fdf_inhibit

      end interface
!
!     New derived types to support blocks
!
      type line_dlist
         character(len=132)                ::  str
         type(line_dlist), pointer         ::  next
         type(line_dlist), pointer         ::  prev
      end type line_dlist

      type block
      private
         type(line_dlist), pointer         ::  mark
         type(line_dlist), pointer         ::  txt
         type(line_dlist), pointer         ::  last
         type(line_dlist), pointer         ::  last_line_returned
      end type block

      interface
         function leqi(s1,s2)
         logical leqi
         character(len=*), intent(in)   :: s1, s2
         end function leqi
      end interface

      CONTAINS

!-------------------------------------------------------------------
      subroutine destroy_bp(bp)
      type(block), pointer       :: bp
      if (associated(bp%txt)) call destroy_dl(bp%txt)
      deallocate(bp)
      end subroutine destroy_bp

!-------------------------------------------------------------------
      recursive subroutine destroy_dl(dlp)
      type(line_dlist), pointer       :: dlp
      if (associated(dlp%next)) call destroy_dl(dlp%next)
      deallocate(dlp)
      end subroutine destroy_dl

!-------------------------------------------------------------------
!
!        To be able to use a generic fdf_block, the two instances
!        (old and new interface) have to be module procedures.
!        Here is fdf_blockf. Note that, to avoid a scope bug in
!        pgf90, we need to use a new name for the f77 routine.
!
         function fdf_blockf(label,unit)
         logical fdf_blockf
         character(len=*), intent(in) :: label
         integer, intent(out) ::  unit

         interface
            function fdf_block_old(label,unit)
            logical fdf_block_old
            character(len=*), intent(in) :: label
            integer, intent(out)  :: unit
            end function fdf_block_old
         end interface

         fdf_blockf = fdf_block_old(label,unit)
         end function fdf_blockf
!
!-------------------------------------------------------------------
!        Fill in block structure
!
         function fdf_blockp(label,bp) result(res)
         use parse
         logical res
         character(len=*), intent(in) :: label
         type(block), pointer         :: bp

         integer unit, ierr
         character(len=132) line
         logical head
         type(line_dlist), pointer         :: dlp
         type(parsed_line), pointer        :: p

         if (associated(bp)) call destroy(bp)
         head = .true.

         res = fdf_blockf(label,unit)
         if (.not.res) return

         allocate(bp)
         nullify(bp%mark)
         nullify(bp%txt)
         loop: DO

           read(unit,fmt='(a)',iostat=ierr) line
           if (ierr .ne. 0) exit loop

           p=>digest(line)
           if (ntokens(p) .ge. 1) then
              if (leqi(tokens(p,1),"%endblock")) exit loop
           endif

           if (head) then
              allocate(bp%txt)
              nullify(bp%txt%prev)
              dlp => bp%txt
              bp%mark=>bp%txt
              bp%last=>bp%txt
              head = .false.
           else
              allocate(dlp%next)
              dlp%next%prev => dlp
              dlp=>dlp%next
           endif

           dlp%str = line
           nullify(dlp%next)
           bp%last => dlp
           call destroy(p)

        enddo loop

        if (.not. associated(bp%txt)) then
           !!! Empty block!!!
           call warn("fdf_blockp: Block is empty...")
           res = .false.
        else
           bp%last_line_returned => bp%txt   
        endif

         end function fdf_blockp
!
!-------------------------------------------------------------------
      subroutine backspace_fdf_block(bp,physical_line)
      type(block), pointer       :: bp
      logical, intent(in), optional :: physical_line

      logical last_physical_line

      last_physical_line = .false.
      if (present(physical_line)) last_physical_line = physical_line

      if (.not. last_physical_line) then
!        Put the mark at the point of the last returned line,
!        as determined in fdf_bline
         bp%mark => bp%last_line_returned
         return
      endif
!
!     Backspace to the previous physical line in the block
!     (i.e., it might be a blank or a comment line)
!
      if (.not. associated(bp%mark)) then  ! We are at the end of block
         bp%mark=>bp%last
      else
         if (.not. associated(bp%mark%prev)) then ! at the beginning
            bp%mark => bp%txt
         else
            bp%mark => bp%mark%prev
         endif
      endif
      end subroutine backspace_fdf_block

!-------------------------------------------------------------------
      subroutine rewind_fdf_block(bp)
      type(block), pointer       :: bp

      if (.not. associated(bp))
     $     call die("rewind: Block not associated")
      if (.not. associated(bp%txt))
     $     call die("rewind: Block text not associated")
      bp%mark=>bp%txt
      bp%last_line_returned=>bp%txt
      end subroutine rewind_fdf_block
!
!
!-------------------------------------------------------------------
!     Get successive non-blank, non-comment lines from block
!     Optionally, if 'any=.true.' is specified, return any line,
!     even if blank or a comment line.
!
      function fdf_bline(bp,line,any) result(res)
      use parse
      logical res
      type(block), pointer       :: bp
      character(len=*), intent(out)  :: line
      logical, intent(in), optional  :: any
      
      type(parsed_line), pointer   :: p
      logical any_line

      any_line = .false.
      if (present(any)) any_line = any

      res = .false.
      loop: do
         if (.not.associated(bp%mark))      return
         line = bp%mark%str
         bp%last_line_returned=>bp%mark
         bp%mark => bp%mark%next
         p => digest(line)
         if ((ntokens(p) .ne. 0)  .or. any_line) exit loop
      enddo loop
      call destroy(p)
      res = .true.
      end function fdf_bline

!-------------------------------------------------------------------
!     Print block
!
      subroutine print_block(bp)
      type(block), pointer       :: bp

      type(line_dlist), pointer       :: p
      character(len=70) :: line

      if (.not. associated(bp)) return
      if (.not. associated(bp%txt)) return
      p=>bp%txt
 5    continue
         if (.not.associated(p)) return
         line = p%str
         write(6,'(a70)') line
         p => p%next
         goto 5
      end subroutine print_block
!
!-------------------------------------------------------------------
!
!     Private copy of die
!
      subroutine die(str)
      character(len=*), optional :: str
      if (present(str)) write(6,'(a)') str
      stop
      end subroutine die
      subroutine warn(str)
      character(len=*) :: str
      write(6,'(2a)') '*WARNING: ', str
      end subroutine warn

      end module fdf





