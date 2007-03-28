c
c Copyright Alberto Garcia, 1996, 1997, 1998
c
c This module implements an interface to the FORTRAN logical unit
c system. Based on code by Richard Maine.
c
c
c Alberto Garcia, December 30, 1996
c Rewritten as a single subroutine 
c with multiple entry points, March 7, 1998
c
c This scheme is actually the closest in spirit to f90 modules, but
c in perfectly legal f77.
c
c---------------------------------------------------------------
c
      subroutine io
c
c     Logical unit management. Units 0 to min_lun-1 are "reserved",
c     since most of the "typical" files (output, etc) use them.
c
c     Logical units min_lun to min_max are managed by this module.
      
      implicit none
c
c----------------------------------------------------------------
c     Module variables
c
      integer stdout, stderr
      integer min_lun, max_lun, nunits
      parameter (min_lun=10, max_lun=99, nunits=max_lun-min_lun+1)
      logical lun_is_free(min_lun:max_lun)

      save stdout, stderr, lun_is_free
c-----------------------------------------------------------------
c
c     Internal and dummy variables
c
      integer i, unit, lun, iostat
      logical used, named, opened
      character filename*50, form*11
c
c-----------------------------------------------------------------
c     Initialization section
c
      data lun_is_free /nunits*.true./
      data stdout, stderr /6,0/
c-----------------------------------------------------------------
c
c     Executable routines
c
c     Simple interfaces to modify standard units
c
      entry io_seterr(unit)
      stderr = unit
      return
      entry io_setout(unit)
      stdout = unit
      return

      entry io_geterr(unit)
      unit = stderr
      return
      entry io_getout(unit)
      unit = stdout
      return
c
c------------------------------------------------------------------     
c
c     Logical unit management
c
      entry io_assign(lun)
c
c     Looks for a free unit and assigns it to lun
c
      do lun= min_lun, max_lun
         if (lun_is_free(lun)) then
            inquire(unit=lun, opened=used, iostat=iostat)
            if (iostat .ne. 0) used = .true.
            lun_is_free(lun) = .false.
            if (.not. used) return
         endif
      enddo
      write(stderr,'(a)') 'No luns available in io_assign'
      stop 'LUN'
c
c===
c
      entry io_reserve(lun)
c
c     Useful to specify that one needs to use a particular unit number
c
c     For example, assume some legacy code expects to work with unit 15:
c
c     call io_reserve(15)   ! this call at the beginning of the program
c     ...
c     open(15,....)
c
      inquire(unit=lun, opened=used, iostat=iostat)
      if (iostat .ne. 0) used = .true.
      if (used) then
         write(stderr,'(a,i3,a)')
     $        'Cannot reserve unit',lun,'. Already connected'
         stop 'LUN'
      endif
      if (lun .ge. min_lun .and. lun .le. max_lun)
     $                      lun_is_free(lun) = .false.

      return
c
c===
c
      entry io_close(lun)
c
c     Use this routine instead of a simple close!!
c
      close(lun)
      if (lun .ge. min_lun .and. lun .le. max_lun)
     $                     lun_is_free(lun) = .true.
      return
c
c===
c
      entry io_status
c
c     Prints a list of the connected logical units and the names of
c     the associated files
c

      write(stdout,'(a)') '******** io_status ********'
      do i = 0, max_lun
         inquire(i,opened=opened,named=named,name=filename,
     $           form=form,iostat=iostat)
         if (iostat .eq. 0) then
            if (opened) then
               if (named) then
                  write(stdout,9000) i, form, filename
               else
                  write(stdout,9000) i, form, 'No name available'
               endif
            endif
         else
            write(stdout,9000) i, 'Iostat error'
         endif
      enddo
      write(stdout,'(a)') '********           ********'

 9000 format(i4,5x,a,5x,a)
      return

      end






