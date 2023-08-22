!!****h* Conquest/biblio *
!!  NAME
!!   biblio_module
!!  PURPOSE
!!   Generate a relevant bibliography based on input flags
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   2019/07/03
!!  MODIFICATION HISTORY
!!   2020/01/03 12:14 dave
!!    Added check in get_ref to catch non-existent references
!!  SOURCE
!!
module biblio

  use global_module,  only: io_lun
  use GenComms,       only: inode, ionode

  implicit none

  character(len=*), parameter :: bibtex_file = "conquest.bib"
  character(len=*), parameter :: ref_fmt = '(4x,a," ",i0,", ",i0," (",i0,")")'
  character(len=*), parameter :: doi_fmt = '(4x,"http://dx.doi.org/",a)'
  character(len=*), parameter :: bib_key_fmt = '("@article{",a,",")'
  character(len=*), parameter :: bib_afield_fmt = '(2x,a," = {{",a,"}},")'
  character(len=*), parameter :: bib_ifield_fmt = '(2x,a," = {",i0,"},")'
  integer, parameter          :: max_refs = 50
  logical                     :: flag_dump_bib

  type type_reference

     character(len=40)   :: key
     character(len=400)  :: authors
     character(len=400)  :: title
     character(len=400)  :: journal
     integer             :: volume
     integer             :: page
     integer             :: year
     character(len=80)   :: doi
     character(len=132)   :: comment

   contains

     procedure, public :: cite_reference
     procedure, public :: write_bib

  end type type_reference

  type type_bibliography

     type(type_reference), allocatable, dimension(:) :: db
     integer             :: nrefs
     logical             :: first

   contains

     procedure, public :: init_bib
     procedure, public :: close_bib
     procedure, public :: add_ref
     procedure, public :: get_ref
     procedure, public :: cite

  end type type_bibliography

contains

  !!****f* biblio/cite_reference *
  !!
  !!  NAME 
  !!   cite_reference
  !!  PURPOSE
  !!   Cite a reference in Conquest_out
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/07/04
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine cite_reference(ref)

    ! passed variables
    class(type_reference) :: ref

    if (inode==ionode) then
       write(io_lun,'(4x,a)') trim(ref%comment)
       write(io_lun,'(4x,a)') trim(ref%title)
       write(io_lun,'(4x,a)') trim(ref%authors)
       write(io_lun,ref_fmt) trim(ref%journal), ref%volume, ref%page, ref%year
       write(io_lun,doi_fmt) trim(ref%doi)
       write(io_lun,*)
    end if

  end subroutine cite_reference
  !!***

  !!****f* biblio/write_bib *
  !!
  !!  NAME 
  !!   write_bib
  !!  PURPOSE
  !!   Generate a BibTeX record as a string
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/07/04
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine write_bib(ref, first)

    use io_module, only: io_assign, io_close

    ! passed variables
    class(type_reference), intent(in) :: ref
    logical, intent(inout)            :: first

    ! local variables
    integer :: lun
    character(len=400)                :: str

   if(inode == ionode) then
    call io_assign(lun) 
    if (first) then
       open(unit=lun,file=bibtex_file,status='replace')
       first = .false.
    else
       open(unit=lun,file=bibtex_file,position='append')
    end if

    write(lun,bib_key_fmt) trim(ref%key)
    write(lun,bib_afield_fmt) "author", trim(ref%authors)
    write(lun,bib_afield_fmt) "title", trim(ref%title)
    write(lun,bib_afield_fmt) "journal", trim(ref%journal)
    write(lun,bib_ifield_fmt) "year", ref%year
    write(lun,bib_ifield_fmt) "volume", ref%volume
    write(lun,bib_ifield_fmt) "pages", ref%page
    write(lun,bib_afield_fmt) "doi", trim(ref%doi)
    write(lun,'("}")')

    call io_close(lun)

   endif !(inode == ionode) 

  end subroutine write_bib
  !!***

  type(type_reference) function get_ref(bib, key)

    use input_module,  only: leqi
    use GenComms,      only: cq_abort

    ! passed variables
    class(type_bibliography), intent(inout) :: bib
    character(*), intent(in)                :: key

    ! local variables
    integer       :: i, done
    character(40) :: k

    k = key

    done = 0
    do i=1,max_refs
       if (leqi(bib%db(i)%key, k)) then
          get_ref = bib%db(i)
          done = 1
          exit
       end if
    end do
    if(done==0) call cq_abort("Unable to find reference key "//key)
  end function get_ref

  subroutine init_bib(bib)

    ! passed variables
    class(type_bibliography), intent(inout) :: bib

    ! local variables
    integer :: i

    allocate(bib%db(max_refs))
    bib%nrefs = 0
    bib%first = .true.

  end subroutine init_bib

  subroutine close_bib(bib)

    ! passed variables
    class(type_bibliography), intent(inout) :: bib

    ! local variables
    integer :: i

    deallocate(bib%db)

  end subroutine close_bib

  !!****f* biblio/add_ref *
  !!
  !!  NAME 
  !!   add_ref
  !!  PURPOSE
  !!   Add a reference
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/07/04
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine add_ref(bib, key, authors, title, journal, volume, page, &
       year, doi, comment)

    use input_module, only: leqi
    use GenComms,     only: cq_abort

    ! passed variables
    class(type_bibliography), intent(inout) :: bib
    character(*), intent(in)  :: key
    character(*), intent(in)  :: authors
    character(*), intent(in)  :: title
    character(*), intent(in)  :: journal
    integer, intent(in)       :: volume
    integer, intent(in)       :: page
    integer, intent(in)       :: year
    character(*), intent(in)  :: doi
    character(*), intent(in)  :: comment

    ! local variables
    type(type_reference), allocatable, target :: reference

    bib%nrefs = bib%nrefs + 1
    if (bib%nrefs > max_refs) &
         call cq_abort("Number of references > max_refs: ", bib%nrefs, max_refs)

    bib%db(bib%nrefs)%key = key
    bib%db(bib%nrefs)%authors = authors
    bib%db(bib%nrefs)%title = title
    bib%db(bib%nrefs)%journal = journal
    bib%db(bib%nrefs)%volume = volume
    bib%db(bib%nrefs)%page = page
    bib%db(bib%nrefs)%year = year
    bib%db(bib%nrefs)%doi = doi
    bib%db(bib%nrefs)%comment = comment

  end subroutine add_ref
  !!***

  !!****f* biblio/cite *
  !!
  !!  NAME 
  !!   cite
  !!  PURPOSE
  !!   Write reference to output file and bibliography
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/07/04
  !!  MODIFICATION HISTORY
  !!   2020/03/24 14:22 dave
  !!    Changed to write reference key only if iprint<2
  !!   2022/09/21 10:58 dave
  !!    Changed to write reference key only if iprint<4
  !!  SOURCE
  !!
  subroutine cite(bib, key, punc, pre)

    use global_module, ONLY: iprint_init
    ! passed variables
    class(type_bibliography), intent(inout) :: bib
    character(*), intent(in)                :: key
    character(2), intent(in), optional      :: punc
    character(*), intent(in), optional      :: pre

    ! local variables
    type(type_reference)                    :: reference
    character(20) :: start

    reference = bib%get_ref(key)
    if(iprint_init>3) then
       call reference%cite_reference
       if (flag_dump_bib) call reference%write_bib(bib%first)
    else
       if(present(pre)) then
          write(start,'(a)') pre
          if(present(punc)) then
             write(io_lun,fmt='(4x,a18,a,a2)',advance='no') adjustl(start), trim(reference%key), punc
          else
             write(io_lun,fmt='(4x,a18,a)',advance='no') adjustl(start), trim(reference%key)
          end if
       else
          if(present(punc)) then
             write(io_lun,fmt='(a,a2)',advance='no') trim(reference%key),punc
          else
             write(io_lun,fmt='(a)',advance='no') trim(reference%key)
          end if
       end if
       call reference%write_bib(bib%first)
    end if

  end subroutine cite
  !!***

end module biblio
