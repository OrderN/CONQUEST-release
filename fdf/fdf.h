c
c Copyright Alberto Garcia, Jose Soler, 1996, 1997, 1998
c---
c     I/O variables for the fdf package. 
c
c     In Fortran 90, all this should go in a module...
c
c
c     ndepth: number of open files (maximum maxdepth)
c     fdf_stack holds their unit numbers.

      integer maxdepth
      parameter (maxdepth=5)
      integer ndepth, fdf_stack(maxdepth)
c
c     Unit numbers for input, output, error notification, and
c     debugging output (the latter active if fdf_debug is true)
c
      integer fdf_in, fdf_out, fdf_err, fdf_log
      common /fdf_io/ fdf_in, fdf_out, fdf_err, fdf_log, 
     $                ndepth, fdf_stack
      logical fdf_debug, fdf_debug2, fdf_started, fdf_donothing
      common /fdf_logicals/ fdf_debug, fdf_debug2, fdf_started,
     $                      fdf_donothing

      save /fdf_io/, /fdf_logicals/
c
c     Line just read and parsing info
c
      character*132 line
      integer maxntokens
      parameter (maxntokens=50)
      integer ntokens
      integer first(maxntokens), last(maxntokens)
      common /fdf_line/ line
      common /fdf_parsing/ ntokens, first, last

      save /fdf_line/, /fdf_parsing/
c---


