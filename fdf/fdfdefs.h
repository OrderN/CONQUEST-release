c
c Copyright Alberto Garcia, Jose Soler, 1996, 1997, 1998
c
c     Declarations for external fdf functions
c
      integer fdf_integer
      real    fdf_single
      real*8  fdf_double, fdf_physical
      logical fdf_boolean, fdf_defined, fdf_enabled, fdf_block
      character*80 fdf_string
c
      external fdf_integer, fdf_block, fdf_single, fdf_double,
     $         fdf_physical, fdf_boolean, fdf_defined, fdf_string,
     $         fdf_enabled
