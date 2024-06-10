module omp_module
  !! This module allows compiling without the OpenMP library
  implicit none

contains
  
  subroutine init_threads(number_of_threads)
    !! Return 0

    integer, intent(out) :: number_of_threads
    
    number_of_threads = 0
    
  end subroutine init_threads
  
end module omp_module
