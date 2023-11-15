module omp_module
  !! The idea of this module is to isolate the dependency on omp_lib to one place
  use omp_lib

  implicit none

contains
  
  subroutine init_threads(number_of_threads)
    !! Return the number of threads given by the omp library call

    integer, intent(out) :: number_of_threads
    
    number_of_threads = omp_get_max_threads()
    !! Outside of a parallel region, get_max_threads returns the number of available threads
    
  end subroutine init_threads
  
end module omp_module
