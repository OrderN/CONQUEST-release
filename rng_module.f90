!!****h* Conquest/rng *
!!  NAME
!!   rng
!!  PURPOSE
!!   A place for all random number generator-related code. Mersenne Twister
!!   code adapted from original C code by Takuji Nishimura and Makoto
!!   Matsumoto, and translated to Fortran by Remi Piatek (see below)
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   2019/05/16
!!  MODIFICATION HISTORY
!!   2019/05/21 zamaan
!!    Changes to intialisation to make rng_normal easier to use. Now one call
!!    returns a single random number, (the other is buffered) but it needs to
!!    be initialised using init_normal
!!  SOURCE
!!  
!-------------------------------------------------------------------------------
!   This is a Fortran translation of the 64-bit version of
!   the Mersenne Twister pseudorandom number generator
!
!   Translated from C-program for MT19937-64 (2004/9/29 version)
!   originally coded by Takuji Nishimura and Makoto Matsumoto
!   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
!
!   Fortran translation by Remi Piatek
!   The University of Copenhagen
!   Department of Economics
!   email: {first}.{last}@econ.ku.dk
!
!-------------------------------------------------------------------------------
!   A C-program for MT19937-64 (2004/9/29 version).
!   Coded by Takuji Nishimura and Makoto Matsumoto.
!
!   This is a 64-bit version of Mersenne Twister pseudorandom number
!   generator.
!
!   Before using, initialize the state by using init_genrand64(seed)  
!   or init_by_array64(init_key, key_length).
!
!   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
!   All rights reserved.                          
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions
!   are met:
!
!     1. Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.
!
!     2. Redistributions in binary form must reproduce the above copyright
!        notice, this list of conditions and the following disclaimer in the
!        documentation and/or other materials provided with the distribution.
!
!     3. The names of its contributors may not be used to endorse or promote 
!        products derived from this software without specific prior written 
!        permission.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!   References:
!   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
!     ACM Transactions on Modeling and 
!     Computer Simulation 10. (2000) 348--357.
!   M. Matsumoto and T. Nishimura,
!     ``Mersenne Twister: a 623-dimensionally equidistributed
!       uniform pseudorandom number generator''
!     ACM Transactions on Modeling and 
!     Computer Simulation 8. (Jan. 1998) 3--30.
!
!   Any feedback is very welcome.
!   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
!   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
!-------------------------------------------------------------------------------
module rng
  
  use datatypes

  implicit none

  integer(wide), parameter, private :: nn       = 312_wide
  integer(wide), parameter, private :: mm       = 156_wide
  integer(wide), parameter, private :: seed_def = 5489_wide
  integer(wide), parameter, private :: matrix_a = -5403634167711393303_wide
  integer(wide), parameter, private :: um       = -2147483648_wide
  ! most significant 33 bits
  integer(wide), parameter, private :: lm       = -2147483647_wide
  ! least significant 31 bits
  real(double), parameter, private      :: pi253_1 = 1._double/(2._double**53 -&
                                                     1._double)
  real(double), parameter, private      :: pi253   = 1._double/(2._double**53)
  real(double), parameter, private      :: pi252   = 1._double/(2._double**52)

  type type_rng

    integer(wide), allocatable  :: seed(:)

    integer(wide) :: mt(nn)     ! array for state vector
    integer       :: mti = nn+1 ! mti==nn+1 means mt(nn) is not initialized

    ! For normal distribution
    logical       :: buffered
    real(double)  :: buffer, sigma, mu

    contains

      procedure, private  :: init_genrand64
      procedure, private  :: init_by_array64
      procedure, private  :: genrand64_int64
      procedure, private  :: genrand64_real1
      procedure, private  :: genrand64_real2
      procedure, private  :: genrand64_real3
      procedure, private  :: get_random_seed
      procedure, public   :: init_rng
      procedure, public   :: init_normal
      procedure, public   :: rng_uniform
      procedure, public   :: rng_normal
      procedure, public   :: rng_integer
  end type type_rng
  !!***

  contains

    !!****f* rng/init_genrand64 *
    !!
    !!  NAME 
    !!   init_genrand64
    !!  PURPOSE
    !!   Initialises mt(nn) with a seed
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/16
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    subroutine init_genrand64(rn, seed)

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)  :: rn
      integer(wide), intent(in)       :: seed

      ! Local variables
      integer                 :: i

      rn%mt(1) = seed
      do i = 1, nn-1
        rn%mt(i+1) = 6364136223846793005_wide * ieor(rn%mt(i), ishft(rn%mt(i), -62)) + i
      end do

      rn%mti = nn

    end subroutine init_genrand64
    !!***

    !!****f* rng/init_by_array64 *
    !!
    !!  NAME 
    !!   init_by_array64
    !!  PURPOSE
    !!   Initialises by array with an array length
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/16
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    subroutine init_by_array64(rn, init_key)

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)  :: rn
      integer(wide), intent(in)       :: init_key(:)

      ! Local variables
      integer(wide), parameter  :: c1 = 3935559000370003845_wide
      integer(wide), parameter  :: c2 = 2862933555777941757_wide
      integer(wide)             :: i, j, k, kk, key_length

      call rn%init_genrand64(19650218_wide)
      key_length = size(init_key)
      i = 1_wide; j = 0_wide
      k = max(nn, key_length)

      do kk = 1, k
        rn%mt(i+1) = ieor(rn%mt(i+1), c1 * ieor(rn%mt(i), ishft(rn%mt(i), -62))) &
                    + init_key(j+1) + j
        i = i+1; j = j+1
        if(i >= nn) then
          rn%mt(1) = rn%mt(nn)
          i = 1
        end if
        if(j >= key_length) j = 0
      end do

      do kk = 1, nn-1
        rn%mt(i+1) = ieor(rn%mt(i+1), c2 * ieor(rn%mt(i), ishft(rn%mt(i), -62))) - i
        i = i+1
        if(i >= nn) then
          rn%mt(1) = rn%mt(nn)
          i = 1
        end if
      end do

      rn%mt(1) = ishft(1_wide, 63)  ! MSB is 1; assuring non-zero initial array

    end subroutine init_by_array64
    !!***

    !!****f* rng/genrand64_int64 *
    !!
    !!  NAME 
    !!   genrand64_int64
    !!  PURPOSE
    !!   Generates random integer on the interval [-2^63,2^63-1]
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/16
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    integer(double) function genrand64_int64(rn)

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)  :: rn

      ! Local variables
      integer(wide) :: mag01(0:1) = (/0_wide, matrix_a/)
      integer(wide) :: x
      integer     :: i

      if(rn%mti >= nn) then ! generate nn words at one time

        ! if init_genrand64() has not been called, a default initial seed is used
        if(rn%mti == nn+1) call rn%init_genrand64(seed_def)

        do i = 1, nn-mm
          x = ior(iand(rn%mt(i),um), iand(rn%mt(i+1), lm))
          rn%mt(i) = ieor(ieor(rn%mt(i+mm), ishft(x, -1)), mag01(iand(x, 1_wide)))
        end do

        do i = nn-mm+1, nn-1
          x = ior(iand(rn%mt(i), um), iand(rn%mt(i+1), lm))
          rn%mt(i) = ieor(ieor(rn%mt(i+mm-nn), ishft(x, -1)), mag01(iand(x, 1_wide)))
        end do

        x = ior(iand(rn%mt(nn), um), iand(rn%mt(1), lm))
        rn%mt(nn) = ieor(ieor(rn%mt(mm), ishft(x, -1)), mag01(iand(x, 1_wide)))

        rn%mti = 0

      end if

      rn%mti = rn%mti + 1
      x = rn%mt(rn%mti)

      x = ieor(x, iand(ishft(x,-29), 6148914691236517205_wide))
      x = ieor(x, iand(ishft(x, 17), 8202884508482404352_wide))
      x = ieor(x, iand(ishft(x, 37),   -2270628950310912_wide))
      x = ieor(x, ishft(x, -43))

      genrand64_int64 = x

    end function genrand64_int64
    !!***

    !!****f* rng/genrand64_real1 *
    !!
    !!  NAME 
    !!   genrand64_real1
    !!  PURPOSE
    !!   Generates a random real number on the interval [0,1]
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/16
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    real(double) function genrand64_real1(rn)

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)  :: rn

        genrand64_real1 = real(ishft(rn%genrand64_int64(), -11), kind=double) * pi253_1

    end function genrand64_real1

    !!****f* rng/genrand64_real2 *
    !!
    !!  NAME 
    !!   genrand64_real2
    !!  PURPOSE
    !!   Generates a random real number on the interval [0,1)
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/16
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    real(double) function genrand64_real2(rn)

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)  :: rn

        genrand64_real2 = real(ishft(rn%genrand64_int64(), -11), kind=double) * pi253

    end function genrand64_real2
    !!***

    !!****f* rng/genrand64_real3 *
    !!
    !!  NAME 
    !!   genrand64_real3
    !!  PURPOSE
    !!   Generates a random real number on the interval (0,1)
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/16
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    real(double) function genrand64_real3(rn)

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)  :: rn

      genrand64_real3 = real(ishft(rn%genrand64_int64(), -12), kind=double)
      genrand64_real3 = (genrand64_real3 + 0.5_double) * pi252

    end function genrand64_real3
    !!***

    !!****f* rng/get_random_seed *
    !!
    !!  NAME 
    !!   get_random_seed
    !!  PURPOSE
    !!   Generates a random real number on the interval (0,1)
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/16
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    subroutine get_random_seed(rn)

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)  :: rn

      ! Local variables
      integer               :: i, un, istat, dt(8), pid, t(2), s, n
      integer(8)            :: count, tms

      call random_seed(size=n)
      allocate(rn%seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         read(un) rn%seed
         close(un)
      else
         ! Fallback to XOR:ing the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         call system_clock(count)
         if (count /= 0) then
            t = transfer(count, t)
         else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                 + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                 + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                 + dt(5) * 60 * 60 * 1000 &
                 + dt(6) * 60 * 1000 + dt(7) * 1000 &
                 + dt(8)
            t = transfer(tms, t)
         end if
         s = ieor(t(1), t(2))
         pid = getpid() + 1099279 ! Add a prime
         s = ieor(s, pid)
         if (n >= 3) then
            rn%seed(1) = t(1) + 36269
            rn%seed(2) = t(2) + 72551
            rn%seed(3) = pid
            if (n > 3) then
               rn%seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
         else
            rn%seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
         end if
      end if

    end subroutine get_random_seed
    !!***

    !!****f* rng/init_rng *
    !!
    !!  NAME 
    !!   init_rng
    !!  PURPOSE
    !!   Conquest wrapper for seed generator and MT rng. Initialise using
    !!   int_seed if provided, otherwise use get_random_seed
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/16
    !!  MODIFICATION HISTORY
    !!   2019/05/21 zamaan
    !!    Now takes rng_seed from global_module, generates a random seed array
    !     if rng_seed < 0 (default = -1)
    !!  SOURCE
    !!
    subroutine init_rng(rn)

      use global_module, only: rng_seed

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)    :: rn

      if (rng_seed > 0)  then
        call init_genrand64(rn, int(rng_seed,wide))
      else
        call rn%get_random_seed
        call rn%init_by_array64(rn%seed)
      end if
    end subroutine init_rng
    !!***

    !!****f* rng/init_normal *
    !!
    !!  NAME 
    !!   init_normal
    !!  PURPOSE
    !!   Initialise generator for normal distribution
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/21
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    subroutine init_normal(rn, sigma, mu)

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)  :: rn
      real(double), intent(in)        :: sigma, mu

      rn%buffered = .false.
      rn%sigma = sigma
      rn%mu = mu

    end subroutine init_normal
    !!***

    !!****f* rng/rng_uniform *
    !!
    !!  NAME 
    !!   rng_uniform
    !!  PURPOSE
    !!   Conquest wrapper for real random number generator.
    !!  USAGE
    !!   interval 1 (default): [0,1], interval 2: [0,1), interval 3: (0,1)
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/16
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    real(double) function rng_uniform(rn, interval)

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)    :: rn
      integer, optional, intent(in)     :: interval

      ! Local variables
      real(double)                      :: num

      if (present(interval)) then
        select case(interval)
        case(1)
          num = rn%genrand64_real1()
        case(2)
          num = rn%genrand64_real2()
        case(3)
          num = rn%genrand64_real3()
        end select
      else
        num = rn%genrand64_real1()
      end if

      rng_uniform = real(num, double)

    end function rng_uniform
    !!***

    !!****f* move_atoms/rng_normal *
    !!  NAME
    !!   rng_normal
    !!  PURPOSE
    !!   Box-Muller transform to generate pairs of numbers drawn from a 
    !!   Gaussian distribution given a uniform distribution, via 
    !!   Marsaglia's polar method (no sine or cosine calls)
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/17
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    real(double) function rng_normal(rn)

      use numbers

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)  :: rn

      ! Local variables
      real(double)              :: x1, x2, y1, y2, w

      if (rn%buffered) then
        rng_normal = rn%buffer
        rn%buffered = .false.
      else
        do
          x1 = rn%rng_uniform()
          x2 = rn%rng_uniform()
          y1 = two*x1 - one
          y2 = two*x2 - one
          w = y1*y1 + y2*y2
          if (w < 1.0) exit
        end do
        w = sqrt(-two*log(w)/w)
        rng_normal = y1*w*rn%sigma + rn%mu
        rn%buffer  = y2*w*rn%sigma + rn%mu
        rn%buffered = .true.
      end if

    end function rng_normal
    !!***

    !!****f* move_atoms/rng_integer *
    !!  NAME
    !!   rng_integer
    !!  PURPOSE
    !!   Generate random numbers in the interval [a,b)
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/05/17
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    integer function rng_integer(rn, a, b)

      use datatypes

      implicit none

      ! Passed variables
      class(type_rng), intent(inout)  :: rn
      integer, intent(in)             :: a, b

      ! Local variables
      integer :: rint

      rint = modulo(rn%genrand64_int64(), (b-a))
      rng_integer = int(rint) + a

    end function  rng_integer
    !!***

end module rng
!!***
