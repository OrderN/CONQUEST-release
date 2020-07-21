! Module for output routines
module output

  use datatypes

  implicit none
  
  real(double), dimension(95) :: atrad

contains

  ! Write OpenDX file for charge or current density
  subroutine write_dx_density(ci)

    use datatypes
    use numbers
    use local, ONLY: current, root_file, nptsx, nptsy, nptsz, grid_x, grid_y, grid_z, gpv, &
         nrptx, nrpty, nrptz, nsampx, nsampy, nsampz
    use dimens, ONLY: volume
    use units, only : BohrToAng
    
    implicit none

    character(len=50), OPTIONAL :: ci
    
    character(len=50) :: filename
    integer :: nx, ny, nz, icount, nrx, nry, nrz

    ! Open file
    if(PRESENT(ci)) then
       filename = trim(ci)//".dx"
    else
       filename = trim(root_file)//".dx"
    end if
    open(unit=17,file=filename)
    write(17,fmt='("#Now give details of grid")')
    write(17,fmt='("object 1 class gridpositions counts ",3i4)') nptsz*nrptz,nptsy*nrpty,nptsx*nrptx
    write(17,fmt='("origin   0.0000    0.0000    0.0000")')
    write(17,fmt='("delta",3f10.4)') zero,   zero,   grid_z*BohrToAng
    write(17,fmt='("delta",3f10.4)') zero,   grid_y*BohrToAng, zero
    write(17,fmt='("delta",3f10.4)') grid_x*BohrToAng, zero,   zero
    write(17,*)
    write(17,fmt='("object 2 class gridconnections counts ",3i4)') nptsz*nrptz,nptsy*nrpty,nptsx*nrptx
    write(17,*)
    write(17,fmt='("object 3 class array type float rank 0 items ",i10," data follows")') nptsx*nptsy*nptsz*nrptx*nrpty*nrptz
    icount = 0
    do nrz=1,nrptz
    do nz=1,nptsz,nsampz
       do nry=1,nrpty
       do ny=1,nptsy,nsampy
          do nrx=1,nrptx
          do nx=1,nptsx,nsampx
             icount = icount + 1
             if(mod(icount,5)==0) then
                write(17,fmt='(f15.8)') current(nx,ny,nz)/gpv
             else
                write(17,fmt='(f15.8)',advance='no') current(nx,ny,nz)/gpv
             endif
          end do
          end do
       end do
       end do
    end do
    end do
    write(17,*)
    write(17,fmt='(" object ",a18," class field")') '"electron density"'
    write(17,fmt='("  component ",a11," 1")') '"positions"'
    write(17,fmt='("  component ",a13," 2")') '"connections"'
    write(17,fmt='("  component ",a6," 3")') '"data"'
    !printf("Total electrons: %20.12lf\n Normalised: %20.12lf\n",n_elec,n_elec/vol);
    !printf("Maximum density: %20.12lf\n",a_max);
    close(unit=17)
  end subroutine write_dx_density

  subroutine write_dx_coords(ci)

    use datatypes
    use numbers
    use local, ONLY: root_file, stm_x_min, stm_y_min, stm_z_min
    use dimens, ONLY: r_super_x, r_super_y, r_super_z, atomicnum
    use species_module, ONLY: n_species
    use global_module, ONLY: ni_in_cell, atom_coord, species_glob
    use units, only : BohrToAng
    
    implicit none

    character(len=50), OPTIONAL :: ci
    
    character(len=50) :: filename
    character(len=2) :: ele
    integer :: i, j, totnabs
    integer, dimension(15,ni_in_cell) :: neigh
    integer, dimension(ni_in_cell) :: nabs
    real(double) :: size1, size2, dist2, d2, dx, dy, dz
    
    call assign_atomic_radii
    ! Open file
    if(PRESENT(ci)) then
       filename = trim(ci)//"XYZ.dx"
    else
       filename = trim(root_file)//"XYZ.dx"
    end if
    open(unit=17,file=filename)
    write(17,fmt='("object 1 class array type float rank 0 items",i7," data follows")') ni_in_cell
    do i=1,ni_in_cell
       write(17,fmt='(f8.4)') atrad(atomicnum(species_glob(i)))
    end do
    write(17,fmt='("attribute ""dep"" string ""positions"" ")')
    write(17,fmt='("object 2 class array type float rank 1 shape 3 items",i7," data follows")') ni_in_cell
    nabs = 0
    neigh = 0
    totnabs = 0
    do i=1,ni_in_cell
       write(17,fmt='(3f12.6)') BohrToAng*(atom_coord(1,i)-stm_x_min),BohrToAng*(atom_coord(2,i)-stm_y_min),&
            BohrToAng*(atom_coord(3,i)-stm_z_min)
       size1 = atrad(atomicnum(species_glob(i)))
       ! Find neighbours
       do j=i+1,ni_in_cell
          size2 = atrad(atomicnum(species_glob(j)))
          dist2 = (size1*size1 + size2*size2)
          dx = BohrToAng*(atom_coord(1,i) - atom_coord(1,j))
          dy = BohrToAng*(atom_coord(2,i) - atom_coord(2,j))
          dz = BohrToAng*(atom_coord(3,i) - atom_coord(3,j))
          d2 = dx*dx + dy*dy + dz*dz
          if(d2<dist2) then
             nabs(i) = nabs(i) + 1
             nabs(j) = nabs(j) + 1
             totnabs = totnabs + 2 ! i->j and j->i
             neigh(nabs(i),i) = j
             neigh(nabs(j),j) = i
          end if
       end do
    end do
    write(17,fmt='("attribute ""dep"" string ""positions"" ")')
    write(17,fmt='("object 3 class array type int rank 1 shape 2 items ",i7," data follows")') totnabs
    do i=1,ni_in_cell
       do j=1,nabs(i)
          write(17,fmt='(2i7)') i-1,neigh(j,i)-1 ! Indexing from 0
       end do
    end do
    write(17,fmt='("attribute ""element type"" string ""lines"" ",/,"attribute ""ref"" string ""positions"" ")')
    write(17,fmt='("#",/,"object ""molecule"" class field")')
    write(17,fmt='("component ""data"" value 1",/,"component ""positions"" value 2")')
    write(17,fmt='("component ""connections"" value 3",/,"attribute ""name"" string ""molecule"" ")')
    write(17,fmt='("#",/,"end")')
    close(unit=17)
  end subroutine write_dx_coords
  
  subroutine write_cube(ci)

    use datatypes
    use numbers
    use local, ONLY: root_file, stm_x_min, stm_y_min, stm_z_min, nsampx, nsampy, nsampz
    use dimens, ONLY: r_super_x, r_super_y, r_super_z, volume, atomicnum
    use local, ONLY: current, root_file, nptsx, nptsy, nptsz, grid_x, grid_y, grid_z, gpv, nrptx, nrpty, nrptz
    use global_module, ONLY: ni_in_cell, atom_coord, species_glob
    use units, ONLY: BohrToAng, AngToBohr
    
    implicit none

    character(len=50), OPTIONAL :: ci
    
    character(len=50) :: filename

    integer :: i, j, icount, nrx, nry, nrz, nx, ny, nz, isym
    logical :: SymSamp
    character(len=3), dimension(7) :: asymm = (/'X  ','Y  ','XY ','Z  ','XZ ','YZ ','XYZ'/)
    real(double) :: shiftx, shifty, shiftz

    ! Check sampling points are symmetric or not
    SymSamp = .true.
    isym = 0
    if (mod(nptsx-1,nsampx).ne.0) isym = 1
    if (mod(nptsy-1,nsampy).ne.0) isym = isym + 2
    if (mod(nptsz-1,nsampz).ne.0) isym = isym + 4
    if (isym.ne.0) then
       SymSamp = .false.
       write(*,*) 'Warning: sampling points not symmetric for ',asymm(isym)
    end if
    ! Open file
    if(PRESENT(ci)) then
       filename = trim(ci)//".cube"
    else
       filename = trim(root_file)//".cube"
    end if
    open(unit=17,file=filename)
    ! Header - improve this later
    write(17,fmt='("Conquest charge")')
    if(SymSamp) then
       write(17,fmt='("Sampling points symmetric")')
    else
       write(17,fmt='("Sampling points asymmetric in directions ",a2)') asymm(isym)
    end if
    ! Atoms, origin location
    write(17,fmt='(i5,3f12.6)') ni_in_cell*nrptx*nrpty*nrptz, zero, zero, zero
    ! Numbers of grid points, grid increment (i.e. spacing)
    write(17,fmt='(i5,3f12.6)') nptsx*nrptx/nsampx,grid_x,zero,zero
    write(17,fmt='(i5,3f12.6)') nptsy*nrpty/nsampy,zero,grid_y,zero
    write(17,fmt='(i5,3f12.6)') nptsz*nrptz/nsampz,zero,zero,grid_z
    ! Atomic coordinates
    do nrx = 1,nrptx
       shiftx = real(nrx-1,double)*r_super_x
       do nry = 1,nrpty
          shifty = real(nry-1,double)*r_super_y
          do nrz = 1,nrptz
             shiftz = real(nrz-1,double)*r_super_z
             do i=1,ni_in_cell
                write(17,fmt='(i5,4f12.6)') atomicnum(species_glob(i)),real(atomicnum(species_glob(i)),double),&
                     atom_coord(1,i)-stm_x_min + shiftx, &
                     atom_coord(2,i)-stm_y_min + shifty, &
                     atom_coord(3,i)-stm_z_min + shiftz
                !write(17,fmt='(i5,4f12.6)') species_glob(i),zero,(atom_coord(j,i)*BohrToAng,j=1,3)
             end do
          end do
       end do
    end do
    ! Charge density, z fastest
    icount = 0
    do nrx=1,nrptx
       do nx=1,nptsx,nsampx
          do nry=1,nrpty
             do ny=1,nptsy,nsampy
                do nrz=1,nrptz
                   do nz=1,nptsz,nsampz
                      icount = icount + 1
                      if(mod(icount,6)==0) then
                         write(17,fmt='(e13.5)') current(nx,ny,nz)/gpv
                      else
                         write(17,fmt='(e13.5)',advance='no') current(nx,ny,nz)/gpv
                      endif
                   end do
                end do
             end do
          end do
       end do
    end do
    close(unit=17)
  end subroutine write_cube

  ! Taken from density_module.f90
  subroutine assign_atomic_radii

    use datatypes
    use numbers, only: zero
    use units, only: AngToBohr

    implicit none

    atrad = zero
    atrad(1)=0.35_double*AngToBohr
    atrad(2)=0.28_double*AngToBohr

    atrad(3)=1.45_double*AngToBohr ! Li
    atrad(4)=1.05_double*AngToBohr
    atrad(5)=0.85_double*AngToBohr
    atrad(6)=0.70_double*AngToBohr
    atrad(7)=0.65_double*AngToBohr
    atrad(8)=0.60_double*AngToBohr
    atrad(9)=0.50_double*AngToBohr ! F
    atrad(10)=0.58_double*AngToBohr

    atrad(11)=1.80_double*AngToBohr ! Na
    atrad(12)=1.50_double*AngToBohr
    atrad(13)=1.25_double*AngToBohr
    atrad(14)=1.10_double*AngToBohr
    atrad(15)=1.00_double*AngToBohr
    atrad(16)=1.00_double*AngToBohr
    atrad(17)=1.00_double*AngToBohr ! Cl
    atrad(18)=1.06_double*AngToBohr

    atrad(19)=2.20_double*AngToBohr ! K
    atrad(20)=1.80_double*AngToBohr
    atrad(21)=1.60_double*AngToBohr
    atrad(22)=1.40_double*AngToBohr
    atrad(23)=1.35_double*AngToBohr
    atrad(24)=1.40_double*AngToBohr
    atrad(25)=1.40_double*AngToBohr
    atrad(26)=1.40_double*AngToBohr ! Fe
    atrad(27)=1.35_double*AngToBohr
    atrad(28)=1.35_double*AngToBohr
    atrad(29)=1.35_double*AngToBohr
    atrad(30)=1.35_double*AngToBohr
    atrad(31)=1.30_double*AngToBohr ! Ga
    atrad(32)=1.25_double*AngToBohr
    atrad(33)=1.15_double*AngToBohr
    atrad(34)=1.15_double*AngToBohr
    atrad(35)=1.15_double*AngToBohr ! Br
    atrad(36)=1.16_double*AngToBohr

    atrad(37)=2.35_double*AngToBohr ! Rb
    atrad(38)=2.00_double*AngToBohr
    atrad(39)=1.80_double*AngToBohr
    atrad(40)=1.55_double*AngToBohr ! Zr
    atrad(41)=1.45_double*AngToBohr
    atrad(42)=1.45_double*AngToBohr
    atrad(43)=1.35_double*AngToBohr
    atrad(44)=1.30_double*AngToBohr ! Ru
    atrad(45)=1.35_double*AngToBohr
    atrad(46)=1.40_double*AngToBohr
    atrad(47)=1.60_double*AngToBohr ! Ag
    atrad(48)=1.55_double*AngToBohr
    atrad(49)=1.55_double*AngToBohr ! In
    atrad(50)=1.45_double*AngToBohr
    atrad(51)=1.45_double*AngToBohr
    atrad(52)=1.40_double*AngToBohr
    atrad(53)=1.40_double*AngToBohr ! I
    atrad(54)=1.40_double*AngToBohr

    atrad(55)=2.60_double*AngToBohr ! Cs
    atrad(56)=2.15_double*AngToBohr
    atrad(57)=1.95_double*AngToBohr ! La
    atrad(58)=1.85_double*AngToBohr
    atrad(59)=1.85_double*AngToBohr
    atrad(60)=1.85_double*AngToBohr
    atrad(61)=1.85_double*AngToBohr
    atrad(62)=1.85_double*AngToBohr
    atrad(63)=1.85_double*AngToBohr ! Eu
    atrad(64)=1.80_double*AngToBohr
    atrad(65)=1.75_double*AngToBohr
    atrad(66)=1.75_double*AngToBohr
    atrad(67)=1.75_double*AngToBohr
    atrad(68)=1.75_double*AngToBohr
    atrad(69)=1.75_double*AngToBohr
    atrad(70)=1.75_double*AngToBohr
    atrad(71)=1.75_double*AngToBohr ! Lu
    atrad(72)=1.55_double*AngToBohr
    atrad(73)=1.45_double*AngToBohr
    atrad(74)=1.35_double*AngToBohr ! W
    atrad(75)=1.35_double*AngToBohr
    atrad(76)=1.30_double*AngToBohr
    atrad(77)=1.35_double*AngToBohr ! Ir
    atrad(78)=1.35_double*AngToBohr
    atrad(79)=1.35_double*AngToBohr
    atrad(80)=1.50_double*AngToBohr ! Hg
    atrad(81)=1.90_double*AngToBohr
    atrad(82)=1.80_double*AngToBohr
    atrad(83)=1.60_double*AngToBohr
    atrad(84)=1.90_double*AngToBohr ! Po
    atrad(85)=1.50_double*AngToBohr
    atrad(86)=1.50_double*AngToBohr
    
    atrad(87)=2.60_double*AngToBohr
    atrad(88)=2.15_double*AngToBohr ! Ra
    atrad(89)=1.95_double*AngToBohr
    atrad(90)=1.80_double*AngToBohr
    atrad(91)=1.80_double*AngToBohr
    atrad(92)=1.75_double*AngToBohr
    atrad(93)=1.75_double*AngToBohr
    atrad(94)=1.75_double*AngToBohr
    atrad(95)=1.75_double*AngToBohr
    return
  end subroutine assign_atomic_radii

  subroutine write_banner

    use datestamp, ONLY: datestr, commentver

    implicit none

    character(len=10) :: today, the_time

    write(*,fmt='(/"CONQUEST charge density conversion and STM image simulation"/)')
    write(*,fmt='("D. R. Bowler (UCL) and A. Nakata (NIMS)")')
    call date_and_time(today, the_time)
    write(*,fmt='(/4x,"This job was run on ",a4,"/",a2,"/",a2," at ",a2,":",a2,/)') &
         today(1:4), today(5:6), today(7:8), the_time(1:2), the_time(3:4)
    write(*,&
          '(/4x,"Code compiled on: ",a,/10x,"Version comment: ",/6x,a//)') &
         datestr, commentver
  end subroutine write_banner

end module output
