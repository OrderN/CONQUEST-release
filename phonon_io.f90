! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module phonon_io
! ------------------------------------------------------------------------------
! Code area XX: spectroscopy
! ------------------------------------------------------------------------------

!!***h* Conquest/phonon_io *
!!  NAME
!!   phonon_io
!!
!!  CREATION DATE
!!   2016/28/06 Z. Raza and L. A. Truflandier
!!
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module phonon_io
  
  use datatypes

  use global_module,             ONLY: area_basis, io_lun, numprocs
  use timer_module,              ONLY: start_timer,stop_timer
  use io_module,                 ONLY: get_file_name
  use input_module,              ONLY: io_assign, io_close
  use GenComms,                  ONLY: inode, ionode  
  use dimens,                    ONLY: r_super_x, r_super_y, r_super_z
  use timer_module,              ONLY: cq_timer
  
  implicit none

  integer :: unit_global_write
  type(cq_timer), save :: tmr_std_ph_write
  
contains
  
  subroutine phonon_global_write()

    use datatypes
    use numbers,                ONLY: zero
    use matrix_data,            ONLY: mat, Hrange, Srange, Xrange, SXrange
    use cover_module,           ONLY: BCS_parts
    use primary_module,         ONLY: bundle
    use global_module,          ONLY: io_lun
    use species_module,         ONLY: species_label
    use support_spec_format,    ONLY: supports_on_atom, flag_one_to_one
    use atomic_density,         ONLY: atomic_density_table   

    use global_module,          ONLY: id_glob ! new
    use group_module,           ONLY: parts   ! new

    implicit none

    ! << Local variables >>
    real(double), dimension(3) :: ip_xyz, nb_xyz
    real(double), dimension(3) :: xyz_ipnb
    real(double)               ::   r_ipnb

    real(double), dimension(3) :: xyz_ghost = 0.0d0
    real(double)               ::   r_ghost = 0.0d0
    real(double)               :: alpha_deg = 90.0d0
    real(double)               :: beta_deg  = 90.0d0
    real(double)               :: gamma_deg = 90.0d0

    integer                    :: ist, npc, nic, ip
    integer                    :: nb_spec, ip_spec
    integer                    :: ip_lab_cell, unit
    character(len=2)           :: nb_name, ip_name
    character(len=23)          :: filename
    character(len=23)          :: filename1, filename2, filename3
    character(len=23)          :: filenamexyz1, filenamexyz2, filenamexyz3 
    character(len=23)          :: filenamecif1, filenamecif2, filenamecif3 

    integer                    :: part, memb, ip_num, iprim, ip_nsup, neigh ! new
    integer                    :: nb_global_part, nb_global_num, nb_nsup    ! new
    integer                    :: i, j, k, n_atoms

    real(double)               :: ip_radi, nb_radi
    integer                    :: unit_write_xyz1, unit_write_xyz2, unit_write_xyz3 
    integer                    :: unit_write_cif1, unit_write_cif2, unit_write_cif3 
    integer                    :: toto
    
    !call start_timer(tmr_std_ph_write)

    call get_file_name('proc_partitions' ,numprocs,inode,filename)
    call get_file_name('proc_atoms_stnd' ,numprocs,inode,filename1)
    call get_file_name('proc_atoms_part' ,numprocs,inode,filename2)
    call get_file_name('proc_atoms_proc' ,numprocs,inode,filename3)

    filenamexyz1 = trim(filename1)//'.xyz'
    filenamexyz2 = trim(filename2)//'.xyz'
    filenamexyz3 = trim(filename3)//'.xyz'

    filenamecif1 = trim(filename1)//'.cif'
    filenamecif2 = trim(filename2)//'.cif'
    filenamecif3 = trim(filename3)//'.cif'


    call io_assign(unit_global_write)
    unit=unit_global_write
    open(unit,file=filename)

    call io_assign(unit_write_xyz1)
    open(unit_write_xyz1,file=filenamexyz1)

    call io_assign(unit_write_xyz2)
    open(unit_write_xyz2,file=filenamexyz2)

    call io_assign(unit_write_xyz3)
    open(unit_write_xyz3,file=filenamexyz3)

    call io_assign(unit_write_cif1)
    open(unit_write_cif1,file=filenamecif1)

    call io_assign(unit_write_cif2)
    open(unit_write_cif2,file=filenamecif2)

    call io_assign(unit_write_cif3)
    open(unit_write_cif3,file=filenamecif3)

    call write_info(unit,inode,10,bundle%groups_on_node,0,0,0,0, &
         0,0,'XX',xyz_ghost,r_ghost,0,0,0,zero,zero,0)

    n_atoms = 0
    ! count atoms on node
    do part = 1,bundle%groups_on_node ! previously np
       n_atoms = n_atoms + bundle%nm_nodgroup(part)
    end do
    write(unit_write_xyz1,'(4X,I8)') n_atoms
    write(unit_write_xyz1,*)
    write(unit_write_xyz2,'(4X,I8)') n_atoms
    write(unit_write_xyz2,*)
    write(unit_write_xyz3,'(4X,I8)') n_atoms
    write(unit_write_xyz3,*)

    write(unit_write_cif1,*) trim('data_Conquest')
    write(unit_write_cif1,*) '_symmetry_space_group_name_H-M P1'
    write(unit_write_cif1,*) '_symmetry_Int_Tables_number    1'
    write(unit_write_cif1,*) '_cell_length_a', r_super_x
    write(unit_write_cif1,*) '_cell_length_b', r_super_y
    write(unit_write_cif1,*) '_cell_length_c', r_super_z
    write(unit_write_cif1,*) '_cell_angle_alpha',  alpha_deg
    write(unit_write_cif1,*) '_cell_angle_beta ',  beta_deg
    write(unit_write_cif1,*) '_cell_angle_gamma',  gamma_deg
    write(unit_write_cif1,*) 'loop_'
    write(unit_write_cif1,*) '_atom_site_label'
    write(unit_write_cif1,*) '_atom_site_type_symbol'
    write(unit_write_cif1,*) '_atom_site_fract_x'
    write(unit_write_cif1,*) '_atom_site_fract_y'
    write(unit_write_cif1,*) '_atom_site_fract_z'


    write(unit_write_cif2,*) 'data_Conquest'
    write(unit_write_cif2,*) '_symmetry_space_group_name_H-M P1'
    write(unit_write_cif2,*) '_symmetry_Int_Tables_number    1'
    write(unit_write_cif2,*) '_cell_length_a', r_super_x
    write(unit_write_cif2,*) '_cell_length_b', r_super_y
    write(unit_write_cif2,*) '_cell_length_c', r_super_z
    write(unit_write_cif2,*) '_cell_angle_alpha',  alpha_deg
    write(unit_write_cif2,*) '_cell_angle_beta ',  beta_deg
    write(unit_write_cif2,*) '_cell_angle_gamma',  gamma_deg
    write(unit_write_cif2,*) 'loop_'
    write(unit_write_cif2,*) '_atom_site_label'
    write(unit_write_cif2,*) '_atom_site_type_symbol'
    write(unit_write_cif2,*) '_atom_site_fract_x'
    write(unit_write_cif2,*) '_atom_site_fract_y'
    write(unit_write_cif2,*) '_atom_site_fract_z'


    write(unit_write_cif3,*) 'data_Conquest'
    write(unit_write_cif3,*) '_symmetry_space_group_name_H-M P1'
    write(unit_write_cif3,*) '_symmetry_Int_Tables_number    1'
    write(unit_write_cif3,*) '_cell_length_a', r_super_x
    write(unit_write_cif3,*) '_cell_length_b', r_super_y
    write(unit_write_cif3,*) '_cell_length_c', r_super_z
    write(unit_write_cif3,*) '_cell_angle_alpha',  alpha_deg
    write(unit_write_cif3,*) '_cell_angle_beta ',  beta_deg
    write(unit_write_cif3,*) '_cell_angle_gamma',  gamma_deg
    write(unit_write_cif3,*) 'loop_'
    write(unit_write_cif3,*) '_atom_site_label'
    write(unit_write_cif3,*) '_atom_site_type_symbol'
    write(unit_write_cif3,*) '_atom_site_fract_x'
    write(unit_write_cif3,*) '_atom_site_fract_y'
    write(unit_write_cif3,*) '_atom_site_fract_z'

    toto = 0
    ! Loop over primary set partition <part>
    part_loop: do part = 1,bundle%groups_on_node ! previously np

       call write_info(unit,inode,20,0,parts%ngnode(parts%inode_beg(inode)+part-1), &
            bundle%nm_nodgroup(part),0,0,0,0,'XX',xyz_ghost,r_ghost,0,0,0,zero,zero,0)

       toto = toto + 1
   
       ! must be members
       check_members: if(bundle%nm_nodgroup(part) > 0) then                           
          ! for each member <memb> of the partition <part>
          iprim = 0
          primary_loop: do memb = 1, bundle%nm_nodgroup(part) ! previously ip             

             ip_num = bundle%nm_nodbeg(part) + memb - 1 ! new
             iprim  = iprim + 1                         ! new
             ip     = bundle%ig_prim(iprim)             ! new

             ip_spec     = bundle%species(ip_num)
             ip_name     = species_label(bundle%species(ip_num))
             ip_radi     = atomic_density_table(ip_spec)%cutoff
             ip_lab_cell = ip

             if (flag_one_to_one) then
                ip_nsup = supports_on_atom(ip_spec)%nsuppfuncs
             else
                ip_nsup = supports_on_atom(ip)%nsuppfuncs
             end if
             !end if

             ip_xyz(1) = bundle%xprim(ip_num)
             ip_xyz(2) = bundle%yprim(ip_num)
             ip_xyz(3) = bundle%zprim(ip_num)             

             call write_info(unit,inode,31,0,0,0,memb,supports_on_atom(ip)%nsuppfuncs, &
                  Xrange,mat(part,Xrange)%n_nab(memb),ip_name,ip_xyz,r_ghost, &
                  ip_spec,ip,ip_lab_cell,ip_radi,zero,0)

             call write_info(unit_write_xyz1,inode,32,0,parts%ngnode(parts%inode_beg(inode)+part-1),0,&
                  memb,supports_on_atom(ip)%nsuppfuncs, &
                  Xrange,mat(part,Xrange)%n_nab(memb),ip_name,ip_xyz,r_ghost, &
                  ip_spec,ip,ip_lab_cell,ip_radi,zero,0)

             call write_info(unit_write_xyz2,inode,33,0,toto,0,&
                  memb,supports_on_atom(ip)%nsuppfuncs, &
                  Xrange,mat(part,Xrange)%n_nab(memb),ip_name,ip_xyz,r_ghost, &
                  ip_spec,ip,ip_lab_cell,ip_radi,zero,0)

             call write_info(unit_write_xyz3,inode,34,0,parts%ngnode(parts%inode_beg(inode)+part-1),0,&
                  memb,supports_on_atom(ip)%nsuppfuncs, &
                  Xrange,mat(part,Xrange)%n_nab(memb),ip_name,ip_xyz,r_ghost, &
                  ip_spec,ip,ip_lab_cell,ip_radi,zero,0)


             call write_info(unit_write_cif1,inode,42,0,parts%ngnode(parts%inode_beg(inode)+part-1),0,&
                  memb,supports_on_atom(ip)%nsuppfuncs, &
                  Xrange,mat(part,Xrange)%n_nab(memb),ip_name,ip_xyz,r_ghost, &
                  ip_spec,ip,ip_lab_cell,ip_radi,zero,0)

             call write_info(unit_write_cif2,inode,43,0,toto,0,&
                  memb,supports_on_atom(ip)%nsuppfuncs, &
                  Xrange,mat(part,Xrange)%n_nab(memb),ip_name,ip_xyz,r_ghost, &
                  ip_spec,ip,ip_lab_cell,ip_radi,zero,0)

             call write_info(unit_write_cif3,inode,44,0,parts%ngnode(parts%inode_beg(inode)+part-1),0,&
                  memb,supports_on_atom(ip)%nsuppfuncs, &
                  Xrange,mat(part,Xrange)%n_nab(memb),ip_name,ip_xyz,r_ghost, &
                  ip_spec,ip,ip_lab_cell,ip_radi,zero,0)

          end do primary_loop
       end if check_members
    end do part_loop
    
    call write_info(unit,inode,0,0,0,0,0,0, &
         0,0,'XX',xyz_ghost,r_ghost,0,0,0,zero,zero,0)
    !
    call io_close(unit)
    !
    call io_close(unit_write_xyz1)
    call io_close(unit_write_xyz2)
    call io_close(unit_write_xyz3)
    !
    call io_close(unit_write_cif1)
    call io_close(unit_write_cif2)
    call io_close(unit_write_cif3)
    !
    !call stop_timer(tmr_std_ph_write,.true.)
    !
    return
  end subroutine phonon_global_write
  
  subroutine write_info(unit,inode,in,npartitions,partition,nmembers,member,nsuppfuncs, &
       range,nneigh,name,xyz,r,spec,atom,labcell,radii,e_hf,dim,matHF,matK)
    

    use datatypes
    use units,                  ONLY: BohrToAng
    use numbers,                ONLY: very_small
    use matrix_data,            ONLY: rcut, mat_name

    implicit none

    ! << Passed variables >>
    real(double), dimension(3), intent(in) :: xyz
    real(double),               intent(in) :: r
    real(double),               intent(in) :: e_hf
    real(double),               intent(in) :: radii

    integer,                    intent(in) :: unit
    integer,                    intent(in) :: inode
    integer,                    intent(in) :: in

    integer,                    intent(in) :: dim
    real(double), optional,     intent(in) :: matHF(dim,dim)
    real(double), optional,     intent(in) :: matK(dim,dim)

    integer,                    intent(in) :: npartitions, partition
    integer,                    intent(in) :: nmembers, member
    integer,                    intent(in) :: nsuppfuncs, nneigh
    integer,                    intent(in) :: atom, spec, labcell
    integer,                    intent(in) :: range

    character(len=2),           intent(in) :: name
    character(len=20)                      :: filename
    
    ! << Local variables >>
    real(double), allocatable, dimension(:,:) :: mat
    real(double), dimension(3) :: xyz_Ang
    real(double)               :: r_Ang

    integer                    :: i, j, k, l    
    CHARACTER(len=2)                :: pte(103)

    DATA pte /  "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
         "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", &
         "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
         "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", &
         "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "St", &
         "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
         "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
         "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
         "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
         "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", &
         "Md", "No", "Lr"/
    
    r_Ang   = r!*BohrtoAng
    xyz_Ang = xyz*BohrtoAng

    !call start_timer(tmr_std_ph_write)

    select case(in)       
    case(1)   
       write(unit,*) 
       write(unit,'(30X,A8,2X,A8,4X,A,11X,A,11X,A,11X,A,8X,A5,X,A7,2X,A3,3X,A5,6X,A3)') &
            '  type  ','lab.cell','x','y','z','r','spec.','neighb.','nsf','radii','exx'
       write(unit,'(A9,X,A,I8,A2,I3,A,6X,2X,A2,X,I8,X,4F12.4,3X,I3,2I8,X,F7.3,F12.8)') &
            '{i\alpha}','{',labcell,'\ ',nsuppfuncs,'}', name, labcell, &
            xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), r, spec, atom, nsuppfuncs, radii,e_hf
       
    case(2)       
       write(unit,'(A9,X,A,I8,A2,I3,A,6X,2X,A2,X,I8,X,4F12.4,3X,I3,2I8,X,F7.3,F14.8)') &
            '{j\beta} ','{',labcell,'\ ',nsuppfuncs,'}', name, labcell, &
            xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), r, spec, atom, nsuppfuncs, radii,e_hf
       write(unit,*) 
       
    case(3)
       write(unit,'(2X,A9,X,A,I8,A2,I3,A,4X,2X,A2,X,I8,X,4F12.4,3X,I3,2I8,X,F7.3)') &
            '{k\gamma}','{',labcell,'\ ',nsuppfuncs,'}', name, labcell, &
            xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), r, spec, atom, nsuppfuncs, radii
       
    case(4)
       write(unit,'(4X,A9,X,A,I8,A2,I3,A,2X,2X,A2,X,I8,X,4F12.4,3X,I3,2I8,X,F7.3)') &
            '{l\delta}','{',labcell,'\ ',nsuppfuncs,'}', name, labcell, &
            xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), r, spec, atom, nsuppfuncs, radii
       
    case(5)
       write(unit,'(8X,A20,X,A,2I3,A,2I3,A,X,F12.8)') 'K{k-\gamma,l-\delta}', &
            '{', npartitions,partition,',',nmembers,member,'}', r       

    case(10)
       open(unit)
       write(unit,'(X,104A)', advance='no') ('=', j=1,104) ; write(unit,'(A)') ' '
       write(unit,'(A15,I8)') '| inode:        ', inode
       write(unit,'(A15,I8)') '| n_partitions: ', npartitions
       write(unit,'(X,104A)', advance='no') ('=', j=1,104) ; write(unit,'(A)') ' '

    case(20)
       write(unit,'(X,104A)', advance='no') ('=', j=1,104) ; write(unit,'(A)') ' '
       write(unit,'(A12,I8,A12,I8)')   ' partition: ', partition,' n_members: ', &
            nmembers
       write(unit,'(X,104A)', advance='no') ('=', j=1,104) ; write(unit,'(A)') ' '

    case(30)   
       write(unit,'(A4,A8,F12.4,X,A)') adjustr(mat_name(range)),'range:  ', &
            rcut(range)*BohrtoAng, 'A'

       !write(unit,'(A12,I8,A12,I8,A16,I8)')'   member:  ', member, ' n_supports:  ', &
       !     nsuppfuncs, ' n_neighbours:  ', nneigh       
       !write(unit,'(A12,6X,A2,4X,3F12.4,X,I3,3I8)') '   member:  ', name, &
       !     xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), spec, atom, labcell,  nsuppfuncs

       write(unit,'(A12,6X,A2,4X,3F12.4,X,I3)') '   member:  ', name, &
            xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), spec
            
       write(unit,'(18X,A,11X,A,11X,A,7X,A11,X,A5,X,A7,X,A8,X,A3,2X,A5)') &
            'x','y','z','distance(A)','spec.','neighb.','lab.cell','nsf','radii'
       !
       !
       !
    case(31)   
       write(unit,'(6X,A2,4X,3F12.4)') name, xyz_Ang(1), xyz_Ang(2), xyz_Ang(3)

    case(32)   
       write(unit,'(6X,A2,4X,3F12.4,4X,I8)') name, xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), partition

    case(33)   
       write(unit,'(6X,A2,4X,3F12.4,4X,I8)') pte(partition), xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), partition

    case(34)   
       write(unit,'(6X,A2,4X,3F12.4,4X,I8)') pte(inode), xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), partition
       !
       !
       !
    case(42) 
       xyz_Ang(1) = xyz_Ang(1)/r_super_x
       xyz_Ang(2) = xyz_Ang(2)/r_super_y
       xyz_Ang(3) = xyz_Ang(3)/r_super_z

       write(unit,'(6X,A2,4X,A2,4X,3F12.4,4X,I8)') name, name, xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), partition

    case(43)   
       xyz_Ang(1) = xyz_Ang(1)/r_super_x
       xyz_Ang(2) = xyz_Ang(2)/r_super_y
       xyz_Ang(3) = xyz_Ang(3)/r_super_z

       write(unit,'(6X,A2,4X,A2,4X,3F12.4,4X,I8)') pte(partition), pte(partition), xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), partition

    case(44)   
       xyz_Ang(1) = xyz_Ang(1)/r_super_x
       xyz_Ang(2) = xyz_Ang(2)/r_super_y
       xyz_Ang(3) = xyz_Ang(3)/r_super_z

       write(unit,'(6X,A2,4X,A2,4X,3F12.4,4X,I8)') pte(inode),  pte(inode), xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), partition
       !
       !
       !
    case(60)
       if (r_Ang > very_small) then
          write(unit,'(6X,A2,4X,4F12.4,X,I3,3I8,F7.3)') name, xyz_Ang(1), xyz_Ang(2), &
               xyz_Ang(3), r_Ang, spec, atom, labcell, nsuppfuncs, radii
       else
          write(unit,'(2X,A3,X,A2,4X,4F12.4,X,I3,3I8,F7.3,X,A3)') '>>>',&
               name, xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), r_Ang, spec, atom, labcell, &
               nsuppfuncs, radii, '<<<'
       end if
    case(70)
       if (r_Ang > very_small) then
          write(unit,'(6X,A2,4X,4F12.4,X,I3,3I8,F16.12)') name, xyz_Ang(1), xyz_Ang(2), &
               xyz_Ang(3), r_Ang, spec, atom, labcell, nsuppfuncs, e_hf
       else
          write(unit,'(2X,A3X,A2,4X,4F12.4,X,I3,3I8,F16.12,X,A3)') '>>>', name, xyz_Ang(1), &
               xyz_Ang(2), xyz_Ang(3), r_Ang, spec, atom, labcell, nsuppfuncs, e_hf, '<<<'
       end if

    case(800)
       write(unit,'(18X,A,11X,A,11X,A,7X,A11,X,A5,X,A7,X,A8,X,A3,2X,A12)') &
            'x','y','z','distance(A)','spec.','neighb.','lab.cell','nsf','Exx(Hartree)'
       
    case(1000)
       write(unit,'(X,104A)', advance='no') ('=', j=1,104) ; write(unit,'(A)') ' '       

    case(0)
       write(unit,'(X,104A)', advance='no') ('=', j=1,104) ; write(unit,'(A)') ' '
       write(unit,'(A15,I8)') '| inode:        ', inode
       write(unit,'(A13)')    '| ...job done'
       write(unit,'(X,104A)', advance='no') ('=', j=1,104) ; write(unit,'(A)') ' '
       close(unit)     

    end select

    call stop_timer(tmr_std_ph_write,.true.)    

    return
  end subroutine write_info  

  subroutine write_matrix(unit,inode,dim,matHF,matK,e_hf2)

    use datatypes
    
    implicit none

    ! << Passed variables >>
    integer,      intent(in) :: unit
    integer,      intent(in) :: inode

    integer,      intent(in) :: dim
    real(double), intent(in) :: matHF(dim,dim)
    real(double), intent(in) :: matK(dim,dim)
    real(double), optional, intent(in) :: e_hf2
    
    ! << Local variables >>
    real(double), allocatable, dimension(:,:) :: mat
    real(double) :: e_hf1

    integer :: i, j, k, l    

    !call start_timer(tmr_std_ph_write)    

    write(unit,'(8X,A17)') 'K matrix elements'
    do i = 1, dim
          write(unit,101) (matK(i,j), j=1,dim)
    end do
    write(unit,*) ' '
    
    write(unit,'(8X,A18)') 'HF matrix elements'
    do i = 1, dim
       write(unit,101) (matHF(i,j), j=1,dim)
    end do
    write(unit,*) ' '
    
    write(unit,'(8X,A20)') 'HF.K matrix elements'
    allocate(mat(dim,dim))
    mat = matmul(matHF,matK)
    do i = 1, dim
       write(unit,101) (mat(i,j), j=1,dim)
    end do
    write(unit,*) ' '
    
    write(unit,'(8X,A20)') 'Tr{HF.K} = e_hf (au)'
    e_hf1 =  sum(mat, mask=reshape(source=(/((i==j,i=1,dim),j=1,dim)/), &
         shape=shape(mat)))
    write(unit,'(4X,2F14.8)') e_hf1, e_hf2
    deallocate(mat)
    write(unit,*) ' '

    call stop_timer(tmr_std_ph_write,.true.)

101 format(4X,100F14.8)
    
    return
  end subroutine write_matrix

  
  subroutine write_tail(unit,inode)  
    
    implicit none

    integer, intent(in) :: unit
    integer, intent(in) :: inode
    integer             :: j
    
    write(unit,'(X,150A)', advance='no') ('=', j=1,150) ; write(unit,'(A)') ' '
    write(unit,'(A15,I8)') '| inode:        ', inode
    write(unit,'(A13)')    '| ...job done'
    write(unit,'(X,150A)', advance='no') ('=', j=1,150) ; write(unit,'(A)') ' '
    
    return
  end subroutine write_tail

  subroutine write_head(unit,inode,npart)  
    
    implicit none

    integer, intent(in) :: unit
    integer, intent(in) :: inode
    integer, intent(in) :: npart
    integer             :: j

    write(unit,'(X,150A)', advance='no') ('=', j=1,150) ; write(unit,'(A)') ' '
    write(unit,'(A15,I8)') '| inode:        ', inode
    write(unit,'(A15,I8)') '| n_partitions: ', npart
    write(unit,'(X,150A)', advance='no') ('=', j=1,150) ; write(unit,'(A)') ' '
    write(unit,'(43X,A8,2X,A8,4X,A,11X,A,11X,A,11X,A,8X,A5,3X,A3,5X,A3,2X,A5,4X,A3)') &
         '  type  ','lab.cell','x','y','z','r','spec.','ind','nsf','radii','exx'

    return
  end subroutine write_head

  subroutine write_part(unit,part,n)  
    
    implicit none

    integer, intent(in) :: unit
    integer, intent(in) :: part
    integer, intent(in) :: n
    integer             :: j

    write(unit,'(X,150A)', advance='no') ('=', j=1,150) ; write(unit,'(A)') ' '
    write(unit,'(A12,I8,A12,I8)') ' partition: ', part,' n_members: ', n
    write(unit,'(X,150A)', advance='no') ('=', j=1,150) ; write(unit,'(A)') ' '
    
    return
  end subroutine write_part


end module phonon_io
