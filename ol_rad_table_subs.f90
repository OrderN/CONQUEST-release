! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module make_rad_tables
! ------------------------------------------------------------------------------
! Area 11: basis functions
! ------------------------------------------------------------------------------

!!****h* Conquest/make_rad_tables *
!!  NAME
!!   make_rad_tables 
!!  PURPOSE
!!   Contains subroutines to create radial tables from paos.
!!  AUTHOR
!!   R. Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
module make_rad_tables
  implicit none

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"
!!***

contains


!!****f* make_rad_tables/make_rad_table_nlpfpao *
!!
!!  NAME 
!!   make_rad_table_nlpfpao
!!  USAGE
!!   make_rad_table_nlpfpao(npts1,l1,delr1,npts2,l2,delr2,table1,table2,count,del_k,k_cut)
!!
!!  PURPOSE
!!  Constructs overlap integral radial tables between NLPF's and PAO functions.
!!
!!
!!  INPUTS
!!   npts1,npts2 : no of points in input tables of NLPF and PAO
!!   l1,delr1,l2,delr2 - ang mom and grid spacings of input tables 1 and 2
!!   del_k - mimimum grid spacing in k-space
!!   k_cut - minimum cut off in k-space
!!   table1,table2 - input tables containing NLPF and PAO info
!!   count - used to index the particular set of overlap integral radial 
!!   tables that are created; count is set by gen_nlpf_supp_tbls
!!  USES
!!   datatypes, ol_int_datatypes, pao_format
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine make_rad_table_nlpfpao(npts1,l1,delr1,npts2,l2,delr2&
       &,table1,table2,count,del_k,k_cut)
    use datatypes
    use numbers, ONLY: zero, twopi
    use ol_int_datatypes !, ONLY : rad_tables_nlpf_pao
    use bessel_integrals !, ONLY : bessloop, maxtwon, complx_fctr
    use cubic_spline_routines !, ONLY : matcharrays
    use GenComms, ONLY: cq_abort

    implicit none

    !code to evaluate radial tables for 
    !overlap integrals of two basis functions
    !29/apr/03 code looks in pretty good shape- RC
    !using this version for a pp_test
    !09/05/03 adding matcharrays subroutine- RC
    !16/05/03 new (matcharrays) grid reprojection subroutines
    !tested and are working nicely.
    !18/09/03 corrected k-space params, added y2 construction too.
    !30/10/03 modifying end of this routine to store nlpf_pao tables 

    real(double), allocatable, dimension(:) :: dummy1,dummy2,dummy1bt,dummy2bt
    real(double), allocatable, dimension(:) :: dummyprod,dumout,ol_out,xin,y2
    real(double), allocatable, dimension(:,:) :: fullradtbl
    real(double), pointer, dimension(:) :: table1
    real(double), pointer, dimension(:) :: table2
        
    real(double) :: rcut1,rcut2,rcut1_fake,rcut2_fake,delr1,delr2
    real(double) :: kcut,deltak,d12,del_k,k_cut,del_r,yp1,ypn
    real(double) :: factor
    integer npts1,npts2,l1,l2,l,npts1_2,npts2_2,lmin,lmax,h,count
    integer i,j,n12,n1_new,n2_new,ke,ld,npts,stat

    deltak = del_k
    kcut = k_cut
    
    call maxtwon(npts1,delr1,npts2,delr2,n12,d12,deltak,kcut)
    !ld+1 is now radial table multiplicity i.e. no. of 'l3' values permitted
    ld = min(l1,l2)
    stat=0
    allocate(fullradtbl(n12/4,ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    allocate(rad_tables_nlpf_pao(count)%rad_tbls(ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating rad_tables in make_rad_table: ",ld+1)
    rcut1 = delr1*(npts1-1)
    rcut2 = delr2*(npts2-1)
    allocate(dummy1(n12),dummy2(n12),dummy1bt(n12/2),dummy2bt(n12/2),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummy in make_rad_table: ",n12)
    allocate(dummyprod(n12/2),dumout(n12/4),ol_out(n12/4),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummpy outs in make_rad_table: ",n12/4)
    rcut1_fake = (n12-1)*d12
    !RC have a problem with k space arrays not having sufficient size here.
    dummy1(1:npts1) = table1(1:npts1)
    dummy2(1:npts2) = table2(1:npts2)
    dummy1(npts1+1:n12) = zero
    dummy2(npts2+1:n12) = zero
    call matcharrays(dummy1,npts1,delr1,dummy2,npts2,delr2,d12,n1_new&
         &,n2_new,n12) 
    !loop to calculate spherical Bessel transforms of basis
    !function 
    call bessloop(dummy1,l1,n1_new,n12,d12,rcut1_fake,dummy1bt)
    call bessloop(dummy2,l2,n2_new,n12,d12,rcut1_fake,dummy2bt)
    !now calculate radial-tables(l) that correspond to the 
    !overlap-integral
    dummyprod=dummy1bt*dummy2bt
    deltak = twopi/(d12+rcut1_fake)
    
    kcut = ((n12/2)-1)*deltak
    lmin = abs(l1-l2)
    lmax = l1+l2
    h=1
    do i=lmin,lmax,2
       call bessloop(dummyprod,i,n12/2,n12/2,deltak,kcut,dumout)
       call complx_fctr(l1,l2,i,factor)
       fullradtbl(1:n12/4,h) = dumout(1:n12/4)*factor 
       h = h+1
    enddo
    npts = 1+((rcut1+rcut2)*kcut/twopi)
    del_r = twopi/(kcut+deltak)
    
    call store_nlpf_pao_tables(npts,del_r,ld,fullradtbl,n12/4,count)
    deallocate(ol_out,dumout,dummyprod,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating dummpy outs in make_rad_table: ",n12/4)
    deallocate(dummy2bt, dummy1bt,dummy2,dummy1,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating rad_tables in make_rad_table: ",ld+1)
    deallocate(fullradtbl,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating fullradtbl in make_rad_table: ",n12/4,ld+1)
  end subroutine make_rad_table_nlpfpao
  !!***
  

!!****f* make_rad_tables/make_rad_table_supp_ke *
!!
!!  NAME 
!!   make_rad_table_supp_ke
!!  USAGE
!!   make_rad_table_supp_ke(npts1,l1,delr1,npts2,l2,delr2,table1,table2,count,del_k,k_cut)
!!
!!  PURPOSE
!!  Constructs overlap integral radial tables <PAO|PAO> & <PAO|KE|PAO>
!!
!!
!!  INPUTS
!!   npts1,npts2 : no of points in input tables of NLPF and PAO
!!   l1,delr1,l2,delr2 - ang mom and grid spacings of input tables 1 and 2
!!   del_k - mimimum grid spacing in k-space
!!   k_cut - minimum cut off in k-space
!!   table1,table2 - input tables containing NLPF and PAO info
!!   count - used to index the particular set of overlap integral radial 
!!   tables that are created; count is set by gen_nlpf_supp_tbls
!!  USES
!!   datatypes, ol_int_datatypes, pao_format
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   09/05/03 adding matcharrays subroutine- RC
!!   16/05/03 new (matcharrays) grid reprojection subroutines
!!    tested and are working nicely.
!!   18/09/03 corrected k-space params, added y2 construction too.
!!   28/10/03 new routine to calculate BOTH pao_pao radial tables
!!    and pao_ke_pao radial tables without option flag
!!
!!  SOURCE
!!  
  subroutine make_rad_table_supp_ke(npts1,l1,delr1,npts2,l2,delr2&
       &,table1,table2,count,del_k,k_cut)

    use datatypes
    use numbers, ONLY: zero, twopi
    use ol_int_datatypes !,ONLY : rad_tables, rad_tables_ke
    use bessel_integrals !,ONLY : maxtwon,bessloop,complx_fctr,multiply_ksq
    use cubic_spline_routines !,ONLY : matcharrays
    use GenComms, ONLY: cq_abort

    implicit none

    real(double), allocatable, dimension(:) :: dummy1,dummy2,dummy1bt,dummy2bt
    real(double), allocatable, dimension(:) :: dummyprod,dummyprod_ke,dumout
    real(double), allocatable, dimension(:) :: ol_out,xin,y2
    real(double), allocatable, dimension(:,:) :: fullradtbl,fullradtbl_ke
    real(double), pointer, dimension(:) :: table1
    real(double), pointer, dimension(:) :: table2
    real(double) :: rcut1,rcut2,rcut1_fake,rcut2_fake,delr1,delr2
    real(double) :: kcut,deltak,d12,del_k,k_cut,del_r,yp1,ypn
    real(double) :: factor
    integer npts1,npts2,l1,l2,l,npts1_2,npts2_2,lmin,lmax,h,count
    integer i,j,n12,n1_new,n2_new,ke,ld,npts,stat

    !to stop incremental corruption of deltak value
    deltak = del_k
    kcut = k_cut
    
    call maxtwon(npts1,delr1,npts2,delr2,n12,d12,deltak,kcut)
    !ld+1 is now radial table multiplicity i.e. no. of 'l3' values permitted
    ld = min(l1,l2)
    !allocate(fullradtbl(ld+1,n12/4),fullradtbl_ke(ld+1,n12/4))
    allocate(fullradtbl(n12/4,ld+1),fullradtbl_ke(n12/4,ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    allocate(rad_tables(count)%rad_tbls(ld+1), rad_tables_ke(count)%rad_tbls(ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating rad_tables in make_rad_table: ",ld+1)
    rcut1 = delr1*(npts1-1)
    rcut2 = delr2*(npts2-1)
    allocate(dummy1(n12),dummy2(n12),dummy1bt(n12/2),dummy2bt(n12/2),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummy in make_rad_table: ",n12)
    allocate(dummyprod(n12/2),dummyprod_ke(n12/2),dumout(n12/4),ol_out(n12/4),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummpy outs in make_rad_table: ",n12/4)
    rcut1_fake = (n12-1)*d12
    !RC have a problem with k space arrays not having sufficient size here.
    dummy1(1:npts1) = table1(1:npts1)
    dummy2(1:npts2) = table2(1:npts2)
    dummy1(npts1+1:n12) = zero
    dummy2(npts2+1:n12) = zero
    call matcharrays(dummy1,npts1,delr1,dummy2,npts2,delr2,d12,n1_new&
         &,n2_new,n12) 
    !loop to calculate spherical Bessel transforms of basis
    !function 
    call bessloop(dummy1,l1,n1_new,n12,d12,rcut1_fake,dummy1bt)
    call bessloop(dummy2,l2,n2_new,n12,d12,rcut1_fake,dummy2bt)
    !now calculate radial-tables(l) that correspond to the 
    !overlap-integral for BOTH pao_pao and pao_ke_pao matrix elements
    dummyprod = dummy1bt*dummy2bt
    dummyprod_ke = dummyprod
    deltak = twopi/(d12+rcut1_fake)
    call multiply_ksq(dummyprod_ke,n12/2,deltak)
    
    kcut = ((n12/2)-1)*deltak
    lmin = abs(l1-l2)
    lmax = l1+l2
    h=1
    !loop to calculate radial tables
    do i=lmin,lmax,2
       call complx_fctr(l1,l2,i,factor)
       
       call bessloop(dummyprod,i,n12/2,n12/2,deltak,kcut,dumout)
       fullradtbl(1:n12/4,h) = dumout(1:n12/4)*factor
       !RC resetting to zero(really just for peace of mind)
       dumout = 0.0_double
       call bessloop(dummyprod_ke,i,n12/2,n12/2,deltak,kcut,dumout)
       fullradtbl_ke(1:n12/4,h) = dumout(1:n12/4)*factor
       h = h+1
    enddo
    npts = 1+((rcut1+rcut2)*kcut/twopi)
    del_r = twopi/(kcut+deltak) 
    
    call store_supp_ke_tables(npts,del_r,ld,fullradtbl,fullradtbl_ke,n12/4,count)
    deallocate(ol_out,dumout,dummyprod_ke,dummyprod,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating dummpy outs in make_rad_table: ",n12/4)
    deallocate(dummy2bt, dummy1bt,dummy2,dummy1,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating rad_tables in make_rad_table: ",ld+1)
    deallocate(fullradtbl_ke,fullradtbl,STAT=stat)    
    if(stat/=0) call cq_abort("Error deallocating fullradtbl in make_rad_table: ",n12/4,ld+1)
  end subroutine make_rad_table_supp_ke
  !!***

!!****f* make_rad_tables/store_supp_ke_tables *
!!
!!  NAME 
!!   store_supp_ke_tables
!!  USAGE
!!   store_supp_ke_tables(npnts,delta_r,num_l,fullradtbl,fullradtbl_ke,count)
!!
!!  PURPOSE
!!  Stores the overlap integral radial tables between <PAO|PAO> & <PAO|KE|PAO>
!!
!!  INPUTS
!!   npnts - no of points per radial table in the set
!!   delta_r - hmm same as before just grid spacing
!!   num_l - multiplicity of radial_tables(l) (= (|l1-l2|->l1+l2),2)
!!   fullradtbl,fullradtbl_ke : all the radial tables for this overlap integral
!!   count : indexing of overlap integral set by gen_rad_tables
!!  USES
!!   datatypes, ol_int_datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  
  subroutine store_supp_ke_tables(npnts,delta_r,num_l,fullradtbl,fullradtbl_ke,nsize,count)
    use datatypes
    use cubic_spline_routines !,ONLY : spline_new
    use ol_int_datatypes !,ONLY : rad_tables, rad_tables_ke
    implicit none
    !routine to copy the radial table information from the dummy structures
    !fullradtbl and fullradtbl_ke into the storage types defined in ol_int_dataypes.
    !module.f90
    integer, intent(in) :: npnts,num_l,count,nsize
    real(double), intent(in) :: delta_r
    real(double), intent(in), dimension(nsize,num_l+1) :: fullradtbl, fullradtbl_ke
    real(double) :: yp1, ypn
    real(double), dimension(npnts) :: xin, y2
    integer :: i

    !defining x axis array for the spline routines
    do i = 1, npnts
       xin(i) = (i-1)*delta_r
    enddo
    
    do i = 1, num_l+1
       allocate(rad_tables(count)%rad_tbls(i)%arr_vals(1:npnts),&
            &rad_tables_ke(count)%rad_tbls(i)%arr_vals(1:npnts))
       allocate(rad_tables(count)%rad_tbls(i)%arr_vals2(1:npnts),&
            &rad_tables_ke(count)%rad_tbls(i)%arr_vals2(1:npnts))

       rad_tables(count)%rad_tbls(i)%arr_vals(1:npnts) = fullradtbl(1:npnts,i)
       rad_tables_ke(count)%rad_tbls(i)%arr_vals(1:npnts) = fullradtbl_ke(1:npnts,i)
       
       !RC calculating and storing tables of 2nd derivatives here
       !pao_pao case first
       y2 = 0.0_double
       yp1 = (fullradtbl(2,i)-fullradtbl(1,i))/delta_r
       ypn = (fullradtbl(npnts,i)-fullradtbl(npnts-1,i))/delta_r
       call spline_new(xin,fullradtbl(1:npnts,i),npnts,yp1,ypn,y2)
       rad_tables(count)%rad_tbls(i)%arr_vals2(1:npnts) = y2(1:npnts)
       !now pao_ke_pao case
       y2 = 0.0_double
       yp1 = (fullradtbl_ke(2,i)-fullradtbl_ke(1,i))/delta_r
       ypn = (fullradtbl_ke(npnts,i)-fullradtbl_ke(npnts-1,i))/delta_r
       call spline_new(xin,fullradtbl_ke(1:npnts,i),npnts,yp1,ypn,y2)
       rad_tables_ke(count)%rad_tbls(i)%arr_vals2(1:npnts) = y2(1:npnts)

       !RC finished storing arrays of 2nd derivatives
       rad_tables(count)%rad_tbls(i)%npnts = npnts
       rad_tables(count)%rad_tbls(i)%del_x = delta_r
       
       rad_tables_ke(count)%rad_tbls(i)%npnts = npnts
       rad_tables_ke(count)%rad_tbls(i)%del_x = delta_r
    enddo

  end subroutine store_supp_ke_tables
  !!***

!!****f* make_rad_tables/store_nlpf_pao_tables *
!!
!!  NAME 
!!   store_nlpf_pao_tables
!!  USAGE
!!   store_nlpf_pao_tables(npnts,delta_r,num_l,fullradtbl,fullradtbl_ke,count)
!!
!!  PURPOSE
!!  Stores the overlap integral radial tables between <NLPF|PAO>
!!
!!  INPUTS
!!   npnts - no of points per radial table in the set
!!   delta_r - hmm same as before just grid spacing
!!   num_l - multiplicity of radial_tables(l) (= (|l1-l2|->l1+l2),2)
!!   fullradtbl,fullradtbl_ke : all the radial tables for this overlap integral
!!   count : indexing of overlap integral set by gen_rad_tables
!!  USES
!!   datatypes, ol_int_datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!

  subroutine store_nlpf_pao_tables(npnts,delta_r,num_l,fullradtbl,nsize,count)
    use datatypes
    use ol_int_datatypes !,ONLY : rad_tables_nlpf_pao
    use cubic_spline_routines !,ONLY : spline_new
    implicit none
    !routine to store the non-local projector function/pao 
    !overlap radial tables in the appropriate derived type.
    integer, intent(in) :: npnts,num_l,count,nsize
    real(double), intent(in) :: delta_r
    real(double), intent(in), dimension(nsize,num_l+1) :: fullradtbl
    real(double) :: yp1, ypn
    real(double), dimension(npnts) :: xin, y2
    integer :: i
    !defining x axis array for the spline routines
    do i = 1, npnts
       xin(i) = (i-1)*delta_r
    enddo
    
    do i = 1, num_l+1
       allocate(rad_tables_nlpf_pao(count)%rad_tbls(i)%arr_vals(1:npnts),&
            &rad_tables_nlpf_pao(count)%rad_tbls(i)%arr_vals2(1:npnts))
      
       rad_tables_nlpf_pao(count)%rad_tbls(i)%arr_vals(1:npnts) = fullradtbl(1:npnts,i)
       
       !RC calculating and storing tables of 2nd derivatives here for
       !nlpf_pao case 
       y2 = 0.0_double
       yp1 = (fullradtbl(2,i)-fullradtbl(1,i))/delta_r
       ypn = (fullradtbl(npnts,i)-fullradtbl(npnts-1,i))/delta_r
       call spline_new(xin,fullradtbl(1:npnts,i),npnts,yp1,ypn,y2)
       rad_tables_nlpf_pao(count)%rad_tbls(i)%arr_vals2(1:npnts) = y2(1:npnts)
       
       !RC finished storing arrays of 2nd derivatives
       rad_tables_nlpf_pao(count)%rad_tbls(i)%npnts = npnts
       rad_tables_nlpf_pao(count)%rad_tbls(i)%del_x = delta_r
       
    enddo
    
  end subroutine store_nlpf_pao_tables
  !!***

!!****f* make_rad_tables/gen_rad_tables *
!!
!!  NAME 
!!   gen_rad_tables
!!  USAGE
!!  gen_rad_tables(nspecies,del_k,kcut,ke_flag)
!!
!!  PURPOSE
!!   Loops over all paos, calling make_rad_table and then storing and indexing
!!   radial tables produced within the rad_tables array. 
!!
!!  INPUTS
!!   nspecies - no of atomic species
!!   del_k - mimimum grid spacing in k-space
!!   kcut - minimum cut off in k-space
!!   ke_flag - if 1 then calculate ke radial tables, otherwise no.
!!  USES
!!   datatypes, ol_int_datatypes, pao_format
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  
  subroutine gen_rad_tables(inode,ionode)

    use datatypes
    use ol_int_datatypes !,ONLY : ol_index,rad_tables,rad_tables_ke
    use pao_format !,ONLY : pao
    use GenComms, ONLY: cq_abort, gcopy, myid, my_barrier

    implicit none

    !routine to loop over the paos and generate the
    !array of radial tables
    integer :: nsp1,nsp2,nz1,nz2,l1,l2,n1,n2,i,j,k,count,ke_flag
    integer :: smallcount,nspecies,lmax,nzmax,inode,ionode,lun
    real(double) :: del1,del2,del
    real(double), pointer, dimension(:) :: fire1,fire2
    
    !RC fixing routine to read in k space parameters
    call get_max_paoparams(lmax,nzmax)
    nspecies = size(pao)
    allocate(ol_index(nspecies,nspecies,1:nzmax,1:nzmax,0:lmax,0:lmax))
    !write(*,*) lmax,nzmax, 'ol_index allocated'
    count = 1
    do i=1,2
       do nsp1 = 1,nspecies
          do nsp2 = 1,nspecies
             do l1 = 0, pao(nsp1)%greatest_angmom
                do l2 = 0, pao(nsp2)%greatest_angmom
                   do nz1 = 1, pao(nsp1)%angmom(l1)%n_zeta_in_angmom
                      do nz2 = 1, pao(nsp2)%angmom(l2)%n_zeta_in_angmom
                         if(i.eq.1) then !allocate storage arrays on first iteration 
                            count = count+1
                         else
                            ol_index(nsp1,nsp2,nz1,nz2,l1,l2) = count
                            ! now calculating radial tables
                            n1 = pao(nsp1)%angmom(l1)%zeta(nz1)%length 
                            n2 = pao(nsp2)%angmom(l2)%zeta(nz2)%length
                            del1 = pao(nsp1)%angmom(l1)%zeta(nz1)%cutoff/&
                                 &(pao(nsp1)%angmom(l1)%zeta(nz1)%length-1)
                            del2 = pao(nsp2)%angmom(l2)%zeta(nz2)%cutoff/&
                              &(pao(nsp2)%angmom(l2)%zeta(nz2)%length-1)
                            
                            allocate(fire1(n1),fire2(n2))
                            !RATHIN switching off for Gaussian testing
                            call unnorm_siesta_tbl(fire1,fire2,n1,n2,&
                                 &pao(nsp1)%angmom(l1)%zeta(nz1)%table,pao(nsp2)%angmom(l2)&
                                 &%zeta(nz2)%table,del1,del2,l1,l2)
                            call make_rad_table_supp_ke(n1,l1,del1,n2,l2,del2,fire1,fire2,&
                                 &count,del_k,kcut)
                            deallocate(fire1,fire2)
                            !RC N.B. WE HAVE NOT DONE THE BELOW FOR KE TABLES!
                            rad_tables(count)%no_of_lvals = min(l1,l2)+1
                            allocate(rad_tables(count)%l_values(min(l1,l2)+1))
                            smallcount = 1
                            do k=abs(l1-l2),l1+l2,2
                               rad_tables(count)%l_values(smallcount)=k
                               smallcount = smallcount+1
                            enddo
                            count = count+1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
       if(i.eq.1) then
          !write(*,*) 'now allocating rad_tables(ke) storage type', ' count =', count
          allocate(rad_tables(count))
          allocate(rad_tables_ke(count))
          count = 1
       else
          continue
       endif
    enddo !corresponding to the i counter
    !all radial tables have now been generated for pao_pao and pao_ke_pao
  end subroutine gen_rad_tables
  !!***
  
!!****f* make_rad_tables/gen_nlpf_supp_tbls *
!!
!!  NAME 
!!   gen_nlpf_supp_tbls
!!  USAGE
!!  gen_nlpf_supp_tbls(inode,ionode)
!!
!!  PURPOSE
!!   Loops over all NLPFs and paos, calling make_rad_table and then storing and indexing
!!   radial tables produced within the rad_tables_nlpf_pao array. 
!!
!!  INPUTS
!!   nspecies - no of atomic species
!!   del_k - mimimum grid spacing in k-space
!!   kcut - minimum cut off in k-space
!!   ke_flag - if 1 then calculate ke radial tables, otherwise no.
!!  USES
!!   datatypes, ol_int_datatypes, pao_format
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!

  subroutine gen_nlpf_supp_tbls(inode,ionode)
    use datatypes
    use pseudo_tm_info, ONLY: pseudo
    use ol_int_datatypes !,ONLY : rad_tables_nlpf_pao,ol_index_nlpf_pao
    use pao_format !,ONLY : pao
    use GenComms, ONLY: cq_abort, gcopy, myid, my_barrier
    implicit none
    !routine to loop over the paos and non-local projector functions
    !and generate the array of radial tables
    integer :: nsp1,nsp2,nz1,nz2,l1,l2,n1,n2,ke_flag,i,j,k,count
    integer :: n_pfnls,i_pfnl
    integer :: smallcount,nspecies,lmax,nzmax,inode,ionode,lun
    real(double) :: del1,del2,del
    real(double), pointer, dimension(:) :: fire1,fire2
    
    !setting allocation bounds for indexing array
    call get_max_pao_nlpfparams(lmax,nzmax)
    nspecies = size(pao)
    !write(*,*) 'allocating indexing array'
    allocate(ol_index_nlpf_pao(nspecies,nspecies,1:nzmax,1:nzmax,0:lmax,0:lmax))
    !write(*,*) lmax,nzmax, 'ol_index_nlpf_pao allocated'
    
    count = 1
    do i=1,2
       do nsp1 = 1,nspecies
          do l1 = 0, pao(nsp1)%greatest_angmom
             !do l2 = 0, pao(nsp2)%greatest_angmom
             do nz1 = 1, pao(nsp1)%angmom(l1)%n_zeta_in_angmom
                do nsp2 = 1,nspecies
                   n_pfnls = pseudo(nsp2)%n_pjnl
                  if(n_pfnls > 0) then
                   do i_pfnl = 1, n_pfnls !1,0 should abort loop
                         l2 = pseudo(nsp2)%pjnl_l(i_pfnl)
                         nz2 = pseudo(nsp2)%pjnl_n(i_pfnl)
                         if(i.eq.1) then !allocate storage arrays on first iteration
                            count = count+1
                         else
                            ol_index_nlpf_pao(nsp1,nsp2,nz1,nz2,l1,l2) = count
                            !write(47+myid,*) nsp1,nsp2,nz1,nz2
                            !write(47+myid,*) l1,l2,count
                            ! now calculating radial tables
                            n1 = pao(nsp1)%angmom(l1)%zeta(nz1)%length 
                            n2 = pseudo(nsp2)%pjnl(i_pfnl)%n
                            del1 = pao(nsp1)%angmom(l1)%zeta(nz1)%cutoff/&
                                 &(pao(nsp1)%angmom(l1)%zeta(nz1)%length-1)
                            del2 = pseudo(nsp2)%pjnl(i_pfnl)%delta
                            
                            allocate(fire1(n1),fire2(n2))
                            call unnorm_siesta_tbl(fire1,fire2,n1,n2,&
                                 &pao(nsp1)%angmom(l1)%zeta(nz1)%table,pseudo(nsp2)%&
                                 &pjnl(i_pfnl)%f,del1,del2,l1,l2)
                            call make_rad_table_nlpfpao(n1,l1,del1,n2,l2,del2,fire1,fire2,count,del_k&
                                 &,kcut)
                            call my_barrier()
                            deallocate(fire1,fire2)
                            rad_tables_nlpf_pao(count)%no_of_lvals = min(l1,l2)+1
                            allocate(rad_tables_nlpf_pao(count)%l_values(min(l1,l2)+1))
                            smallcount = 1
                            do k=abs(l1-l2),l1+l2,2
                            rad_tables_nlpf_pao(count)%l_values(smallcount)=k
                            smallcount = smallcount+1
                            enddo
                            count = count+1
                         endif
                   enddo
                  endif ! (n_pfnls > 0) 
                enddo
             enddo
          enddo
       enddo
       if(i.eq.1) then
          !write(*,*) 'now allocating rad_tables storage type', 'count =', count
          allocate(rad_tables_nlpf_pao(count))
          count = 1
       else
          continue
       endif
    enddo !corresponding to the i counter
    
  end subroutine gen_nlpf_supp_tbls
  !!***

!!****f* make_rad_tables/get_max_paoparams *
!!
!!  NAME 
!!   get_max_paoparams
!!  USAGE
!!   get_max_paoparams(lmax,nzmax)
!!  PURPOSE
!!   Searches through all paos and returns maximum l value and 
!!   maximum zeta value.
!!  INPUTS
!!   Pao's stored in pao_format
!! 
!!  USES
!!   datatypes, pao_format
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!

subroutine get_max_paoparams(lmax,nzmax)
  use datatypes
  use pao_format !,ONLY : pao
  implicit none
  !code to analyse pao information and retrieve 
  !the maximum angular momentum and zeta value associated with
  !any of the paos.
  integer, intent(out) :: lmax,nzmax
  integer :: n_sp,n_zeta,nt,l,nspecies_tot
  
  !loop to figure out lmax and nzmax
  lmax = 0 !initialising
  nzmax = 0 !initialising
        
  nspecies_tot = size(pao)
  do n_sp = 1,nspecies_tot
     if(lmax.lt.pao(n_sp)%greatest_angmom) then
        lmax = pao(n_sp)%greatest_angmom
     else
        continue
     endif
     do l = 0, pao(n_sp)%greatest_angmom
        if(nzmax.lt.pao(n_sp)%angmom(l)%n_zeta_in_angmom) then
           nzmax = pao(n_sp)%angmom(l)%n_zeta_in_angmom
        else
           continue
        endif
     enddo
  enddo
  
end subroutine get_max_paoparams
!!***

!!****f* make_rad_tables/get_max_pao_nlpfparams *
!!
!!  NAME 
!!   get_max_pao_nlpfparams
!!  USAGE
!!   get_max_pao_nlpfparams(lmax,nzmax)
!!  PURPOSE
!!   Searches through all paos and returns maximum l value and 
!!   maximum zeta value.
!!  INPUTS
!!   Pao's stored in pao_format
!! 
!!  USES
!!   datatypes, pao_format, pseudo_tm_info
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!

subroutine get_max_pao_nlpfparams(lmax,nzmax)
  use datatypes
  use pao_format !, ONLY ; pao
  use pseudo_tm_info !,ONLY : pseudo
  implicit none
  !code to analyse pao AND nlpf information and retrieve 
  !the maximum angular momentum and zeta value associated with
  !any of the paos/nlpfs
  integer, intent(out) :: lmax,nzmax
  integer :: n_sp,n_zeta,nt,l,nspecies_tot
  
  !loop to figure out lmax and nzmax
  lmax = 0 !initialising
  nzmax = 0 !initialising
        
  nspecies_tot = size(pao)
  do n_sp = 1,nspecies_tot
     if(lmax.lt.pao(n_sp)%greatest_angmom) then
        lmax = pao(n_sp)%greatest_angmom
     else
        continue
     endif
     if(lmax.lt.pseudo(n_sp)%n_pjnl) then
        lmax = pseudo(n_sp)%n_pjnl
     else
        continue
     endif    
    

     do l = 0, pao(n_sp)%greatest_angmom
        if(nzmax.lt.max(pao(n_sp)%angmom(l)%n_zeta_in_angmom&
             &,pseudo(n_sp)%n_pjnl)) then
           nzmax = max(pao(n_sp)%angmom(l)%n_zeta_in_angmom,pseudo(n_sp)%n_pjnl)
        else
           continue
        endif
     enddo
     
  enddo
  
end subroutine get_max_pao_nlpfparams
!!***

!!****f* make_rad_tables/unnorm_siesta_tbl *
!!
!!  NAME 
!!   unnorm_siesta_tbl
!!  USAGE
!!   unnorm_siesta_tbl(dummy1,dummy2,n1,n2,intable1,intable2,d1,d2,l1,l2)
!!  PURPOSE
!!   Takes away the 1/r**l normalisation from the SIESTA input PAOs
!!
!!  INPUTS
!!   dummy1,dummy2 : UNnormalized pao radial tables
!!   n1,n2 : n points in tables
!!   intable1,intable2 : SIESTA input PAO tables
!!   l1,l2 : angular momentum values of input radial tables
!!   d1,d2 : grid spacings of input radial tables
!!  USES
!!   datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!

subroutine unnorm_siesta_tbl(dummy1,dummy2,n1,n2,intable1,intable2,d1,d2,l1,l2)
  use datatypes
  !routine to multiply siesta input radial tables by r**l
  implicit none
  integer, intent(in) :: n1,l1,n2,l2
  real(double), intent(in) :: intable1(n1), intable2(n2)
  real(double), intent(inout) :: dummy1(n1), dummy2(n2)
  real(double), intent(in) :: d1,d2
  
  integer :: i
  real(double) :: r

  if(l1.gt.0) then
     do i = 1, n1
        r = (i-1)*d1
        dummy1(i) = intable1(i)*(r**(l1))
     enddo
  else
     dummy1 = intable1 !r**l = 1
  endif
  
  if(l2.gt.0) then
     do i = 1, n2
        r = (i-1)*d2
        dummy2(i) = intable2(i)*(r**l2)
     enddo
  else
     dummy2 = intable2 !r**l = 1
  endif
  
end subroutine unnorm_siesta_tbl
!!***

!!****f* make_rad_tables/get_support_pao_rep *
!!
!!  NAME 
!!   get_support_pao_rep
!!  USAGE
!!   get_support_pao_rep(inode,ionode)
!!  PURPOSE
!!   Gets the  coefficients (if read in) for representation of the
!!   support functions in terms of paos / else generates them randomly
!!  INPUTS
!!   Pao's stored in pao_format
!! 
!!  USES
!!   datatypes, pao_format, pseudo_tm_info
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   2006/06/13 08:01 dave
!!    Changing allocation of support function coefficients
!!   2006/06/13 16:58 dave
!!    Changing coefficient defaults, adding check for basis set sanity
!!  SOURCE
!!
  subroutine get_support_pao_rep(inode,ionode)
    use datatypes
    use numbers
    use global_module, ONLY: ni_in_cell, species_glob, iprint_basis
    use primary_module, ONLY: bundle
    use GenComms, ONLY : gcopy, my_barrier, cq_abort
    use species_module, ONLY : nsf_species, npao_species
    use support_spec_format, ONLY : supports_on_atom, flag_paos_atoms_in_cell, mx_pao_coeff_atoms, &
         allocate_supp_coeff_array, associate_supp_coeff_array, support_gradient, support_elec_gradient, &
         coefficient_array,grad_coeff_array, elec_grad_coeff_array, coeff_array_size, &
         read_option, symmetry_breaking, support_pao_file
    use pao_format, ONLY: pao

    implicit none
    !code to allocate latest support_pao data structure
    integer, intent(in) :: inode, ionode
    integer :: lmax, l, m, acz, i, total_size, n_sup, species_i
    integer :: acz_thisl, count, lun, ios, j, k, n, myflag, nacz, stat, idum, count_pao
    real(double) :: coeff, sum, tmp
    logical :: warn_flag
    !integer, external :: time
    
    warn_flag = .false.
    if(inode==ionode.AND.iprint_basis>1) &
         write(*,*) 'Flag for all pao coeffs on all procs is ',flag_paos_atoms_in_cell
    if(flag_paos_atoms_in_cell) then
       mx_pao_coeff_atoms = ni_in_cell
    else
       mx_pao_coeff_atoms = bundle%n_prim
       ! This is TEMPORARY !
       call cq_abort("Hard Failure in get_support_pao_rep: MUST store ALL PAO coefficients.")
    end if
    allocate(supports_on_atom(mx_pao_coeff_atoms))
    allocate(support_gradient(mx_pao_coeff_atoms))
    allocate(support_elec_gradient(mx_pao_coeff_atoms))
    total_size = 0
    do i = 1, mx_pao_coeff_atoms
       if(flag_paos_atoms_in_cell) then
          species_i = species_glob(i)
       else
          species_i = bundle%species(i)
       end if
       n_sup = nsf_species(species_i)
       if(iprint_basis>2) write(*,*) 'atom, supp: ',i,species_i,n_sup
       supports_on_atom(i)%nsuppfuncs = n_sup
       allocate(supports_on_atom(i)%supp_func(n_sup))
       support_gradient(i)%nsuppfuncs = n_sup
       allocate(support_gradient(i)%supp_func(n_sup))
       support_elec_gradient(i)%nsuppfuncs = n_sup
       allocate(support_elec_gradient(i)%supp_func(n_sup))
       if(iprint_basis>2) write(*,*) 'atom, lmax: ',i,pao(species_i)%greatest_angmom
       supports_on_atom(i)%lmax = pao(species_i)%greatest_angmom
       allocate(supports_on_atom(i)%naczs(0:supports_on_atom(i)%lmax))
       count = 0
       do l = 0, supports_on_atom(i)%lmax
          supports_on_atom(i)%naczs(l) = pao(species_i)%angmom(l)%n_zeta_in_angmom
          if(supports_on_atom(i)%naczs(l)>0) count = count + (2*l+1) ! Accumulate no. of ang. mom. components
       enddo
       ! If number of support functions is less than total number of ang. mom. components (ignoring 
       ! for now multiple zetas) then there is a formal problem with basis set: we require the user
       ! to set an additional flag to assert that this is really desired
       if(count>supports_on_atom(i)%nsuppfuncs.AND..NOT.warn_flag) then 
          if(.NOT.symmetry_breaking.OR..NOT.read_option) then
             if(inode==ionode) then
                write(*,fmt='("You have a major problem with your basis set.")')
                write(*,fmt='("There are less support functions than the minimal angular momentum")')
                write(*,fmt='("components.  Either increase number of support functions on species ",i4)') &
                     species_glob(i)
                write(*,fmt='("to",i4," or set flag BasisSet.SymmetryBreaking to T")') &
                     count
                write(*,fmt='("You must also specify the basis set coefficients.")')
                write(*,fmt='("Use read_pao_coeffs T and support_pao_file <filename>.")')
             end if
             call cq_abort("Basis set error for species ",species_glob(i))
          else if(symmetry_breaking.AND.read_option) then
             if(.NOT.warn_flag) then
                write(*,fmt='("You have a major problem with your basis set.")')
                write(*,fmt='("There are ",i4," support functions and",i4," angular momentum")') &
                     supports_on_atom(i)%nsuppfuncs,count                     
                write(*,fmt='("components.  But as BasisSet.SymmetryBreaking is set T we will continue.")')
                warn_flag = .true.
             end if
          end if
       end if
       count = pao(species_i)%count
       if(inode==ionode.AND.iprint_basis>2) write(*,*) 'PAO basis size: ',count
       do k=1,supports_on_atom(i)%nsuppfuncs
          supports_on_atom(i)%supp_func(k)%ncoeffs = count
          support_gradient(i)%supp_func(k)%ncoeffs = count
          support_elec_gradient(i)%supp_func(k)%ncoeffs = count
          total_size = total_size + count
       end do
    end do
    call allocate_supp_coeff_array(total_size)
    coeff_array_size = total_size
    call associate_supp_coeff_array(supports_on_atom,mx_pao_coeff_atoms,coefficient_array,total_size)
    call associate_supp_coeff_array(support_gradient,mx_pao_coeff_atoms,grad_coeff_array,total_size)
    call associate_supp_coeff_array(support_elec_gradient,mx_pao_coeff_atoms,elec_grad_coeff_array,total_size)
    if(read_option) then
       if(inode == ionode) then
          call io_assign(lun)
          open(unit=lun,file=support_pao_file,status='old',iostat=ios)
          if(ios/=0) call cq_abort("Error opening file for PAOs ")
       endif
       !-------------------------------------------------------!
       !now read in the support function PAO coefficients      !
       !                   atom by atom                        !
       !-------------------------------------------------------!
       if(flag_paos_atoms_in_cell) then ! We can read EVERYTHING in and broadcast
          do i = 1, mx_pao_coeff_atoms
             if(inode==ionode) then
                if(iprint_basis>2) write(*,*) 'Atom: ',i
                read(lun,*) n_sup
                if(supports_on_atom(i)%nsuppfuncs/=n_sup) &
                     call cq_abort("n_sup mismatch in PAO reading: ",supports_on_atom(i)%nsuppfuncs,n_sup)
                if(iprint_basis>2) write(*,*) 'N_sup: ',n_sup
                read(lun,*) lmax
                if(supports_on_atom(i)%lmax/=lmax) &
                     call cq_abort("lmax mismatch in PAO reading: ",supports_on_atom(i)%lmax,lmax)
                if(iprint_basis>2) write(*,*) 'lmax: ',lmax
                count = 0
                do l = 0, supports_on_atom(i)%lmax
                   read(lun,*) nacz
                   if(supports_on_atom(i)%naczs(l)/=nacz) &
                        call cq_abort("zeta mismatch in PAO reading: ",supports_on_atom(i)%naczs(l),nacz)
                   count = count + (2*l+1)*supports_on_atom(i)%naczs(l)
                   if(iprint_basis>2) write(*,*) 'zeta: ',l,nacz,count
                enddo
                do k = 1, supports_on_atom(i)%nsuppfuncs
                   if(count/=supports_on_atom(i)%supp_func(k)%ncoeffs) &
                        call cq_abort("coeff size mismatch in PAO reading: ",count,supports_on_atom(i)%supp_func(k)%ncoeffs)
                end do
                !end if
                do k = 1, supports_on_atom(i)%nsuppfuncs
                   !   if(inode==ionode) then
                   !-------------------------------------------------------!
                   !          read and store the coefficients              !     
                   !-------------------------------------------------------!
                   count = 1
                   sum = zero
                   do l = 0, supports_on_atom(i)%lmax
                      do acz = 1, supports_on_atom(i)%naczs(l)
                         do m = -l,l
                            read(lun,*) supports_on_atom(i)%supp_func(k)%coefficients(count)
                            sum = sum + supports_on_atom(i)%supp_func(k)%coefficients(count)* &
                                 supports_on_atom(i)%supp_func(k)%coefficients(count)
                            if(iprint_basis>2) write(*,*) 'Coeff: ',count,supports_on_atom(i)%supp_func(k)%coefficients(count)
                            count = count+1
                         enddo ! m=-l,l
                      enddo ! acz=1,supports_on_atom(i)%naczetas(l)
                   enddo ! do l=0,lmax
                   sum = sqrt(sum)
                   supports_on_atom(i)%supp_func(k)%coefficients = supports_on_atom(i)%supp_func(k)%coefficients/sum
                enddo ! k=1,nsuppfuncs
             endif ! inode==ionode
          enddo ! mx_pao_coeff_atoms
          call gcopy(coefficient_array, total_size)
       else
          call cq_abort("This option not implemented")
          do i = 1, ni_in_cell ! We need ALL atoms read in
             if(inode==ionode) then
                if(iprint_basis>2) write(*,*) 'Atom: ',i
                read(lun,*) n_sup
                if(supports_on_atom(i)%nsuppfuncs/=n_sup) &
                     call cq_abort("n_sup mismatch in PAO reading: ",supports_on_atom(i)%nsuppfuncs,n_sup)
                if(iprint_basis>2) write(*,*) 'N_sup: ',n_sup
                read(lun,*) lmax
                if(supports_on_atom(i)%lmax/=lmax) &
                     call cq_abort("lmax mismatch in PAO reading: ",supports_on_atom(i)%lmax,lmax)
                if(iprint_basis>2) write(*,*) 'lmax: ',lmax
                count = 0
                do l = 0, supports_on_atom(i)%lmax
                   read(lun,*) nacz
                   if(supports_on_atom(i)%naczs(l)/=nacz) &
                        call cq_abort("zeta mismatch in PAO reading: ",supports_on_atom(i)%naczs(l),nacz)
                   count = count + (2*l+1)*supports_on_atom(i)%naczs(l)
                   if(iprint_basis>2) write(*,*) 'zeta: ',l,nacz,count
                enddo
                do k = 1, supports_on_atom(i)%nsuppfuncs
                   if(count/=supports_on_atom(i)%supp_func(k)%ncoeffs) &
                        call cq_abort("coeff size mismatch in PAO reading: ",count,supports_on_atom(i)%supp_func(k)%ncoeffs)
                end do
             end if
             if(inode==ionode) then
                do k = 1, supports_on_atom(i)%nsuppfuncs
                   !-------------------------------------------------------!
                   !          read and store the coefficients              !     
                   !-------------------------------------------------------!
                   count = 1
                   sum = zero
                   do l = 0, supports_on_atom(i)%lmax
                      do acz = 1, supports_on_atom(i)%naczs(l)
                         do m = -l,l
                            read(lun,*) supports_on_atom(i)%supp_func(k)%coefficients(count)
                            sum = sum + supports_on_atom(i)%supp_func(k)%coefficients(count)* &
                                 supports_on_atom(i)%supp_func(k)%coefficients(count)
                            if(iprint_basis>2) write(*,*) 'Coeff: ',count,supports_on_atom(i)%supp_func(k)%coefficients(count)
                            count = count+1
                         enddo ! m=-l,l
                      enddo ! acz=1,supports_on_atom(i)%naczetas(l)
                   enddo ! do l=0,lmax
                   sum = sqrt(sum)
                   supports_on_atom(i)%supp_func(k)%coefficients = supports_on_atom(i)%supp_func(k)%coefficients/sum
                enddo ! k=1,nsuppfuncs
             end if
             ! Here we need to determine WHICH processor this belongs to, and send it
             ! We'll need a send on ionode and recv on owning proc, and then a barrier.
          enddo ! mx_pao_coeff_atoms
       end if
       if(inode==ionode) call io_close(lun)
    else ! Initialise from PAO structures
       do i = 1, mx_pao_coeff_atoms
          if(flag_paos_atoms_in_cell) then
             species_i = species_glob(i)
          else
             species_i = bundle%species(i)
          end if
          if(nsf_species(species_i)==npao_species(species_i)) then ! We have one-to-one
             do k = 1, supports_on_atom(i)%nsuppfuncs
                count = 1
                do l = 0, supports_on_atom(i)%lmax
                   do acz = 1, supports_on_atom(i)%naczs(l)
                      do m = -l,l
                         if(k==count) then
                            supports_on_atom(i)%supp_func(k)%coefficients(count) = one
                         else
                            supports_on_atom(i)%supp_func(k)%coefficients(count) = zero
                         end if
                         count = count+1
                      enddo ! m=-l,l
                   enddo ! acz=1,supports_on_atom(i)%naczetas(l)
                enddo ! do l=0,lmax
             enddo ! k=1,nsuppfuncs
          else
             do k = 1, supports_on_atom(i)%nsuppfuncs
                count_pao = 1 ! which PAO 
                sum = zero
                do l = 0, supports_on_atom(i)%lmax
                   do acz = 1, supports_on_atom(i)%naczs(l)
                      ! Count indexes which angular momentum channel we're in
                      count = 1
                      if(l>0.AND.supports_on_atom(i)%naczs(l)>0) then
                         do m=0,l-1
                            if(supports_on_atom(i)%naczs(m)>0) count = count + 2*m+1
                         end do
                      end if
                      if(inode==ionode.AND.iprint_basis>2) write(*,*) 'Counter: ',count,count_pao,k
                      do m = -l,l
                         if(k==count) then
                            if(acz==1) then
                               supports_on_atom(i)%supp_func(k)%coefficients(count_pao) = one
                            else
                               supports_on_atom(i)%supp_func(k)%coefficients(count_pao) = 0.1_double
                            end if
                         else
                            supports_on_atom(i)%supp_func(k)%coefficients(count_pao) = zero
                         end if
                         sum = sum + supports_on_atom(i)%supp_func(k)%coefficients(count_pao)* &
                              supports_on_atom(i)%supp_func(k)%coefficients(count_pao)
                         count = count + 1
                         count_pao = count_pao + 1
                      enddo ! m=-l,l
                   enddo ! acz=1,supports_on_atom(i)%naczetas(l)
                enddo ! do l=0,lmax
                sum = sqrt(sum)
                supports_on_atom(i)%supp_func(k)%coefficients = supports_on_atom(i)%supp_func(k)%coefficients/sum
             enddo ! k=1,nsuppfuncs
          end if ! Is this one-to-one PAOs-to-SFs
       enddo ! i=1,mx_pao_coeff_atoms
       call writeout_support_functions(inode,ionode)
    end if
    call my_barrier
    
  end subroutine get_support_pao_rep
!!***
  
  subroutine writeout_support_functions(inode,ionode)
    use datatypes
    use support_spec_format, ONLY : supports_on_atom, support_pao_file
    use global_module, ONLY: ni_in_cell !this is our total no of atoms
    use GenComms, ONLY : gcopy, my_barrier
    use species_module, ONLY : n_species
    implicit none
    !code to write out latest support_pao data structure
    integer, intent(in) :: inode, ionode
    integer :: i,j,k,l,m,n,myflag,aczeta,count,lun,ios
    
    !write(*,*) "Entering writeout ",inode
    if(inode == ionode) then
       call io_assign(lun)
       open(unit=lun,file=support_pao_file,iostat=ios)
       do i = 1, ni_in_cell
          !write(*,*) "Atom ",i,inode
          write(lun,*) supports_on_atom(i)%nsuppfuncs
          write(lun,*) supports_on_atom(i)%lmax
          do l = 0, supports_on_atom(i)%lmax
             write(lun,*) supports_on_atom(i)%naczs(l)
          enddo
          do j = 1, supports_on_atom(i)%nsuppfuncs
             count = 1
             do l = 0, supports_on_atom(i)%lmax
                do aczeta = 1, supports_on_atom(i)%naczs(l)
                   do m = -l,l
                      write(lun,*) supports_on_atom(i)%supp_func(j)%coefficients(count)
                      count = count + 1
                   enddo ! m
                enddo ! acz
             enddo ! l
          enddo ! nsuppfuncs
       enddo ! ni_in_cell
    end if
    call my_barrier()
    if(inode==ionode) call io_close(lun)
!    write(*,*) "Finished writeout ",inode
                
  end subroutine writeout_support_functions

! =====================================================================
!   sbrt ran2: generates uniform random numbers in the
!   interval [0,1]. Taken from Numerical Recipes, 1st edition,
!   page 197.
! ---------------------------------------------------------------------
  subroutine ran2(x,idum)
                                                                                
    use datatypes
                                                                                
    implicit none
                                                                                
    integer, parameter :: m=714025,ia=1366,ic=150889
    real(double), parameter :: rm=1.0_double/m
                                                                                
    integer,save :: ir(97) !!   I add save statement 13/1/2000 TM
    integer,save :: iff=0
    integer,save :: j,iy   !!   I add save statement 13/1/2000 TM
    integer :: idum
    real(double) :: x
 
!    data iff/0/
 
    if((idum.lt.0).or.(iff.eq.0)) then
      iff=1
      idum=mod(ic-idum,m)
      do j=1,97
        idum=mod(ia*idum+ic,m)
        ir(j)=idum
      enddo
      idum=mod(ia*idum+ic,m)
      iy=idum
    endif
!       write(*,11) idum,iy,iff,m,j,1+(97*iy)/m
    j=1+(97*iy)/m
 11 format('idum,iy,iff,m,old j,new j in ran2= ',6i12)
    if((j.gt.97).or.(j.lt.1)) pause
    iy=ir(j)
    x=iy*rm
    idum=mod(ia*idum+ic,m)
    ir(j)=idum
    return
  end subroutine ran2

end module make_rad_tables
!!***


