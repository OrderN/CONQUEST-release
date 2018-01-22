! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module make_rad_tables
! ------------------------------------------------------------------------------
! Code area 11: basis functions
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
!!   2008/02/06 08:29 dave
!!    Changed for output to file not stdout
!!   2008/02/10 ast
!!    Added timers
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2017/02/21 16:00 nakata
!!    commented out get_support_pao_rep and writeout_support_functions
!!    which are no longer used
!!   2017/03/08 15:00 nakata
!!    removed get_support_pao_rep, writeout_support_functions and ran2
!!    which are no longer used
!!  SOURCE
module make_rad_tables

  use global_module,          only: io_lun, area_basis
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_basis, tmr_std_allocation

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
!!   09:29, 27/11/2007 drb 
!!    Bug fix: zero dummy1bt and dummy2bt before use
!!   2008/02/10 ast
!!    Added timers
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
    real(double), allocatable, dimension(:) :: dummyprod,dumout,ol_out
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
    call start_timer(tmr_std_allocation)
    allocate(fullradtbl(n12/4,ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    allocate(rad_tables_nlpf_pao(count)%rad_tbls(ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating rad_tables in make_rad_table: ",ld+1)
    call stop_timer(tmr_std_allocation)
    rcut1 = delr1*(npts1-1)
    rcut2 = delr2*(npts2-1)
    call start_timer(tmr_std_allocation)
    allocate(dummy1(n12),dummy2(n12),dummy1bt(n12/2),dummy2bt(n12/2),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummy in make_rad_table: ",n12)
    allocate(dummyprod(n12/2),dumout(n12/4),ol_out(n12/4),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummpy outs in make_rad_table: ",n12/4)
    call stop_timer(tmr_std_allocation)
    rcut1_fake = (n12-1)*d12
    !RC have a problem with k space arrays not having sufficient size here.
    dummy1(1:npts1) = table1(1:npts1)
    dummy2(1:npts2) = table2(1:npts2)
    dummy1(npts1+1:n12) = zero
    dummy2(npts2+1:n12) = zero
    dummy1bt = zero
    dummy2bt = zero
    call matcharrays(dummy1,npts1,delr1,dummy2,npts2,delr2,d12,n1_new&
         &,n2_new,n12) 
    !loop to calculate spherical Bessel transforms of basis
    !function 
    call bessloop(dummy1,l1,n1_new,n12,d12,rcut1_fake,dummy1bt,1)
    call bessloop(dummy2,l2,n2_new,n12,d12,rcut1_fake,dummy2bt,1)
    !now calculate radial-tables(l) that correspond to the 
    !overlap-integral
    dummyprod=dummy1bt*dummy2bt
    deltak = twopi/(d12+rcut1_fake)
    
    kcut = ((n12/2)-1)*deltak
    lmin = abs(l1-l2)
    lmax = l1+l2
    h=1
    do i=lmin,lmax,2
       call bessloop(dummyprod,i,n12/2,n12/2,deltak,kcut,dumout,0)
       call complx_fctr(l1,l2,i,factor)
       fullradtbl(1:n12/4,h) = dumout(1:n12/4)*factor 
       h = h+1
    enddo
    npts = 1+((rcut1+rcut2)*kcut/twopi)
    del_r = twopi/(kcut+deltak)
    
    call store_nlpf_pao_tables(npts,del_r,ld,fullradtbl,n12/4,count)
    call start_timer(tmr_std_allocation)
    deallocate(ol_out,dumout,dummyprod,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating dummpy outs in make_rad_table: ",n12/4)
    deallocate(dummy2bt, dummy1bt,dummy2,dummy1,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating rad_tables in make_rad_table: ",ld+1)
    deallocate(fullradtbl,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    call stop_timer(tmr_std_allocation)
  end subroutine make_rad_table_nlpfpao
  !!***
  
  ! Neutral atom Projector functions
  ! This calculates tables for elements like < phi_i | VNA_j | phi_j >
  subroutine make_rad_table_paoNApao(npts1,l1,delr1,npts2,l2,delr2&
       &,table1,table2,count,del_k,k_cut)
    use datatypes
    use numbers, ONLY: zero, twopi
    use ol_int_datatypes , ONLY : rad_tables_paoNApao
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
    !30/10/03 modifying end of this routine to store pao_na_pao tables 

    real(double), allocatable, dimension(:) :: dummy1,dummy2,dummy1bt,dummy2bt
    real(double), allocatable, dimension(:) :: dummyprod,dumout,ol_out
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
    call start_timer(tmr_std_allocation)
    allocate(fullradtbl(n12/4,ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    allocate(rad_tables_paoNApao(count)%rad_tbls(ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating rad_tables in make_rad_table: ",ld+1)
    call stop_timer(tmr_std_allocation)
    rcut1 = delr1*(npts1-1)
    rcut2 = delr2*(npts2-1)
    call start_timer(tmr_std_allocation)
    allocate(dummy1(n12),dummy2(n12),dummy1bt(n12/2),dummy2bt(n12/2),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummy in make_rad_table: ",n12)
    allocate(dummyprod(n12/2),dumout(n12/4),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummpy outs in make_rad_table: ",n12/4)
    call stop_timer(tmr_std_allocation)
    rcut1_fake = (n12-1)*d12
    !RC have a problem with k space arrays not having sufficient size here.
    dummy1(1:npts1) = table1(1:npts1)
    dummy2(1:npts2) = table2(1:npts2)
    dummy1(npts1+1:n12) = zero
    dummy2(npts2+1:n12) = zero
    dummy1bt = zero
    dummy2bt = zero
    call matcharrays(dummy1,npts1,delr1,dummy2,npts2,delr2,d12,n1_new&
         &,n2_new,n12) 
    !loop to calculate spherical Bessel transforms of basis
    !function 
    call bessloop(dummy1,l1,n1_new,n12,d12,rcut1_fake,dummy1bt,1)
    call bessloop(dummy2,l2,n2_new,n12,d12,rcut1_fake,dummy2bt,1)
    !now calculate radial-tables(l) that correspond to the 
    !overlap-integral
    dummyprod=dummy1bt*dummy2bt
    deltak = twopi/(d12+rcut1_fake)
    
    kcut = ((n12/2)-1)*deltak
    lmin = abs(l1-l2)
    lmax = l1+l2
    h=1
    do i=lmin,lmax,2
       call bessloop(dummyprod,i,n12/2,n12/2,deltak,kcut,dumout,0)
       call complx_fctr(l1,l2,i,factor)
       fullradtbl(1:n12/4,h) = dumout(1:n12/4)*factor 
       h = h+1
    enddo
    npts = 1+((rcut1+rcut2)*kcut/twopi)
    del_r = twopi/(kcut+deltak)
    
    call store_table(npts,del_r,ld,fullradtbl,n12/4,count,rad_tables_paoNApao)
    call start_timer(tmr_std_allocation)
    deallocate(dumout,dummyprod,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating dummpy outs in make_rad_table: ",n12/4)
    deallocate(dummy2bt, dummy1bt,dummy2,dummy1,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating rad_tables in make_rad_table: ",ld+1)
    deallocate(fullradtbl,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    call stop_timer(tmr_std_allocation)

  end subroutine make_rad_table_paoNApao

  ! This calculates tables for elements like < phi_i | VNA_j | phi_i >
  subroutine make_rad_table_paopaoNA(npts1,l1,delr1,npts2,l2,delr2&
       &,table1,table2,count,del_k,k_cut)
    use datatypes
    use numbers, ONLY: zero, twopi
    use ol_int_datatypes , ONLY : rad_tables_paopaoNA!, rad_tables_paodpaoNA
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
    !30/10/03 modifying end of this routine to store pao_na_pao tables 

    real(double), allocatable, dimension(:) :: dummy1,dummy2,dummy1bt,dummy2bt
    real(double), allocatable, dimension(:) :: dummyprod,dumout,ol_out
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
    call start_timer(tmr_std_allocation)
    allocate(fullradtbl(n12/4,ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    allocate(rad_tables_paopaoNA(count)%rad_tbls(ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating rad_tables in make_rad_table: ",ld+1)
    call stop_timer(tmr_std_allocation)
    rcut1 = delr1*(npts1-1)
    rcut2 = delr2*(npts2-1)
    call start_timer(tmr_std_allocation)
    allocate(dummy1(n12),dummy2(n12),dummy1bt(n12/2),dummy2bt(n12/2),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummy in make_rad_table: ",n12)
    allocate(dummyprod(n12/2),dumout(n12/4),ol_out(n12/4),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummpy outs in make_rad_table: ",n12/4)
    call stop_timer(tmr_std_allocation)
    rcut1_fake = (n12-1)*d12
    !RC have a problem with k space arrays not having sufficient size here.
    dummy1(1:npts1) = table1(1:npts1)
    dummy2(1:npts2) = table2(1:npts2)
    dummy1(npts1+1:n12) = zero
    dummy2(npts2+1:n12) = zero
    dummy1bt = zero
    dummy2bt = zero
    call matcharrays(dummy1,npts1,delr1,dummy2,npts2,delr2,d12,n1_new&
         &,n2_new,n12)
    ! NB This is different to normal routines because we've got both spherical harmonics
    ! in the bra, and VNA in the ket, so we have to transform table1 for each l value below
    !loop to calculate spherical Bessel transforms of basis
    !function 
    !call bessloop(dummy1,l1,n1_new,n12,d12,rcut1_fake,dummy1bt,1)
    !call bessloop(dummy2,l2,n2_new,n12,d12,rcut1_fake,dummy2bt,1)
    !dummy1bt = dummy1bt * dummy2bt ! This is now phi_ia x phi_ib
    !dummy2bt = zero
    !dummy2 = zero
    !dummy2(1:npts2) = table3(1:npts2)
    !call matcharrays(dummy1,npts1,delr1,dummy2,npts2,delr2,d12,n1_new&
    !     &,n2_new,n12) 
    call bessloop(dummy2,0,n2_new,n12,d12,rcut1_fake,dummy2bt,1) ! l=0 for VNA
    !now calculate radial-tables(l) that correspond to the 
    !overlap-integral
    dummyprod=dummy1bt*dummy2bt
    deltak = twopi/(d12+rcut1_fake)
    
    kcut = ((n12/2)-1)*deltak
    lmin = abs(l1-l2)
    lmax = l1+l2
    h=1
    do i=lmin,lmax,2
       dummy1bt = zero
       call bessloop(dummy1,i,n1_new,n12,d12,rcut1_fake,dummy1bt,1)
       dummyprod=dummy1bt*dummy2bt
       call bessloop(dummyprod,i,n12/2,n12/2,deltak,kcut,dumout,0)
       !call complx_fctr2(l1,l2,i,factor) ! I *think* that the i factors cancel
       fullradtbl(1:n12/4,h) = dumout(1:n12/4)!*factor 
       h = h+1
    enddo
    npts = 1+((rcut1+rcut2)*kcut/twopi)
    del_r = twopi/(kcut+deltak)
    
    call store_table(npts,del_r,ld,fullradtbl,n12/4,count,rad_tables_paopaoNA)
    call start_timer(tmr_std_allocation)
    deallocate(ol_out,dumout,dummyprod,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating dummpy outs in make_rad_table: ",n12/4)
    deallocate(dummy2bt, dummy1bt,dummy2,dummy1,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating rad_tables in make_rad_table: ",ld+1)
    deallocate(fullradtbl,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    call stop_timer(tmr_std_allocation)

  end subroutine make_rad_table_paopaoNA
  
  ! Neutral atom Projector functions
  subroutine make_rad_table_napfpao(npts1,l1,delr1,npts2,l2,delr2&
       &,table1,table2,count,del_k,k_cut)
    use datatypes
    use numbers, ONLY: zero, twopi
    use ol_int_datatypes !, ONLY : rad_tables_napf_pao
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
    !30/10/03 modifying end of this routine to store napf_pao tables 

    real(double), allocatable, dimension(:) :: dummy1,dummy2,dummy1bt,dummy2bt
    real(double), allocatable, dimension(:) :: dummyprod,dumout,ol_out
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
    call start_timer(tmr_std_allocation)
    allocate(fullradtbl(n12/4,ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    allocate(rad_tables_napf_pao(count)%rad_tbls(ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating rad_tables in make_rad_table: ",ld+1)
    call stop_timer(tmr_std_allocation)
    rcut1 = delr1*(npts1-1)
    rcut2 = delr2*(npts2-1)
    call start_timer(tmr_std_allocation)
    allocate(dummy1(n12),dummy2(n12),dummy1bt(n12/2),dummy2bt(n12/2),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummy in make_rad_table: ",n12)
    allocate(dummyprod(n12/2),dumout(n12/4),ol_out(n12/4),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummpy outs in make_rad_table: ",n12/4)
    call stop_timer(tmr_std_allocation)
    rcut1_fake = (n12-1)*d12
    !RC have a problem with k space arrays not having sufficient size here.
    dummy1(1:npts1) = table1(1:npts1)
    dummy2(1:npts2) = table2(1:npts2)
    dummy1(npts1+1:n12) = zero
    dummy2(npts2+1:n12) = zero
    dummy1bt = zero
    dummy2bt = zero
    call matcharrays(dummy1,npts1,delr1,dummy2,npts2,delr2,d12,n1_new&
         &,n2_new,n12) 
    !loop to calculate spherical Bessel transforms of basis
    !function 
    call bessloop(dummy1,l1,n1_new,n12,d12,rcut1_fake,dummy1bt,1)
    call bessloop(dummy2,l2,n2_new,n12,d12,rcut1_fake,dummy2bt,1)
    !now calculate radial-tables(l) that correspond to the 
    !overlap-integral
    dummyprod=dummy1bt*dummy2bt
    deltak = twopi/(d12+rcut1_fake)
    
    kcut = ((n12/2)-1)*deltak
    lmin = abs(l1-l2)
    lmax = l1+l2
    h=1
    do i=lmin,lmax,2
       call bessloop(dummyprod,i,n12/2,n12/2,deltak,kcut,dumout,0)
       call complx_fctr(l1,l2,i,factor)
       fullradtbl(1:n12/4,h) = dumout(1:n12/4)*factor 
       h = h+1
    enddo
    npts = 1+((rcut1+rcut2)*kcut/twopi)
    del_r = twopi/(kcut+deltak)
    
    call store_napf_pao_tables(npts,del_r,ld,fullradtbl,n12/4,count)
    call start_timer(tmr_std_allocation)
    deallocate(ol_out,dumout,dummyprod,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating dummpy outs in make_rad_table: ",n12/4)
    deallocate(dummy2bt, dummy1bt,dummy2,dummy1,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating rad_tables in make_rad_table: ",ld+1)
    deallocate(fullradtbl,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    call stop_timer(tmr_std_allocation)

  end subroutine make_rad_table_napfpao
  
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
!!   09:29, 27/11/2007 drb 
!!    Bug fix: zero dummy1bt and dummy2bt before use
!!   2008/02/10 ast
!!    Added timers
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
    real(double), allocatable, dimension(:) :: ol_out
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
    call start_timer(tmr_std_allocation)
    allocate(fullradtbl(n12/4,ld+1),fullradtbl_ke(n12/4,ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    allocate(rad_tables(count)%rad_tbls(ld+1), rad_tables_ke(count)%rad_tbls(ld+1),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating rad_tables in make_rad_table: ",ld+1)
    call stop_timer(tmr_std_allocation)
    rcut1 = delr1*(npts1-1)
    rcut2 = delr2*(npts2-1)
    call start_timer(tmr_std_allocation)
    allocate(dummy1(n12),dummy2(n12),dummy1bt(n12/2),dummy2bt(n12/2),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummy in make_rad_table: ",n12)
    allocate(dummyprod(n12/2),dummyprod_ke(n12/2),dumout(n12/4),ol_out(n12/4),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating dummpy outs in make_rad_table: ",n12/4)
    call stop_timer(tmr_std_allocation)
    rcut1_fake = (n12-1)*d12
    !RC have a problem with k space arrays not having sufficient size here.
    dummy1(1:npts1) = table1(1:npts1)
    dummy2(1:npts2) = table2(1:npts2)
    dummy1(npts1+1:n12) = zero
    dummy2(npts2+1:n12) = zero
    dummy1bt = zero
    dummy2bt = zero
    call matcharrays(dummy1,npts1,delr1,dummy2,npts2,delr2,d12,n1_new&
         &,n2_new,n12) 
    !loop to calculate spherical Bessel transforms of basis
    !function 
    call bessloop(dummy1,l1,n1_new,n12,d12,rcut1_fake,dummy1bt,1)
    call bessloop(dummy2,l2,n2_new,n12,d12,rcut1_fake,dummy2bt,1)
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
       
       call bessloop(dummyprod,i,n12/2,n12/2,deltak,kcut,dumout,0)
       fullradtbl(1:n12/4,h) = dumout(1:n12/4)*factor
       !RC resetting to zero(really just for peace of mind)
       dumout = 0.0_double
       call bessloop(dummyprod_ke,i,n12/2,n12/2,deltak,kcut,dumout,0)
       fullradtbl_ke(1:n12/4,h) = dumout(1:n12/4)*factor
       h = h+1
    enddo
    npts = 1+((rcut1+rcut2)*kcut/twopi)
    del_r = twopi/(kcut+deltak) 
    
    call store_supp_ke_tables(npts,del_r,ld,fullradtbl,fullradtbl_ke,n12/4,count)
    call start_timer(tmr_std_allocation)
    deallocate(ol_out,dumout,dummyprod_ke,dummyprod,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating dummpy outs in make_rad_table: ",n12/4)
    deallocate(dummy2bt, dummy1bt,dummy2,dummy1,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating rad_tables in make_rad_table: ",ld+1)
    deallocate(fullradtbl_ke,fullradtbl,STAT=stat)    
    if(stat/=0) call cq_abort("Error deallocating fullradtbl in make_rad_table: ",n12/4,ld+1)
    call stop_timer(tmr_std_allocation)
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
!!   2008/02/10 ast
!!    Added timers
!!  SOURCE
!!
  subroutine store_supp_ke_tables(npnts, delta_r, num_l, fullradtbl, &
                                  fullradtbl_ke, nsize, count)
    use datatypes
    use cubic_spline_routines !,ONLY : spline_new
    use ol_int_datatypes !,ONLY : rad_tables, rad_tables_ke
    use GenComms,       only: cq_abort
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none
    !routine to copy the radial table information from the dummy structures
    !fullradtbl and fullradtbl_ke into the storage types defined in ol_int_dataypes.
    !module.f90
    integer, intent(in) :: npnts,num_l,count,nsize
    real(double), intent(in) :: delta_r
    real(double), intent(in), dimension(nsize,num_l+1) :: fullradtbl, fullradtbl_ke
    real(double) :: yp1, ypn
    real(double), dimension(:), allocatable :: xin, y2
    integer :: i, stat

    allocate(xin(npnts), y2(npnts), STAT=stat)
    if (stat /= 0) call cq_abort("store_supp_ke_tables: Error alloc mem: ", npnts)
    call reg_alloc_mem(area_basis, 2*npnts, type_dbl)

    !defining x axis array for the spline routines
    do i = 1, npnts
       xin(i) = (i-1)*delta_r
    enddo
    
    do i = 1, num_l+1
       call start_timer(tmr_std_allocation)
       allocate(rad_tables(count)%rad_tbls(i)%arr_vals(1:npnts),&
            &rad_tables_ke(count)%rad_tbls(i)%arr_vals(1:npnts))
       allocate(rad_tables(count)%rad_tbls(i)%arr_vals2(1:npnts),&
            &rad_tables_ke(count)%rad_tbls(i)%arr_vals2(1:npnts))
       call stop_timer(tmr_std_allocation)

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

    deallocate(xin, y2, STAT=stat)
    if (stat /= 0) call cq_abort("store_supp_ke_tables: Error dealloc mem")
    call reg_dealloc_mem(area_basis, 2*npnts, type_dbl)

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
!!   2008/02/10 ast
!!    Added timers
!!  SOURCE
!!

  subroutine store_nlpf_pao_tables(npnts,delta_r,num_l,fullradtbl,nsize,count)
    use datatypes
    use ol_int_datatypes !,ONLY : rad_tables_nlpf_pao
    use cubic_spline_routines !,ONLY : spline_new
    use GenComms,       only: cq_abort
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    implicit none
    !routine to store the non-local projector function/pao 
    !overlap radial tables in the appropriate derived type.
    integer, intent(in) :: npnts,num_l,count,nsize
    real(double), intent(in) :: delta_r
    real(double), intent(in), dimension(nsize,num_l+1) :: fullradtbl
    real(double) :: yp1, ypn
    real(double), dimension(:), allocatable :: xin, y2
    integer :: i, stat
    
    allocate(xin(npnts), y2(npnts), STAT=stat)
    if (stat /= 0) call cq_abort("store_nlpf_pao_tables: Error alloc mem: ", npnts)
    call reg_alloc_mem(area_basis, 2*npnts, type_dbl)

    !defining x axis array for the spline routines
    do i = 1, npnts
       xin(i) = (i-1)*delta_r
    enddo
    
    do i = 1, num_l+1
       call start_timer(tmr_std_allocation)
       allocate(rad_tables_nlpf_pao(count)%rad_tbls(i)%arr_vals(1:npnts),&
            &rad_tables_nlpf_pao(count)%rad_tbls(i)%arr_vals2(1:npnts))
       call stop_timer(tmr_std_allocation)
      
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
    
    deallocate(xin, y2, STAT=stat)
    if (stat /= 0) call cq_abort("store_nlpf_pao_tables: Error dealloc mem")
    call reg_dealloc_mem(area_basis, 2*npnts, type_dbl)

  end subroutine store_nlpf_pao_tables
  !!***

  subroutine store_table(npnts,delta_r,num_l,fullradtbl,nsize,count,rad_table)
    use datatypes
    use ol_int_datatypes, ONLY : ol_integral
    use cubic_spline_routines !,ONLY : spline_new
    use GenComms,       only: cq_abort
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    implicit none
    !routine to store the non-local projector function/pao 
    !overlap radial tables in the appropriate derived type.
    integer, intent(in) :: npnts,num_l,count,nsize
    real(double), intent(in) :: delta_r
    real(double), intent(in), dimension(nsize,num_l+1) :: fullradtbl
    type(ol_integral), dimension(:) :: rad_table
    
    real(double) :: yp1, ypn
    real(double), dimension(:), allocatable :: xin, y2
    integer :: i, stat
    
    allocate(xin(npnts), y2(npnts), STAT=stat)
    if (stat /= 0) call cq_abort("store_table: Error alloc mem: ", npnts)
    call reg_alloc_mem(area_basis, 2*npnts, type_dbl)

    !defining x axis array for the spline routines
    do i = 1, npnts
       xin(i) = (i-1)*delta_r
    enddo
    
    do i = 1, num_l+1
       call start_timer(tmr_std_allocation)
       allocate(rad_table(count)%rad_tbls(i)%arr_vals(1:npnts), rad_table(count)%rad_tbls(i)%arr_vals2(1:npnts))
       call stop_timer(tmr_std_allocation)
      
       rad_table(count)%rad_tbls(i)%arr_vals(1:npnts) = fullradtbl(1:npnts,i)
       
       !RC calculating and storing tables of 2nd derivatives here for
       !pao-NA-pao case 
       y2 = 0.0_double
       yp1 = (fullradtbl(2,i)-fullradtbl(1,i))/delta_r
       ypn = (fullradtbl(npnts,i)-fullradtbl(npnts-1,i))/delta_r
       call spline_new(xin,fullradtbl(1:npnts,i),npnts,yp1,ypn,y2)
       rad_table(count)%rad_tbls(i)%arr_vals2(1:npnts) = y2(1:npnts)
       
       !RC finished storing arrays of 2nd derivatives
       rad_table(count)%rad_tbls(i)%npnts = npnts
       rad_table(count)%rad_tbls(i)%del_x = delta_r
       
    enddo
    
    deallocate(xin, y2, STAT=stat)
    if (stat /= 0) call cq_abort("store_table: Error dealloc mem")
    call reg_dealloc_mem(area_basis, 2*npnts, type_dbl)

  end subroutine store_table
  
  ! Neutral atom Projector functions
  subroutine store_napf_pao_tables(npnts,delta_r,num_l,fullradtbl,nsize,count)
    use datatypes
    use ol_int_datatypes !,ONLY : rad_tables_napf_pao
    use cubic_spline_routines !,ONLY : spline_new
    use GenComms,       only: cq_abort
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    implicit none
    !routine to store the non-local projector function/pao 
    !overlap radial tables in the appropriate derived type.
    integer, intent(in) :: npnts,num_l,count,nsize
    real(double), intent(in) :: delta_r
    real(double), intent(in), dimension(nsize,num_l+1) :: fullradtbl
    real(double) :: yp1, ypn
    real(double), dimension(:), allocatable :: xin, y2
    integer :: i, stat
    
    allocate(xin(npnts), y2(npnts), STAT=stat)
    if (stat /= 0) call cq_abort("store_napf_pao_tables: Error alloc mem: ", npnts)
    call reg_alloc_mem(area_basis, 2*npnts, type_dbl)

    !defining x axis array for the spline routines
    do i = 1, npnts
       xin(i) = (i-1)*delta_r
    enddo
    
    do i = 1, num_l+1
       call start_timer(tmr_std_allocation)
       allocate(rad_tables_napf_pao(count)%rad_tbls(i)%arr_vals(1:npnts),&
            &rad_tables_napf_pao(count)%rad_tbls(i)%arr_vals2(1:npnts))
       call stop_timer(tmr_std_allocation)
      
       rad_tables_napf_pao(count)%rad_tbls(i)%arr_vals(1:npnts) = fullradtbl(1:npnts,i)
       
       !RC calculating and storing tables of 2nd derivatives here for
       !napf_pao case 
       y2 = 0.0_double
       yp1 = (fullradtbl(2,i)-fullradtbl(1,i))/delta_r
       ypn = (fullradtbl(npnts,i)-fullradtbl(npnts-1,i))/delta_r
       call spline_new(xin,fullradtbl(1:npnts,i),npnts,yp1,ypn,y2)
       rad_tables_napf_pao(count)%rad_tbls(i)%arr_vals2(1:npnts) = y2(1:npnts)
       
       !RC finished storing arrays of 2nd derivatives
       rad_tables_napf_pao(count)%rad_tbls(i)%npnts = npnts
       rad_tables_napf_pao(count)%rad_tbls(i)%del_x = delta_r
       
    enddo
    
    deallocate(xin, y2, STAT=stat)
    if (stat /= 0) call cq_abort("store_napf_pao_tables: Error dealloc mem")
    call reg_dealloc_mem(area_basis, 2*npnts, type_dbl)

  end subroutine store_napf_pao_tables

  
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
!!   2008/02/10 ast
!!    Added timers
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
    
    call start_timer(tmr_std_basis)
    !RC fixing routine to read in k space parameters
    call get_max_paoparams(lmax,nzmax)
    nspecies = size(pao)
    call start_timer(tmr_std_allocation)
    allocate(ol_index(nspecies,nspecies,1:nzmax,1:nzmax,0:lmax,0:lmax))
    call stop_timer(tmr_std_allocation)
    !write(io_lun,*) lmax,nzmax, 'ol_index allocated'
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
                            
                            call start_timer(tmr_std_allocation)
                            allocate(fire1(n1),fire2(n2))
                            call stop_timer(tmr_std_allocation)
                            !RATHIN switching off for Gaussian testing
                            call unnorm_siesta_tbl(fire1,fire2,n1,n2,&
                                 &pao(nsp1)%angmom(l1)%zeta(nz1)%table,pao(nsp2)%angmom(l2)&
                                 &%zeta(nz2)%table,del1,del2,l1,l2)
                            call make_rad_table_supp_ke(n1,l1,del1,n2,l2,del2,fire1,fire2,&
                                 &count,del_k,kcut)
                            call start_timer(tmr_std_allocation)
                            deallocate(fire1,fire2)
                            call stop_timer(tmr_std_allocation)
                            !RC N.B. WE HAVE NOT DONE THE BELOW FOR KE TABLES!
                            rad_tables(count)%no_of_lvals = min(l1,l2)+1
                            call start_timer(tmr_std_allocation)
                            allocate(rad_tables(count)%l_values(min(l1,l2)+1))
                            call stop_timer(tmr_std_allocation)
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
          !write(io_lun,*) 'now allocating rad_tables(ke) storage type', ' count =', count
          call start_timer(tmr_std_allocation)
          allocate(rad_tables(count))
          allocate(rad_tables_ke(count))
          call stop_timer(tmr_std_allocation)
          count = 1
       else
          continue
       endif
    enddo !corresponding to the i counter
    !all radial tables have now been generated for pao_pao and pao_ke_pao
    call stop_timer(tmr_std_basis)
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
!!   2008/02/10 ast
!!    Added timers
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
    
    call start_timer(tmr_std_basis)
    !setting allocation bounds for indexing array
    call get_max_pao_nlpfparams(lmax,nzmax)
    nspecies = size(pao)
    !write(io_lun,*) 'allocating indexing array'
    call start_timer(tmr_std_allocation)
    allocate(ol_index_nlpf_pao(nspecies,nspecies,1:nzmax,1:nzmax,0:lmax,0:lmax))
    call stop_timer(tmr_std_allocation)
    !write(io_lun,*) lmax,nzmax, 'ol_index_nlpf_pao allocated'
    
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
                            
                            call start_timer(tmr_std_allocation)
                            allocate(fire1(n1),fire2(n2))
                            call stop_timer(tmr_std_allocation)
                            call unnorm_siesta_tbl(fire1,fire2,n1,n2,&
                                 &pao(nsp1)%angmom(l1)%zeta(nz1)%table,pseudo(nsp2)%&
                                 &pjnl(i_pfnl)%f,del1,del2,l1,l2)
                            call make_rad_table_nlpfpao(n1,l1,del1,n2,l2,del2,fire1,fire2,count,del_k&
                                 &,kcut)
                            call my_barrier()
                            call start_timer(tmr_std_allocation)
                            deallocate(fire1,fire2)
                            call stop_timer(tmr_std_allocation)
                            rad_tables_nlpf_pao(count)%no_of_lvals = min(l1,l2)+1
                            call start_timer(tmr_std_allocation)
                            allocate(rad_tables_nlpf_pao(count)%l_values(min(l1,l2)+1))
                            call stop_timer(tmr_std_allocation)
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
          !write(io_lun,*) 'now allocating rad_tables storage type', 'count =', count
          call start_timer(tmr_std_allocation)
          allocate(rad_tables_nlpf_pao(count))
          call stop_timer(tmr_std_allocation)
          count = 1
       else
          continue
       endif
    enddo !corresponding to the i counter
    call stop_timer(tmr_std_basis)
    
  end subroutine gen_nlpf_supp_tbls
  !!***

  ! Neutral atom Projector functions
  subroutine gen_paoNApao_tbls(inode,ionode)
    use datatypes
    use pseudo_tm_info, ONLY: pseudo
    use ol_int_datatypes ,ONLY : rad_tables_paoNApao, rad_tables_paopaoNA, ol_index_paopao!, &
         !rad_tables_paodNApao, rad_tables_paoNAdpao, rad_tables_paodpaoNA
    use pao_format !,ONLY : pao
    use GenComms, ONLY: cq_abort, gcopy, myid, my_barrier
    implicit none
    !routine to loop over the paos and non-local projector functions
    !and generate the array of radial tables
    integer :: nsp1,nsp2,nz1,nz2,l1,l2,n1,n2,ke_flag,i,j,k,count, count2
    integer :: n_pfnas,i_pfna
    integer :: smallcount,nspecies,lmax,nzmax,inode,ionode,lun
    real(double) :: del1,del2,del
    real(double), pointer, dimension(:) :: fire1,fire2,fire3,fire4
    
    call start_timer(tmr_std_basis)
    !setting allocation bounds for indexing array
    call get_max_paoparams(lmax,nzmax)
    nspecies = size(pao)
    call start_timer(tmr_std_allocation)
    allocate(ol_index_paopao(nspecies,nspecies,1:nzmax,1:nzmax,0:lmax,0:lmax))
    nspecies = size(pao)
    !write(io_lun,*) 'allocating indexing array'
    
    count = 1
    do i=1,2
       do nsp1 = 1,nspecies
          do nsp2 = 1,nspecies
             do l1 = 0, pao(nsp1)%greatest_angmom
                do l2 = 0, paoVNA(nsp2)%greatest_angmom
                   do nz1 = 1, pao(nsp1)%angmom(l1)%n_zeta_in_angmom
                      do nz2 = 1, paoVNA(nsp2)%angmom(l2)%n_zeta_in_angmom
                         if(i.eq.1) then !allocate storage arrays on first iteration
                            count = count+1
                         else
                            n1 = pao(nsp1)%angmom(l1)%zeta(nz1)%length 
                            n2 = paoVNA(nsp2)%angmom(l2)%zeta(nz2)%length!pseudo(nsp2)%pjna(count2)%n
                            del1 = pao(nsp1)%angmom(l1)%zeta(nz1)%delta
                            del2 = paoVNA(nsp2)%angmom(l2)%zeta(nz2)%delta!pseudo(nsp2)%pjna(count2)%delta
                            call start_timer(tmr_std_allocation)
                            allocate(fire1(n1),fire2(n2),fire3(n2),fire4(n2))
                            call stop_timer(tmr_std_allocation)
                            call unnorm_siesta_tbl(fire1,fire2,n1,n2,&
                                 pao(nsp1)%angmom(l1)%zeta(nz1)%table,&
                                 paoVNA(nsp2)%angmom(l2)%zeta(nz2)%table,del1,del2,l1,l2)
                            call make_rad_table_paoNApao(n1,l1,del1,n2,l2,del2,fire1,fire2,&
                                 count,del_k,kcut)
                            call my_barrier()
                            call start_timer(tmr_std_allocation)
                            deallocate(fire1,fire2,fire3,fire4)
                            call stop_timer(tmr_std_allocation)
                            rad_tables_paoNApao(count)%no_of_lvals = min(l1,l2)+1
                            call start_timer(tmr_std_allocation)
                            allocate(rad_tables_paoNApao(count)%l_values(min(l1,l2)+1))
                            call stop_timer(tmr_std_allocation)
                            smallcount = 1
                            do k=abs(l1-l2),l1+l2,2
                               rad_tables_paoNApao(count)%l_values(smallcount)=k
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
       if(i==1) then
          !write(io_lun,*) 'now allocating rad_tables storage type', 'count =', count
          call start_timer(tmr_std_allocation)
          allocate(rad_tables_paoNApao(count))
          call stop_timer(tmr_std_allocation)
          count = 1
       endif
    enddo !corresponding to the i counter
    ! Now find < ia | V_j | ib >
    count = 1
    do i=1,2
       do nsp1 = 1,nspecies ! i
          do nsp2 = 1,nspecies ! j
             do l1 = 0, paopao(nsp1)%greatest_angmom ! ia
                do l2 = 0, paopao(nsp1)%greatest_angmom ! ib ! Deliberate use of nsp1 ! 
                   do nz1 = 1, paopao(nsp1)%angmom(l1)%n_zeta_in_angmom ! ia
                      do nz2 = 1, paopao(nsp1)%angmom(l2)%n_zeta_in_angmom ! ib ! Deliberate use of nsp1 ! 
                         if(i.eq.1) then !allocate storage arrays on first iteration
                            count = count+1
                         else
                            ol_index_paopao(nsp1,nsp2,nz1,nz2,l1,l2) = count
                            n1 = paopao(nsp1)%angmom(l1)%zeta(nz1)%length 
                            n2 = pseudo(nsp2)%vna%n
                            del1 = paopao(nsp1)%angmom(l1)%zeta(nz1)%delta
                            del2 = pseudo(nsp2)%vna%delta
                            call start_timer(tmr_std_allocation)
                            allocate(fire1(n1),fire2(n1),fire3(n1)) ! Again deliberate but n1=n2
                            call stop_timer(tmr_std_allocation)
                            call unnorm_siesta_tbl(fire1,fire2,n1,n1,&
                                 paopao(nsp1)%angmom(l1)%zeta(nz1)%table,&
                                 paopao(nsp1)%angmom(l2)%zeta(nz2)%table,del1,del1,l1,l2)
                            fire1(:) = fire1(:)*fire2(:)
                            deallocate(fire2)
                            allocate(fire2(n2))
                            fire2(:) = pseudo(nsp2)%vna%f(:)
                            call make_rad_table_paopaoNA(n1,l1,del1,n2,l2,del2,fire1,fire2,&!fire3,&
                                 count,del_k ,kcut)
                            call my_barrier()
                            call start_timer(tmr_std_allocation)
                            deallocate(fire1,fire2,fire3)
                            call stop_timer(tmr_std_allocation)
                            rad_tables_paopaoNA(count)%no_of_lvals = min(l1,l2)+1
                            call start_timer(tmr_std_allocation)
                            allocate(rad_tables_paopaoNA(count)%l_values(min(l1,l2)+1))
                            call stop_timer(tmr_std_allocation)
                            smallcount = 1
                            do k=abs(l1-l2),l1+l2,2
                               rad_tables_paopaoNA(count)%l_values(smallcount)=k
                               smallcount = smallcount+1
                            enddo
                            count = count+1
                         end if
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
       if(i==1) then
          !write(io_lun,*) 'now allocating rad_tables storage type', 'count =', count
          call start_timer(tmr_std_allocation)
          allocate(rad_tables_paopaoNA(count))
          call stop_timer(tmr_std_allocation)
          count = 1
       endif
    enddo !corresponding to the i counter
    call stop_timer(tmr_std_basis)

  end subroutine gen_paoNApao_tbls

  !!****f* make_rad_tables/gen_napf_supp_tbls *
  !!
  !!  NAME 
  !!   gen_napf_supp_tbls
  !!  USAGE
  !!  gen_napf_supp_tbls(inode,ionode)
  !!
  !!  PURPOSE
  !!   Loops over all NAPFs and paos, calling make_rad_table and then storing and indexing
  !!   radial tables produced within the rad_tables_napf_pao array. 
  !!
  !!   Essentially identical to gen_nlpf_supp_tbls  
  !!  INPUTS
  !!
  !!  USES
  !!   datatypes, ol_int_datatypes, pao_format
  !!  AUTHOR
  !!   D. R. Bowler (following R Choudhury and NW (Mizuho))
  !!  CREATION DATE
  !!   2017/06/01
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  subroutine gen_napf_supp_tbls(inode,ionode)
    use datatypes
    use pseudo_tm_info, ONLY: pseudo
    use ol_int_datatypes !,ONLY : rad_tables_nlpf_pao,ol_index_nlpf_pao
    use pao_format !,ONLY : pao
    use GenComms, ONLY: cq_abort, gcopy, myid, my_barrier
    implicit none
    !routine to loop over the paos and non-local projector functions
    !and generate the array of radial tables
    integer :: nsp1,nsp2,nz1,nz2,l1,l2,n1,n2,ke_flag,i,j,k,count
    integer :: n_pfnas,i_pfna
    integer :: smallcount,nspecies,lmax,nzmax,inode,ionode,lun
    real(double) :: del1,del2,del
    real(double), pointer, dimension(:) :: fire1,fire2

    call start_timer(tmr_std_basis)
    !setting allocation bounds for indexing array
    call get_max_pao_napfparams(lmax,nzmax)
    nspecies = size(pao)
    !write(io_lun,*) 'allocating indexing array'
    call start_timer(tmr_std_allocation)
    allocate(ol_index_napf_pao(nspecies,nspecies,1:nzmax,1:nzmax,0:lmax,0:lmax))
    call stop_timer(tmr_std_allocation)
    !write(io_lun,*) lmax,nzmax, 'ol_index_napf_pao allocated'

    count = 1
    do i=1,2
       do nsp1 = 1,nspecies
          do l1 = 0, pao(nsp1)%greatest_angmom
             !do l2 = 0, pao(nsp2)%greatest_angmom
             do nz1 = 1, pao(nsp1)%angmom(l1)%n_zeta_in_angmom
                do nsp2 = 1,nspecies
                   n_pfnas = pseudo(nsp2)%n_pjna
                   if(n_pfnas > 0) then
                      do i_pfna = 1, n_pfnas !1,0 should abort loop
                         l2 = pseudo(nsp2)%pjna_l(i_pfna)
                         nz2 = pseudo(nsp2)%pjna_n(i_pfna)
                         if(i.eq.1) then !allocate storage arrays on first iteration
                            count = count+1
                         else
                            ol_index_napf_pao(nsp1,nsp2,nz1,nz2,l1,l2) = count
                            !write(47+myid,*) nsp1,nsp2,nz1,nz2
                            !write(47+myid,*) l1,l2,count
                            ! now calculating radial tables
                            n1 = pao(nsp1)%angmom(l1)%zeta(nz1)%length 
                            n2 = pseudo(nsp2)%pjna(i_pfna)%n
                            del1 = pao(nsp1)%angmom(l1)%zeta(nz1)%delta
                            del2 = pseudo(nsp2)%pjna(i_pfna)%delta

                            call start_timer(tmr_std_allocation)
                            allocate(fire1(n1),fire2(n2))
                            call stop_timer(tmr_std_allocation)
                            call unnorm_siesta_tbl(fire1,fire2,n1,n2,&
                                 &pao(nsp1)%angmom(l1)%zeta(nz1)%table,pseudo(nsp2)%&
                                 &pjna(i_pfna)%f,del1,del2,l1,l2)
                            ! NB this line restores the NA projector - we don't need the r**l scaling
                            fire2 = pseudo(nsp2)%pjna(i_pfna)%f
                            call make_rad_table_napfpao(n1,l1,del1,n2,l2,del2,fire1,fire2,count,del_k&
                                 &,kcut)
                            call my_barrier()
                            call start_timer(tmr_std_allocation)
                            deallocate(fire1,fire2)
                            call stop_timer(tmr_std_allocation)
                            rad_tables_napf_pao(count)%no_of_lvals = min(l1,l2)+1
                            call start_timer(tmr_std_allocation)
                            allocate(rad_tables_napf_pao(count)%l_values(min(l1,l2)+1))
                            call stop_timer(tmr_std_allocation)
                            smallcount = 1
                            do k=abs(l1-l2),l1+l2,2
                               rad_tables_napf_pao(count)%l_values(smallcount)=k
                               smallcount = smallcount+1
                            enddo
                            count = count+1
                         endif
                      enddo
                   endif ! (n_pfnas > 0) 
                enddo
             enddo
          enddo
       enddo
       if(i.eq.1) then
          !write(io_lun,*) 'now allocating rad_tables storage type', 'count =', count
          call start_timer(tmr_std_allocation)
          allocate(rad_tables_napf_pao(count))
          call stop_timer(tmr_std_allocation)
          count = 1
       else
          continue
       endif
    enddo !corresponding to the i counter
    call stop_timer(tmr_std_basis)

  end subroutine gen_napf_supp_tbls
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
!!   2015/06/10 16:01 dave
!!    Fixed problems defining maximum l and n_zeta
!!   2017/02/24 15:36 dave
!!    Changed to check max l and zeta for both paos and non-local pseudopotential projectors
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
     if(lmax < pao(n_sp)%greatest_angmom) lmax = pao(n_sp)%greatest_angmom
     if(lmax < pseudo(n_sp)%lmax) lmax = pseudo(n_sp)%lmax

     do l = 0, pao(n_sp)%greatest_angmom
        if(nzmax < pao(n_sp)%angmom(l)%n_zeta_in_angmom) nzmax = pao(n_sp)%angmom(l)%n_zeta_in_angmom
     enddo
     do l = 1,pseudo(n_sp)%n_pjnl
        if(nzmax<pseudo(n_sp)%pjnl_n(l)) nzmax = pseudo(n_sp)%pjnl_n(l)
     end do
     
  enddo
  
end subroutine get_max_pao_nlpfparams
!!***

! Neutral atom Projector functions
subroutine get_max_pao_napfparams(lmax,nzmax)
  use datatypes
  use pao_format !, ONLY ; pao
  use pseudo_tm_info !,ONLY : pseudo
  use pseudopotential_common, ONLY: maxL_neutral_atom_projector
  implicit none
  !code to analyse pao AND nlpf information and retrieve 
  !the maximum angular momentum and zeta value associated with
  !any of the paos/nlpfs
  integer, intent(out) :: lmax,nzmax
  integer :: n_sp,n_zeta,nt,l,nspecies_tot

  !loop to figure out lmax and nzmax
  lmax = maxL_neutral_atom_projector!0 !initialising
  nzmax = 0 !initialising

  nspecies_tot = size(pao)
  do n_sp = 1,nspecies_tot
     if(lmax < pao(n_sp)%greatest_angmom) lmax = pao(n_sp)%greatest_angmom
     !if(lmax < pseudo(n_sp)%pjna_l) lmax = pseudo(n_sp)%pjna_l

     do l = 0, pao(n_sp)%greatest_angmom
        if(nzmax<pao(n_sp)%angmom(l)%n_zeta_in_angmom) nzmax = pao(n_sp)%angmom(l)%n_zeta_in_angmom
        !if(nzmax < max(pao(n_sp)%angmom(l)%n_zeta_in_angmom,!pseudo(n_sp)%n_pjna)) &
        !     nzmax = max(pao(n_sp)%angmom(l)%n_zeta_in_angmom,!pseudo(n_sp)%n_pjna)
     enddo
     do nt = 1,pseudo(n_sp)%n_pjna
        if(nzmax<pseudo(n_sp)%pjna_n(nt)) nzmax = pseudo(n_sp)%pjna_n(nt)
     end do
  enddo

end subroutine get_max_pao_napfparams

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

end module make_rad_tables
