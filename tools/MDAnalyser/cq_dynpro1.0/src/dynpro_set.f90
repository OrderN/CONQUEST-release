MODULE dynpro_set

  USE kinds,         ONLY : DP 
  USE io_files,      ONLY : tmp_dir, prefix, iunpun, xmlpun, outdir

  IMPLICIT NONE
  SAVE
  
CONTAINS
  !==---------------------------------------------------------------------=
  SUBROUTINE dynpro_read_input()

    USE dynpro_inp

    IMPLICIT NONE
    
    NAMELIST /inputdynpro/ prefix, outdir, fileout, ndr, time_step, ct_atom, ldebug, &
         lcenter, lau_unit, lgrace, lsample, lrdf, lsphere, lsub, lvacf, lpsd, lnmr, &
         prefix_cqdata, lcq_coordforce, lcq_mdframes

    NAMELIST /inputsample/ sp_start, sp_end, sp_n, iprint 
    
    NAMELIST /inputrdf/ rd_atom, rmin, rmax, rint, nr
    
    NAMELIST /inputsphere/ sph_atom, sph_atn, sph_rcut, sph_n, sph_m, &
         sph_q, langle_spef, ldihed_spef, langle, angle_atom1, angle_atom2, ldist, ldist_spef, &
         dist_atom, sph_p
    
    !NAMELIST /inputvacf/  m_cof, win, n_seg, loverlap, temp, lspe, &
    !     qcfactor, lspe, ang_grid, ang_start, ang_stop, ang_test

    NAMELIST /inputvacf/ ang_grid, ang_start, ang_stop, ang_test

    NAMELIST /inputpsd/ m_cof, win, n_seg, loverlap, temp, qcfactor, lspe
    
    NAMELIST /inputdip/ ldip

    NAMELIST /inputnmr/ spin, quad, prefix_nmr, lrandom, lstat, verbose_nmr, &
         cor_function, cor_order

    NAMELIST /inputspec/ ityp_read, atm_read
    
    CALL input_from_file( )
  
    READ(5, inputdynpro, iostat=ios)

    
    IF (lsample) THEN
       REWIND(stdin)
       print*, 'read SAMPLE parameters'
       READ( 5, inputsample, iostat=ios)
    END IF

    IF (lrdf) THEN
       REWIND(stdin)
       print*, 'read RDF parameters'
       READ( 5, inputrdf, iostat=ios)
    END IF

    IF (lsphere) THEN
       REWIND(stdin)
       print*, 'read SPHERE parameters'
       READ( 5, inputsphere, iostat=ios)
    END IF

    IF (lvacf) THEN
       REWIND(stdin)     
       print*, 'read VACF parameters'
       READ( 5, inputvacf, iostat=ios)
    END IF

    IF (lpsd) THEN
       REWIND(stdin)
       print*, 'read PSD parameters'
       READ( 5, inputpsd, iostat=ios)
    END IF

    !print*, 'read_psd'
    print*, 'M_COF   ', m_cof
    !print*, 'N_SEG   ', n_seg
    !print*, 'win     ', win    
    !print*, 'loverlap', loverlap
    !print*, 'temp', temp
    !print*, 'lspe', lspe
    
    IF (lnmr) THEN
       REWIND(stdin)
       print*, 'read NMR parameters'
       READ( 5, inputnmr, iostat=ios)
    END IF

    IF (ldip) THEN
       REWIND(stdin)
       print*, 'read DIP parameters'
       READ( 5, inputdip, iostat=ios)
    END IF

    ! set input file names
    filecel = TRIM(outdir) // TRIM(prefix) // '.cel'
    filepos = TRIM(outdir) // TRIM(prefix) // '.pos'
    filefor = TRIM(outdir) // TRIM(prefix) // '.for'
    filevel = TRIM(outdir) // TRIM(prefix) // '.vel'
    fileevp = TRIM(outdir) // TRIM(prefix) // '.evp'

    if ( lcq_coordforce ) then
       filecqd = TRIM(outdir) // TRIM(prefix_cqdata) // '.dat'
    else
       filecqd = TRIM(outdir) // TRIM(prefix_cqdata)
    end if
    
!    fileeig = TRIM(outdir) // TRIM(prefix_nmr) // TRIM(prefix) // '.eig'
!    filevec = TRIM(outdir) // TRIM(prefix_nmr) // TRIM(prefix) // '.vec'

    fileeig =  TRIM(prefix_nmr)  // '.eig'
    filevec =  TRIM(prefix_nmr)  // '.vec'

    filedip = TRIM(outdir) // TRIM(prefix) // '.dip'
    filespe = TRIM(outdir) // TRIM(prefix) // '.spe'
    
    ! set output file names

    filemol = TRIM(fileout) // '.xyz'
    filesph = TRIM(fileout) // '.sph'
    filerdf = TRIM(fileout) // '.rdf'
    filesort= TRIM(fileout) // '.sor'
    fileclas= TRIM(fileout) // '.class'
    filedbg = TRIM(fileout) // '.dbg'
    filevac = TRIM(fileout) // '.vac'
    filepsd = TRIM(fileout) // '.psd'

   
  ! filetens    = TRIM(outdir) // TRIM(prefix_nmr) // TRIM(prefix) // '.tens'
  ! fileran     = TRIM(outdir) // TRIM(prefix_nmr) // TRIM(prefix) // '.rand'
  ! fileR       = TRIM(outdir) // TRIM(prefix_nmr) // TRIM(prefix) // '.R'
  ! fileRiso    = TRIM(outdir) // TRIM(prefix_nmr) // TRIM(prefix) // '.Riso'
  ! fileRant    = TRIM(outdir) // TRIM(prefix_nmr) // TRIM(prefix) // '.Rant'
  ! fileRsym    = TRIM(outdir) // TRIM(prefix_nmr) // TRIM(prefix) // '.Rsym'
  ! fileRsph    = TRIM(outdir) // TRIM(prefix_nmr) // TRIM(prefix) // '.R_kq'

  ! filetoto    = TRIM(outdir) // TRIM(prefix_nmr) // TRIM(prefix) // '.toto'


    filetens    = TRIM(prefix_nmr) // '.tens'
    fileran     = TRIM(prefix_nmr) // '.rand'
    fileR       = TRIM(prefix_nmr) // '.R'
    fileRiso    = TRIM(prefix_nmr) // '.Riso'
    fileRant    = TRIM(prefix_nmr) // '.Rant'
    fileRsym    = TRIM(prefix_nmr) // '.Rsym'
    fileRsph    = TRIM(prefix_nmr) // '.R_kq'
    filetoto    = TRIM(prefix_nmr) // '.toto'



    filesph_tmp = TRIM(fileout) // '.sph_m'
    filesph_pc  = TRIM(fileout) // '.sph_p'
    
    filedip_tcf = TRIM(fileout) // '.diptcf'
    
    filesub_out = 'sub.struct_perm.lk'
    filesub_in  = 'sub.struct_perm.in'
    
    
    IF (lgrace) THEN
       filener = TRIM(fileout) // '.nrj.agr'
    ELSE
       filener = TRIM(fileout) // '.nrj.dat'
    END IF


    !CALL input_from_file( )
    
  END SUBROUTINE dynpro_read_input

  
!!$  SUBROUTINE dynpro_read_xml()
!!$
!!$    USE iotk_module
!!$    USE xml_io_base
!!$    USE dynpro_inp
!!$    USE dynpro_mod,    ONLY : write_xml, check_file
!!$    
!!$    IMPLICIT NONE
!!$
!!$    !INTEGER, INTENT(inout) :: n_xml
!!$    INTEGER                :: i, j, k, l, m, tmp_xml
!!$
!!$    filepp = restart_dir( outdir, ndr )
!!$    !
!!$    filepp = TRIM( filepp ) // '/' // TRIM(xmlpun)
!!$    !
!!$    CALL iotk_open_read( xmlunit, file = TRIM( filepp ), BINARY = .FALSE., ROOT = attr )
!!$    CALL iotk_scan_begin( xmlunit, "STATUS" )
!!$    CALL iotk_scan_empty( xmlunit, "STEP", attr )
!!$    CALL iotk_scan_attr( attr, "ITERATION", iteration_xml )
!!$    CALL iotk_scan_dat( xmlunit, "TIME", time_xml )                     
!!$    CALL iotk_scan_end( xmlunit, "STATUS" )
!!$    CALL iotk_scan_begin( xmlunit, "IONS", FOUND = found )
!!$    IF( .NOT. found ) THEN
!!$       CALL errore( ' dynpro ', ' IONS not found in data-file.xml ', 1 )
!!$    END IF
!!$    CALL iotk_scan_dat( xmlunit, "NUMBER_OF_ATOMS", nat_xml )
!!$    CALL iotk_scan_dat( xmlunit, "NUMBER_OF_SPECIES", nsp_xml )
!!$    
!!$    ALLOCATE(ityp_xml(nat_xml))
!!$    ALLOCATE(label_xml(nat_xml))
!!$    ALLOCATE(atm_xml(nsp_xml))
!!$
!!$    ityp_xml = 0 ; label_xml = 'XX' ; atm_xml = 'XX'
!!$    
!!$    DO i = 1, nsp_xml
!!$       CALL iotk_scan_begin( xmlunit, "SPECIE" // TRIM( iotk_index( i ) ), FOUND = found )
!!$       IF( .NOT. found ) THEN
!!$          CALL errore( ' cppp ', "SPECIE" // TRIM( iotk_index( i ) ) // ' not found in data-file.xml ', 1 )
!!$       END IF
!!$       CALL iotk_scan_dat( xmlunit, "ATOM_TYPE", atm_xml(i) )!
!!$       CALL iotk_scan_end( xmlunit, "SPECIE" // TRIM( iotk_index( i ) ) )
!!$    END DO
!!$    
!!$    DO i = 1, nat_xml
!!$       CALL iotk_scan_empty( xmlunit, "ATOM" // TRIM( iotk_index( i ) ), attr )
!!$       CALL iotk_scan_attr( attr, "SPECIES", label_xml(i) )
!!$       CALL iotk_scan_attr( attr, "INDEX",   ityp_xml(i)  )
!!$    END DO
!!$    
!!$    CALL iotk_scan_end( xmlunit, "IONS" )
!!$    CALL iotk_scan_begin( xmlunit, "PLANE_WAVES" )
!!$    CALL iotk_scan_empty( xmlunit, "FFT_GRID", attr )    
!!$    CALL iotk_scan_attr( attr, "nr1", nr1 )
!!$    CALL iotk_scan_attr( attr, "nr2", nr2 )
!!$    CALL iotk_scan_attr( attr, "nr3", nr3 )
!!$    CALL iotk_scan_end( xmlunit, "PLANE_WAVES" )   
!!$    
!!$    ALLOCATE(stau_xml( 3, nat_xml))
!!$    ALLOCATE(svel_xml( 3, nat_xml))
!!$    ALLOCATE(force_xml(3, nat_xml))
!!$    
!!$    stau_xml = 0.0d0 ; svel_xml = 0.0d0 ; force_xml = 0.0d0 
!!$    
!!$    CALL iotk_scan_begin( xmlunit, "TIMESTEPS", attr )
!!$    CALL iotk_scan_begin( xmlunit, "STEP0" )
!!$    CALL iotk_scan_begin( xmlunit, "IONS_POSITIONS" )
!!$    CALL iotk_scan_dat( xmlunit, "stau", stau_xml )
!!$    CALL iotk_scan_dat( xmlunit, "svel", svel_xml )
!!$    CALL iotk_scan_dat( xmlunit, "force", force_xml )
!!$    CALL iotk_scan_end( xmlunit, "IONS_POSITIONS" )
!!$    CALL iotk_scan_begin( xmlunit, "CELL_PARAMETERS" )
!!$    CALL iotk_scan_dat( xmlunit, "ht", ht_xml )
!!$    CALL iotk_scan_end( xmlunit, "CELL_PARAMETERS" )
!!$    CALL iotk_scan_end( xmlunit, "STEP0" )
!!$    CALL iotk_scan_end( xmlunit, "TIMESTEPS" )
!!$    !  ispin = 1
!!$    CALL iotk_close_read( xmlunit )
!!$    
!!$    CALL check_file(filecel,3      ,1,i)
!!$    CALL check_file(filepos,nat_xml,1,j)
!!$    CALL check_file(filevel,nat_xml,1,k)
!!$    CALL check_file(filefor,nat_xml,1,l)
!!$    CALL check_file(fileevp,1      ,0,m)
!!$    
!!$    tmp_xml = MODULO((i+j+k+l+m), 5)
!!$
!!$    !IF (tmp_xml /= 0) THEN
!!$    !   PRINT *,"Error -DynPro read- check the lengh of the files"
!!$    !   STOP
!!$    !END IF
!!$
!!$    n_xml = i
!!$
!!$    CALL write_xml(time_xml, iteration_xml, nat_xml, nsp_xml)
!!$
!!$  END  SUBROUTINE dynpro_read_xml

END MODULE dynpro_set
