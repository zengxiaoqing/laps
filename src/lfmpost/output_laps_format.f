!dis   
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis    
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis   

  SUBROUTINE output_laps_format ( ht, u3, v3, w3, om, t3, sh, rh3, lwc, ice, &
                                  rai, sno, pic, ref, pty, u, v, w, t, td, &
                                  rh, lcb, lct, msl, p, ps, lil, tpw, r01, rto, &
                                  s01, sto, th, the, pbe, nbe, lcv, cce, lmt, &
                                  lmr, llr, spt, lhe, li, hi, vis, terdot, &
                                  lwout,swout,shflux,lhflux,pblhgt,ground_t,&
                                  press_levels, &
                                  lfmprd_dir, laps_data_root, domnum, &
                                  laps_reftime, laps_valtime, nx, ny, nz, &
                                  realtime )

    ! Creates the LAPS *.fua and *.fsf files for model output.  The names of the
    ! input variables correspond to their netCDF names found in the fua.cdl
    ! and fsf.cdl files.
    !
    ! Use of this routine requires a valid MM5_DATA_ROOT that contains
    ! a subdirectory in MM5_DATA_ROOT/mm5prd/dxx/fua (and fsf), where the
    ! xx is the two digit domain number.  The fsf.cdl and fua.cdl files
    ! in LAPS_DATA_ROOT/cdl should also have dimensions equal to nx/ny.
    ! 
    ! Basically, this routine assumes the model output grids are identical
    ! (projection, resolution, dimensions, etc.) to the domain defined
    ! in LAPS_DATA_ROOT/static/nest7grid.parms.
    !
    !
    !
    ! History
    ! =======
    ! Initial version:  Brent Shaw, NOAA/FSL 5 Jan 01
    ! Modified to output files to domain specific directories in
    ! MM5_DATA_ROOT instead of LAPS_DATA_ROOT.  2 Nov 01

    ! Note that the only variable that is supported in the netCDF files but
    ! not produced by this routine is the fire index.  

    IMPLICIT NONE
   
    ! Argument declarations
    INTEGER, INTENT(IN)             :: nx
    INTEGER, INTENT(IN)             :: ny
    INTEGER, INTENT(IN)             :: nz
    REAL, INTENT(IN)                :: ht    ( nx , ny , nz )
    REAL, INTENT(IN)                :: u3    ( nx , ny , nz )
    REAL, INTENT(IN)                :: v3    ( nx , ny , nz )
    REAL, INTENT(IN)                :: w3    ( nx , ny , nz )
    REAL, INTENT(IN)                :: om    ( nx , ny , nz )
    REAL, INTENT(IN)                :: t3    ( nx , ny , nz )
    REAL, INTENT(IN)                :: sh    ( nx , ny , nz )
    REAL, INTENT(IN)                :: rh3   ( nx , ny , nz )
    REAL, INTENT(IN)                :: lwc   ( nx , ny , nz )
    REAL, INTENT(IN)                :: ice   ( nx , ny , nz )
    REAL, INTENT(IN)                :: rai   ( nx , ny , nz )
    REAL, INTENT(IN)                :: sno   ( nx , ny , nz )
    REAL, INTENT(IN)                :: pic   ( nx , ny , nz )
    REAL, INTENT(IN)                :: ref   ( nx , ny , nz )
    REAL, INTENT(IN)                :: pty   ( nx , ny , nz )
    REAL, INTENT(IN)                :: u     ( nx , ny )
    REAL, INTENT(IN)                :: v     ( nx , ny ) 
    REAL, INTENT(IN)                :: w     ( nx , ny )
    REAL, INTENT(IN)                :: t     ( nx , ny )
    REAL, INTENT(IN)                :: td    ( nx , ny )
    REAL, INTENT(IN)                :: rh    ( nx , ny )
    REAL, INTENT(IN)                :: lcb   ( nx , ny )
    REAL, INTENT(IN)                :: lct   ( nx , ny )
    REAL, INTENT(IN)                :: msl   ( nx , ny )
    REAL, INTENT(IN)                :: p     ( nx , ny )
    REAL, INTENT(IN)                :: ps    ( nx , ny )
    REAL, INTENT(IN)                :: lil   ( nx , ny )
    REAL, INTENT(IN)                :: tpw   ( nx , ny )
    REAL, INTENT(IN)                :: r01   ( nx , ny )
    REAL, INTENT(IN)                :: rto   ( nx , ny )
    REAL, INTENT(IN)                :: s01   ( nx , ny )
    REAL, INTENT(IN)                :: sto   ( nx , ny )
    REAL, INTENT(IN)                :: th    ( nx , ny )
    REAL, INTENT(IN)                :: the   ( nx , ny )
    REAL, INTENT(IN)                :: pbe   ( nx , ny )
    REAL, INTENT(IN)                :: nbe   ( nx , ny )
    REAL, INTENT(IN)                :: lcv   ( nx , ny )
    REAL, INTENT(IN)                :: cce   ( nx , ny )
    REAL, INTENT(IN)                :: lmt   ( nx , ny )
    REAL, INTENT(IN)                :: lmr   ( nx , ny )
    REAL, INTENT(IN)                :: llr   ( nx , ny )
    REAL, INTENT(IN)                :: spt   ( nx , ny )
    REAL, INTENT(IN)                :: lhe   ( nx , ny )
    REAL, INTENT(IN)                :: li    ( nx , ny )
    REAL, INTENT(IN)                :: hi    ( nx , ny )
    REAL, INTENT(IN)                :: vis   ( nx , ny )
    REAL, INTENT(IN)                :: terdot( nx , ny )
    REAL, INTENT(IN)                :: lwout ( nx , ny )
    REAL, INTENT(IN)                :: swout ( nx , ny )
    REAL, INTENT(IN)                :: shflux( nx , ny )
    REAL, INTENT(IN)                :: lhflux( nx , ny )
    REAL, INTENT(IN)                :: pblhgt( nx , ny )
    REAL, INTENT(IN)                :: ground_t( nx , ny )
    REAL, INTENT(IN)                :: press_levels (nz)
    CHARACTER(LEN=*),INTENT(IN)     :: lfmprd_dir 
    CHARACTER(LEN=*),INTENT(IN)     :: laps_data_root
    INTEGER,INTENT(IN)              :: domnum
    INTEGER, INTENT(IN)             :: laps_reftime
    INTEGER, INTENT(IN)             :: laps_valtime
    LOGICAL, INTENT(IN)             :: realtime

    ! Locals
    CHARACTER(LEN=2)             :: domnum_str
    INTEGER, PARAMETER           :: nvar3d = 15 ! Equals # of 3d arrays above!
    INTEGER, PARAMETER           :: nvar2d = 38 ! # of 2d arrays above!
    REAL, ALLOCATABLE               :: laps_data ( : , : , : )
    INTEGER, ALLOCATABLE            :: levels ( : )
    CHARACTER(LEN=3),ALLOCATABLE    :: varname (: )
    CHARACTER(LEN=10),ALLOCATABLE   :: varunits( : )
    CHARACTER(LEN=4), ALLOCATABLE   :: varlvltype ( : )
    CHARACTER(LEN=132),ALLOCATABLE  :: varcomment ( : )
    INTEGER                         :: startind
    INTEGER                         :: stopind
    INTEGER                         :: istatus
    CHARACTER(LEN=255)              :: output_dir
    REAL, ALLOCATABLE               :: cdl_levels(:)

    INTEGER                         :: nz_one
    REAL                            :: cdl_levels_one
    CHARACTER(LEN=9)                :: asctime, gtime
    CHARACTER(LEN=4)                :: fcst_hhmm
    CHARACTER(LEN=256)              :: output_file
    CHARACTER(LEN=256)              :: donefile
    INTEGER                         :: fnlen,extlen
   

    WRITE(domnum_str, '(I2.2)') domnum
    extlen = 3
    nz_one = 1
    cdl_levels_one = 0.
    ALLOCATE (cdl_levels(nz))
    cdl_levels = press_levels(nz:1:-1)
    ! Lets make the fua file first (contains 3d variables).  First, allocate
    ! the big array and all of the metadata arrays

    PRINT *, 'Outputting LAPS format (fua/fsf) files...'
    ALLOCATE (laps_data ( nx , ny , nz * nvar3d ) )
    ALLOCATE (varname   ( nvar3d * nz ) )
    ALLOCATE (varunits  ( nvar3d * nz ) )
    ALLOCATE (varlvltype( nvar3d * nz ) )
    ALLOCATE (varcomment( nvar3d * nz ) )
    ALLOCATE (levels    ( nvar3d * nz ) )
    
    ! Initialize varcomment
    varcomment(:) = '                                                  ' //&
                    '                                                  ' // &
                    '                                '                  
    startind = 1
    stopind = nz

    ! Heights in meters
    laps_data(:,:,startind:stopind) = ht
    varname(startind:stopind) = 'HT '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast geopotential height in meters'

    ! U-wind in m/s
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = u3
    varname(startind:stopind) = 'U3 '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast U wind component in m/s'
    
    ! V-wind in m/s
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = v3
    varname(startind:stopind) = 'V3 '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast V wind component in m/s'

    ! W in m/s
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = w3
    varname(startind:stopind) = 'W3 '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast W wind component in m/s'

    ! Omega in Pa/s
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = om
    varname(startind:stopind) = 'OM '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast W wind component in Pa/s'

    ! Temperature in K
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = t3
    varname(startind:stopind) = 'T3 '   
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast temperature in K'

    ! Specific humidity in kg/kg
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = sh
    varname(startind:stopind) = 'SH '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast specific humidity in kg/kg'
 
    ! Relative humidity wrt liquid in %
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = rh3
    varname(startind:stopind) = 'RH3'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast relative humidity wrt liquid in %'

    ! Cloud liquid content in kg/m2
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = lwc
    varname(startind:stopind) = 'LWC'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast cloud liquid content in kg/m2'

    ! Cloud ice content in kg/m2
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = ice
    varname(startind:stopind) = 'ICE'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast cloud ice content in kg/m2'

    ! Rain water content in kg/m2
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = rai
    varname(startind:stopind) = 'RAI'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast rain content in kg/m2'

    ! Snow content in kg/m2
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = sno
    varname(startind:stopind) = 'SNO'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast snow content in kg/m2'

    ! Graupel content in kg/m2
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = pic
    varname(startind:stopind) = 'PIC'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast graupel content in kg/m2'

    ! Reflectivity in dBZ
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = ref
    varname(startind:stopind) = 'REF'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast radar reflectivity in dBZ'
 
    ! Coded precipitation type
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = pty
    varname(startind:stopind) = 'PTY'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Forecast precip type in coded values'

    ! Write out the 3D stuff using LAPS library routine
    output_dir = TRIM(lfmprd_dir) // '/d' // domnum_str // '/fua/'

    ! Build the output file name so we can create a "donefile" if 
    ! running in realtime mode

    CALL make_fnam_lp(laps_reftime, gtime, istatus)
    CALL make_fcst_time(laps_valtime, laps_reftime, fcst_hhmm, istatus)
    CALL cv_i4tim_asc_lp(laps_valtime,asctime,istatus)
    CALL cvt_fname_v3(output_dir,gtime,fcst_hhmm,'fua',extlen,output_file, &
         fnlen,istatus)

    PRINT *, 'Writing 3d fields to ', TRIM(output_file)
    CALL write_laps_lvls(laps_reftime, laps_valtime, output_dir, 'fua', &
                        nx,ny,nz*nvar3d,nz*nvar3d,varname,levels,&
                        varlvltype, varunits, varcomment, nz, cdl_levels,&
                        laps_data, istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, 'Error writing LAPS 3D (fua) file.'
    ENDIF
    PRINT *, 'Done writing 3d data.'
    IF ((realtime) .AND. (istatus .EQ. 1)) THEN
      donefile = TRIM(output_file) // '.done'
      OPEN(77,FILE=donefile,STATUS='UNKNOWN')
      CLOSE(77)
    ENDIF
    DEALLOCATE (laps_data)
    DEALLOCATE (varname)
    DEALLOCATE (levels)
    DEALLOCATE (varlvltype)
    DEALLOCATE (varunits)
    DEALLOCATE (varcomment)
    ! Do 2D variables
    ALLOCATE (laps_data ( nx , ny , nvar2d ) )
    ALLOCATE (varname   ( nvar2d ) )
    ALLOCATE (varlvltype( nvar2d ) ) 
    ALLOCATE (varunits  ( nvar2d ) )
    ALLOCATE (varcomment( nvar2d ) )
    ALLOCATE (levels    ( nvar2d ) )

    ! Initialize varcomment

    varcomment(:) = '                                                  '// &
                    '                                                  '// &
                    '                                '

    levels(:) = 0.
    startind = 1 
    laps_data(:,:,startind) = u
    varname(startind) = 'USF'
    varcomment(startind) = 'Forecast U wind component in m/s'
    
    startind = startind + 1
    laps_data(:,:,startind) = v
    varname(startind) = 'VSF'
    varcomment(startind) = 'Forecast V wind component in m/s'

    startind = startind + 1
    laps_data(:,:,startind) = w
    varname(startind) = 'WSF'
    varcomment(startind) = 'Forecast W wind component in m/s'

    startind = startind + 1
    laps_data(:,:,startind) = t
    varname(startind) = 'TSF'
    varcomment(startind) = 'Forecast temperature in K'

    startind = startind + 1
    laps_data(:,:,startind) = td
    varname(startind) = 'DSF'
    varcomment(startind) = 'Forecast dewpoint temperature in K'

    startind = startind + 1
    laps_data(:,:,startind) = rh
    varname(startind) = 'RH '
    varcomment(startind) = 'Forecast relative humidity wrt liquid in %'

    startind = startind + 1
    laps_data(:,:,startind) = lcb
    varname(startind) = 'LCB'
    varcomment(startind) = 'Forecast cloud base in m above sea-level'

    startind = startind + 1
    laps_data(:,:,startind) = lct
    varname(startind) = 'LCT'
    varcomment(startind) = 'Forecast cloud top in m above sea-level'

    startind = startind + 1
    laps_data(:,:,startind) = msl
    varname(startind) = 'SLP'
    varcomment(startind) = 'Forecast mean sea-level pressure in Pa'

    startind = startind + 1
    laps_data(:,:,startind) = p
    varname(startind) = 'P  '
    varcomment(startind) = 'Forecast reduced pressure in Pa'

    startind = startind + 1
    laps_data(:,:,startind) = ps
    varname(startind) = 'PSF'
    varcomment(startind) = 'Forecast surface pressure in Pa'

    startind = startind + 1
    laps_data(:,:,startind) = lil
    varname(startind) = 'LIL'
    varcomment(startind) = 'Forecast integrated liquid water depth in m'

    startind = startind + 1
    laps_data(:,:,startind) = tpw
    varname(startind) = 'TPW'
    varcomment(startind) = 'Forecast total precipitable water depth in m'

    startind = startind + 1
    laps_data(:,:,startind) = r01
    varname(startind) = 'R01'
    varcomment(startind) = 'Forecast incremental liquid precip in m'

    startind = startind + 1
    laps_data(:,:,startind) = rto
    varname(startind) = 'RTO'
    varcomment(startind) = 'Forecast total accum liquid precip in m'

    startind = startind + 1
    laps_data(:,:,startind) = s01
    varname(startind) = 'S01'
    varcomment(startind) = 'Forecast incremental snow depth in m'

    startind = startind + 1
    laps_data(:,:,startind) = sto
    varname(startind) = 'STO'
    varcomment(startind) = 'Forecast total accum snow depth in m'

    startind = startind + 1
    laps_data(:,:,startind) = th
    varname(startind) = 'TH '
    varcomment(startind) = 'Forecast potential temperature in K'

    startind = startind + 1
    laps_data(:,:,startind) = the
    varname(startind) = 'THE'
    varcomment(startind) = 'Forecast equivalent potential temperature in K'

    startind = startind + 1
    laps_data(:,:,startind) = pbe
    varname(startind) = 'PBE'
    varcomment(startind) = 'Forecast CAPE in J/kg'

    startind = startind + 1
    laps_data(:,:,startind) = nbe
    varname(startind) = 'NBE'
    varcomment(startind) = 'Forecast CIN in J/kg'

    startind = startind + 1
    laps_data(:,:,startind) = lcv
    varname(startind) = 'LCV'
    varcomment(startind) = 'Forecast cloud fraction.  Note the model ' //&
                           'currently only uses 0 or 1!'
    startind = startind + 1
    laps_data(:,:,startind) = cce
    varname(startind) = 'CCE'
    varcomment(startind) = 'Forecast cloud ceiling in m AGL'


    startind = startind + 1
    laps_data(:,:,startind) = lmt
    varname(startind) = 'LMT'
    varcomment(startind) = 'Forecast radar echo tops in m above sea-level'

    startind = startind + 1
    laps_data(:,:,startind) = lmr
    varname(startind) = 'LMR'
    varcomment(startind) = 'Forecast column max reflectivity in dBZ'

    startind = startind + 1
    laps_data(:,:,startind) = llr
    varname(startind) = 'LLR'
    varcomment(startind) = 'Forecast low-level reflectivity in dBZ'

    startind = startind + 1
    laps_data(:,:,startind) = spt
    varname(startind) = 'SPT'
    varcomment(startind) = 'Forecast precip type using coded values'

    startind = startind + 1
    laps_data(:,:,startind) = lhe
    varname(startind) = 'LHE'
    varcomment(startind) = 'Forecast storm-relative helicity in m2/s2'

    startind = startind + 1
    laps_data(:,:,startind) = li
    varname(startind) = 'LI '
    varcomment(startind) = 'Forecast lifted index in K'

    startind = startind + 1
    laps_data(:,:,startind) = hi
    varname(startind) = 'HI '
    varcomment(startind) = 'Forecast heat index in K'

    startind = startind + 1
    laps_data(:,:,startind) = vis
    varname(startind) = 'VIS'
    varcomment(startind) = 'Forecast visibility in m'

    startind = startind + 1
    laps_data(:,:,startind) = terdot
    varname(startind) = 'TER'
    varcomment(startind) = 'Model background terrain'

    startind = startind + 1
    laps_data(:,:,startind) = lwout
    varname(startind) = 'LWO'
    varcomment(startind) = 'Forecast Outgoing LW Radiation'

    startind = startind + 1
    laps_data(:,:,startind) = swout  
    varname(startind) = 'SWO'
    varcomment(startind) = 'Forecast Outgoing SW Radiation'

    startind = startind + 1
    laps_data(:,:,startind) = shflux 
    varname(startind) = 'SHF'
    varcomment(startind) = 'Forecast Sensible Heat Flux'

    startind = startind + 1
    laps_data(:,:,startind) = lhflux 
    varname(startind) = 'LHF'
    varcomment(startind) = 'Forecast Latent Heat Flux'

    startind = startind + 1
    laps_data(:,:,startind) = pblhgt 
    varname(startind) = 'BLH'
    varcomment(startind) = 'Forecast Boundary Layer Height'

    startind = startind + 1
    laps_data(:,:,startind) = ground_t
    varname(startind) = 'TGD'
    varcomment(startind) = 'Forecast Ground Temperature'

    output_dir = TRIM(lfmprd_dir) // '/d' // domnum_str // '/fsf/' 

    CALL cvt_fname_v3(output_dir,gtime,fcst_hhmm,'fsf',extlen,output_file, &
         fnlen,istatus)

    PRINT *, 'Writing 2d fields to ', TRIM(output_file)
    
    CALL write_laps_lvls(laps_reftime, laps_valtime, TRIM(output_dir), 'fsf', &
                         nx, ny, nvar2d, nvar2d, varname, levels, varlvltype, &
                         varunits, varcomment, nz_one, cdl_levels_one,  &  
                         laps_data, istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, 'Error writing LAPS 2d (fsf) file.'
    ENDIF
    PRINT *, 'Done writing 2d data.'
    IF ( (realtime) .AND. (istatus .EQ. 1)) THEN
      donefile = TRIM(output_file) // '.done'
      OPEN(77,FILE=donefile, STATUS='UNKNOWN')
      CLOSE(77) 
    ENDIF
    DEALLOCATE (laps_data)
    DEALLOCATE (varname)
    DEALLOCATE (levels)
    DEALLOCATE (varlvltype)
    DEALLOCATE (varunits)
    DEALLOCATE (varcomment)
    DEALLOCATE (cdl_levels)
    RETURN

  END SUBROUTINE output_laps_format


    
