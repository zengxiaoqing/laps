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
                                  rai, sno, pic, tke, ref, pty, u, v, w,t,td, &
                                  rh, lcb, lct, msl, p, ps, lil, tpw, r01, rto, &
                                  s01, sto, th, the, pbe, nbe, lcv, cce, lmt, &
                                  lmr, llr, spt, lhe, li, hi, vis, terdot, &
                                  lwout,swout,shflux,lhflux,pblhgt,ground_t,&
                                  upb, vpb, vnt,ham,hah,fwi, &
                                  press_levels, &
                                  lfmprd_dir, laps_data_root, domnum, &
                                  laps_reftime, laps_valtime, nx, ny, nz, &
                                  realtime,write_to_lapsdir,model_name, &
                                  make_donefile )

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
    REAL, INTENT(IN)                :: tke   ( nx , ny , nz )
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
    REAL, INTENT(IN)                :: upb ( nx , ny )
    REAL, INTENT(IN)                :: vpb ( nx , ny )
    REAL, INTENT(IN)                :: vnt   ( nx , ny )
    REAL, INTENT(IN)                :: ham   ( nx , ny )
    REAL, INTENT(IN)                :: hah   ( nx , ny )
    REAL, INTENT(IN)                :: fwi   ( nx , ny )

    REAL, INTENT(IN)                :: press_levels (nz)
    CHARACTER(LEN=*),INTENT(IN)     :: lfmprd_dir 
    CHARACTER(LEN=*),INTENT(IN)     :: laps_data_root
    INTEGER,INTENT(IN)              :: domnum
    INTEGER, INTENT(IN)             :: laps_reftime
    INTEGER, INTENT(IN)             :: laps_valtime
    LOGICAL, INTENT(IN)             :: realtime
    LOGICAL, INTENT(IN)             :: write_to_lapsdir
    CHARACTER(LEN=32), INTENT(IN)   :: model_name
    LOGICAL, INTENT(IN)             :: make_donefile
    ! Locals
    CHARACTER(LEN=2)             :: domnum_str
    INTEGER, PARAMETER           :: nvar3d = 16 ! Equals # of 3d arrays above!
    INTEGER, PARAMETER           :: nvar2d = 44 ! # of 2d arrays above!
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
    CHARACTER(LEN=255)              :: cdl_dir
    REAL, ALLOCATABLE               :: cdl_levels(:)

    INTEGER                         :: nz_one
    REAL                            :: cdl_levels_one
    CHARACTER(LEN=9)                :: asctime, gtime
    CHARACTER(LEN=4)                :: fcst_hhmm
    CHARACTER(LEN=256)              :: output_file
    CHARACTER(LEN=256)              :: donefile
    INTEGER                         :: fnlen,extlen
    INTEGER                         :: doneunit 

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
    varcomment(startind:stopind) = 'Geopotential Height                   '

    ! U-wind in m/s
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = u3
    varname(startind:stopind) = 'U3 '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'U-component Wind                '
    
    ! V-wind in m/s
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = v3
    varname(startind:stopind) = 'V3 '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'V-component wind                '

    ! W in m/s
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = w3
    varname(startind:stopind) = 'W3 '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Vertical Velocity               '

    ! Omega in Pa/s
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = om
    varname(startind:stopind) = 'OM '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Pressure Vertical Velocity       '

    ! Temperature in K
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = t3
    varname(startind:stopind) = 'T3 '   
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Temperature              '

    ! Specific humidity in kg/kg
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = sh
    varname(startind:stopind) = 'SH '
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Specific Humidity                  '
 
    ! Relative humidity wrt liquid in %
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = rh3
    varname(startind:stopind) = 'RH3'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Relative Humidity                         '

    ! Cloud liquid content in kg/m2
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = lwc
    varname(startind:stopind) = 'LWC'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Cloud Liquid Water                    '

    ! Cloud ice content in kg/m2
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = ice
    varname(startind:stopind) = 'ICE'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Cloud Ice                          '

    ! Rain water content in kg/m2
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = rai
    varname(startind:stopind) = 'RAI'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Rain Concentration            '

    ! Snow content in kg/m2
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = sno
    varname(startind:stopind) = 'SNO'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Snow Concentration            '

    ! Graupel content in kg/m2
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = pic
    varname(startind:stopind) = 'PIC'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Graupel Concentration            '

    ! Reflectivity in dBZ
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = ref
    varname(startind:stopind) = 'REF'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Sim. Radar Reflectivity           '
 
    ! Coded precipitation type
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = pty
    varname(startind:stopind) = 'PTY'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Precip. Type                        '

    ! Turbulent Kinetic Energy
    startind=stopind+1
    stopind = startind + nz - 1
    laps_data(:,:,startind:stopind) = tke
    varname(startind:stopind) = 'TKE'
    levels(startind:stopind) = NINT(press_levels)
    varcomment(startind:stopind) = 'Turbulent Kinetic Energy            '

    ! Write out the 3D stuff using LAPS library routine
    IF (.NOT. write_to_lapsdir) THEN
      output_dir = TRIM(lfmprd_dir) // '/d' // domnum_str // '/fua/'
    ELSE
      output_dir = TRIM(laps_data_root) // '/lapsprd/fua/' // &
                   TRIM(model_name) // '/'
    ENDIF
    cdl_dir = TRIM(laps_data_root) // '/cdl/'

    ! Build the output file name so we can create a "donefile" if 
    ! running in realtime mode

    CALL make_fnam_lp(laps_reftime, gtime, istatus)
    CALL make_fcst_time(laps_valtime, laps_reftime, fcst_hhmm, istatus)
    CALL cv_i4tim_asc_lp(laps_valtime,asctime,istatus)
    CALL cvt_fname_v3(output_dir,gtime,fcst_hhmm,'fua',extlen,output_file, &
         fnlen,istatus)

    PRINT *, 'Writing 3d fields to ', TRIM(output_file)
    CALL write_laps_lfm(laps_reftime, laps_valtime, output_dir, cdl_dir, &
                        'fua', &
                        nx,ny,nz*nvar3d,nz*nvar3d,varname,levels,&
                        varlvltype, varunits, varcomment, nz, cdl_levels,&
                        laps_data, istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, 'Error writing LAPS 3D (fua) file.'
    ENDIF
    PRINT *, 'Done writing 3d data.'
    IF ((realtime) .AND. (istatus .EQ. 1) .AND. (make_donefile)) THEN
      donefile = TRIM(output_file) // '.done'
      CALL get_file_unit(doneunit)
      OPEN(UNIT=doneunit,FILE=donefile,STATUS='UNKNOWN')
      CLOSE(doneunit)
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
    varcomment(startind) = 'Sfc U-Component Wind            '
    
    startind = startind + 1
    laps_data(:,:,startind) = v
    varname(startind) = 'VSF'
    varcomment(startind) = 'Sfc V-Component Wind            '

    startind = startind + 1
    laps_data(:,:,startind) = w
    varname(startind) = 'WSF'
    varcomment(startind) = 'Sfc Vertical Velocity           '

    startind = startind + 1
    laps_data(:,:,startind) = t
    varname(startind) = 'TSF'
    varcomment(startind) = 'Sfc Temperature          '

    startind = startind + 1
    laps_data(:,:,startind) = td
    varname(startind) = 'DSF'
    varcomment(startind) = 'Sfc Dewpoint Temperature          '

    startind = startind + 1
    laps_data(:,:,startind) = rh
    varname(startind) = 'RH '
    varcomment(startind) = 'Sfc Relative Humidity                     '

    startind = startind + 1
    laps_data(:,:,startind) = lcb
    varname(startind) = 'LCB'
    varcomment(startind) = 'Cloud Base ASL                          '

    startind = startind + 1
    laps_data(:,:,startind) = lct
    varname(startind) = 'LCT'
    varcomment(startind) = 'Cloud Top ASL                          '

    startind = startind + 1
    laps_data(:,:,startind) = msl
    varname(startind) = 'SLP'
    varcomment(startind) = 'Sea-level Pressure                    '

    startind = startind + 1
    laps_data(:,:,startind) = p
    varname(startind) = 'P  '
    varcomment(startind) = 'Reduced Pressure               '

    startind = startind + 1
    laps_data(:,:,startind) = ps
    varname(startind) = 'PSF'
    varcomment(startind) = 'Surface Pressure               '

    startind = startind + 1
    laps_data(:,:,startind) = lil
    varname(startind) = 'LIL'
    varcomment(startind) = 'Integrated Liquid Water                     '

    startind = startind + 1
    laps_data(:,:,startind) = tpw
    varname(startind) = 'TPW'
    varcomment(startind) = 'Total Precipitable Water                     '

    startind = startind + 1
    laps_data(:,:,startind) = r01
    varname(startind) = 'R01'
    varcomment(startind) = 'Incremental Tot. Liq. Precip           '

    startind = startind + 1
    laps_data(:,:,startind) = rto
    varname(startind) = 'RTO'
    varcomment(startind) = 'Run-total Liq. Precip Accum            '

    startind = startind + 1
    laps_data(:,:,startind) = s01
    varname(startind) = 'S01'
    varcomment(startind) = 'Incremental Snow Depth              '

    startind = startind + 1
    laps_data(:,:,startind) = sto
    varname(startind) = 'STO'
    varcomment(startind) = 'Run-total Snow Accum                '

    startind = startind + 1
    laps_data(:,:,startind) = th
    varname(startind) = 'TH '
    varcomment(startind) = 'Sfc Potential Temperature          '

    startind = startind + 1
    laps_data(:,:,startind) = the
    varname(startind) = 'THE'
    varcomment(startind) = 'Sfc Equiv. Potential Temperature              '

    startind = startind + 1
    laps_data(:,:,startind) = pbe
    varname(startind) = 'PBE'
    varcomment(startind) = 'CAPE                 '

    startind = startind + 1
    laps_data(:,:,startind) = nbe
    varname(startind) = 'NBE'
    varcomment(startind) = 'CIN                 '

    startind = startind + 1
    laps_data(:,:,startind) = lcv
    varname(startind) = 'LCV'
    varcomment(startind) = 'Cloud Fraction                          ' 

    startind = startind + 1
    laps_data(:,:,startind) = cce
    varname(startind) = 'CCE'
    varcomment(startind) = 'Cloud Ceiling AGL              '


    startind = startind + 1
    laps_data(:,:,startind) = lmt
    varname(startind) = 'LMT'
    varcomment(startind) = 'Sim. Radar Echo Tops                         '

    startind = startind + 1
    laps_data(:,:,startind) = lmr
    varname(startind) = 'LMR'
    varcomment(startind) = 'Sim. Composite Reflectivity            '

    startind = startind + 1
    laps_data(:,:,startind) = llr
    varname(startind) = 'LLR'
    varcomment(startind) = 'Sim. Sfc. Reflectivity                '

    startind = startind + 1
    laps_data(:,:,startind) = spt
    varname(startind) = 'SPT'
    varcomment(startind) = 'Sfc Precip. Type                       '

    startind = startind + 1
    laps_data(:,:,startind) = lhe
    varname(startind) = 'LHE'
    varcomment(startind) = 'Storm Relative Helicity                  '

    startind = startind + 1
    laps_data(:,:,startind) = li
    varname(startind) = 'LI '
    varcomment(startind) = 'Lifted Index              '

    startind = startind + 1
    laps_data(:,:,startind) = hi
    varname(startind) = 'HI '
    varcomment(startind) = 'Heat Index              '

    startind = startind + 1
    laps_data(:,:,startind) = vis
    varname(startind) = 'VIS'
    varcomment(startind) = 'Sfc. Visibility              '

    startind = startind + 1
    laps_data(:,:,startind) = terdot
    varname(startind) = 'TER'
    varcomment(startind) = 'Model Terrain           '

    startind = startind + 1
    laps_data(:,:,startind) = lwout
    varname(startind) = 'LWO'
    varcomment(startind) = 'Outgoing LW Radiation         '

    startind = startind + 1
    laps_data(:,:,startind) = swout  
    varname(startind) = 'SWO'
    varcomment(startind) = 'Outgoing SW Radiation         '

    startind = startind + 1
    laps_data(:,:,startind) = shflux 
    varname(startind) = 'SHF'
    varcomment(startind) = 'Sensible Heat Flux            '

    startind = startind + 1
    laps_data(:,:,startind) = lhflux 
    varname(startind) = 'LHF'
    varcomment(startind) = 'Latent Heat Flux          '

    startind = startind + 1
    laps_data(:,:,startind) = pblhgt 
    varname(startind) = 'BLH'
    varcomment(startind) = 'Boundary Layer Depth         '

    startind = startind + 1
    laps_data(:,:,startind) = ground_t
    varname(startind) = 'TGD'
    varcomment(startind) = 'Ground Temperature        '

    startind = startind + 1
    laps_data(:,:,startind) = upb
    varname(startind) = 'UPB'
    varcomment(startind) = 'U-component Wind in PBL    '

    startind = startind + 1
    laps_data(:,:,startind) = vpb
    varname(startind) = 'VPB'
    varcomment(startind) = 'V-component Wind in PBL    '

    startind = startind + 1
    laps_data(:,:,startind) = vnt
    varname(startind) = 'VNT'
    varcomment(startind) = 'Ventilation Index          '

    startind = startind + 1
    laps_data(:,:,startind) = ham
    varname(startind) = 'HAM'
    varcomment(startind) = 'Mid-Level Haines Index     '

    startind = startind + 1
    laps_data(:,:,startind) = hah
    varname(startind) = 'HAH'
    varcomment(startind) = 'High-Level Haines Index    '

    startind = startind + 1
    laps_data(:,:,startind) = fwi
    varname(startind) = 'FWI'
    varcomment(startind) = 'Fosberg Fire Wx Index     '

    IF (.NOT. write_to_lapsdir) THEN
      output_dir = TRIM(lfmprd_dir) // '/d' // domnum_str // '/fsf/' 
    ELSE
      output_dir = TRIM(laps_data_root) // '/lapsprd/fsf/' // &
                   TRIM(model_name) // '/'
    ENDIF

    CALL cvt_fname_v3(output_dir,gtime,fcst_hhmm,'fsf',extlen,output_file, &
         fnlen,istatus)

    PRINT *, 'Writing 2d fields to ', TRIM(output_file)
    
    CALL write_laps_lfm(laps_reftime, laps_valtime, TRIM(output_dir), cdl_dir, &
                         'fsf', &
                         nx, ny, nvar2d, nvar2d, varname, levels, varlvltype, &
                         varunits, varcomment, nz_one, cdl_levels_one,  &  
                         laps_data, istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, 'Error writing LAPS 2d (fsf) file.'
    ENDIF
    PRINT *, 'Done writing 2d data.'
    IF ( (realtime) .AND. (istatus .EQ. 1).AND.(make_donefile)) THEN
      donefile = TRIM(output_file) // '.done'
      CALL get_file_unit(doneunit)
      OPEN(UNIT=doneunit,FILE=donefile, STATUS='UNKNOWN')
      CLOSE(doneunit) 
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


    
