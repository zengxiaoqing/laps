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
!dis 

MODULE vis5d

! Module to contain Fortran routines needed by the model post processor
! to create Vis5D output.  This module also depends on the v5d.c 
! Vis5D API.  

  USE setup
  USE time_utils
  IMPLICIT NONE
  PRIVATE
  include 'v5df.h'
  INTEGER                               :: v5dtimes(MAXTIMES)
  INTEGER                               :: v5ddates(MAXTIMES)  
  INTEGER                               :: num2dvar
  INTEGER                               :: num3dvar
  INTEGER                               :: totalvars
  INTEGER                               :: ret      
  CHARACTER(LEN=10)                     :: varname(MAXVARS) 
  CHARACTER(LEN=255)                    :: v5dfile 
  PUBLIC v5dinit
  PUBLIC v5dout
  PUBLIC v5dend

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE v5dinit(compression)

    ! Initializes the Vis5D file.  Only argument provided is the compression
    ! factor.  The rest of the items are input from the setup module to
    ! call various other routines.

    IMPLICIT NONE

    INCLUDE 'v5df.h'
      
    ! Declare Arguments
    INTEGER, INTENT(IN)                   :: compression

    ! Local Variables
    INTEGER                               :: yyyy, yy, mm, dd
    INTEGER                               :: hh, min, doy
    INTEGER                               :: i,j,k,t
    INTEGER, PARAMETER                    :: vertical_code = 3
    REAL                                  :: levels(MAXLEVELS)
    INTEGER                               :: nlevels(MAXVARS)
    INTEGER                               :: varindex
    INTEGER                               :: v5d_proj
    REAL                                  :: h,msf,pole_lat
    REAL, PARAMETER                       :: rad_per_deg = .0174533
    REAL                                  :: v5d_proj_args(100)
    CHARACTER(LEN=9)                      :: file_date_string
    INTEGER                               :: v5d_ntimes
    CHARACTER(LEN=3)                      :: domnum_str
    !!!!! Begin Code !!!!!
    
    ! Check compression argument
    IF ( (compression .NE. 1).AND. &
         (compression .NE. 2).AND. &
         (compression .NE. 4) ) THEN
      PRINT '(A,I2)', 'V5DINIT: Invalid compression: ', compression
      PRINT '(A)', '  Must be 1, 2, or 4 bytes per point.'
      STOP
    ENDIF
    
    ! Build the date/time values for vis5d
    v5d_ntimes = 0
    DO t = 1,num_times_to_proc,time_index_inc
      v5d_ntimes = v5d_ntimes+1
      READ(times_to_proc(t),FMT='(I4.4)') yyyy
      READ(times_to_proc(t),FMT='(2X,I2.2)') yy
      READ(times_to_proc(t),FMT='(5X,I2.2)') mm
      READ(times_to_proc(t),FMT='(8X,I2.2)') dd
      READ(times_to_proc(t),FMT='(11X,I2.2)') hh
      READ(times_to_proc(t),FMT='(14x,I2.2)') min
      doy = compute_day_of_year(yyyy,mm,dd)
      v5ddates(v5d_ntimes)= yy *1000 + doy
      v5dtimes(v5d_ntimes) = hh*10000 + min*100
      PRINT '(A,I5.5,1x,I6.6)', 'Vis5D Date/Time: ', v5ddates(v5d_ntimes), &
         v5dtimes(v5d_ntimes)
    ENDDO

    ! Build the vis5d output file name

    WRITE(file_date_string, FMT='(I5.5,I4.4)') v5ddates(1), v5dtimes(1)/100
    WRITE(domnum_str, '("d",I2.2)') domain_num
    v5dfile = TRIM(lfmprd_dir) // '/' // domnum_str // &
      '/v5d/' // file_date_string // '.v5d'
    
    ! Build the levels array (descending pressures in mb)
    levels(1:kprs) = prslvl/100.

    ! Build the projection arguments for vis5d

    IF (TRIM(projection).EQ.'LAMBERT CONFORMAL') THEN
      v5d_proj = 2
      v5d_proj_args(1) = proj%truelat2
      v5d_proj_args(2) = proj%truelat1
      v5d_proj_args(3) = float(ny)-proj%polej
      v5d_proj_args(4) = proj%polei
      v5d_proj_args(5) = -proj%stdlon
      v5d_proj_args(6) = proj%dx * 0.001
    ELSE IF (TRIM(projection).EQ.'POLAR STEREOGRAPHIC')THEN
      v5d_proj = 3
      IF (proj%truelat1 .GE. 0.) THEN
        pole_lat = 90.
        h = 1.
      ELSE
        pole_lat=  -90.
        h = -1.
      ENDIF
      v5d_proj_args(1) = pole_lat
      v5d_proj_args(2) = -proj%stdlon
      v5d_proj_args(3) = float(ny)-proj%polej
      v5d_proj_args(4) = proj%polei
     
      ! Compute the grid spacing at the pole using
      ! map scale factor
      msf = (1.+h*SIN(proj%truelat1*rad_per_deg)) / &
            (1.+h*SIN(pole_lat*rad_per_deg))
      v5d_proj_args(5) = proj%dx * 0.001 / msf
      print *, 'dx at pole = ', v5d_proj_args(5)
    ELSE IF (TRIM(projection).EQ.'MERCATOR') THEN
      v5d_proj = 1
      v5d_proj_args(1) = corner_lats(2)  ! north boundary
      v5d_proj_args(2) = -corner_lons(2)  ! West boundary
      ! Compute increment between rows in degrees
      v5d_proj_args(3) = (corner_lats(2)-corner_lats(1))/(FLOAT(ny)-1.)
      ! Compute increment between columns in degrees
      IF (corner_lons(3) .GT. corner_lons(2)) THEN 
        v5d_proj_args(4) = (corner_lons(3)-corner_lons(2))/(FLOAT(nx)-1.)
      ELSE
        v5d_proj_args(4) = (corner_lons(3)-(corner_lons(2)-360.))/(FLOAT(nx)-1.)
      ENDIF
    ELSE
     PRINT '(2A)', 'V5DINIT: Unsupported map projection:', projection
     STOP
    ENDIF
    
!   CALL getproj(proj_tmp,proj_cent_lon, proj_cent_lat, grid_spacing*0.001, &
!           truelat1,truelat2,nx,ny,v5d_proj,v5d_proj_args, corner_lats(2), &
!             corner_lats(1), corner_lons(2), corner_lons(3))

    ! Set up the vertical levels
    levels(1:kprs) = prslvl(:)*0.01

    ! Set up the varnames and number of levels.  Do the 3d variables first.
    ! These consist of all of the fua variables plus vorticity and 3D
    ! dewpoint, which will be derived by the calling routine.

    varindex = 1
    ! Height (fua/ht)
    varname(varindex) = 'HGT       '
    nlevels(varindex) = kprs
    varindex = varindex + 1
                                          
    ! Temperature (fua/t3)
    varname(varindex) = 'T         '
    nlevels(varindex) = kprs
    varindex = varindex + 1

    ! Dewpoint (not in fua)
    varname(varindex) = 'TD        '
    nlevels(varindex) = kprs
    varindex = varindex + 1  

    ! RH (fua/rh3)
    varname(varindex) = 'RH        '
    nlevels(varindex) = kprs
    varindex = varindex + 1 

    ! U (fua/u3)
    varname(varindex) = 'U         '
    nlevels(varindex) = kprs
    varindex = varindex + 1 

    ! V (fua/v3)
    varname(varindex) = 'V         '
    nlevels(varindex) = kprs
    varindex = varindex + 1  

    ! W (fua/w3)
    varname(varindex) = 'W         '
    nlevels(varindex) = kprs
    varindex = varindex + 1 

    ! Omega  (fua/om)
    varname(varindex) = 'OM        '        
    nlevels(varindex) = kprs
    varindex = varindex + 1  

    ! Absolute Vorticity (Not in LAPS files)
    varname(varindex) = 'AVORT     '
    nlevels(varindex) = kprs
    varindex = varindex + 1

    ! Specific Humidity (fua/sh)
    varname(varindex) = 'SH        '
    nlevels(varindex) = kprs
    varindex = varindex + 1    
                               
    ! Cloud liquid water (fua/lwc)
    varname(varindex) = 'LWC       '
    nlevels(varindex) = kprs
    varindex = varindex + 1    

    ! Cloud Ice (fua/ice)
    varname(varindex) = 'ICE       '
    nlevels(varindex) = kprs
    varindex = varindex + 1    

    ! Rain (fua/rai)
    varname(varindex) = 'RAI       '
    nlevels(varindex) = kprs
    varindex = varindex + 1    

    ! Snow (fua/rai)
    varname(varindex) = 'SNO       '
    nlevels(varindex) = kprs
    varindex = varindex + 1    

    ! Precip Ice/graupel (fua/pic)
    varname(varindex) = 'PIC       '
    nlevels(varindex) = kprs
    varindex = varindex + 1    
 
    ! Total Condensate
    varname(varindex) = 'COND      '
    nlevels(varindex) = kprs
    varindex = varindex + 1

    ! Radar Reflectivity (fua/ref)
    varname(varindex) = 'REF       '
    nlevels(varindex) = kprs
    varindex = varindex + 1    

    ! Precip Type Code (fua/pty)
    varname(varindex) = 'PTY       '
    nlevels(varindex) = kprs
    varindex = varindex + 1

    ! Turbulent kinetic energy (fua/tke)
    varname(varindex) = 'TKE       '
    nlevels(varindex) = kprs
    varindex = varindex + 1

    ! Print out number of 3D variables as a sanity check
    num3dvar = varindex - 1
    PRINT '(A,I3)', 'V5DINIT: Number of 3D variables to output: ', num3dvar

    ! Set up 2D variables for output 

    ! 1000-500 mb thickness (not in fua/fsf)
    varname(varindex) = 'THICK     '
    nlevels(varindex) = 1
    varindex = varindex + 1
    ! Surface Temperature (fsf/t)
    varname(varindex) = 'T_SFC     '
    nlevels(varindex) = 1
    varindex = varindex + 1
   
    ! Surface dewpoint (fsf/td)
    varname(varindex) = 'TD_SFC    '
    nlevels(varindex) = 1
    varindex = varindex + 1

   ! Surface RH (fsf/rh)
    varname(varindex) = 'RH_SFC    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Surface
    ! Surface u (fsf/u)
    varname(varindex) = 'U_SFC     '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Surface v (fsf/v)
    varname(varindex) = 'V_SFC     '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Surface w (fsf/w)
    varname(varindex) = 'W_SFC     '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Mean SLP (fsf/msl)
    varname(varindex) = 'MSLP      '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Surface pressure (fsf/ps)
    varname(varindex) = 'P_SFC     '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Cloud bases (fsf/lcb)
    varname(varindex) = 'CLDBAS    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Cloud tops (fsf/lct)
    varname(varindex) = 'CLDTOP    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Cloud cover (fsf/lcv)
    varname(varindex) = 'CLDCOV    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Cloud ceiling (fsf/cce)
    varname(varindex) = 'CLDCEI    '
    nlevels(varindex) = 1
    varindex = varindex + 1
    
    ! Integrated liquid water depth (fsf/lil)
    varname(varindex) = 'INTLIQ    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Total Precip. Water (fsf/tpw)
    varname(varindex) = 'TPW       '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Incremental precip accum (fsf/r01)
    varname(varindex) = 'PCPINC    '
    nlevels(varindex) = 1
    varindex = varindex + 1
 
   ! Sim. total precip accum (fsf/rto)
    varname(varindex) = 'PCPTOT    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Incremental convective precip accum (not in LAPS output)
    varname(varindex) = 'CONINC    '
    nlevels(varindex) = 1
    varindex = varindex + 1

   ! Sim. total convective precip accum (not in LAPS output)
    varname(varindex) = 'CONTOT    '
    nlevels(varindex) = 1
    varindex = varindex + 1

 
   ! Incremental snow accum (fsf/s01)
    varname(varindex) = 'SNOINC    '
    nlevels(varindex) = 1
    varindex = varindex + 1
 
   ! Sim. total precip accum (fsf/sto)
    varname(varindex) = 'SNOTOT    '
    nlevels(varindex) = 1
    varindex = varindex + 1
    
    ! Sfc precip type code (fsf/spt)
    varname(varindex) = 'SFCPCP    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Sfc potential temp (fsf_th)
    varname(varindex) = 'THETA     '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Sfc equive pot temp (fsf/the)
    varname(varindex) = 'THETAE    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Sfc CAPE (fsf/pbe)
    varname(varindex) = 'CAPE      '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Sfc CIN (fsf/nbe)
    varname(varindex) = 'CIN       '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Sfc lifted index (fsf/li)
    varname(varindex) = 'LI        '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Storm rel. helicity (fsf/lhe)
    varname(varindex) = 'SRHEL     '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Low-level reflectivity (fsf/llr)
    varname(varindex) = 'LLREF'
    nlevels(varindex) = 1
    varindex = varindex + 1
 
    ! Max reflectivity (fsf/lmr)
    varname(varindex) = 'MAXREF    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Echo tops (fsf/lmt)
    varname(varindex) = 'ECHTOP    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Heat index (fsf/hi)
    varname(varindex) = 'HTIDX     '
    nlevels(varindex) = 1
    varindex = varindex + 1
 
    ! Visibility (fsf/vis)
    varname(varindex) = 'VIS       '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Model terrain (setup)
    varname(varindex) = 'TOPO      '
    nlevels(varindex) = 1
    varindex = varindex + 1    

    ! Snow Cover 
    varname(varindex) = 'SNOWCOV   '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Outgoing LW Radiation
    varname(varindex) = 'LWOUT     '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Outgoing SW Radiaion
    varname(varindex) = 'SWOUT     '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Incoming LW Radiation
    varname(varindex) = 'LWDOWN    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Incoming SW Radiation
    varname(varindex) = 'SWDOWN    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Sensible Heat Flux
    varname(varindex) = 'SHFLUX    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Latent Heat Flux
    varname(varindex) = 'LHFLUX    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! PBL Height
    varname(varindex) = 'PBLHGT    '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Mean U wind in PBL
    varname(varindex) = 'U_PBL     '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Mean V wind in PBL
    varname(varindex) = 'V_PBL     '
    nlevels(varindex) = 1
    varindex = varindex + 1

    ! Ground Temp
    varname(varindex) = 'T_GND     '
    nlevels(varindex) = 1
    varindex = varindex + 1
 
    ! Sanity check print of two-d variables
    num2dvar = varindex - 1 - num3dvar
    PRINT '(A,I3)', 'V5DINIT: Number of 2d variables to output: ', num2dvar
    totalvars = num2dvar + num3dvar

    ! Initialize the file with v5dcreate
    ret = v5dcreate(v5dfile,v5d_ntimes, &
                    totalvars, ny,nx, nlevels,&
                    varname, v5dtimes, v5ddates, compression,v5d_proj, &
                    v5d_proj_args, vertical_code, levels)
    IF (ret .eq. 0) THEN 
      PRINT '(A)', '---Error opening Vis5D output file: v5dfile' 
      CALL ABORT
    ENDIF
    PRINT '(A)', 'Vis5D header written.'
    RETURN
  END SUBROUTINE v5dinit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE v5dout(timestep,zprs,tprs,tdprs,rhprs,uprs,vprs,wprs,omprs, &
             abs_vort,shprs,cldliqmr_prs,cldicemr_prs,rainmr_prs, &
                  snowmr_prs,graupelmr_prs,refl_prs,pcptype_prs,tkeprs, &
                  thick_10_5,tsfc,tdsfc,rhsfc,usfc,vsfc,wsfc, &
                  pmsl,psfc,cldbase,cldtop,cldamt,ceiling, &
                  intliqwater,totpcpwater,pcp_inc,pcp_tot, &
                  con_pcp_inc, con_pcp_tot, &
                  snow_inc,snow_tot,pcptype_sfc,thetasfc,thetaesfc, &
                  cape,cin,liftedind,srhel,refl_sfc,max_refl,echo_tops, &
                  heatind,visibility,snowcover,lwout,swout,lwdown,swdown, &
                  shflux,lhflux, &
                  pblhgt,upbl, vpbl, ground_t)     

  ! Subroutine to output the above variables into the already created
  ! vis5d file.  Note that the variables must be output in the same order
  ! (i.e., the same index) as the names were defined in v5dinit.  

    IMPLICIT NONE
    include 'v5df.h'

    ! Declare the arguments
    INTEGER                 :: timestep
    REAL                    :: zprs(nx,ny,kprs)
    REAL                    :: tprs(nx,ny,kprs)
    REAL                    :: tdprs(nx,ny,kprs)
    REAL                    :: rhprs(nx,ny,kprs)
    REAL                    :: uprs(nx,ny,kprs)
    REAL                    :: vprs(nx,ny,kprs)
    REAL                    :: wprs(nx,ny,kprs)
    REAL                    :: omprs(nx,ny,kprs)
    REAL                    :: abs_vort(nx,ny,kprs)
    REAL                    :: shprs(nx,ny,kprs)
    REAL                    :: cldliqmr_prs(nx,ny,kprs)
    REAL                    :: cldicemr_prs(nx,ny,kprs)
    REAL                    :: rainmr_prs(nx,ny,kprs)
    REAL                    :: snowmr_prs(nx,ny,kprs)
    REAL                    :: graupelmr_prs(nx,ny,kprs)
    REAL                    :: refl_prs(nx,ny,kprs)
    REAL                    :: pcptype_prs(nx,ny,kprs)
    REAL                    :: tkeprs(nx,ny,kprs)
    REAL                    :: thick_10_5 (nx,ny)
    REAL                    :: tsfc(nx,ny)
    REAL                    :: tdsfc(nx,ny)
    REAL                    :: rhsfc(nx,ny)
    REAL                    :: usfc(nx,ny)
    REAL                    :: vsfc(nx,ny)
    REAL                    :: wsfc(nx,ny)
    REAL                    :: pmsl(nx,ny)
    REAL                    :: psfc(nx,ny)
    REAL                    :: cldbase(nx,ny)
    REAL                    :: cldtop(nx,ny)
    REAL                    :: cldamt(nx,ny)
    REAL                    :: ceiling(nx,ny)
    REAL                    :: intliqwater(nx,ny)
    REAL                    :: totpcpwater(nx,ny)
    REAL                    :: pcp_inc(nx,ny)
    REAL                    :: pcp_tot(nx,ny)
    REAL                    :: con_pcp_inc(nx,ny)
    REAL                    :: con_pcp_tot(nx,ny)

    REAL                    :: snow_inc(nx,ny)
    REAL                    :: snow_tot(nx,ny)
    REAL                    :: pcptype_sfc(nx,ny)
    REAL                    :: thetasfc(nx,ny)
    REAL                    :: thetaesfc(nx,ny)
    REAL                    :: cape(nx,ny)
    REAL                    :: cin(nx,ny)
    REAL                    :: liftedind(nx,ny)
    REAL                    :: srhel(nx,ny)
    REAL                    :: refl_sfc(nx,ny)
    REAL                    :: max_refl(nx,ny)
    REAL                    :: echo_tops(nx,ny)
    REAL                    :: heatind(nx,ny)
    REAL                    :: visibility(nx,ny)
    REAL                    :: snowcover(nx,ny)
    REAL                    :: lwout(nx,ny)
    REAL                    :: swout(nx,ny)
    REAL                    :: lwdown(nx,ny)
    REAL                    :: swdown(nx,ny) 
    REAL                    :: shflux(nx,ny)
    REAL                    :: lhflux(nx,ny)
    REAL                    :: pblhgt(nx,ny)
    REAL                    :: upbl(nx,ny)
    REAL                    :: vpbl(nx,ny)
    REAL                    :: ground_t(nx,ny)

    ! Local variables
    REAL, ALLOCATABLE       :: data2d(:,:)
    REAL, ALLOCATABLE       :: data3d(:,:,:)
    REAL, ALLOCATABLE       :: condmr(:,:,:)
    INTEGER                 :: varindex
    INTEGER                 :: ret
   
    varindex = 1

    ! Do the 3D variables first
    ALLOCATE(data3d(nx,ny,kprs))

    CALL post2v5d(zprs,nx,ny,kprs,data3d)
    data3d = data3d*0.1 
    PRINT '(2A)', 'V5DOUT: Writing zprs*0.1 as varname: ', varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1

    CALL post2v5d(tprs,nx,ny,kprs,data3d)
    data3d = data3d - 273.15
    PRINT '(2A)', 'V5DOUT: Writing tprs-273.15 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1  

    CALL post2v5d(tdprs,nx,ny,kprs,data3d)
    data3d = data3d - 273.15
    PRINT '(2A)', 'V5DOUT: Writing tdprs-273.15 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1     

    CALL post2v5d(rhprs,nx,ny,kprs,data3d)
    PRINT '(2A)', 'V5DOUT: Writing rhprs as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1   

    CALL post2v5d(uprs,nx,ny,kprs,data3d)
    PRINT '(2A)', 'V5DOUT: Writing uprs as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1  

    CALL post2v5d(vprs,nx,ny,kprs,data3d)
    PRINT '(2A)', 'V5DOUT: Writing vprs as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1  

    CALL post2v5d(wprs,nx,ny,kprs,data3d)
    PRINT '(2A)', 'V5DOUT: Writing wprs as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1   

    CALL post2v5d(omprs,nx,ny,kprs,data3d)
    PRINT '(2A)', 'V5DOUT: Writing omprs as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1  

    CALL post2v5d(abs_vort,nx,ny,kprs,data3d)
    data3d = data3d * 100000
    PRINT '(2A)', 'V5DOUT: Writing abs_vort*10000 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1  

    CALL post2v5d(shprs,nx,ny,kprs,data3d)
    data3d = data3d * 1000.
    PRINT '(2A)', 'V5DOUT: Writing shprs*1000 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1  

    CALL post2v5d(cldliqmr_prs,nx,ny,kprs,data3d)
    data3d = data3d * 1000.
    PRINT '(2A)', 'V5DOUT: Writing cldliqmr_prs*1000 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1  

    CALL post2v5d(cldicemr_prs,nx,ny,kprs,data3d)
    data3d = data3d * 1000.
    PRINT '(2A)', 'V5DOUT: Writing cldicemr_prs*1000 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1 

    CALL post2v5d(rainmr_prs,nx,ny,kprs,data3d)
    data3d = data3d * 1000.
    PRINT '(2A)', 'V5DOUT: Writing rainmr_prs*1000 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1   

    CALL post2v5d(snowmr_prs,nx,ny,kprs,data3d)
    data3d = data3d * 1000.
    PRINT '(2A)', 'V5DOUT: Writing snowmr_prs*1000 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1   

    CALL post2v5d(graupelmr_prs,nx,ny,kprs,data3d)
    data3d = data3d * 1000.
    PRINT '(2A)', 'V5DOUT: Writing graupelmr_prs*1000 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1   

    ! Generate total condensate field
    ALLOCATE(condmr(nx,ny,kprs))
    condmr = cldliqmr_prs + cldicemr_prs + rainmr_prs + snowmr_prs &
             + graupelmr_prs
    CALL post2v5d(condmr,nx,ny,kprs,data3d)
    data3d = data3d * 1000.
    PRINT '(2A)', 'V5DOUT: Writing condmr*1000 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    DEALLOCATE(condmr)
    varindex = varindex + 1

    CALL post2v5d(refl_prs,nx,ny,kprs,data3d)
    WHERE(data3d .LT. 0.) data3d = 0.
    PRINT '(2A)', 'V5DOUT: Writing refl_prs as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1   

    CALL post2v5d(pcptype_prs,nx,ny,kprs,data3d)
    PRINT '(2A)', 'V5DOUT: Writing pcptype_prs as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1  

    CALL post2v5d(tkeprs,nx,ny,kprs,data3d)
    PRINT '(2A)', 'V5DOUT: Writing tkeprs as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data3d)
    varindex = varindex + 1

    DEALLOCATE(data3d)
    ! Now do the 2-d variables 
    ALLOCATE(data2d(nx,ny))

    CALL post2v5d(thick_10_5,nx,ny,1,data2d)
    data2d = data2d * 0.1
    PRINT '(2A)', 'V5DOUT: Writing thick_10_5*0.1 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1 

    CALL post2v5d(tsfc,nx,ny,1,data2d)
    data2d = data2d - 273.15
    PRINT '(2A)', 'V5DOUT: Writing tsfc-273.15 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1    

    CALL post2v5d(tdsfc,nx,ny,1,data2d)
    data2d = data2d - 273.15
    PRINT '(2A)', 'V5DOUT: Writing tdsfc-273.15 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1     

    CALL post2v5d(rhsfc,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing rhsfc as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1     

    CALL post2v5d(usfc,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing usfc as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1     

    CALL post2v5d(vsfc,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing vsfc as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(wsfc,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing wsfc as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(pmsl,nx,ny,1,data2d)
    data2d = data2d * 0.01
    PRINT '(2A)', 'V5DOUT: Writing pmsl*0.01 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(psfc,nx,ny,1,data2d)
    data2d = data2d * 0.01
    PRINT '(2A)', 'V5DOUT: Writing psfc*0.01 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1     

    CALL post2v5d(cldbase,nx,ny,1,data2d)
    data2d = data2d * 0.1
    PRINT '(2A)', 'V5DOUT: Writing cldbase*0.1 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(cldtop,nx,ny,1,data2d)
    data2d = data2d * 0.1
    PRINT '(2A)', 'V5DOUT: Writing cldtop*0.1 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(cldamt,nx,ny,1,data2d)
    data2d = data2d*100.
    PRINT '(2A)', 'V5DOUT: Writing cldamt*100 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1    

    CALL post2v5d(ceiling,nx,ny,1,data2d)
    data2d = data2d * 0.1 
    PRINT '(2A)', 'V5DOUT: Writing ceiling*0.1 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(intliqwater,nx,ny,1,data2d)
    data2d = data2d
    PRINT '(2A)', 'V5DOUT: Writing intliqwater (mm) as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1 

    CALL post2v5d(totpcpwater,nx,ny,1,data2d)
    data2d = data2d
    PRINT '(2A)', 'V5DOUT: Writing totpcpwater (mm) as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(pcp_inc,nx,ny,1,data2d)
    data2d = data2d * 39.37
    PRINT '(2A)', 'V5DOUT: Writing pcp_inc*39.37 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(pcp_tot,nx,ny,1,data2d)
    data2d = data2d * 39.37 
    PRINT '(2A)', 'V5DOUT: Writing pcp_tot*39.37 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(con_pcp_inc,nx,ny,1,data2d)
    data2d = data2d * 39.37
    PRINT '(2A)', 'V5DOUT: Writing con_pcp_inc*39.37 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(con_pcp_tot,nx,ny,1,data2d)
    data2d = data2d * 39.37
    PRINT '(2A)', 'V5DOUT: Writing con_pcp_tot*39.37 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(snow_inc,nx,ny,1,data2d)
    data2d = data2d * 39.37 
    PRINT '(2A)', 'V5DOUT: Writing snow_inc*39.37 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(snow_tot,nx,ny,1,data2d)
    data2d = data2d * 39.37 
    PRINT '(2A)', 'V5DOUT: Writing snow_tot*39.37 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1   

    CALL post2v5d(pcptype_sfc,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing pcptype_sfc as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(thetasfc,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing thetasfc as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(thetaesfc,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing thetaesfc as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1  

    CALL post2v5d(cape,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing cape as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1 

    CALL post2v5d(cin,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing cin as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1    

    CALL post2v5d(liftedind,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing liftedind as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1   

    CALL post2v5d(srhel,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing srhel as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1    

    CALL post2v5d(refl_sfc,nx,ny,1,data2d)
    WHERE(refl_sfc .lt. 0) refl_sfc = 0.
    PRINT '(2A)', 'V5DOUT: Writing refl_sfc as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1    

    CALL post2v5d(max_refl,nx,ny,1,data2d)
    WHERE(data2d .lt. 0) data2d = 0.
    PRINT '(2A)', 'V5DOUT: Writing max_refl as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1       

    CALL post2v5d(echo_tops,nx,ny,1,data2d)
    data2d = data2d*0.1
    PRINT '(2A)', 'V5DOUT: Writing echo_tops*0.1 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1       

    CALL post2v5d(heatind,nx,ny,1,data2d)
    data2d = data2d - 273.15
    PRINT '(2A)', 'V5DOUT: Writing heatind - 273.15 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(visibility,nx,ny,1,data2d)
    data2d = data2d *0.001
    PRINT '(2A)', 'V5DOUT: Writing visibility*0.001 as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1   

    CALL post2v5d(terdot,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing terdot as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1    

    CALL post2v5d(snowcover,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing snowcover as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(lwout,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing lwout as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(swout,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing swout as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

   CALL post2v5d(lwdown,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing lwdown as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(swdown,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing swsdown as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(shflux,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing shflux as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(lhflux,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing lhflux as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(pblhgt,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing pblhgt as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(upbl,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing upbl as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(vpbl,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing vpbl as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    CALL post2v5d(ground_t,nx,ny,1,data2d)
    PRINT '(2A)', 'V5DOUT: Writing ground_t as varname: ',&
              varname(varindex)
    ret = v5dwrite(timestep,varindex,data2d)
    varindex = varindex + 1

    DEALLOCATE(data2d)
    RETURN 
  END SUBROUTINE v5dout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE v5dend

    IMPLICIT NONE
    include 'v5df.h'
    CHARACTER(LEN=255) donefile
    INTEGER :: flagunit
    ret = v5dclose()
    IF (ret .eq. 0) THEN
      PRINT '(A)', 'Problem closing the Vis5D data file!'
    ENDIF
    IF (realtime) THEN
       donefile = TRIM(v5dfile) // '.done'
       CALL get_file_unit(flagunit)
       OPEN(UNIT=flagunit,FILE=donefile,STATUS='unknown')
       CLOSE (flagunit)
    ENDIF
      
  END SUBROUTINE v5dend
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE POST2V5D(POSTDATA, IMAX, JMAX, KMAX, V5DDATA)

!-----------------------------------------------------------------------
!  TRANSFORM DATA IN POST SPACE (I INCREASING EAST, J INCREASING NORTH, 
!    K INCREASING UPWARD)  INTO VIS5D SPACE (I INCREASING SOUTHWARD, 
!    J INCREASING EAST, K INCREASING UPWARD) 
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER                      :: I
      INTEGER                      :: IMAX
      INTEGER                      :: J
      INTEGER                      :: JMAX
      INTEGER                      :: K
      INTEGER                      :: JMAXP1
      INTEGER                      :: KMAX
      REAL                         :: POSTDATA   (IMAX,JMAX,KMAX)
      REAL                         :: V5DDATA    (JMAX,IMAX,KMAX)

!-----------------------------------------------------------------------
!  
!-----------------------------------------------------------------------

      JMAXP1 = JMAX + 1

      DO K = 1,KMAX
        DO J = 1,JMAX
          DO I = 1,IMAX

            V5DDATA(J,I,K) = POSTDATA(I, JMAXP1-J,K)

          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE POST2V5D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE vis5d
