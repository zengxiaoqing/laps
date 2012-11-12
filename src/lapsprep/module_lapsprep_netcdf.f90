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


MODULE lapsprep_netcdf

! PURPOSE
! =======
! Module to contain the various output routines needed for lapsprep.
!
! SUBROUTINES CONTAINED
! =====================
! output_netcdf_format  - Used to support RAMS 4.x initializations
!
! REMARKS
! =======

!
! HISTORY
! =======
! 28 Nov 2000 -- Original -- Brent Shaw
! 14 Nov 2001 -- Modified -- John Snook

  USE setup
  USE laps_static
  USE date_pack

  PRIVATE
  PUBLIC output_netcdf_format
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE output_netcdf_format(pr,ht,tp,mr,uw,vw,ww,slp,spr,lwc,ice,rai,sno,pic,rh,snocov,tskin)

  ! Subroutine to output data in netcdf format.  The skin temp from
  ! LAPS (lsx/tgd) is used for SST.
 
  ! Note that the P array has a bogus 2001 mb value as the last entry
  ! to designate the surface for MM5 applications.  The surface values
  ! of the state variables are contained in the last layer of their
  ! respective 3D arrays. 

  ! The u/v winds are grid-relative, and are rotated to make them true.

  ! P levels must be written as integer values

  ! RH must be written as a fraction 

  ! 3d Variables are stored top down, and are flipped to bottom up.

  IMPLICIT NONE

  include 'netcdf.inc'

  REAL                   :: pr(z3+1)      !Pressure levels (mb)
  REAL                   :: uw(x,y,z3+1)  !U-component of wind wrt grid (m/s)
  REAL                   :: vw(x,y,z3+1)  !V-component of wind wrt grid (m/s) 
  REAL                   :: ww(x,y,z3+1)  !W-component of wind (m/s) 
  REAL                   :: tp(x,y,z3+1)  !Temperature (K)
  REAL                   :: ht(x,y,z3+1)  !Geopotential Height (m)
  REAL                   :: mr(x,y,z3+1)  !Mixing ratio (kg/kg)
  REAL                   :: slp(x,y)      !MSL Pressure (Pa)
  REAL                   :: spr(x,y)      !Surface Pressure (Pa)
  REAL                   :: lwc(x,y,z3)   !Cloud water mr (kg/kg)
  REAL                   :: ice(x,y,z3)   !Cloud ice mr (kg/kg)
  REAL                   :: rai(x,y,z3)   !Precip rain mr (kg/kg)
  REAL                   :: sno(x,y,z3)   !Precip snow mr (kg/kg)
  REAL                   :: pic(x,y,z3)   !Precip ice mr (kg/kg)
  REAL                   :: rh(x,y,z3+1)  !Relative Humidity (%)
  REAL                   :: snocov(x,y)   !Snow cover (fract)
  REAL                   :: tskin(x,y)    !Skin Temp (K)

  ! Local Variables

  real, allocatable              :: pr3(:)
  real, allocatable              :: ht3(:,:,:)
  real, allocatable              :: tp3(:,:,:)
  real, allocatable              :: mr3(:,:,:)
  real, allocatable              :: rh3(:,:,:)
  real, allocatable              :: uw3(:,:,:)
  real, allocatable              :: vw3(:,:,:)
  real, allocatable              :: ww3(:,:,:)
  real, allocatable              :: sht(:,:)
  real, allocatable              :: stp(:,:)
  real, allocatable              :: srh(:,:)
  real, allocatable              :: smr(:,:)
  real, allocatable              :: suw(:,:)
  real, allocatable              :: svw(:,:)
  real, allocatable              :: lwcf(:,:,:)
  real, allocatable              :: icef(:,:,:)
  real, allocatable              :: raif(:,:,:)
  real, allocatable              :: snof(:,:,:)
  real, allocatable              :: picf(:,:,:)
  real, allocatable              :: scvf(:,:)
  real, allocatable              :: tskf(:,:)

  INTEGER                        :: yyyyddd, valid_mm, valid_dd
  INTEGER                        :: icode,ncid,nid,idimid(3),start,count
  INTEGER                        :: i,j,k
  INTEGER*2                      :: short
  REAL, ALLOCATABLE              :: ut(:,:,:)
  REAL, ALLOCATABLE              :: vt(:,:,:)
  REAL                           :: lat0
  REAL*8                         :: reftime
  CHARACTER (LEN=256)            :: output_file_name
  CHARACTER (LEN=128)            :: agridtype
  CHARACTER (LEN=15)             :: date_string

  ! Separate surface data from upper air data.

  allocate(pr3(z3))
  pr3 = pr(1:z3)

  allocate(ht3(x,y,z3))
  ht3 = ht(1:x,1:y,1:z3)

  allocate(tp3(x,y,z3))
  tp3 = tp(1:x,1:y,1:z3)

  allocate(mr3(x,y,z3))
  mr3 = mr(1:x,1:y,1:z3)

  allocate(rh3(x,y,z3))
  rh3 = rh(1:x,1:y,1:z3)

  allocate(uw3(x,y,z3))
  uw3 = uw(1:x,1:y,1:z3)

  allocate(vw3(x,y,z3))
  vw3 = vw(1:x,1:y,1:z3)

  allocate(ww3(x,y,z3))
  ww3 = ww(1:x,1:y,1:z3)

  allocate(sht(x,y))
  sht = ht(1:x,1:y,z3+1)

  allocate(stp(x,y))
  stp = tp(1:x,1:y,z3+1)

  allocate(smr(x,y))
  smr = mr(1:x,1:y,z3+1)

  allocate(srh(x,y))
  srh = rh(1:x,1:y,z3+1)

  allocate(suw(x,y))
  suw = uw(1:x,1:y,z3+1)

  allocate(svw(x,y))
  svw = vw(1:x,1:y,z3+1)

  allocate(scvf(x,y))
  scvf = snocov(1:x,1:y)

  allocate(tskf(x,y))
  tskf = tskin(1:x,1:y)

  ! Flip 3d arrays.

  call flip_array(1,1,z3,pr3)
  call flip_array(x,y,z3,ht3)
  call flip_array(x,y,z3,tp3)
  call flip_array(x,y,z3,mr3)
  call flip_array(x,y,z3,rh3)
  call flip_array(x,y,z3,uw3)
  call flip_array(x,y,z3,vw3)

  if (hotstart) then
    allocate(lwcf(x,y,z3))
    lwcf = lwc
    call flip_array(x,y,z3,lwcf)

    allocate(icef(x,y,z3))
    icef = ice
    call flip_array(x,y,z3,icef)

    allocate(raif(x,y,z3))
    raif = rai
    call flip_array(x,y,z3,raif)

    allocate(snof(x,y,z3))
    snof = sno
    call flip_array(x,y,z3,snof)

    allocate(picf(x,y,z3))
    picf = pic
    call flip_array(x,y,z3,picf)
    call flip_array(x,y,z3,ww3)
  endif

  ! Build the output file name
 
  output_prefix = TRIM(laps_data_root)// '/lapsprd/lapsprep/cdf/LAPS'
  yyyyddd = valid_yyyy*1000 + valid_jjj
  CALL wrf_date_to_ymd(yyyyddd, valid_yyyy, valid_mm, valid_dd) 
  WRITE(date_string,'(I4.4,"-",I2.2,"-",I2.2,"-",I2.2,I2.2)') &
          valid_yyyy, valid_mm, valid_dd, valid_hh, valid_min
  output_file_name = TRIM(output_prefix) // ':' // date_string

  !  Create LAPS model netcdf file.

  icode=nf_create(trim(output_file_name),nf_clobber,ncid)
  if (icode .ne. 0) then
     print *,'Could not open output file: ',trim(output_file_name)
     stop
  endif

  !  Define netcdf grid dimensions.
 
  icode=nf_def_dim(ncid,'nx',x,1)
  icode=nf_def_dim(ncid,'ny',y,2)
  icode=nf_def_dim(ncid,'nz',z3,3)
  icode=nf_def_dim(ncid,'len',128,4)
  icode=nf_def_dim(ncid,'pt',1,5)

  !  Define static and grid variables.

  icode=nf_def_var(ncid,'grid_type',nf_char,1,4,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',15,'grid projection')
  icode=nf_def_var(ncid,'Nx',nf_short,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',18,'number of x points')
  icode=nf_def_var(ncid,'Ny',nf_short,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',18,'number of y points')
  icode=nf_def_var(ncid,'Np',nf_short,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',25,'number of pressure levels')
  icode=nf_def_var(ncid,'Dx',nf_real,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',16,'x grid increment')
  icode=nf_put_att_text(ncid,nid,'units',10,'kilometers')
  icode=nf_def_var(ncid,'Dy',nf_real,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',16,'y grid increment')
  icode=nf_put_att_text(ncid,nid,'units',10,'kilometers')
  icode=nf_def_var(ncid,'Lat0',nf_real,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',13,'pole latitude')
  icode=nf_put_att_text(ncid,nid,'units',13,'degrees north')
  icode=nf_def_var(ncid,'Lat1',nf_real,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',14,'first latitude')
  icode=nf_put_att_text(ncid,nid,'units',13,'degrees north')
  icode=nf_def_var(ncid,'Lat2',nf_real,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',15,'second latitude')
  icode=nf_put_att_text(ncid,nid,'units',13,'degrees north')
  icode=nf_def_var(ncid,'Lon0',nf_real,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',14,'pole longitude')
  icode=nf_put_att_text(ncid,nid,'units',12,'degrees east')
  icode=nf_def_var(ncid,'SwLat',nf_real,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',18,'Southwest latitude')
  icode=nf_put_att_text(ncid,nid,'units',13,'degrees north')
  icode=nf_def_var(ncid,'SwLon',nf_real,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',19,'Southwest longitude')
  icode=nf_put_att_text(ncid,nid,'units',12,'degrees east')
  icode=nf_def_var(ncid,'NeLat',nf_real,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',18,'Northeast latitude')
  icode=nf_put_att_text(ncid,nid,'units',13,'degrees north')
  icode=nf_def_var(ncid,'NeLon',nf_real,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',19,'Northeast longitude')
  icode=nf_put_att_text(ncid,nid,'units',12,'degrees east')
  icode=nf_def_var(ncid,'reftime',nf_double,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',14,'reference time')
  icode=nf_put_att_text(ncid,nid,'units',35,'seconds since (1970-1-1 00:00:00.0)')
  icode=nf_def_var(ncid,'valtime',nf_double,1,5,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',10,'valid time')
  icode=nf_put_att_text(ncid,nid,'units',35,'seconds since (1970-1-1 00:00:00.0)')
  icode=nf_def_var(ncid,'level',nf_real,1,3,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',13,'level of data')
  icode=nf_put_att_text(ncid,nid,'units',2,'mb')
  icode=nf_def_var(ncid,'pr',nf_real,1,3,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',8,'pressure')
  icode=nf_put_att_text(ncid,nid,'units',7,'Pascals')

  !  Define surface variables.

  idimid(1)=1
  idimid(2)=2
  icode=nf_def_var(ncid,'sht',nf_real,2,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',14,'surface height')
  icode=nf_put_att_text(ncid,nid,'units',6,'meters')
  icode=nf_def_var(ncid,'spr',nf_real,2,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',16,'surface pressure')
  icode=nf_put_att_text(ncid,nid,'units',7,'pascals')
  icode=nf_def_var(ncid,'slp',nf_real,2,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',23,'mean sea-level pressure')
  icode=nf_put_att_text(ncid,nid,'units',7,'pascals')
  icode=nf_def_var(ncid,'stp',nf_real,2,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',19,'surface temperature')
  icode=nf_put_att_text(ncid,nid,'units',6,'kelvin')
  icode=nf_def_var(ncid,'smr',nf_real,2,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',20,'surface mixing ratio')
  icode=nf_put_att_text(ncid,nid,'units',5,'kg/kg')
  icode=nf_def_var(ncid,'srh',nf_real,2,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',25,'surface relative humidity')
  icode=nf_put_att_text(ncid,nid,'units',7,'percent')
  icode=nf_def_var(ncid,'suw',nf_real,2,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',14,'surface u-wind')
  icode=nf_put_att_text(ncid,nid,'units',14,'meters/seconds')
  icode=nf_def_var(ncid,'svw',nf_real,2,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',14,'surface v-wind')
  icode=nf_put_att_text(ncid,nid,'units',14,'meters/seconds')
  icode=nf_def_var(ncid,'scv',nf_real,2,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',10,'snow cover')
  icode=nf_put_att_text(ncid,nid,'units',8,'fraction')
  icode=nf_def_var(ncid,'tsk',nf_real,2,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',16,'skin temperature')
  icode=nf_put_att_text(ncid,nid,'units',6,'kelvin')

  !  Define upper-air variables.

  idimid(3)=3
  icode=nf_def_var(ncid,'ht',nf_real,3,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',6,'height')
  icode=nf_put_att_text(ncid,nid,'units',6,'meters')
  icode=nf_def_var(ncid,'tp',nf_real,3,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',11,'temperature')
  icode=nf_put_att_text(ncid,nid,'units',6,'kelvin')
  icode=nf_def_var(ncid,'mr',nf_real,3,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',12,'mixing ratio')
  icode=nf_put_att_text(ncid,nid,'units',5,'kg/kg')
  icode=nf_def_var(ncid,'rh',nf_real,3,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',17,'relative humidity')
  icode=nf_put_att_text(ncid,nid,'units',7,'percent')
  icode=nf_def_var(ncid,'uw',nf_real,3,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',19,'u-component of wind')
  icode=nf_put_att_text(ncid,nid,'units',13,'meters/second')
  icode=nf_def_var(ncid,'vw',nf_real,3,idimid,nid)
  icode=nf_put_att_text(ncid,nid,'long_name',19,'v-component of wind')
  icode=nf_put_att_text(ncid,nid,'units',13,'meters/second')

  if (hotstart) then
     icode=nf_def_var(ncid,'lwc',nf_real,3,idimid,nid)
     icode=nf_put_att_text(ncid,nid,'long_name',31  &
          ,'cloud liquid water mixing ratio')
     icode=nf_put_att_text(ncid,nid,'units',5,'kg/kg')
     icode=nf_def_var(ncid,'ice',nf_real,3,idimid,nid)
     icode=nf_put_att_text(ncid,nid,'long_name',22  &
          ,'cloud ice mixing ratio')
     icode=nf_put_att_text(ncid,nid,'units',5,'kg/kg')
     icode=nf_def_var(ncid,'sno',nf_real,3,idimid,nid)
     icode=nf_put_att_text(ncid,nid,'long_name',31  &
          ,'precipitating snow mixing ratio')
     icode=nf_put_att_text(ncid,nid,'units',5,'kg/kg')
     icode=nf_def_var(ncid,'rai',nf_real,3,idimid,nid)
     icode=nf_put_att_text(ncid,nid,'long_name',31  &
          ,'precipitating rain mixing ratio')
     icode=nf_put_att_text(ncid,nid,'units',5,'kg/kg')
     icode=nf_def_var(ncid,'pic',nf_real,3,idimid,nid)
     icode=nf_put_att_text(ncid,nid,'long_name',30  &
          ,'precipitating ice mixing ratio')
     icode=nf_put_att_text(ncid,nid,'units',5,'kg/kg')
     icode=nf_def_var(ncid,'ww',nf_real,3,idimid,nid)
     icode=nf_put_att_text(ncid,nid,'long_name',19  &
          ,'w-component of wind')
     icode=nf_put_att_text(ncid,nid,'units',13,'meters/second')
  endif

  icode=nf_enddef(ncid)

  !  Create reftime.

  call adate_to_i4time(laps_file_time,reftime)

  !  Fill static and grid variables.

  IF      ( grid_type(1:8)  .EQ. 'latlon'                   ) THEN
    agridtype='Lat-Lon'
  ELSE IF ((grid_type(1:24) .EQ. 'secant lambert conformal').OR. &
           (grid_type(1:28) .EQ. 'tangential lambert conformal')) THEN
    agridtype='Lambert-Conformal'
  ELSE IF ( grid_type(1:19) .EQ. 'polar stereographic'      ) THEN
    agridtype='Polar-Stereographic'
  ELSE
    PRINT '(A,A,A)','netcdf unsupported map projection: ', &
        TRIM(grid_type),'.  I quit.'
    STOP 'unsupported projection'
  END IF

  start=1
  count=len_trim(agridtype)
  icode=nf_inq_varid(ncid,'grid_type',nid)
  icode=nf_put_vara_text(ncid,nid,start,count,agridtype)
  short=x
  icode=nf_inq_varid(ncid,'Nx',nid)
  icode=nf_put_var_int2(ncid,nid,short)
  short=y
  icode=nf_inq_varid(ncid,'Ny',nid)
  icode=nf_put_var_int2(ncid,nid,short)
  short=z3
  icode=nf_inq_varid(ncid,'Np',nid)
  icode=nf_put_var_int2(ncid,nid,short)
  icode=nf_inq_varid(ncid,'Dx',nid)
  icode=nf_put_var_real(ncid,nid,dx*1000.)
  icode=nf_inq_varid(ncid,'Dy',nid)
  icode=nf_put_var_real(ncid,nid,dy*1000.)
  lat0=90.
  icode=nf_inq_varid(ncid,'Lat0',nid)
  icode=nf_put_var_real(ncid,nid,lat0)
  icode=nf_inq_varid(ncid,'Lat1',nid)
  icode=nf_put_var_real(ncid,nid,latin1)
  icode=nf_inq_varid(ncid,'Lat2',nid)
  icode=nf_put_var_real(ncid,nid,latin2)
  icode=nf_inq_varid(ncid,'Lon0',nid)
  icode=nf_put_var_real(ncid,nid,lov)
  icode=nf_inq_varid(ncid,'SwLat',nid)
  icode=nf_put_var_real(ncid,nid,la1)
  icode=nf_inq_varid(ncid,'SwLon',nid)
  icode=nf_put_var_real(ncid,nid,lo1)
  icode=nf_inq_varid(ncid,'NeLat',nid)
  icode=nf_put_var_real(ncid,nid,la2)
  icode=nf_inq_varid(ncid,'NeLon',nid)
  icode=nf_put_var_real(ncid,nid,lo2)
  icode=nf_inq_varid(ncid,'reftime',nid)
  icode=nf_put_var_double(ncid,nid,reftime)
  icode=nf_inq_varid(ncid,'valtime',nid)
  icode=nf_put_var_double(ncid,nid,reftime)
  icode=nf_inq_varid(ncid,'level',nid)
  icode=nf_put_var_real(ncid,nid,pr3)
  pr3=pr3*100.
  icode=nf_inq_varid(ncid,'pr',nid)
  icode=nf_put_var_real(ncid,nid,pr3)

  ! Rotate the surface u and v winds to true if polar grid projection.

  ALLOCATE (ut(x,y,1)) ! Array for true u-winds
  ALLOCATE (vt(x,y,1)) ! Array for true v-winds

  IF ((agridtype(1:3) .EQ. 'Pol').OR.(agridtype(1:3) .EQ. 'Lam')) THEN
    DO j=1,y
    DO i=1,x
      CALL uvgrid_to_uvtrue(suw(i,j),svw(i,j), &
                            ut(i,j,1),vt(i,j,1),&
                            lons(i,j))
    ENDDO 
    ENDDO
  ELSE
    ut(:,:,1) = suw
    vt(:,:,1) = svw
  ENDIF

  !  Fill surface variables.
 
  icode=nf_inq_varid(ncid,'sht',nid)
  icode=nf_put_var_real(ncid,nid,sht)
  icode=nf_inq_varid(ncid,'spr',nid)
  icode=nf_put_var_real(ncid,nid,spr)
  icode=nf_inq_varid(ncid,'slp',nid)
  icode=nf_put_var_real(ncid,nid,slp)
  icode=nf_inq_varid(ncid,'stp',nid)
  icode=nf_put_var_real(ncid,nid,stp)
  icode=nf_inq_varid(ncid,'smr',nid)
  icode=nf_put_var_real(ncid,nid,smr)
  icode=nf_inq_varid(ncid,'srh',nid)
  icode=nf_put_var_real(ncid,nid,srh)
  icode=nf_inq_varid(ncid,'suw',nid)
  icode=nf_put_var_real(ncid,nid,ut)
  icode=nf_inq_varid(ncid,'svw',nid)
  icode=nf_put_var_real(ncid,nid,vt)
  icode=nf_inq_varid(ncid,'scv',nid)
  icode=nf_put_var_real(ncid,nid,scvf)
  icode=nf_inq_varid(ncid,'tsk',nid)
  icode=nf_put_var_real(ncid,nid,tskf)

  DEALLOCATE (ut)
  DEALLOCATE (vt)

  ! Rotate the u and v winds to true if polar grid projection.

  ALLOCATE (ut(x,y,z3)) ! Array for true u-winds
  ALLOCATE (vt(x,y,z3)) ! Array for true v-winds

  IF ((agridtype(1:3) .EQ. 'Pol').OR.(agridtype(1:3) .EQ. 'Lam')) THEN
    level_loop: DO k = 1,z3
    rotate_winds_j: DO j=1,y
    rotate_winds_i: DO i=1,x
      CALL uvgrid_to_uvtrue(uw3(i,j,k),vw3(i,j,k), &
                            ut(i,j,k),vt(i,j,k),&
                            lons(i,j))
    ENDDO rotate_winds_i
    ENDDO rotate_winds_j
    ENDDO level_loop
  ELSE
    ut = uw3
    vt = vw3
  ENDIF

  !  Fill upper air variables.

  icode=nf_inq_varid(ncid,'ht',nid)
  icode=nf_put_var_real(ncid,nid,ht3)
  icode=nf_inq_varid(ncid,'tp',nid)
  icode=nf_put_var_real(ncid,nid,tp3)
  icode=nf_inq_varid(ncid,'mr',nid)
  icode=nf_put_var_real(ncid,nid,mr3)
  icode=nf_inq_varid(ncid,'rh',nid)
  icode=nf_put_var_real(ncid,nid,rh3)
  icode=nf_inq_varid(ncid,'uw',nid)
  icode=nf_put_var_real(ncid,nid,ut)
  icode=nf_inq_varid(ncid,'vw',nid)
  icode=nf_put_var_real(ncid,nid,vt)

  if (hotstart) then
    icode=nf_inq_varid(ncid,'lwc',nid)
    icode=nf_put_var_real(ncid,nid,lwcf)
    icode=nf_inq_varid(ncid,'ice',nid)
    icode=nf_put_var_real(ncid,nid,icef)
    icode=nf_inq_varid(ncid,'rai',nid)
    icode=nf_put_var_real(ncid,nid,raif)
    icode=nf_inq_varid(ncid,'sno',nid)
    icode=nf_put_var_real(ncid,nid,snof)
    icode=nf_inq_varid(ncid,'pic',nid)
    icode=nf_put_var_real(ncid,nid,picf)
    icode=nf_inq_varid(ncid,'ww',nid)
    icode=nf_put_var_real(ncid,nid,ww3)
    DEALLOCATE(lwcf)
    deallocate(icef)
    deallocate(raif)
    deallocate(snof)
    deallocate(picf)
  endif
     
  DEALLOCATE (ut)
  DEALLOCATE (vt)
  deallocate (pr3)
  deallocate (ht3)
  deallocate (tp3)
  deallocate (mr3)
  deallocate (rh3)
  deallocate (uw3)
  deallocate (vw3)
  deallocate (ww3)
  deallocate (sht)
  deallocate (stp)
  deallocate (smr)
  deallocate (srh)
  deallocate (suw)
  deallocate (svw)
  deallocate (scvf)
  deallocate (tskf)

  ! Close the netcdf file.

  icode=nf_close(ncid)
  
  RETURN
  END SUBROUTINE output_netcdf_format

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE flip_array(nx,ny,nz,a)

  implicit none

  integer nx,ny,nz,i,j,k,kk

  real a(nx,ny,nz),tmp(nx,ny,nz/2)

  kk=nz+1
  do k=1,nz/2
     kk=kk-1
     do j=1,ny
     do i=1,nx
        tmp(i,j,k)=a(i,j,k)
        a(i,j,k)=a(i,j,kk)
     enddo
     enddo
  enddo

  kk=nz/2+1
  do k=(nz+1)/2+1,nz
     kk=kk-1
     do j=1,ny
     do i=1,nx
        a(i,j,k)=tmp(i,j,kk)
     enddo
     enddo
  enddo

  RETURN
  END SUBROUTINE flip_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE adate_to_i4time(adate,i4time)

  implicit none

  real*8 i4time

  integer*4 iyear,iday,ihour,imin,lp

  character*9 adate

  read(adate(1:2),'(i2)') iyear
  read(adate(3:5),'(i3)') iday
  read(adate(6:7),'(i2)') ihour
  read(adate(8:9),'(i2)') imin

  !  Valid for years 1960-2060.

  if (iyear .lt. 60) iyear = iyear + 100

  lp = (iyear + 3 - 60) / 4

  i4time = (iyear-60)  * 31536000  &
         + (iday-1+lp) * 86400     &
         + ihour       * 3600      &
         + imin        * 60

  RETURN
  END SUBROUTINE adate_to_i4time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
END MODULE lapsprep_netcdf
