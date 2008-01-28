!dis    Forecast Systems Laboratory
!dis    NOAA/OAR/ERL/FSL
!dis    325 Broadway
!dis    Boulder, CO     80303
!dis
!dis    Forecast Research Division
!dis    Local Analysis and Prediction Branch
!dis    LAPS
!dis
!dis    This software and its documentation are in the public domain and
!dis    are furnished "as is."  The United States government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  They assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis    technical support to users.
!dis
!dis    Permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  All modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  If significant modifications or enhancements
!dis    are made to this software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

SUBROUTINE GSI_BkObs

!==========================================================
!  This routine converts LAPS background and observations
!  to GSI format files.
!
!	Background  --> wrf_inout;
!	Observation --> bufr files, e.g., prepqc.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	3-2006
!	Modified: YUANFU XIE	10-2007 adding height to wrf.
!==========================================================

  USE LAPS_Parm

  IMPLICIT NONE

  ! Generate wrf_inout file: background.
  CALL GSI_Bkg(n(1),n(2),n(3),lat,lon, &
	      dxy,u_wind3d,v_wind3d)

  ! Generate bufr files for observation: prepqc.laps
  CALL GSI_Obs(i4time,asctime,lat,lon,n(1),n(2),n(3), &
	      nobs_point,obs_point,n_tobs,obs_temp,maxtobs)
  !CALL GSI_Obs_test
  

END SUBROUTINE GSI_BkObs



SUBROUTINE GSI_Bkg(imax,jmax,kmax,xlat,xlong, &
                  grid_spacing,u_laps_bkg,v_laps_bkg)

!==========================================================
!  This routine converts LAPS background into GSI wrf_inout
!  on the GSI mass coordinate.
!
!  HISTORY: MAR. 2006 by YUANFU XIE.
!==========================================================

  USE LAPS_Parm

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: imax,jmax,kmax  	! Dimensions
  REAL, INTENT(IN) :: xlat(imax,jmax),xlong(imax,jmax) 
  REAL, INTENT(IN) :: grid_spacing
  REAL, INTENT(IN) :: u_laps_bkg(imax,jmax,kmax)   ! u bkg
  REAL, INTENT(IN) :: v_laps_bkg(imax,jmax,kmax)   ! v bkg

  ! Local variables:
  REAL, PARAMETER :: cp=1004.0, rc=287.0,t0=300.0 !273.15
  CHARACTER varname*3, fnm*9, hr*2, mins*2, jday*5, filename*150
  INTEGER :: i4time_sys,namelen
  INTEGER :: istatus,i,j,k
  REAL :: t_mass_bkg(imax,jmax,kmax)
  REAL :: sh_mass_bkg(imax,jmax,kmax)
  REAL :: geo_mass_bkg(imax,jmax,kmax)	! Geopotential
  REAL :: u_mass_bkg(imax,jmax,kmax),v_mass_bkg(imax,jmax,kmax)
  REAL :: dam(imax,jmax),pdam(imax,jmax)
  REAL :: znw(kmax),znu(kmax-1),mapfac_m(imax,jmax)

  ! Times:
  CHARACTER*19 :: times

  ! System time:
  CALL get_systime_all(i4time_sys,fnm,hr,mins,asctime,jday,istatus)
  IF (i4time .NE. i4time_sys) then
    PRINT*,'gsibkg: error: reading background in wrong time'
    STOP
  ENDIF

  ! Dry air mass in column (base state):
  CALL dryairmass(dam,pdam,imax,jmax,kmax,pressr1d, &
                  height3d,pressr3d,temptr3d)

  ! As default, a uniform mass vertical grid is used:
  ! eta = (p_d(k)-p_d(top))/(p_d(sfc)-p_d(top)):
  ! DO k=1,kmax
  !  znw(k) = 1.0-float(k-1)/float(kmax-1)
  ! ENDDO
  ! DO k=1,kmax-1
  !   znu(k) = 1.0-(float(k-1)+0.5)/float(kmax-1)
  ! ENDDO
  DO k=1,kmax
    znw(k) = (pressr1d(k)-pressr1d(kmax))/(pressr1d(1)-pressr1d(kmax))
  ENDDO
  DO k=1,kmax-1
    znu(k) = (0.5*(pressr1d(k)+pressr1d(k+1))-pressr1d(kmax))/(pressr1d(1)-pressr1d(kmax))
  ENDDO

  ! Map factor: use LAPS routine: get_sigma(lat,lon,fac,istatus)
  ! mapfac_m = 1.0		! Test
  DO j=1,jmax
    DO i=1,imax
      CALL get_sigma(xlat(i,j),xlong(i,j),mapfac_m(i,j),istatus)
    ENDDO
  ENDDO

  ! T background:
  ! Convert T to perturbation potential temperature (theta-t0):
  DO k=1,kmax
    temptr3d(1:imax,1:jmax,k) = temptr3d(1:imax,1:jmax,k)* &
                               (100000.0/pressr1d(k))**(rc/cp)-t0
  ENDDO
  CALL laps2mass(temptr3d,imax,jmax,kmax,pressr1d,dam,znw,4,0,t_mass_bkg)

  ! QVAPOR background:
  CALL laps2mass(sphumd3d,imax,jmax,kmax,pressr1d,dam,znw,4,0,sh_mass_bkg)

  ! U background:
  CALL laps2mass(u_laps_bkg,imax,jmax,kmax,pressr1d,dam,znw,4,0,u_mass_bkg)

  ! V background:
  CALL laps2mass(v_laps_bkg,imax,jmax,kmax,pressr1d,dam,znw,4,0,v_mass_bkg)

  ! Height background:
  CALL laps2mass(height3d,imax,jmax,kmax,pressr1d,dam,znw,4,0,geo_mass_bkg)
  ! Geopotential:
  geo_mass_bkg = 9.80665*geo_mass_bkg

  ! Write out state variable for post-processing:
  CALL get_directory('log',filename,namelen)
  filename = filename(1:namelen)//'fort.12'
  open(unit=12,file=filename(1:namelen+7),form='unformatted')
  WRITE(12) znw,pressr1d,dam,znu
  close(12)

  times(1:4) = asctime(8:11)
  times(5:5) = '-'
  SELECT CASE (asctime(4:6))
    CASE ('JAN')
	times(6:7) = '01'
    CASE ('FEB')
	times(6:7) = '02'
    CASE ('MAR')
	times(6:7) = '03'
    CASE ('APR')
	times(6:7) = '04'
    CASE ('MAY')
	times(6:7) = '05'
    CASE ('JUN')
	times(6:7) = '06'
    CASE ('JUL')
	times(6:7) = '07'
    CASE ('AUG')
	times(6:7) = '08'
    CASE ('SEP')
	times(6:7) = '09'
    CASE ('OCT')
	times(6:7) = '10'
    CASE ('NOV')
	times(6:7) = '11'
    CASE ('DEC')
	times(6:7) = '12'
    CASE default
	PRINT*,'GSIb: error: invalid month: ',asctime(4:6)
	STOP
  END SELECT
  times(8:8) = '-'
  times(9:10) = asctime(1:2)
  times(11:11) = '_'
  times(12:13) = asctime(13:14)
  times(14:14) = ':'
  times(15:16) = asctime(15:16)
  times(17:19) = ':00'

  ! Write the variables into a wrf_inout netcdf file:
  CALL wrfbkgout(times,imax,jmax,kmax,pressr1d(kmax), &
     		 znu,znw,grid_spacing,mapfac_m,xlat, &
     		 xlong,dam,pdam,t_mass_bkg,geo_mass_bkg, &
     		 sh_mass_bkg,u_mass_bkg,v_mass_bkg,topo)

END SUBROUTINE GSI_Bkg


SUBROUTINE GSI_Obs(i4time,asctime,lat,lon,imax,jmax,kmax, &
     		  nobs,obs_point,n_tobs,obs_temp,maxtobs)

!==========================================================
!  This routine generates a bufr format data file for GSI.
!
!  HISTORY: MAR. 2006 by YUANFU XIE.
!==========================================================

  ! GSIobs reads LAPS ingest observation data (obs_point)
  ! and converts it into prepbufr format saved in prepqc.laps
  ! so that GSI 3DVAR can use.

  IMPLICIT NONE

  INCLUDE 'barnesob.inc'

  CHARACTER*16, intent(in) :: asctime
  INTEGER, intent(in) :: i4time,nobs,imax,jmax,kmax
  REAL, intent(in) :: lat(imax,jmax),lon(imax,jmax)
  TYPE (barnesob) :: obs_point(*)
  INTEGER, INTENT(IN) :: n_tobs,maxtobs
  REAL, INTENT(IN) :: obs_temp(maxtobs,12)

  INTEGER, parameter :: MXMN = 8
  INTEGER, parameter :: MXLV = 255
  REAL*8, parameter :: missing = 10.0e10
  REAL :: rlat,rlon,ztopsa,pres_3d(imax,jmax,kmax),p
  REAL*8 :: r8arr ( MXMN, MXLV ),rval

  INTEGER, parameter :: MXBF = 16000
  INTEGER :: ibfmsg ( MXBF/4 )

  INTEGER :: yearofreport,mnthofreport,daysofreport
  INTEGER :: hourofreport,istatus
  INTEGER :: idate,nlv,nlvst,jj,libf,ierr,iobs,zero

  CHARACTER	:: cval*8
  EQUIVALENCE ( cval, rval )

  ! Static access:
  CHARACTER :: dirstc*256,dir*200
  INTEGER :: dirlen

  ! Ingest obs:
  PRINT*,'Total wind obs into GSI: ',nobs
  PRINT*,'Total temp obs into GSI: ',n_tobs

  ! Read background total pressure 3d field:
  CALL get_pres_3d(i4time,imax,jmax,kmax,pres_3d,istatus)

  ! Open the BUFR messages file.
  OPEN( UNIT = 11, FILE = 'prepqc.laps', form='unformatted' )

  ! Open the BUFR tables file.
  CALL get_directory('static',dirstc,dirlen)
  dir = dirstc(1:dirlen)//'prepobs_prep.bufrtable'
  OPEN( UNIT = 12, FILE = dir(1:dirlen+23) )

  ! Associate the tables file with the messages file, and 
  ! identify the latter to the BUFRLIB software.

  CALL OPENBF  ( 11, 'OUT', 12 )

  ! Report date:
  yearofreport = 2005
  mnthofreport = 11
  daysofreport = 21
  hourofreport = 15
  yearofreport = 2006
  mnthofreport = 2
  daysofreport = 7
  hourofreport = 14
  zero = ichar('0')
  yearofreport = (ichar(asctime(8:8))-zero)*1000+ &
     		 (ichar(asctime(9:9))-zero)*100+ &
     		 (ichar(asctime(10:10))-zero)*10+ &
     		 (ichar(asctime(11:11))-zero)
  SELECT CASE (asctime(4:6))
    CASE ('JAN')
	mnthofreport = 1
    CASE ('FEB')
	mnthofreport = 2
    CASE ('MAR')
	mnthofreport = 3
    CASE ('APR')
	mnthofreport = 4
    CASE ('MAY')
	mnthofreport = 5
    CASE ('JUN')
	mnthofreport = 6
    CASE ('JUL')
	mnthofreport = 7
    CASE ('AUG')
	mnthofreport = 8
    CASE ('SEP')
	mnthofreport = 9
    CASE ('OCT')
	mnthofreport = 10
    CASE ('NOV')
	mnthofreport = 11
    CASE ('DEC')
	mnthofreport = 12
    CASE default
	PRINT*,'wrfbkgout: error: invalid month: ',asctime(4:6)
	STOP
  END SELECT
  daysofreport = (ichar(asctime(1:1))-zero)*10+ &
     		  ichar(asctime(2:2))-zero
  hourofreport = (ichar(asctime(13:13))-zero)*10+ &
     		  ichar(asctime(14:14))-zero

  ! For every wind observation data:
  DO iobs=1,nobs

    idate = ( ( yearofreport ) * 1000000 ) + 	&
     	    ( ( mnthofreport ) * 10000 )  +	&
            ( ( daysofreport ) * 100 ) +	&
     	    (   hourofreport       ) 

    ! Open a rawinsonde BUFR message in order to store the new
    ! data subset (i.e. report).

    CALL OPENMB  ( 11, 'ADPUPA', idate )

    ! Store the report date-time within the data subset.

    r8arr (1,1) = ( yearofreport )
    r8arr (2,1) = ( mnthofreport )
    r8arr (3,1) = ( daysofreport )
    r8arr (4,1) = ( hourofreport )

    CALL UFBSEQ  ( 11, r8arr, MXMN, 1, nlv, 'UARTM' )

    !   Store the station identification information within the
    !   data subset.
    !   cval = ( station ID, e.g. '72403', 'DBBH', etc.)
    cval = '72403'

    ! Convert to lat/lon:
    CALL rlapsgrid_to_latlon(obs_point(iobs)%ri,obs_point(iobs)%rj, &
                             lat,lon,imax,jmax,rlat,rlon,ierr)

    r8arr (1,1) = rval
    r8arr (3,1) = rlat		!( station latitude )
    r8arr (2,1) = rlon		!( station longitude )
    r8arr (4,1) = (obs_point(iobs)%i4time-i4time)/3600.0 !( obs time )
    r8arr (5,1) = 210.0		!( prepbufr report type )
    r8arr (6,1) = obs_point(iobs)%elev!( station elevation )

    CALL UFBINT  ( 11, r8arr, MXMN, 1, nlv, &
    		  'SID XOB YOB DHR TYP ELV ')

    ! Store the level data within the data subset.

    ! For LAPS data ingest, treat obs individually for now:
    nlvst = 1 	!( number of data levels to be stored)

    DO jj = 1, nlvst

      ! Use trilinear_laps.f to interpolate pressure value at ri,rj,rk:
      CALL trilinear_laps(obs_point(iobs)%ri, &
    			  obs_point(iobs)%rj, &
     			  obs_point(iobs)%rk,imax,jmax,kmax, &
     			  pres_3d,p)
      r8arr(1,jj) = p/100.0	! (in mb)

      ! r8arr (1,jj) = ztopsa(obs_point(iobs)%elev)

      !899.0 !missing !( pressure, in mb, for level jj )
      r8arr (2,jj) = missing !( SPECIFIC HUMIDITY OBSERVATION)
      r8arr (3,jj) = missing !( temperature in C for level jj )
      r8arr (4,jj) = missing !( height, in M, for level jj )
      r8arr (5,jj) = obs_point(iobs)%value(1) !( U M/S, for level jj )
      r8arr (6,jj) = obs_point(iobs)%value(2) !( V M/S, for level jj )
      r8arr (7,jj) = missing  !( Precipitable water mm, for level jj )
      r8arr (8,jj) = missing

    ENDDO

    CALL UFBINT ( 11, r8arr, MXMN, nlvst, nlv, &
     		'POB QOB TOB ZOB UOB VOB PWO CAT' )
    r8arr = 0.0
    CALL UFBINT ( 11, r8arr, MXMN, nlvst, nlv, &
                'PQM QQM TQM ZQM WQM NUL PWQ' )

    ! ( store any other available values in a similar manner )
    ! Once all data values have been stored for this data subset,
    ! we are now ready to store the data subset into the message.

    CALL WRITSA  ( 11, ibfmsg, libf )

  ENDDO

  ! For every temperature observation data:
  DO iobs=1,n_tobs

    idate = ( ( yearofreport ) * 1000000 ) + 	&
     	    ( ( mnthofreport ) * 10000 )  +	&
            ( ( daysofreport ) * 100 ) +	&
     	    (   hourofreport       ) 

    ! Open a rawinsonde BUFR message in order to store the new
    ! data subset (i.e. report).

    CALL OPENMB  ( 11, 'ADPUPA', idate )

    ! Store the report date-time within the data subset.

    r8arr (1,1) = ( yearofreport )
    r8arr (2,1) = ( mnthofreport )
    r8arr (3,1) = ( daysofreport )
    r8arr (4,1) = ( hourofreport )

    CALL UFBSEQ  ( 11, r8arr, MXMN, 1, nlv, 'UARTM' )

    !   Store the station identification information within the
    !   data subset.
    !   cval = ( station ID, e.g. '72403', 'DBBH', etc.)
    cval = '72403'

    ! Convert to lat/lon:
    CALL rlapsgrid_to_latlon(obs_temp(iobs,1),obs_temp(iobs,2), &
                             lat,lon,imax,jmax,rlat,rlon,ierr)

    r8arr (1,1) = rval
    r8arr (3,1) = rlat		!( station latitude )
    r8arr (2,1) = rlon		!( station longitude )
    r8arr (4,1) = 0.0 		!( obs time; Enhance this later)
    r8arr (5,1) = 132.0		!( prepbufr report type )
    r8arr (6,1) = missing !( station elevation )

    CALL UFBINT  ( 11, r8arr, MXMN, 1, nlv, &
    		  'SID XOB YOB DHR TYP ELV ')

    ! Store the level data within the data subset.

    ! For LAPS data ingest, treat obs individually for now:
    nlvst = 1 	!( number of data levels to be stored)

    DO jj = 1, nlvst

      ! Use trilinear_laps.f to interpolate pressure value at ri,rj,rk:
      CALL trilinear_laps(obs_temp(iobs,1), &
    			  obs_temp(iobs,2), &
     			  obs_temp(iobs,3),imax,jmax,kmax, &
     			  pres_3d,p)
      r8arr(1,jj) = p/100.0	! (in mb)

      ! r8arr (1,jj) = ztopsa(obs_point(iobs)%elev)

      !899.0 !missing !( pressure, in mb, for level jj )
      r8arr (2,jj) = missing !( SPECIFIC HUMIDITY OBSERVATION)
      r8arr (3,jj) = obs_temp(iobs,7) !( temperature in C for level jj )
      r8arr (4,jj) = missing !( height, in M, for level jj )
      r8arr (5,jj) = missing !( U M/S, for level jj )
      r8arr (6,jj) = missing !( V M/S, for level jj )
      r8arr (7,jj) = missing  !( Precipitable water mm, for level jj )
      r8arr (8,jj) = missing

    ENDDO

    CALL UFBINT ( 11, r8arr, MXMN, nlvst, nlv, &
     		'POB QOB TOB ZOB UOB VOB PWO CAT' )
    r8arr = 0.0
    CALL UFBINT ( 11, r8arr, MXMN, nlvst, nlv, &
                'PQM QQM TQM ZQM WQM NUL PWQ' )

    ! ( store any other available values in a similar manner )
    ! Once all data values have been stored for this data subset,
    ! we are now ready to store the data subset into the message.

    CALL WRITSA  ( 11, ibfmsg, libf )

  ENDDO
	    
  ! Forcibly flush the last BUFR message, if any, from the
  ! BUFRLIB software.

  CALL WRITSA  ( -11, ibfmsg, libf )
  CALL CLOSBF  ( 11 )

  PRINT*,'Wind data converted: ',nobs
  PRINT*,'Temp data converted: ',n_tobs

END SUBROUTINE GSI_Obs
