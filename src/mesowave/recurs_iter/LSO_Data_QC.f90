SUBROUTINE LSO_Data

!==========================================================
!  This routine reads the LSO from LAPS.
!
!  HISTORY: MAY. 2004 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  REAL, PARAMETER :: temp0 = 273.16
  REAL :: badsfc

  INTEGER :: maxtime,mintime,hour,minute,seconds,midnight
  INTEGER :: imx,imm,mem_error,i,j,nsts,iobs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: idsts
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: naccn
  REAL    :: vmx,vmm
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: numob
  REAL,    ALLOCATABLE, DIMENSION(:,:,:) :: stobs

  ! Get the value for bad surface data:
  CALL get_sfc_badflag(badsfc,istatus)

  CALL get_laps_info(nx,ny,stanlat,stanlat2,stanlon, &
                     laps_cycle_time,badflag,maxstations,maproj)

  PRINT*,'START to allocate space'

  ALLOCATE(lat(nx,ny), lon(nx,ny), ldf(nx,ny),STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     print*,'LSO_Data_QC: cannot allocate space for lat/lon/ldf'
     STOP
  ENDIF
  ALLOCATE (olaps(nvlaps,maxstations*ncycles), &
            olat(maxstations*ncycles),olon(maxstations*ncycles), &
            otime(maxstations*ncycles),wght(maxstations*ncycles), &
            STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     print*,'LSO_Data_QC: cannot allocate space for olaps,olat,olon,otime,wght'
     STOP
  ENDIF

  ! Space for background fields:
  ALLOCATE(bkgd(nx,ny,ncycles,nvlaps),STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     print*,'LSO_Data: cannot allocate space for background fields'
     STOP
  ENDIF

  PRINT*,'Space allocated'
	
  ext_s = 'nest7grid'
  var_s = 'LAT'
  CALL rd_laps_static(dir_s,ext_s,nx,ny,1,var_s,units, &
                      comment,lat ,grid_spacingy,istatus)
  ! Debug: 
  print*,'Gridspace Y: ',grid_spacingy

  var_s = 'LON'
  CALL rd_laps_static(dir_s,ext_s,nx,ny,1,var_s,units, &
                      comment,lon ,grid_spacingx,istatus)
  ! Debug: 
  print*,'Gridspace X: ',grid_spacingx

  var_s = 'LDF'
  CALL rd_laps_static(dir_s,ext_s,nx,ny,1,var_s,units, &
                      comment,ldf ,grid_spacingx,istatus)
  DO j=1,ny
     DO i=1,nx
	IF (ldf(i,j) .GT. 1.0) ldf(i,j) = 1.0
	IF (ldf(i,j) .LT. 0.0) ldf(i,j) = 0.0
     ENDDO
  ENDDO

  ! clear unused portion of the array:
  olaps = badflag
  otime = -1.0   ! to search the maximum hour
  bkgd = 0.0
  CALL lso_reader_meso(maxstations,nvlaps,ncycles, &
       laps_cycle_time,badflag,olaps,wght,olat,olon, &
       otime,istarttime,nx,ny,bkgd)

  ! Debug: PRINT*,'Max: ',maxstations,laps_cycle_time
  ! Debug: PRINT*,'OBS: ',olaps(1:6,1),lat(1,1),lon(1,1)
  ! Debug: PRINT*,'OBS: ',olaps(1:6,2),lat(nx,ny),lon(nx,ny)

  ! Convert John's data to iterative recursive filter data:
  dm(1,1:2) = 1.0-nfic
  dm(2,1) = FLOAT(nx)+nfic
  dm(2,2) = FLOAT(ny)+nfic

  ! Midnight:
  midnight = 0
  IF (MAXVAL(otime) .GE. 2300.00) midnight = 1

  ! Space for QC analysis: 
  ALLOCATE(idsts(ncycles*maxstations*nvlaps), STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     PRINT*,'LSO_Data_QC: error in allocating space for idsts'
     STOP
  ENDIF

  nobs = 0
  maxtime = -100
  mintime = 100000000
  DO ilaps=1,maxstations*ncycles
     DO ivar=1,nvlaps
	IF ((olaps(ivar,ilaps) .NE. badflag) .AND. &
            (olaps(ivar,ilaps) .NE. badsfc)) THEN

 	   ! Valid obs:
	   nobs = nobs+1
	   IF (nobs .GT. mobs) THEN
	      PRINT*,'LSO_Data: too many obs!'
	      STOP
	   ENDIF

	   ! Separating hours and minutes:
	   hour = otime(ilaps)/100
	   minute = otime(ilaps)-hour*100

	   ! Handle data cross midnight:
 	   IF ((midnight .EQ. 1) .AND. (hour .LE. 2)) hour = 24+hour

	   ! Debug: PRINT*,'oll: ',ivar,olaps(ivar,ilaps), &
           ! olat(ilaps),olon(ilaps),otime(ilaps),hour,minute

	   seconds = hour*3600.0+minute*60.0

	   IF (maxtime .LT. seconds) maxtime = seconds
	   IF (mintime .GT. seconds) mintime = seconds

           vid(nobs) = ivar
	   o(1,nobs) = olaps(ivar,ilaps)

	   ! Convert F to K:
	   IF ((ivar .EQ. 1) .OR. (ivar .EQ. 5)) &
	      o(1,nobs) = (o(1,nobs)-32.0)*5.0/9.0+temp0

           ! Convert Altimeter to pascal:
           IF ((ivar .EQ. 4) .OR. (ivar .EQ. 6)) &
              o(1,nobs) = o(1,nobs)*100.0

           CALL latlon_to_rlapsgrid(olat(ilaps),olon(ilaps),lat,lon, &
                                 nx,ny,o(2,nobs),o(3,nobs),istatus)
	   ! Debug: print*,'Grid loc: ',o(2,nobs),o(3,nobs),nobs
           o(4,nobs) = seconds

	   ! Weight:
	   w(nobs) = 1.0
	ENDIF
     ENDDO
  ENDDO
  PRINT*,'MAX TIME: ',maxtime/3600,MOD(maxtime,3600)/60
  PRINT*,'MIN TIME: ',mintime/3600,MOD(mintime,3600)/60
  PRINT*,'ISTARTTIME: ',istarttime,MOD(istarttime,86400)/3600.0
  IF ((maxtime-mintime)/3600 .GT. 3) THEN
     PRINT*,'LSO_Data_QC: it is hard coded to run'
     PRINT*,' the analysis less than 3 hour!'
     STOP
  ENDIF

  PRINT*,'Total number of OBS: ',nobs

  ! QC:
  ! 1. Find out number of different stations and identical ones:
  nsts = 0
  DO iobs=1,nobs
     DO i=1,iobs-1
        IF ((ABS(o(2,i)-o(2,iobs)) .LT. 1.0e-3) .AND. &
            (ABS(o(3,i)-o(3,iobs)) .LT. 1.0e-3)) THEN
           idsts(iobs) = idsts(i)
           GOTO 11
        ENDIF
     ENDDO
     nsts = nsts+1
     idsts(iobs) = nsts
11   CONTINUE
  ENDDO
  PRINT*,'Number of obs stations: ',nsts
  
  ALLOCATE(stobs(ncycles*15,nvlaps,nsts), STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     PRINT*,'LSO_Data_QC: error to allocate space for stobs'
     STOP
  ENDIF
  ALLOCATE(numob(ncycles*15,nvlaps,nsts), STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     PRINT*,'LSO_Data_QC: error to allocate space for numob'
     STOP
  ENDIF
  ALLOCATE(naccn(nvlaps,nsts), STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     PRINT*,'LSO_Data_QC: error to allocate space for naccn'
     STOP
  ENDIF
  naccn = 0
  DO iobs=1,nobs
     naccn(vid(iobs),idsts(iobs)) = naccn(vid(iobs),idsts(iobs))+1
     numob(naccn(vid(iobs),idsts(iobs)),vid(iobs),idsts(iobs)) = iobs
     stobs(naccn(vid(iobs),idsts(iobs)),vid(iobs),idsts(iobs)) = o(1,iobs)
  ENDDO

  PRINT*,'naccn(1,1) = ',naccn(1,1),stobs(1:naccn(1,1),1,1)
  DO iobs=1,naccn(1,1)
     PRINT*,'First: ',o(2:4,numob(iobs,1,1)),olat(numob(iobs,1,1)), &
                                             olon(numob(iobs,1,1))
  ENDDO

  ! Time interval:
  dm(1,3) = MOD(istarttime,86400)  	! Second of the time
  dm(2,3) = dm(1,3)+(ncycles-1)*laps_cycle_time

  ! Spacings: x, y and t
  d(1:3) = (dm(2,1:3)-dm(1,1:3))/FLOAT(n(1:3)-1)

  PRINT*,'Analysis start TIME: ',INT(dm(1,3))/3600,MOD(INT(dm(1,3)),3600)/60
  PRINT*,'Analysis endng TIME: ',INT(dm(2,3))/3600,MOD(INT(dm(2,3)),3600)/60

  ! Check dimension for variables:
  IF (nvlaps .GT. mv) THEN
     PRINT*,'ReadObsn: Too much variables'
     STOP
  ENDIF
  ! Check if number of obs is fit into the array:
  IF (nobs .GT. mobs) THEN
     PRINT*,'ReadObsn: Too many obs'
     STOP
  ENDIF
  IF (nobs .LE. 0) THEN
     PRINT*,'ReadObsn: no obs to analyze'
     STOP
  ENDIF

  PRINT*,'Spacing: ',d

  ! WRITE surface.dat for testing:
  ! OPEN(unit=10,file='surface.dat',form='formatted')

  ! Domain info:
  PRINT*,'DOMAIN: ', dm(1,1),dm(2,1),dm(1,2),dm(2,2),dm(1,3),dm(2,3)
  vmx = -1000.0
  vmm = 1000.0
  DO ilaps=1,nobs

     IF (vid(ilaps) .EQ. 2) THEN
	IF (vmx .LT. o(1,ilaps)) THEN
	   vmx = o(1,ilaps)
	   imx = ilaps
        ENDIF
        IF (vmm .GT. o(1,ilaps)) THEN
	   vmm = o(1,ilaps)
	   imm = ilaps
        ENDIF
     ENDIF
     ! WRITE(10,1) vid(ilaps),o(1:4,ilaps),w(ilaps)

  ENDDO
1 FORMAT(i2,5e14.6)
  ! CLOSE(10)
  PRINT*,'MAX/MIN U obs: ',vmx,vmm,imx,imm,w(imm)

  PRINT*,'Time interval: ',MINVAL(o(4,1:nobs)),MAXVAL(o(4,1:nobs))

  ! Deallocate the memory:
  DEALLOCATE(lat,lon,STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     PRINT*,'LSO_Data_QC: error in deallocating lat/lon'
     STOP
  ENDIF
  DEALLOCATE(olaps,olat,olon,otime,wght,STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     print*,'LSO_Data_QC: error in deallocating olaps,olat,olon,otime,wght'
     STOP
  ENDIF
  
  DEALLOCATE(stobs, STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     PRINT*,'LSO_Data_QC: error to deallocate space for stobs'
     STOP
  ENDIF
  DEALLOCATE(numob, STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     PRINT*,'LSO_Data_QC: error to deallocate space for numob'
     STOP
  ENDIF
  DEALLOCATE(naccn, STAT=mem_error)
  IF (mem_error .NE. 0) THEN
     PRINT*,'LSO_Data_QC: error to deallocate space for naccn'
     STOP
  ENDIF
  DEALLOCATE(idsts,STAT=mem_error)

END SUBROUTINE LSO_Data
