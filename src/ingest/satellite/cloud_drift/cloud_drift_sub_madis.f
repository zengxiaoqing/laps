      subroutine get_cloud_drift_madis 
     ~           (i4time_sys, i4_window, filename, istatus)     

      IMPLICIT NONE
      include 'netcdf.inc'
      INTEGER, INTENT(IN)  :: i4time_sys, i4_window
      character*(*),INTENT(IN) ::  filename
      INTEGER, INTENT(OUT) :: istatus
      INTEGER, PARAMETER  :: lun_cdw = 11

      INTEGER :: i,nobs, ncid,nf_vid,nf_status
      REAL  r_missing_data
      REAL, ALLOCATABLE :: obLat(:),obLon(:),
     ~                     pressure(:),windDir(:), windSpd(:)

      INTEGER,PARAMETER::double=SELECTED_REAL_KIND(p=13,r=200) 
      REAL(kind=double), ALLOCATABLE :: validTime(:)
      CHARACTER(LEN=1),ALLOCATABLE ::windDirDD(:),windSpdDD(:)
      INTEGER, PARAMETER :: unix2i4 = 315619200
      CHARACTER(LEN=9)  :: a9timeObs
      INTEGER :: obtime, i4dif 
      INTEGER :: nkept,nreject
      call get_r_missing_data(r_missing_data,istatus)
      if ( istatus .ne. 1 )  then
         write (6,*) 'Error getting r_missing_data'
         return
      endif
      print *, "Opening ", TRIM(filename)
      ! Open the netcdf file
      nf_status = NF_OPEN(filename,NF_NOWRITE,ncid)
      IF (nf_status .NE. NF_NOERR) THEN
        print *, "Problem opening madis netcdf file."
        print *, "Filename:",trim(filename)
        istatus = 0
        return
      endif
      ! Get the number of records
      nf_status = NF_INQ_DIMID(ncid,'recNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Problem getting recNum'
        print *, ' in get_cloud_drift_madis'
        istatus = 0
        return
      endif
      nf_status = NF_INQ_DIMLEN(ncid,nf_vid,nobs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
        istatus = 0
        return
      endif
C
      print *, " MADIS Cloud Drift: nobs = ", nobs

      print *, " Allocating arrays..."
      ALLOCATE(obLat(nobs))
      ALLOCATE(obLon(nobs))
      ALLOCATE(validTime(nobs))
      ALLOCATE(pressure(nobs))
      ALLOCATE(windDir(nobs))
      ALLOCATE(windDirDD(nobs))
      ALLOCATE(windSpd(nobs))
      ALLOCATE(windSpdDD(nobs))
 
      print *, "Reading data" 
      ! Read in the data
      nf_status = NF_INQ_VARID(ncid,"obLat",nf_vid)
      IF (nf_status .NE. NF_NOERR) THEN
        print *," Problem getting obLat"
        istatus = 0
        GOTO 900
      ENDIF
      nf_status = NF_GET_VAR_REAL(ncid,nf_vid,obLat)

      nf_status = NF_INQ_VARID(ncid,"obLon",nf_vid)
      IF (nf_status .NE. NF_NOERR) THEN
        print *," Problem getting obLon"
        istatus = 0
        GOTO 900
      ENDIF
      nf_status = NF_GET_VAR_REAL(ncid,nf_vid,obLon)

      nf_status = NF_INQ_VARID(ncid,"validTime",nf_vid)
      IF (nf_status .NE. NF_NOERR) THEN
        print *," Problem getting validTime"
        istatus = 0
        GOTO 900
      ENDIF
      nf_status = NF_GET_VAR_DOUBLE(ncid,nf_vid,validTime)

      nf_status = NF_INQ_VARID(ncid,"pressure",nf_vid)
      IF (nf_status .NE. NF_NOERR) THEN
        print *," Problem getting pressure"
        istatus = 0
        GOTO 900
      ENDIF
      nf_status = NF_GET_VAR_REAL(ncid,nf_vid,pressure)

      nf_status = NF_INQ_VARID(ncid,"windDir",nf_vid)
      IF (nf_status .NE. NF_NOERR) THEN
        print *," Problem getting windDir"
        istatus = 0
        GOTO 900
      ENDIF
      nf_status = NF_GET_VAR_REAL(ncid,nf_vid,windDir)

      nf_status = NF_INQ_VARID(ncid,"windSpd",nf_vid)
      IF (nf_status .NE. NF_NOERR) THEN
        print *," Problem getting windSpd"
        istatus = 0
        GOTO 900
      ENDIF
      nf_status = NF_GET_VAR_REAL(ncid,nf_vid,windSpd)

      nf_status = NF_INQ_VARID(ncid,"windDirDD",nf_vid)
      IF (nf_status .NE. NF_NOERR) THEN
        print *," Problem getting windDirDD"
        istatus = 0
        GOTO 900
      ENDIF
      nf_status = NF_GET_VAR_TEXT(ncid,nf_vid,windDirDD)

      nf_status = NF_INQ_VARID(ncid,"windSpdDD",nf_vid)
      IF (nf_status .NE. NF_NOERR) THEN
        print *," Problem getting windSpdDD"
        istatus = 0
        GOTO 900
      ENDIF
      nf_status = NF_GET_VAR_TEXT(ncid,nf_vid,windSpdDD)

 
      call open_ext(lun_cdw,i4time_sys,'cdw',istatus)
      nreject =0
      nkept  = 0
      do i= 1,nobs
C       If this ob passes the quality flag
         obtime = NINT(validTime(i))+unix2i4
         i4dif = i4time_sys - obtime
         IF (ABS(i4dif) .LE. i4_window) THEN
           IF (windSpdDD(i) .NE. "C" .OR. 
     ~         windSpd(i) .LT. 0. .OR.
     ~         windSpd(i) .GT. 125.) THEN
             windSpd(i) = r_missing_data
           ENDIF
           IF (windDirDD(i) .NE. "C" .OR.
     ~         windDir(i) .LT. 0. .OR.
     ~         windDir(i) .GT. 360.) THEN
             windDir(i) = r_missing_data
           ENDIF
           CALL c_time2fname(NINT(validTime(i)),a9timeObs)

           write (lun_cdw,'(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)') 
     ~                 obLat(i), obLon(i), pressure(i), 
     ~                 windDir(i), windSpd(i), a9timeObs

           nkept = nkept+1
         ELSE
           nreject = nreject +1
         ENDIF
      enddo
      print *, "Total Obs/# kept/#rejected:",nobs,nkept,nreject
      istatus = 1
900   DEALLOCATE(obLat)
      DEALLOCATE(obLon)
      DEALLOCATE(validTime)
      DEALLOCATE(pressure)
      DEALLOCATE(windDir)
      DEALLOCATE(windDirDD)
      DEALLOCATE(windSpd)
      DEALLOCATE(windSpdDD)
      return
      end
