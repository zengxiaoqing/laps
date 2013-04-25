! Searches a directory for WRF v2.1 ARW files in netCDF format
!  (1 output time per file) that either match or surround (temporally)
!  a desired LAPS i4time
!
!  If it finds a perfect match, it returns nfiles=1 with the filename.
!  Else, if there are two bounding files, nfiles=2
!  Otherwise, nfiles=0

SUBROUTINE get_acceptable_wrf(bgpath,i4time_needed,nfiles, &
           i4times_found,filenames)

  USE time_utils

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: bgpath
  INTEGER,INTENT(IN)            :: i4time_needed
  INTEGER,INTENT(OUT)           :: nfiles
  INTEGER,INTENT(OUT)           :: i4times_found(2)
  CHARACTER(LEN=256)            :: filenames(3) ! plus 1 dim. by Wei-Ting (130312) to put previous time
  CHARACTER(LEN=256)            :: thisfile
  CHARACTER(LEN=256)            :: all_files(20000)

  ! Local vars
  INTEGER  :: delta_t, delta_t_1, delta_t_2,total_files,i, istatus
  INTEGER, PARAMETER  :: max_files = 500
  INTEGER, PARAMETER  :: dtmax = 86400
  INTEGER, PARAMETER  :: dtmin = -86400
  CHARACTER(LEN=11),PARAMETER   :: filter = "wrfout_d01_"
  CHARACTER(LEN=24)             :: wrftime
  INTEGER                       :: i4time
  INTEGER                       :: fname_len

  delta_t_1 = dtmax
  delta_t_2 = dtmin  
  CALL get_file_names(bgpath,total_files,all_files,max_files,istatus)

  DO i = 1, total_files 

    thisfile = all_files(i)
    fname_len = LEN_TRIM(all_files(i))
    IF (thisfile(fname_len-29:fname_len-19) .EQ. filter) THEN
      wrftime = thisfile(fname_len-18:fname_len) // ".0000" 
      CALL mm5_to_lapstime(wrftime,i4time)
      IF (i4time .EQ. i4time_needed) THEN
         print *, "File matches time needed: ",TRIM(thisfile)
         nfiles = 1
         i4times_found(1) = i4time
         filenames(1)     = thisfile
         IF (i .GT. 1) THEN
            filenames(2)  = all_files(i-1) ! Modified by Wei-Ting (130326) to get previous time
         ENDIF
         return
      ELSE
       ! delta_t = i4time_needed - i4time
        delta_t = i4time - i4time_needed
        IF (delta_t .GT. 0) THEN
          IF (delta_t .LT. delta_t_1) THEN
            delta_t_1 = delta_t
            i4times_found(1) = i4time
            filenames(1) = thisfile
          ENDIF
        ELSE
          IF (delta_t .GT. delta_t_2) THEN
            delta_t_2 = delta_t
            i4times_found(2) = i4time
            filenames(2) = thisfile
            IF (i .GT. 1) THEN
               filenames(3) = all_files(i-1) ! added by Wei-Ting (130326) to get previous time
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  
  ! If we are here, then we did not find an exact
  ! match.  We need to check if we got appropriate
  ! bounding files
 
  IF ( (i4times_found(1)-i4time_needed .LT. dtmax).AND. & ! add -i4time_needed by Wei-Ting (130326)
       (i4times_found(2)-i4time_needed .GT. dtmin) ) THEN ! add -i4time_needed by Wei-Ting (130326)
     nfiles = 2
  ELSE
     nfiles = 0
  ENDIF
  RETURN
END SUBROUTINE get_acceptable_wrf
