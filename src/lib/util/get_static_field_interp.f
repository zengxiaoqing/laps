      SUBROUTINE get_static_field_interp(ctype,i4time,nx,ny,data
     1,istatus)
   
! Returns a time-interpolated (valid for time) 2D array of static climatological
! albedo using the monthly values in the static file.  The monthly values are 
! valid on the 15th day of each month.  This routine only interpolates to the
! nearest day and does not account for leap years, but this should not be any
! big deal.  
!
! J.Smart 4-02: Subroutine taken from WRFSI software (module_wrfsi_static.F)
!               and modified for LAPS use to get albedo.

      IMPLICIT NONE
      CHARACTER(LEN=9)    :: ctime9
      CHARACTER           :: ctype*(*)

      INTEGER             :: nx,ny
      REAL                :: data(nx,ny)
    
      INTEGER             :: midmonth_day(12)
      INTEGER             :: valid_day
      INTEGER             :: yyyyjjj
      INTEGER             :: istatus
      REAL                :: sss
      INTEGER             :: i4time
      INTEGER             :: m, d1, d2, m1, m2
      CHARACTER(LEN=3)    :: var_2d
      REAL, ALLOCATABLE   :: data1(:,:),data2(:,:)
      REAL                :: w1, w2
      INTEGER             :: INT_FILE(9)
      INTEGER             :: I
      INTEGER             :: JDAY
!     INTEGER             :: NYEAR,JDAY,NHOUR,MIN,MONTH,NDAY

! midmonth_day is the julian day of the year corresponding to the 15th day
! of each month for a standard (non-leap) year

c     istatus = 0

      DATA midmonth_day / 15, 43, 74, 105, 135, 166, 196,
     &  227, 258, 288, 319, 349 /

! Convert date string into integer yyyyjjj and sss
      call make_fnam_lp(i4time,ctime9,istatus)
      if(istatus.ne.1)then
         print*,'Error returned: make_fnam_lp'
         return
      endif
      read(ctime9(3:5),'(i3)')valid_day

      PRINT*,'Time-interp monthly static field to day: ',valid_day

    ! Find bounding months
      IF ((valid_day .LT. midmonth_day(1)) .OR.
     &    (valid_day .GT. midmonth_day(12))) THEN
    ! December and January are bounding months
           d1 = midmonth_day(12)
           d2 = midmonth_day(1)
           m1 = 12
           m2 = 1
      ELSE
        find_bounds: DO m = 1, 11
           d1 = midmonth_day(m)
           d2 = midmonth_day(m+1)
           IF (valid_day .EQ. d1) THEN
                d2 = d1
                m1 = m
                m2 = m1
                EXIT find_bounds
           ELSE IF (valid_day .EQ. d2) THEN
                d1 = d2
                m1 = m + 1
                m2 = m1
                EXIT find_bounds
           ELSE IF ((valid_day .GT. d1).AND.(valid_day .LT. d2)) THEN
                m1 = m
                m2 = m + 1
                EXIT find_bounds
           ENDIF
        ENDDO find_bounds
      ENDIF

! If d1 = d2, then we don't need any interpolation, just get that month's 
! data values
      IF ( d1 .EQ. d2) THEN
           IF(ctype.eq.'albedo')THEN
              WRITE(var_2d, '("a",I2.2)') m1
           ELSEIF(ctype.eq.'green')THEN
              WRITE(var_2d, '("g",I2.2)') m1
           ENDIF
           call read_static_grid(nx,ny,var_2d,data,istatus)
           if(istatus.ne.1)then
c             print*,' Error reading LAPS static: ',var_2d
              return
           endif
      ELSE
! We need to get the two months of bounding data and time interpolate
         ALLOCATE(data1 (nx,ny))
         ALLOCATE(data2 (nx,ny))
         IF(ctype.eq.'albedo')THEN
            WRITE(var_2d, '("a",I2.2)') m1
         ELSEIF(ctype.eq.'green')THEN
            WRITE(var_2d, '("g",I2.2)') m1
         ENDIF
         call read_static_grid(nx,ny,var_2d,data1,istatus)
         if(istatus .ne. 1)then
c           print*,' Error reading LAPS static: ',var_2d
            return
         endif
         IF(ctype.eq.'albedo')THEN
            WRITE(var_2d, '("a",I2.2)') m2
         ELSEIF(ctype.eq.'green')THEN
            WRITE(var_2d, '("g",I2.2)') m2
         ENDIF
         call read_static_grid(nx,ny,var_2d,data2,istatus)
         if(istatus .ne. 1)then
c           print*,' Error reading LAPS static: ',var_2d
            return
         endif

!       Compute weights
        IF (d2 .GT. d1) THEN
            w1 = ( FLOAT(d2) - FLOAT(valid_day) ) / FLOAT(d2-d1)
        ELSE ! We must be between Dec 15 and Jan 15
            IF (valid_day .LT. midmonth_day(1)) THEN ! We are in January
                w1 = ( FLOAT(d2) - FLOAT(valid_day) ) / 31.
            ELSE ! We are in December
                w1 = ( 366. - FLOAT(valid_day) +
     &                        FLOAT(midmonth_day(1))) / 31.
            ENDIF
        ENDIF
        w2 = 1. - w1
        data = w1*data1 + w2*data2
        DEALLOCATE(data1)
        DEALLOCATE(data2)
      ENDIF

c     istatus=1
      RETURN
      END SUBROUTINE get_static_field_interp
