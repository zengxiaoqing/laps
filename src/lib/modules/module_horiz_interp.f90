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

MODULE horiz_interp

  ! Contains routines for doing horizontal interpolation
  IMPLICIT NONE
  PRIVATE

  ! Interp method codes for standard interpolation
  INTEGER, PUBLIC, PARAMETER          :: METHOD_NEAREST = 0
  INTEGER, PUBLIC, PARAMETER          :: METHOD_LINEAR  = 1
  INTEGER, PUBLIC, PARAMETER          :: METHOD_HIGHER = 2
  
  ! Data type codes for masked interpolation
  INTEGER, PUBLIC, PARAMETER          :: FLAG_FIELD = 0
  INTEGER, PUBLIC, PARAMETER          :: CATEGORY_FIELD = 1
  INTEGER, PUBLIC, PARAMETER          :: VALUE_FIELD = 2
  PUBLIC interpolate_standard
  PUBLIC interpolate_masked_val

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE interpolate_standard(nx_in, ny_in, data_in, &
                                  nx_out, ny_out, xloc, yloc, &
                                  method, data_out)
    !
    ! Performs standard interpolation (i.e., non-masked fields)
    IMPLICIT NONE

    ! Dimensions and grid we are interpolating from
    INTEGER, INTENT(IN)        :: nx_in
    INTEGER, INTENT(IN)        :: ny_in
    REAL, INTENT(IN)           :: data_in(nx_in,ny_in)

    ! Dimensions of output grid
    INTEGER, INTENT(IN)        :: nx_out
    INTEGER, INTENT(IN)        :: ny_out

    ! Pre-computed real i/j indices of input grid
    ! for each point in output grid
    REAL, INTENT(IN)           :: xloc(nx_out,ny_out)
    REAL, INTENT(IN)           :: yloc(nx_out,ny_out)

    ! Code indicating interp method (see module top)
    INTEGER, INTENT(IN)        :: method
   
    ! Output interpolated data
    REAL, INTENT(OUT)          :: data_out(nx_out,ny_out)
   
    ! Locals
    INTEGER                    :: i
    INTEGER                    :: j
    INTEGER                    :: ilo, jlo
    REAL                       :: dlon, dlat
    REAL                       :: deltalat, deltalon

    SELECT CASE(method)
      CASE(METHOD_NEAREST)
        DO j = 1, ny_out
          DO i = 1, nx_out
            ilo = NINT(xloc(i,j))
            jlo = NINT(yloc(i,j))
            data_out(i,j) = data_in( ilo, jlo )
          ENDDO
        ENDDO

      CASE(METHOD_LINEAR)
        DO j = 1, ny_out
          DO i = 1, nx_out
            ilo = MIN(FLOOR(xloc(i,j)), nx_in-1)
            jlo = MIN(FLOOR(yloc(i,j)), ny_in-1)
            dlon = xloc(i,j) - ilo
            dlat = yloc(i,j) - jlo
            deltalat = 1.
            deltalon = 1.
            data_out(i,j) = ((deltalon-dlon)*((deltalat-dlat)*data_in(ilo,jlo) &
                          + dlat*data_in(ilo,jlo+1))   &
                          + dlon*((deltalat-dlat)*data_in(ilo+1,jlo)     &
                          + dlat*data_in(ilo+1,jlo+1))) &
                          / ( deltalat * deltalon ) 
          ENDDO
        ENDDO
      CASE(METHOD_HIGHER) 
        DO j = 1, ny_out
          DO i = 1, nx_out
            data_out(i,j) = BINT( xloc(i,j), yloc(i,j), data_in, nx_in, ny_in, 0)
          ENDDO
        ENDDO
      CASE DEFAULT
        PRINT '(A,I2)', 'Uknown interpolation method code: ', method
        STOP 'INTERPOLATE_STANDARD'
    END SELECT
    RETURN
  END SUBROUTINE interpolate_standard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE interpolate_masked_val(nx_src, ny_src, lmask_src, data_src, &
                                 nx_out, ny_out, lmask_out, data_out, &
                                 isrcr, jsrcr, make_srcmask, &
                                 min_val, max_val, def_val, val_mask, &
                                 method)

    IMPLICIT none
    
    INTEGER, INTENT(IN)           :: nx_src
    INTEGER, INTENT(IN)           :: ny_src
    REAL, INTENT(INOUT)           :: lmask_src(nx_src,ny_src)
    REAL, INTENT(IN)              :: data_src(nx_src,ny_src)
    INTEGER, INTENT(IN)           :: nx_out, ny_out
    REAL, INTENT(IN)              :: lmask_out(nx_out,ny_out)
    REAL, INTENT(INOUT)             :: data_out(nx_out,ny_out)
    REAL, INTENT(IN)              :: isrcr(nx_out,ny_out)
    REAL, INTENT(IN)              :: jsrcr(nx_out,ny_out)
    LOGICAL, INTENT(IN)           :: make_srcmask
    REAL, INTENT(IN)              :: min_val, max_val, def_val, val_mask
    INTEGER, INTENT(IN)           :: method

    INTEGER                       :: i,j,k,l,ilo,jlo, ii, jj, kk, ll
    INTEGER                       :: search_rad
    INTEGER                       :: search_rad_max
    REAL, PARAMETER               :: bad_point = -99999.
    LOGICAL                       :: bad_points_flag
    LOGICAL                       :: fixed_bad
    INTEGER                       :: total_points
    INTEGER                       :: points_val_mask
    INTEGER                       :: points_skipped
    INTEGER                       :: num_to_fix
    INTEGER                       :: fixed_with_out
    INTEGER                       :: fixed_with_src
    INTEGER                       :: fixed_with_def
    INTEGER                       :: close_pts
    INTEGER                       :: ic,jc,iic,jjc
    REAL                          :: deltagx, deltagy
    REAL                          :: dix, djy
    REAL                          :: data_close(4)
    REAL                          :: distance(4)
    REAL                          :: distx,disty
    REAL                          :: sum_dist
    REAL                          :: a,b,c,d,e,f,g,h
    REAL                          :: stl(4,4)
    REAL                          :: valb
    REAL                          :: srcmask(nx_src,ny_src)

    ! Can we use the source data land mask that was input, or do we need to make it
    ! from the min_val/max_val arguments?  

    IF (make_srcmask) THEN
      PRINT '(A)', 'INTERPOLATE_MASKED_VAL: Making source data land mask...'
    
      ! This code has the flexibility to handled either water data points (lmask =0)
      ! or land points (lmask = 1). So if we are making the landmask using valid 
      ! range, but we do not know anything about the variable other than the 
      ! valid mask value, we need to figure out which points are 0 and which should be
      ! 1.  To do this, initialize the array to the invalid value, which we determine
      ! from the valid value.
      IF (val_mask .EQ. 1) THEN
        lmask_src(:,:) = 0
      ELSE
        lmask_src(:,:) = 1
      ENDIF
      
      ! Now figure out where to set the mask using a WHERE statement
      WHERE((data_src .GE. min_val).AND.(data_src .LE. max_val)) lmask_src = val_mask
    ELSE
      PRINT '(A)', 'INTERPOLATE_MASKED: Using source landmask field.'
    ENDIF

   
    bad_points_flag = .false.
    
    ! Initialize counters
    total_points = 0
    points_val_mask = 0
    points_skipped = 0
    num_to_fix = 0
    fixed_with_out = 0
    fixed_with_src = 0
    fixed_with_def = 0

    ! Select interpolation method.  Putting the case statement out here
    ! increases the amount of replicated code but should be more efficient
    ! than checking this condition at every grid point.

    SELECT CASE(method)
  
      CASE(METHOD_NEAREST)
        ! Use nearest neigbor
        PRINT '(A)', 'Masked interpolation using nearest neighbor value...'
        out_j_loop_1: DO j = 1, ny_out
          out_i_loop_1: DO i = 1, nx_out

            total_points = total_points + 1
            ! We only need to process this point if the lmask_out is equal
            ! to val_mask.  For example, if one is processing soil parameters,
            ! which are only valid at land points (lmask = 1), then the user
            ! passes in 1. for the val_mask.  Any output point that is water
            ! (lmask = 0) will then be skipped.  During this loop, if we have
            ! points that cannot be properly assigned, we will mark them as bad
            ! and set the bad_points_lag.

            IF (lmask_out(i,j) .EQ. val_mask) THEN
              ! Process this point 
              points_val_mask = points_val_mask + 1
              ilo = NINT(isrcr(i,j))
              jlo = NINT(jsrcr(i,j))
        
              ! See if this point can be used
              IF (lmask_src(ilo,jlo).EQ. lmask_out(i,j)) THEN
                data_out(i,j) = data_src(ilo,jlo)
              ELSE
                data_out(i,j) = bad_point
                num_to_fix = num_to_fix + 1
                bad_points_flag = .true.
              ENDIF 
            ELSE
              ! The output grid does not require a value for this point
              ! But do not zero out in case this is a field begin
              ! done twice (once for water and once for land, e.g.
              ! SKINTEMP
              !data_out(i,j) = 0.
              points_skipped = points_skipped + 1
            ENDIF
          ENDDO out_i_loop_1
        ENDDO out_j_loop_1
      
      CASE (METHOD_LINEAR)
        ! Use a 4-point interpolation
        PRINT '(A)', 'Masked interpolation using 4-pt linear interpolation...'
        deltagx = 1.
        deltagy = 1.
        out_j_loop_2: DO j = 1, ny_out
          out_i_loop_2: DO i = 1, nx_out

            total_points = total_points + 1
            ! We only need to process this point if the lmask_out is equal
            ! to val_mask.  For example, if one is processing soil parameters,
            ! which are only valid at land points (lmask = 1), then the user
            ! passes in 1. for the val_mask.  Any output point that is water
            ! (lmask = 0) will then be skipped.  During this loop, if we have
            ! points that cannot be properly assigned, we will mark them as bad
            ! and set the bad_points_lag.

            IF (lmask_out(i,j) .EQ. val_mask) THEN
              ! Process this point
              points_val_mask = points_val_mask + 1
              ilo = MIN(FLOOR(isrcr(i,j)),nx_src-1)
              jlo = MIN(FLOOR(jsrcr(i,j)),ny_src-1)
              dix = isrcr(i,j) - FLOAT(ilo)
              djy = jsrcr(i,j) - FLOAT(jlo)
 
              ! Loop around the four surrounding points
              ! and count up the number of points we can use based
              ! on common mask value
              close_pts = 0
              sum_dist = 0.
              outer_four_j: DO jc = 0,1
                outer_four_i: DO ic = 0,1
                  iic = ilo + ic
                  jjc = jlo + jc
                  IF (lmask_src(iic,jjc).EQ. lmask_out(i,j)) THEN
                    close_pts = close_pts + 1
                    data_close(close_pts) = data_src(iic,jjc)
                       
                    ! Compute distance to this valid point
                    ! in grid units and add to sum of distances (actually,
                    ! we are doing a weight, which is inversely proportional
                    ! to distance)
                    IF (ic .EQ. 0) THEN
                      distx = deltagx - dix
                    ELSE
                      distx =  dix
                    ENDIF
                    IF (jc .EQ. 0) THEN
                      disty = deltagy - djy
                    ELSE
                      disty = djy
                    ENDIF  
                    distance(close_pts) = SQRT(distx**2+disty**2)
                    sum_dist = sum_dist + distance(close_pts)
                  ENDIF 
                ENDDO outer_four_i
              ENDDO outer_four_j
 
              ! Did we find at least one point in the surrounding four 
              ! that was usable?

              IF (close_pts .GT. 0) THEN
               
                ! If we have all four points, then do bilinear interpolation
                IF (close_pts .EQ. 4) THEN
                   data_out(i,j) = ((deltagx - dix)*((deltagy-djy)*data_close(1) &
                                 + djy*data_close(3)) &
                                 + dix*((deltagy-djy)*data_close(2) &
                                 + djy*data_close(4))) &
                                 / (deltagx * deltagy)
                ELSE IF ((close_pts .GT. 1).AND.(close_pts .LT. 4)) THEN

                  ! Simple distance-weighted average by computing
                  ! the sum of all distances to each point and using
                  ! each individual distance divided by the total
                  ! distance as the weighting

                  data_out(i,j) = 0.
                  DO k = 1, close_pts
                    data_out(i,j) = data_out(i,j) + &
                                    (distance(k)/sum_dist) * data_close(k)
                  ENDDO
                ELSE
                  ! Set output value = to one point we found
                  data_out(i,j) = data_close(1)
                ENDIF
              ELSE
                bad_points_flag = .true.
                data_out(i,j) = bad_point
                num_to_fix = num_to_fix + 1   
              ENDIF 
            ELSE
              ! The output grid does not require a value for this point
              ! But do not zero out in case this is a field begin                             ! done twice (once for water and once for land, e.g.                            ! SKINTEMP
              !data_out(i,j) = 0.
              points_skipped = points_skipped + 1
            ENDIF
          ENDDO out_i_loop_2
        ENDDO out_j_loop_2                                        

      CASE (METHOD_HIGHER)
        ! 16-point interpolation 
        out_j_loop_3: DO j = 1, ny_out
          out_i_loop_3: DO i = 1, nx_out

            total_points = total_points + 1
            ! We only need to process this point if the lmask_out is equal
            ! to val_mask.  For example, if one is processing soil parameters,
            ! which are only valid at land points (lmask = 1), then the user
            ! passes in 1. for the val_mask.  Any output point that is water
            ! (lmask = 0) will then be skipped.  During this loop, if we have
            ! points that cannot be properly assigned, we will mark them as bad
            ! and set the bad_points_lag.

            IF (lmask_out(i,j) .EQ. val_mask) THEN
              ! Process this point
              points_val_mask = points_val_mask + 1   
 
              ! Do a 4x4 loop around in the input data around the output
              ! point to get our 16 points of influence.  Only use those
              ! that are appropriately masked
              valb = 0.
              close_pts = 0
              ilo = INT(isrcr(i,j)+0.00001)
              jlo = INT(jsrcr(i,j)+0.00001)
              dix = isrcr(i,j) - ilo
              djy = jsrcr(i,j) - jlo
              IF ( (ABS(dix).GT.0.0001).OR.(ABS(djy).GT.0.0001) ) THEN
                ! Do the interpolation loop
                stl(:,:) = 0.
                loop_16_1: DO k = 1,4
                  kk = ilo + k - 2
                  IF ((kk .LT. 1).OR.(kk .GT. nx_src)) CYCLE loop_16_1
                  loop_16_2: DO l = 1, 4
                    ll = jlo + l - 2
                    IF ((ll .LT. 1).OR.(ll .GT. ny_src)) CYCLE loop_16_2

                    ! Check land mask at this source point
                    IF (lmask_src(kk,ll).NE. val_mask) CYCLE loop_16_2
                    ! If we are here, then mask tests passed
                    stl(k,l) = data_src(kk,ll) 
                    IF ( (stl(k,l) .EQ. 0.).AND.(min_val.LE.0.).AND. &
                                                (max_val.GE.0.) ) THEN
                      stl = 1.E-5
                    ENDIF
                    close_pts = close_pts + 1
                  ENDDO loop_16_2
                ENDDO loop_16_1
  
                ! Did we find any valid points?

                IF ( (close_pts .GT. 0).AND. ( &
                  (stl(2,2).GT.0.).AND.(stl(2,3).GT.0.).AND. &
                  (stl(3,2).GT.0.).AND.(stl(3,3).GT.0.)  ) ) THEN
                  a = oned(dix,stl(1,1),stl(2,1),stl(3,1),stl(4,1))
                  b = oned(dix,stl(1,2),stl(2,2),stl(3,2),stl(4,2))
                  c = oned(dix,stl(1,3),stl(2,3),stl(3,3),stl(4,3))
                  d = oned(dix,stl(1,4),stl(2,4),stl(3,4),stl(4,4))
                  valb = oned(djy,a,b,c,d)
                  IF (close_pts .NE. 16) THEN
                    e = oned(djy,stl(1,1),stl(1,2),stl(1,3),stl(1,4))
                    f = oned(djy,stl(2,1),stl(2,2),stl(2,3),stl(2,4))
                    g = oned(djy,stl(3,1),stl(3,2),stl(3,3),stl(3,4))
                    h = oned(djy,stl(4,1),stl(4,2),stl(4,3),stl(4,4))
                    valb = (valb+oned(dix,e,f,g,h)) * 0.5
                  ENDIF
                  data_out(i,j) = valb
            
                ELSE
                  bad_points_flag = .true.
                  data_out(i,j) = bad_point
                  num_to_fix = num_to_fix + 1
                ENDIF
              ELSE
                ! We are right on a source point, so try to use it
                IF (lmask_src(ilo,jlo).EQ.val_mask) THEN
                  data_out(i,j) = data_src(ilo,jlo)
                ELSE
                  bad_points_flag = .true.
                  data_out(i,j) = bad_point
                  num_to_fix = num_to_fix + 1 
                ENDIF
              ENDIF
            ELSE
              ! The output grid does not require a value for this point
              !data_out(i,j) = 0.
              points_skipped = points_skipped + 1      
            ENDIF
          ENDDO out_i_loop_3
        ENDDO out_j_loop_3        

    END SELECT

     ! Do we need to correct bad points?
  
    IF (bad_points_flag) THEN
     
      search_rad_max = 10
      fix_bad_j: DO j = 1, ny_out
        fix_bad_i: DO i = 1, nx_out
 
          IF (data_out(i,j).NE. bad_point) CYCLE fix_bad_i
          
          ! First, search for nearest non-bad point in the output domain
          ! which is usually higher resolution. 
          fixed_bad = .false.
          search_out_loop: DO search_rad = 1, search_rad_max
            search_out_j: DO ll = -(search_rad-1), (search_rad-1),1
              jj = j + ll
              IF ((jj .LT. 1).OR.(jj .GT. ny_out)) CYCLE search_out_j
              search_out_i: DO kk = -(search_rad), search_rad, 1
                 ii = i + kk
                 IF ((ii .LT. 1).OR.(ii .GT. nx_out)) CYCLE search_out_j
                 IF ((data_out(ii,jj).NE.bad_point).AND. &
                    (lmask_out(ii,jj) .EQ. val_mask) ) THEN
                  data_out(i,j) = data_out(ii,jj)
                  fixed_bad = .true.
                  fixed_with_out = fixed_with_out + 1
                  EXIT search_out_loop
                ENDIF
              ENDDO search_out_i
            ENDDO search_out_j
          ENDDO search_out_loop

          ! Did we fix the point?  If not, then do same search on src data.
          IF (.NOT. fixed_bad) THEN
            search_rad_max = 10
            search_src_loop: DO search_rad = 1, search_rad_max
              search_src_j: DO ll = -(search_rad-1), (search_rad-1),1
                jj = NINT(jsrcr(i,j)) + ll
                IF ((jj .LT. 1).OR.(jj .GT. ny_src)) CYCLE search_src_j
                search_src_i: DO kk = -(search_rad), search_rad, 1
                   ii = NINT(isrcr(i,j)) + kk
                   IF ((ii .LT. 1).OR.(ii .GT. nx_src)) CYCLE search_src_j
                   IF (lmask_src(ii,jj).EQ.val_mask) THEN
                     data_out(i,j) = data_src(ii,jj)
                     fixed_bad = .true.
                     fixed_with_src = fixed_with_src + 1
                     EXIT search_src_loop
                   ENDIF
                ENDDO search_src_i
              ENDDO search_src_j
            ENDDO search_src_loop
          ENDIF
          ! Now is the point fixed?  If not, we have to use a default value.
          IF (.NOT.fixed_bad) THEN
            fixed_with_def = fixed_with_def + 1
            data_out(i,j) = def_val
            PRINT '(A,F10.3,A,2I5)', 'INTERPOLATE_MASKED: Bogus value of ', def_val, &
                ' used at point ', i, j
          ENDIF
         
        ENDDO fix_bad_i
      ENDDO fix_bad_j
    ENDIF
    PRINT '(A)',     '----------------------------------------'
    PRINT '(A)',     'MASKED INTERPOLATION SUMMARY: '
    PRINT '(A,I10)', '  TOTAL POINTS IN GRID:       ', total_points
    PRINT '(A,I10)', '  POINTS NEEDING VALUES:      ', points_val_mask
    PRINT '(A,I10)', '  POINTS NOT REQUIRED:        ', points_skipped
    PRINT '(A,I10)', '  POINTS NEEDING FIX:         ', num_to_fix
    PRINT '(A,I10)', '  POINTS FIXED WITH OUT GRID: ', fixed_with_out
    PRINT '(A,I10)', '  POINTS FIXED WITH SRC GRID: ', fixed_with_src
    PRINT '(A,I10)', '  POINTS FIXED WITH DEF VAL:  ', fixed_with_def
    PRINT '(A)',     '----------------------------------------'
    RETURN
  END SUBROUTINE interpolate_masked_val        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   FUNCTION bint(xx,yy,list,iii,jjj,ibint)
   
      !  16 Point interpolation
   
      IMPLICIT NONE
   
      REAL :: xx , yy
      INTEGER :: ibint , iii, jjj
      REAL list(iii,jjj),stl(4,4)
   
      INTEGER :: ib , jb, n , i , j , k , kk , l , ll
      REAL :: bint , x , y , a , b , c , d , e , f , g , h
      ib=iii-ibint
      jb=jjj-ibint
      bint = 0.0
      n = 0
      i = INT(xx+0.00001)
      j = INT(yy+0.00001)
      x = xx - i
      y = yy-j
     
      IF ( ( ABS(x).GT.0.0001 ) .OR. ( abs(y).gt.0.0001 ) ) THEN
         stl(:,:)=1.E-10      
         loop_1 : DO k = 1,4
            kk = i + k - 2
            IF ( ( kk .LT. 1) .OR. ( kk .GT. ib ) ) THEN
               CYCLE loop_1
            END IF
            loop_2 : DO l = 1,4
               stl(k,l) = 0.
               ll = j + l - 2
               IF ( ( ll .GT. jb ) .OR. ( ll .LT. 1 ) ) THEN
                  CYCLE loop_2
               END IF
               stl(k,l) = list(kk,ll)
               n = n + 1
               IF ( stl(k,l) .EQ. 0. ) THEN
                  stl(k,l) = 1.E-20
               END IF
            END DO loop_2
         END DO loop_1
         a = oned(x,stl(1,1),stl(2,1),stl(3,1),stl(4,1))
         b = oned(x,stl(1,2),stl(2,2),stl(3,2),stl(4,2))
         c = oned(x,stl(1,3),stl(2,3),stl(3,3),stl(4,3))
         d = oned(x,stl(1,4),stl(2,4),stl(3,4),stl(4,4))
         bint = oned(y,a,b,c,d)
   
         IF(n.NE.16) THEN
            e = oned(y,stl(1,1),stl(1,2),stl(1,3),stl(1,4))
            f = oned(y,stl(2,1),stl(2,2),stl(2,3),stl(2,4))
            g = oned(y,stl(3,1),stl(3,2),stl(3,3),stl(3,4))
            h = oned(y,stl(4,1),stl(4,2),stl(4,3),stl(4,4))
            bint = (bint+oned(x,e,f,g,h)) * 0.5
         END IF
   
      ELSE
         bint = list(i,j)
      END IF
   END FUNCTION bint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION oned(x,a,b,c,d) 
   
      IMPLICIT NONE
   
      REAL :: x,a,b,c,d,oned
   
      oned = 0.                
   
      IF      ( x .EQ. 0. ) THEN
         oned = b      
      ELSE IF ( x .EQ. 1. ) THEN
         oned = c      
      END IF
   
      IF(b*c.NE.0.) THEN
         IF ( a*d .EQ. 0. ) THEN
            IF      ( ( a .EQ. 0 ) .AND. ( d .EQ. 0 ) ) THEN
               oned = b*(1.0-x)+c*x                                        
            ELSE IF ( a .NE. 0. ) THEN
               oned = b+x*(0.5*(c-a)+x*(0.5*(c+a)-b))            
            ELSE IF ( d .NE. 0. ) THEN
               oned = c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c)) 
            END IF
         ELSE
            oned = (1.0-x)*(b+x*(0.5*(c-a)+x*(0.5*(c+a)-b)))+ &
                   x*(c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c)))                                   
         END IF
      END IF
   
   END FUNCTION oned 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                     
END MODULE horiz_interp
