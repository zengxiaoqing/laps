  SUBROUTINE grib_ua_vars(table_version,center_id,subcenter_id, &
             process_id,laps_reftime,laps_valtime,period_sec,igds,nx,ny,kprs, &
             plvlmb,zprs,uprs,vprs,wprs,omprs,tprs,shprs,rhprs,&
             cldliqmr_prs,cldicemr_prs,rainmr_prs,snowmr_prs,graupelmr_prs, &
             pcptype_prs,refl_prs,tkeprs,funit,startb,nbytes)

    USE grib
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: table_version
    INTEGER,INTENT(IN)          :: center_id
    INTEGER,INTENT(IN)          :: subcenter_id
    INTEGER,INTENT(IN)          :: process_id
    INTEGER,INTENT(IN)          :: laps_reftime
    INTEGER,INTENT(IN)          :: laps_valtime
    INTEGER,INTENT(IN)          :: period_sec
    INTEGER,INTENT(IN)          :: igds(18)
    INTEGER,INTENT(IN)          :: nx
    INTEGER,INTENT(IN)          :: ny
    INTEGER,INTENT(IN)          :: kprs
    REAL,INTENT(IN)             :: plvlmb(kprs)
    REAL,INTENT(IN)             :: zprs(nx,ny,kprs)
    REAL,INTENT(IN)             :: uprs(nx,ny,kprs)
    REAL,INTENT(IN)             :: vprs(nx,ny,kprs)
    REAL,INTENT(IN)             :: wprs(nx,ny,kprs)
    REAL,INTENT(IN)             :: omprs(nx,ny,kprs)
    REAL,INTENT(IN)             :: tprs(nx,ny,kprs)
    REAL,INTENT(IN)             :: shprs(nx,ny,kprs)
    REAL,INTENT(IN)             :: rhprs(nx,ny,kprs)
    REAL,INTENT(IN)             :: cldliqmr_prs(nx,ny,kprs)
    REAL,INTENT(IN)             :: cldicemr_prs(nx,ny,kprs)
    REAL,INTENT(IN)             :: rainmr_prs(nx,ny,kprs)
    REAL,INTENT(IN)             :: snowmr_prs(nx,ny,kprs)
    REAL,INTENT(IN)             :: graupelmr_prs(nx,ny,kprs)
    REAL,INTENT(IN)             :: pcptype_prs(nx,ny,kprs)
    REAL,INTENT(IN)             :: refl_prs(nx,ny,kprs)
    REAL,INTENT(IN)             :: tkeprs(nx,ny,kprs)
    INTEGER,INTENT(IN)          :: funit
    INTEGER,INTENT(IN)          :: startb
    INTEGER,INTENT(OUT)         :: nbytes

    INTEGER                     :: i,j,k
    INTEGER                     :: lvl
    INTEGER                     :: itype
    INTEGER                     :: istatus
    INTEGER                     :: id(27)
    INTEGER                     :: param
    INTEGER                     :: leveltype
    INTEGER                     :: level1
    INTEGER                     :: level2
    INTEGER                     :: yyyyr
    INTEGER                     :: mmr
    INTEGER                     :: ddr
    INTEGER                     :: hhr
    INTEGER                     :: minr
    INTEGER                     :: timeunit
    INTEGER                     :: timerange
    INTEGER                     :: timeperiod1
    INTEGER                     :: timeperiod2
    INTEGER                     :: scalep10
    CHARACTER(LEN=24)           :: atime
    REAL                        :: fld(nx*ny)
    CHARACTER(LEN=3)            :: amonth
    CHARACTER(LEN=3)            :: amonths(12)
    INTEGER                     :: fcsttime_now
    INTEGER                     :: fcsttime_prev
    INTEGER                     :: itot
     INTEGER                    :: startbyte
    DATA amonths /'JAN','FEB','MAR','APR','MAY','JUN', &
                  'JUL','AUG','SEP','OCT','NOV','DEC'/
    ! Compute year, month, day of month, hour, and minute from laps_reftime

    CALL cv_i4tim_asc_lp(laps_reftime,atime,istatus)
    READ(atime,'(I2.2,x,A3,x,I4.4,x,I2.2,x,I2.2)') ddr,amonth,yyyyr, &
       hhr,minr
    DO i = 1, 12
      IF (amonth .eq. amonths(i)) THEN
         mmr = i
         EXIT
      ENDIF
    ENDDO

    ! Determine appropriate timeunit
      
    IF ( MOD(period_sec,3600) .EQ. 0) THEN
      ! Time unit shoud be hours
      timeunit = 1
      fcsttime_now = (laps_valtime-laps_reftime)/3600
      IF (fcsttime_now .GT. 0) THEN
        fcsttime_prev = fcsttime_now - (period_sec/3600)
      ELSE
        fcsttime_prev = 0
      ENDIF
    ELSE
      ! Time unit in minutes
      timeunit = 0
      fcsttime_now = (laps_valtime-laps_reftime)/60
      IF (fcsttime_now .GT. 0) THEN
        fcsttime_prev = fcsttime_now - (period_sec/60)
      ELSE
        fcsttime_prev = 0
      ENDIF
    ENDIF
    
    ! Grib up each variable at each level...

    nbytes = startb - 1
    startbyte = startb + nbytes

    levelloop: DO k = 1, kprs
      lvl = NINT(plvlmb(k))

      IF (lvl .LE. 1000) THEN
      ! Geopotential Height
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = zprs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  7
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 0
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte = nbytes+1
      ! Temperature          
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = tprs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  11
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 2
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte = nbytes+1
      ! Specific Humidity    
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = shprs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  51
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 8
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1 
      ! Relative Humidity
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = rhprs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  52
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 1
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1
      ! U wind component  
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = uprs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  33
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 1
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1
      ! V wind component  
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = vprs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  34
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 1
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1
      ! W wind component  
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = wprs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  40
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 3
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1
      ! Omega
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = omprs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  39
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 3
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1
      ! Cloud liquid
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = cldliqmr_prs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  153
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 6
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1
      ! Cloud ice
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = cldicemr_prs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  178
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 6
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1
      ! Rain          
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = rainmr_prs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  170
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 6
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1
      ! Snow          
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = snowmr_prs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  171
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 6
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1

      ! Graupel
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = graupelmr_prs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  179
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 6
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1

      ! Precip Type Code
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = pcptype_prs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  136
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 0
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1    

      ! Radar reflectivity
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = refl_prs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  128
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 0
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1  

      ! Turbulent kinetic energy
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = tkeprs(i,j+1,k)
        ENDDO
      ENDDO
      itype = 0
      param =  158
      leveltype = 100
      level1 = lvl
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 = 3
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1 

      ENDIF ! Only levels 1000mb and higher
    ENDDO levelloop
    RETURN
  END SUBROUTINE grib_ua_vars
