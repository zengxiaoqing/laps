  SUBROUTINE grib_sfc_vars(table_version,center_id,subcenter_id,&
     process_id,laps_reftime,laps_valtime,period_sec,igds, &
     nx,ny,tsfc,tdsfc,rhsfc,usfc,vsfc, &
     wsfc,pmsl,psfc,totpcpwater,pcp_inc,pcp_tot,snow_inc, snow_tot,&
     thetasfc,thetaesfc,cape,cin,srhel,liftedind,terdot,lwout, &
     swout,lwdown, swdown, shflux,lhflux,pblhgt,ground_t, &
     clwmrsfc,icemrsfc,rainmrsfc,snowmrsfc,graupmrsfc, &
     cldamt,cldbase,cldtop, &
     visibility,ceiling,echo_tops,max_refl,refl_sfc, &
     pcptype_sfc,funit,startb,nbytes)

    USE grib
    IMPLICIT NONE

    INTEGER,INTENT(IN)          :: table_version
    INTEGER,INTENT(IN)          :: center_id
    INTEGER,INTENT(IN)          :: process_id
    INTEGER,INTENT(IN)          :: subcenter_id
    INTEGER,INTENT(IN)          :: laps_reftime
    INTEGER,INTENT(IN)          :: laps_valtime
    INTEGER,INTENT(IN)          :: period_sec
    INTEGER,INTENT(IN)          :: igds(18)
    INTEGER,INTENT(IN)          :: nx
    INTEGER,INTENT(IN)          :: ny
    REAL,INTENT(IN)             :: tsfc(nx,ny)
    REAL,INTENT(IN)             :: tdsfc(nx,ny)
    REAL,INTENT(IN)             :: rhsfc(nx,ny)
    REAL,INTENT(IN)             :: usfc(nx,ny)
    REAL,INTENT(IN)             :: vsfc(nx,ny)
    REAL,INTENT(IN)             :: wsfc(nx,ny)
    REAL,INTENT(IN)             :: pmsl(nx,ny)
    REAL,INTENT(IN)             :: psfc(nx,ny)
    REAL,INTENT(IN)             :: totpcpwater(nx,ny)
    REAL,INTENT(IN)             :: pcp_inc(nx,ny)
    REAL,INTENT(IN)             :: pcp_tot(nx,ny)
    REAL,INTENT(IN)             :: snow_inc(nx,ny)
    REAL,INTENT(IN)             :: snow_tot(nx,ny)
    REAL,INTENT(IN)             :: thetasfc(nx,ny)
    REAL,INTENT(IN)             :: thetaesfc(nx,ny)
    REAL,INTENT(IN)             :: cape(nx,ny)
    REAL,INTENT(IN)             :: cin(nx,ny)
    REAL,INTENT(IN)             :: srhel(nx,ny)
    REAL,INTENT(IN)             :: liftedind(nx,ny)
    REAL,INTENT(IN)             :: terdot(nx,ny)
    REAL,INTENT(IN)             :: lwout(nx,ny)
    REAL,INTENT(IN)             :: swout(nx,ny)
    REAL,INTENT(IN)             :: lwdown(nx,ny)
    REAL,INTENT(IN)             :: swdown(nx,ny)
    REAL,INTENT(IN)             :: shflux(nx,ny)
    REAL,INTENT(IN)             :: lhflux(nx,ny)
    REAL,INTENT(IN)             :: pblhgt(nx,ny)
    REAL,INTENT(IN)             :: ground_t(nx,ny)
    REAL,INTENT(IN)             :: clwmrsfc(nx,ny)
    REAL,INTENT(IN)             :: icemrsfc(nx,ny)
    REAL,INTENT(IN)             :: rainmrsfc(nx,ny)
    REAL,INTENT(IN)             :: snowmrsfc(nx,ny)
    REAL,INTENT(IN)             :: graupmrsfc(nx,ny)
    REAL,INTENT(IN)             :: cldamt(nx,ny)
    REAL,INTENT(IN)             :: cldbase(nx,ny)
    REAL,INTENT(IN)             :: cldtop(nx,ny)
    REAL,INTENT(IN)             :: visibility(nx,ny)
    REAL,INTENT(IN)             :: ceiling(nx,ny)
    REAL,INTENT(IN)             :: echo_tops(nx,ny)
    REAL,INTENT(IN)             :: max_refl(nx,ny)
    REAL,INTENT(IN)             :: refl_sfc(nx,ny) 
    REAL,INTENT(IN)             :: pcptype_sfc(nx,ny)

    INTEGER,INTENT(IN)          :: funit
    INTEGER,INTENT(IN)          :: startb
    INTEGER,INTENT(OUT)         :: nbytes

    INTEGER                     :: i,j
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
    INTEGER                     :: startbyte
    REAL, PARAMETER             :: rmissing = 1.e36
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
    
    ! Grib up each variable...
    nbytes = startb - 1
    startbyte = nbytes+startb
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! Surface Temperature
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = tsfc(i,j+1)
      ENDDO
    ENDDO
    print *,'Gribbing Tsfc Min/Max = ',minval(fld),maxval(fld)
    itype = 0
    param = 11
    leveltype = 105
    level1 = 2
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
    startbyte = nbytes + 1 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Surface dewpoint
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = tdsfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing TdSfc Min/Max = ',minval(fld),maxval(fld)
    itype = 0
    param = 17
    leveltype = 105
    level1 = 2
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

    !!!!!!!!!!!!!!!!!!!!!!!
    ! Surface RH
    !!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = rhsfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Sfc Rh Min/Max = ',minval(fld),maxval(fld)
    itype = 0
    param = 52
    leveltype = 105
    level1 = 2
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Surface U wind  
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = usfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Sfc U Min/Max = ', minval(fld),maxval(fld)
    itype = 0
    param = 33
    leveltype = 105
    level1 = 10
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
   
    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Surface V wind  
    !!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = vsfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Sfc V Min/Max = ', minval(fld),maxval(fld)
    itype = 0
    param = 34
    leveltype = 105
    level1 = 10
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
    startbyte = nbytes+1
   
    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Surface W wind  
    !!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = wsfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Sfc W Min/Max = ', minval(fld),maxval(fld)
    itype = 0
    param = 40
    leveltype = 105
    level1 = 10
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 = 4
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)    
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Sea Level Pressure
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = pmsl(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing MSLP Min/MAx = ',minval(fld),maxval(fld)
    itype = 0
    param =  2
    leveltype = 102
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Surface Pressure
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = psfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Sfc P Min/Max = ', minval(fld),maxval(fld)
    itype = 0
    param =  1
    leveltype = 1  
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! total preciptable water
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = totpcpwater(i,j+1)*1000.
      ENDDO
    ENDDO
    print *, 'Gribbing total PW Min/Max = ',minval(fld),maxval(fld)
    itype = 0
    param =  54
    leveltype = 200
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Potential Temperature
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = thetasfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Theta Min/Max =', minval(fld),maxval(fld)
    itype = 0
    param =  13
    leveltype = 105
    level1 = 2
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Equiv Potential Temperature
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = thetaesfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing ThetaE Min/Max = ',minval(fld),maxval(Fld)
    itype = 0
    param =  14
    leveltype = 105
    level1 = 2
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! CAPE                        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = cape(i,j+1)
      ENDDO
    ENDDO
    print *,'Gribbing CAPE Min/MAx = ', minval(fld),maxval(fld)
    itype = 0
    param =  157
    leveltype = 1
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! CIN                         
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = cin(i,j+1)
      ENDDO
    ENDDO
    print *,'Gribbing CIN Min/Max =',minval(fld),maxval(fld)
    itype = 0
    param =  156
    leveltype = 1
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Storm Relative helicity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = srhel(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing SRH Min/MAx =', minval(fld),maxval(fld)
    itype = 0
    param =  190
    leveltype = 1
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! lifted index 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = liftedind(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing LI Min/Max = ',minval(fld),maxval(fld)
    itype = 0
    param =  131
    leveltype = 1
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Ground T 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = ground_t(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Tgd Min/Max = ',minval(fld),maxval(fld)
    itype = 0
    param =  11
    leveltype = 1
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Cloud Liquid     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = clwmrsfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing CLW Min/Max = ',minval(fld),maxval(fld)
    itype = 0
    param =  153
    leveltype = 1
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 = 4
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Cloud Ice   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = icemrsfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Ice Min/Max = ',minval(fld),maxval(fld)
    itype = 0
    param =  178
    leveltype = 1
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 = 4
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Rain MR     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = rainmrsfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Rain Min/Max = ', minval(fld),maxval(fld)
    itype = 0
    param =  170
    leveltype = 1
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 = 4
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Snow MR     
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = snowmrsfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing snow Min/Max = ',minval(fld),maxval(fld)
    itype = 0
    param =  171
    leveltype = 1
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 = 4
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Graupel MR     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = graupmrsfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Graupel Min/Max = ',minval(fld),maxval(Fld)
    itype = 0
    param =  179
    leveltype = 1
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 = 4
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Total Cloud Cover
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = cldamt(i,j+1)*100.
      ENDDO
    ENDDO
    print *, 'Gribbing LCV = ',minval(fld),maxval(Fld)
    itype = 0
    param =  71  
    leveltype = 1   
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! cloud base       
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = cldbase(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing cloud base Min/MAx = ', minval(fld),maxval(fld)
    itype = 0
    param =  138
    leveltype = 102
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! cloud top        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = cldtop(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Cloud Top Min/Max = ', minval(fld),maxval(Fld)
    itype = 0
    param =  139
    leveltype = 102
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! cloud ceiling 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = ceiling(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Cloud Ceiling Min/Max = ', minval(fld),maxval(Fld)
    itype = 0
    param =  137
    leveltype = 1
    level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Surface Radar Reflectivity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = refl_sfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Sfc Refl Min/Max = ', minval(fld),maxval(fld)
    itype = 0
    param =  128
    leveltype = 1
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 =  0
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1   

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Column Max Radar Reflectivity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = max_refl(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Col Max Refl Min/Max = ', minval(fld),maxval(fld)
    itype = 0
    param =  129
    leveltype = 1
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 =  0
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Max Echo Tops
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = echo_tops(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Max Echo Top Min/Max = ', minval(fld),maxval(fld)
    itype = 0
    param =  130
    leveltype = 102
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 =  0
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1                            
                                                                                   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Surface Precip Type
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = pcptype_sfc(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Ptype Sfc Min/Max = ', minval(fld),maxval(fld)
    itype = 0
    param =  136
    leveltype = 1
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 =  0
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1  
           
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Visbility 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = visibility(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing Vis Min/Max = ', minval(fld),maxval(fld)
    itype = 0
    param =  20
    leveltype = 1
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 = -3
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Some fields only present in true forecast fields
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (fcsttime_now .GT. 0) THEN
    
      ! Incremental Precip
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = pcp_inc(i,j+1)*1000.
        ENDDO
      ENDDO
      print *, 'Gribbing PCPINC Min/MAx =',minval(fld),maxval(fld)
      itype = 0
      param =  61
      leveltype = 1
      level1 = 0
      level2 = 0
      timerange = 4
      timeperiod1 = fcsttime_prev
      timeperiod2 = fcsttime_now
      scalep10 =  1
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Total Accum Precip
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = pcp_tot(i,j+1)*1000.
        ENDDO
      ENDDO
      print *,'Gribbing TotPCP Min/Max = ',minval(fld),maxval(fld)
      itype = 0
      param =  142
      leveltype = 1
      level1 = 0
      level2 = 0
      timerange = 4
      timeperiod1 = 0
      timeperiod2 = fcsttime_now
      scalep10 =  1
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Incremental Snow   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = snow_inc(i,j+1)
        ENDDO
      ENDDO
      print *, 'Gribbing SNOWINC Min/MAx = ',minval(fld),maxval(Fld)
      itype = 0
      param =  66  
      leveltype = 1
      level1 = 0
      level2 = 0
      timerange = 4
      timeperiod1 = fcsttime_prev
      timeperiod2 = fcsttime_now
      scalep10 =  4
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Total Accum  Snow   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = snow_tot(i,j+1)
        ENDDO
      ENDDO
      print *, 'Gribbing SNOW TOT Min/MAx =',minval(fld),maxval(fld)
      itype = 0
      param =  141
      leveltype = 1
      level1 = 0
      level2 = 0
      timerange = 4
      timeperiod1 = 0
      timeperiod2 = fcsttime_now
      scalep10 =  4
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Outgoing longwave radiation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (MAXVAL(lwout) .LT. rmissing) THEN
        DO j = 0,ny-1
          DO i = 1,nx
            fld(j*nx + i) = lwout(i,j+1)
          ENDDO
        ENDDO
        print *, 'GRibbing LWOUT Min/Max = ',minval(fld),maxval(fld)
        itype = 0
        param =  212
        leveltype = 8
        level1 = 0
        level2 = 0
        timerange = 0
        timeperiod1 = fcsttime_now
        timeperiod2 = 0
        scalep10 =  0
        CALL make_id(table_version,center_id,subcenter_id,process_id, &
                     param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
        CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
        nbytes = nbytes + itot
        startbyte=nbytes+1
      ENDIF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Outgoing shortwave radiation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (MAXVAL(swout) .LT. rmissing) THEN
        DO j = 0,ny-1
          DO i = 1,nx
            fld(j*nx + i) = swout(i,j+1)
          ENDDO
        ENDDO
        print *, 'Gribbing SWOUT Min/Max = ',minval(fld),maxval(fld)
        itype = 0
        param =  211
        leveltype = 8
        level1 = 0
        level2 = 0
        timerange = 0
        timeperiod1 = fcsttime_now
        timeperiod2 = 0
        scalep10 =  0
        CALL make_id(table_version,center_id,subcenter_id,process_id, &
                     param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
        CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
        nbytes = nbytes + itot
        startbyte=nbytes+1
      ENDIF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Incoming shortwave radiation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (MAXVAL(swdown) .LT. rmissing) THEN
        DO j = 0,ny-1
          DO i = 1,nx
            fld(j*nx + i) = swdown(i,j+1)
          ENDDO
        ENDDO
        print *, 'Gribbing SWDOWN Min/Max = ',minval(fld),maxval(fld)
        itype = 0
        param =  111
        leveltype = 1
        level1 = 0
        level2 = 0
        timerange = 0
        timeperiod1 = fcsttime_now
        timeperiod2 = 0
        scalep10 =  0
        CALL make_id(table_version,center_id,subcenter_id,process_id, &
                     param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
        CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
        nbytes = nbytes + itot
        startbyte=nbytes+1
      ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Incoming longwave radiation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (MAXVAL(lwdown) .LT. rmissing) THEN
        DO j = 0,ny-1
          DO i = 1,nx
            fld(j*nx + i) = lwdown(i,j+1)
          ENDDO
        ENDDO
        print *, 'Gribbing LWDOWN Min/Max = ',minval(fld),maxval(fld)
        itype = 0
        param =  112
        leveltype = 1
        level1 = 0
        level2 = 0
        timerange = 0
        timeperiod1 = fcsttime_now
        timeperiod2 = 0
        scalep10 =  0
        CALL make_id(table_version,center_id,subcenter_id,process_id, &
                     param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
        CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
        nbytes = nbytes + itot
        startbyte=nbytes+1
      ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Sensible Heat Flux
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (MAXVAL(shflux) .LT. rmissing) THEN
        DO j = 0,ny-1
          DO i = 1,nx
            fld(j*nx + i) = shflux(i,j+1)
          ENDDO
        ENDDO
        print *, 'Gribbing SHF Min/Max =',minval(fld),maxval(fld)
        itype = 0
        param =  122
        leveltype = 1
        level1 = 0
        level2 = 0
        timerange = 0
        timeperiod1 = fcsttime_now
        timeperiod2 = 0
        scalep10 =  0
        CALL make_id(table_version,center_id,subcenter_id,process_id, &
                     param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
        CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
        nbytes = nbytes + itot
        startbyte=nbytes+1
      ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Latent Heat Flux
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (MAXVAL(lhflux) .LT. rmissing) THEN
        DO j = 0,ny-1
          DO i = 1,nx
            fld(j*nx + i) = lhflux(i,j+1)
          ENDDO
        ENDDO
        print *,'Gribbing LHF Min/Max = ',minval(fld),maxval(fld)
        itype = 0
        param =  121
        leveltype = 1
        level1 = 0
        level2 = 0
        timerange = 0
        timeperiod1 = fcsttime_now
        timeperiod2 = 0
        scalep10 =  0
        CALL make_id(table_version,center_id,subcenter_id,process_id, &
                     param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                   hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                   scalep10,id)
        CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
        nbytes = nbytes + itot
        startbyte=nbytes+1
      ENDIF 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  PBL Height
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO j = 0,ny-1
        DO i = 1,nx
          fld(j*nx + i) = pblhgt(i,j+1)
        ENDDO
      ENDDO
      print *, 'Gribbing PBLHGT Min/MAx =',minval(fld),maxval(fld)
      itype = 0
      param =  221
      leveltype = 1
      level1 = 0
      level2 = 0
      timerange = 0
      timeperiod1 = fcsttime_now
      timeperiod2 = 0
      scalep10 =  0
      CALL make_id(table_version,center_id,subcenter_id,process_id, &
                   param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
      CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes = nbytes + itot
      startbyte=nbytes+1
    ENDIF

    ! Topography
    DO j = 0,ny-1
      DO i = 1,nx
        fld(j*nx + i) = terdot(i,j+1)
      ENDDO
    ENDDO
    print *, 'Gribbing TOPO Min/Max = ',minval(fld),maxval(Fld)
    itype = 0
    param =   7
    leveltype = 1
    level1 = 0
    level2 = 0
    timerange = 0
    timeperiod1 = fcsttime_now
    timeperiod2 = 0
    scalep10 =  0
    CALL make_id(table_version,center_id,subcenter_id,process_id, &
                 param,leveltype,level1,level2,yyyyr,mmr,ddr, &
                 hhr,minr,timeunit,timerange,timeperiod1,timeperiod2, &
                 scalep10,id)
    CALL write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
    nbytes = nbytes + itot
    startbyte=nbytes+1
    RETURN
  END SUBROUTINE grib_sfc_vars
        
