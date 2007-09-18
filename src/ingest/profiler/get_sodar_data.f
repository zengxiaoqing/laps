      subroutine get_sodar_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*170 filename

      integer level, maxStaticIds, nInventoryBins, recNum,nf_fid,
     +     nf_vid, nf_status
C
C  Open netcdf File for reading
C
      nf_status=NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),filename
        istatus=0
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of level
C
      nf_status=NF_INQ_DIMID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim level'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,level)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim level'
      endif
C
C Get size of maxStaticIds
C
      nf_status=NF_INQ_DIMID(nf_fid,'maxStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxStaticIds'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,maxStaticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxStaticIds'
      endif
C
C Get size of nInventoryBins
C
      nf_status=NF_INQ_DIMID(nf_fid,'nInventoryBins',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nInventoryBins'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,nInventoryBins)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nInventoryBins'
      endif
C
C Get size of recNum
C
      nf_status=NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      call read_sodar_data(nf_fid, level, maxStaticIds,
     +     nInventoryBins, recNum, i4time_sys, ilaps_cycle_time,
     +     NX_L, NY_L, i4time_earliest, i4time_latest, lun_out,
     +     istatus)

      return
      end
C
C
      subroutine read_sodar_data(nf_fid, level, maxStaticIds,
     +     nInventoryBins, recNum, i4time_sys, ilaps_cycle_time,
     +     NX_L, NY_L, i4time_earliest, i4time_latest, lun_out,
     +     istatus)


      include 'netcdf.inc'
      integer level, maxStaticIds, nInventoryBins, recNum,nf_fid,
     +     nf_vid, nf_status
      integer assetId(recNum), firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, invTime(recNum),
     +     inventory(maxStaticIds), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     nStaticIds, prevRecord(recNum), uQcFlag( level, recNum),
     +     uStdQcFlag( level, recNum), ugQcFlag( level, recNum),
     +     vQcFlag( level, recNum), vStdQcFlag( level, recNum),
     +     vgQcFlag( level, recNum), wQcFlag( level, recNum),
     +     wStdQcFlag( level, recNum), wdQcFlag( level, recNum),
     +     windDir( level, recNum), wsQcFlag( level, recNum)
      real elevation(recNum), latitude(recNum), levels( level,
     +     recNum), longitude(recNum), uComponent( level, recNum),
     +     uGustComponent( level, recNum), uStdDevComponent( level,
     +     recNum), vComponent( level, recNum), vGustComponent(
     +     level, recNum), vStdDevComponent( level, recNum),
     +     wComponent( level, recNum), wStdDevComponent( level,
     +     recNum), windSpeed( level, recNum)
      double precision observationTime(recNum), receiptTime(recNum),
     +     reportTime(recNum)
      character*6 providerId(recNum)
      character*24 dataProvider(recNum)
      character*30 staticIds(maxStaticIds)

!     Declarations for 'write_pro' call
      integer iwmostanum(recNum)
      character a9time_ob_r(recNum)*9
      character c8_obstype*8
      real height_m(level)
      real dir_deg(level)
      real spd_mps(level)

      logical l_closest_time, l_closest_time_i
      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting r_missing_data'
          return
      endif
      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a,
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_domain_perimeter'
          return
      endif

      call read_sodar_netcdf(nf_fid, level, maxStaticIds, 
     +     nInventoryBins, recNum, assetId, firstInBin, 
     +     firstOverflow, globalInventory, invTime, inventory, 
     +     isOverflow, lastInBin, lastRecord, nStaticIds, prevRecord, 
     +     uQcFlag, uStdQcFlag, ugQcFlag, vQcFlag, vStdQcFlag, 
     +     vgQcFlag, wQcFlag, wStdQcFlag, wdQcFlag, windDir, 
     +     wsQcFlag, elevation, latitude, levels, longitude, 
     +     uComponent, uGustComponent, uStdDevComponent, vComponent, 
     +     vGustComponent, vStdDevComponent, wComponent, 
     +     wStdDevComponent, windSpeed, observationTime, receiptTime, 
     +     reportTime, dataProvider, providerId, staticIds)
C
C The netcdf variables are filled - your pro write call may go here
C
!     Initial loop through obs to get times and stanums
      do iob = 1,recNum
          iwmostanum(iob) = assetId(iob)
          if(abs(observationTime(iob)) .le. 1e10)then
              i4time_ob = idint(observationTime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

      c8_obstype = 'SODAR   '

      do iob = 1,recNum
          call convert_array(levels(:,iob),height_m,level
     1                      ,'none',r_missing_data,istatus)

          call addcon_miss(height_m,elevation(iob),height_m,level,1)

          call convert_array_i2r(windDir(:,iob),dir_deg,level
     1                      ,'none',r_missing_data,istatus)
          call apply_qc_rsa(wdQcFlag(:,iob),dir_deg,level)
          call convert_array(dir_deg,dir_deg,level
     1                      ,'none',r_missing_data,istatus)

          call convert_array(windSpeed(:,iob),spd_mps,level
     1                      ,'none',r_missing_data,istatus)
          call apply_qc_rsa(wsQcFlag(:,iob),spd_mps,level)
          call convert_array(spd_mps,spd_mps,level
     1                      ,'none',r_missing_data,istatus)

          l_closest_time = l_closest_time_i(iwmostanum,a9time_ob_r
     1                                     ,recNum,iob,i4time_sys
     1                                     ,istatus)

          if(l_closest_time)then
!             call 'write_pro' for a single profile
              call open_ext(lun_out,i4time_sys,'pro',istatus)

              call write_pro(lun_out
     +                      ,1,level,1
     +                      ,assetId(iob)
     +                      ,latitude(iob),longitude(iob),elevation(iob)
     +                      ,providerId(iob)
     +                      ,a9time_ob_r(iob),c8_obstype
     +                      ,level
     +                      ,height_m
     +                      ,dir_deg
     +                      ,spd_mps
     +                      ,istatus)
          endif ! valid profile

      enddo ! iob
      return
      end
C
C  Subroutine to read the file 
C
      subroutine read_sodar_netcdf(nf_fid, level, maxStaticIds, 
     +     nInventoryBins, recNum, assetId, firstInBin, 
     +     firstOverflow, globalInventory, invTime, inventory, 
     +     isOverflow, lastInBin, lastRecord, nStaticIds, prevRecord, 
     +     uQcFlag, uStdQcFlag, ugQcFlag, vQcFlag, vStdQcFlag, 
     +     vgQcFlag, wQcFlag, wStdQcFlag, wdQcFlag, windDir, 
     +     wsQcFlag, elevation, latitude, levels, longitude, 
     +     uComponent, uGustComponent, uStdDevComponent, vComponent, 
     +     vGustComponent, vStdDevComponent, wComponent, 
     +     wStdDevComponent, windSpeed, observationTime, receiptTime, 
     +     reportTime, dataProvider, providerId, staticIds)
C
      include 'netcdf.inc'
      integer level, maxStaticIds, nInventoryBins, recNum,nf_fid, 
     +     nf_vid, nf_status
      integer assetId(recNum), firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, invTime(recNum),
     +     inventory(maxStaticIds), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     nStaticIds, prevRecord(recNum), uQcFlag( level, recNum),
     +     uStdQcFlag( level, recNum), ugQcFlag( level, recNum),
     +     vQcFlag( level, recNum), vStdQcFlag( level, recNum),
     +     vgQcFlag( level, recNum), wQcFlag( level, recNum),
     +     wStdQcFlag( level, recNum), wdQcFlag( level, recNum),
     +     windDir( level, recNum), wsQcFlag( level, recNum)
      real elevation(recNum), latitude(recNum), levels( level,
     +     recNum), longitude(recNum), uComponent( level, recNum),
     +     uGustComponent( level, recNum), uStdDevComponent( level,
     +     recNum), vComponent( level, recNum), vGustComponent(
     +     level, recNum), vStdDevComponent( level, recNum),
     +     wComponent( level, recNum), wStdDevComponent( level,
     +     recNum), windSpeed( level, recNum)
      double precision observationTime(recNum), receiptTime(recNum),
     +     reportTime(recNum)
      character*6 providerId(recNum)
      character*24 dataProvider(recNum)
      character*30 staticIds(maxStaticIds)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      elevation    "Elevation above MSL"
C
      nf_status=NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
C
C     Variable        NETCDF Long Name
C      latitude     "Station latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
C
C     Variable        NETCDF Long Name
C      levels       "Instrument Level, height above station"
C
      nf_status=NF_INQ_VARID(nf_fid,'levels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var levels'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,levels)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var levels'
      endif
C
C     Variable        NETCDF Long Name
C      longitude    "Station longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
C
C     Variable        NETCDF Long Name
C      uComponent   "u (eastward) component"
C
      nf_status=NF_INQ_VARID(nf_fid,'uComponent',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uComponent'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,uComponent)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uComponent'
      endif
C
C     Variable        NETCDF Long Name
C      uGustComponent"Gust u (eastward) component"
C
      nf_status=NF_INQ_VARID(nf_fid,'uGustComponent',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uGustComponent'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,uGustComponent)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uGustComponent'
      endif
C
C     Variable        NETCDF Long Name
C      uStdDevComponent"Std Dev u (eastward) component"
C
      nf_status=NF_INQ_VARID(nf_fid,'uStdDevComponent',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uStdDevComponent'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,uStdDevComponent)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uStdDevComponent'
      endif
C
C     Variable        NETCDF Long Name
C      vComponent   "v (northward) component"
C
      nf_status=NF_INQ_VARID(nf_fid,'vComponent',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vComponent'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,vComponent)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vComponent'
      endif
C
C     Variable        NETCDF Long Name
C      vGustComponent"Gust v (northward) component"
C
      nf_status=NF_INQ_VARID(nf_fid,'vGustComponent',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vGustComponent'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,vGustComponent)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vGustComponent'
      endif
C
C     Variable        NETCDF Long Name
C      vStdDevComponent"Std Dev v (northward) component"
C
      nf_status=NF_INQ_VARID(nf_fid,'vStdDevComponent',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vStdDevComponent'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,vStdDevComponent)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vStdDevComponent'
      endif
C
C     Variable        NETCDF Long Name
C      wComponent   "w (upward) component"
C
      nf_status=NF_INQ_VARID(nf_fid,'wComponent',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wComponent'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,wComponent)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wComponent'
      endif
C
C     Variable        NETCDF Long Name
C      wStdDevComponent"Std Dev w (upward) component"
C
      nf_status=NF_INQ_VARID(nf_fid,'wStdDevComponent',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wStdDevComponent'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,wStdDevComponent)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wStdDevComponent'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeed    "Wind Speed (scalar)"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      assetId      "RSA Asset Identifier"
C
      nf_status=NF_INQ_VARID(nf_fid,'assetId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var assetId'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,assetId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var assetId'
      endif
C
C     Variable        NETCDF Long Name
C      firstInBin   
C
      nf_status=NF_INQ_VARID(nf_fid,'firstInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstInBin'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,firstInBin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstInBin'
      endif
C
C     Variable        NETCDF Long Name
C      firstOverflow
C
      nf_status=NF_INQ_VARID(nf_fid,'firstOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstOverflow'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,firstOverflow)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstOverflow'
      endif
C
C     Variable        NETCDF Long Name
C      globalInventory
C
      nf_status=NF_INQ_VARID(nf_fid,'globalInventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var globalInventory'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,globalInventory)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var globalInventory'
      endif
C
C     Variable        NETCDF Long Name
C      invTime      
C
      nf_status=NF_INQ_VARID(nf_fid,'invTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var invTime'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,invTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var invTime'
      endif
C
C     Variable        NETCDF Long Name
C      inventory    
C
      nf_status=NF_INQ_VARID(nf_fid,'inventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var inventory'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,inventory)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var inventory'
      endif
C
C     Variable        NETCDF Long Name
C      isOverflow   
C
      nf_status=NF_INQ_VARID(nf_fid,'isOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isOverflow'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,isOverflow)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isOverflow'
      endif
C
C     Variable        NETCDF Long Name
C      lastInBin    
C
      nf_status=NF_INQ_VARID(nf_fid,'lastInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastInBin'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,lastInBin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastInBin'
      endif
C
C     Variable        NETCDF Long Name
C      lastRecord   
C
      nf_status=NF_INQ_VARID(nf_fid,'lastRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastRecord'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,lastRecord)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastRecord'
      endif
C
C     Variable        NETCDF Long Name
C      nStaticIds   
C
      nf_status=NF_INQ_VARID(nf_fid,'nStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nStaticIds'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,nStaticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nStaticIds'
      endif
C
C     Variable        NETCDF Long Name
C      prevRecord   
C
      nf_status=NF_INQ_VARID(nf_fid,'prevRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var prevRecord'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,prevRecord)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var prevRecord'
      endif
C
C     Variable        NETCDF Long Name
C      uQcFlag      "RSA mini-Sodar wind u-component quality control flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'uQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uQcFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,uQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      uStdQcFlag   "RSA mini-Sodar wind standard deviation u-component quality control flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'uStdQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uStdQcFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,uStdQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uStdQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      ugQcFlag     "RSA mini-Sodar wind gust u-component quality control flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'ugQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ugQcFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,ugQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ugQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      vQcFlag      "RSA mini-Sodar wind v-component quality control flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'vQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vQcFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,vQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      vStdQcFlag   "RSA mini-Sodar wind standard deviation v-component quality control flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'vStdQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vStdQcFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,vStdQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vStdQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      vgQcFlag     "RSA mini-Sodar wind gust v-component quality control flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'vgQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vgQcFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,vgQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vgQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      wQcFlag      "RSA mini-Sodar wind w-component quality control flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'wQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wQcFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,wQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      wStdQcFlag   "RSA mini-Sodar wind standard deviation w-component quality control flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'wStdQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wStdQcFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,wStdQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wStdQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      wdQcFlag     "RSA Wind direction quality control flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'wdQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wdQcFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,wdQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wdQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      windDir      "Wind Direction (scalar)"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
C
C     Variable        NETCDF Long Name
C      wsQcFlag     "RSA Wind speed quality control flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'wsQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wsQcFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,wsQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wsQcFlag'
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      observationTime"Time of observation"
C
      nf_status=NF_INQ_VARID(nf_fid,'observationTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,observationTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      endif
C
C     Variable        NETCDF Long Name
C      receiptTime  "File time stamp (time file was received)"
C
      nf_status=NF_INQ_VARID(nf_fid,'receiptTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var receiptTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,receiptTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var receiptTime'
      endif
C
C     Variable        NETCDF Long Name
C      reportTime   "Time of observation"
C
      nf_status=NF_INQ_VARID(nf_fid,'reportTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reportTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportTime'
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      dataProvider "Name of organization responsible for delivering the data"
C
      nf_status=NF_INQ_VARID(nf_fid,'dataProvider',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,dataProvider)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
      endif
C
C     Variable        NETCDF Long Name
C      providerId   "Alphanumeric station name"
C
      nf_status=NF_INQ_VARID(nf_fid,'providerId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var providerId'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,providerId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var providerId'
      endif
C
C     Variable        NETCDF Long Name
C      staticIds    
C
      nf_status=NF_INQ_VARID(nf_fid,'staticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staticIds'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,staticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staticIds'
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
