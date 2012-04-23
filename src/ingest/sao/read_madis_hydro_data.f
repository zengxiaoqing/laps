C
C  Subroutine to read the file "MADIS - Hydrological Surface" 
C
      subroutine read_madis_hydro_netcdf(nf_fid, ICcheckNum, 
     +     QCcheckNum, maxStaticIds, nInventoryBins, recNum, 
     +     elevation, latitude, longitude, precip12hr, precip12hrQCD, 
     +     precip1hr, precip1hrQCD, precip24hr, precip24hrQCD, 
     +     precip3hr, precip3hrQCD, precip5min, precip5minQCD, 
     +     precip6hr, precip6hrQCD, precipAccum, precipAccumQCD, 
     +     riverFlow, riverStage, ICT, QCT, dataProvider, 
     +     handbook5Id, homeWFO, precip12hrDD, precip1hrDD, 
     +     precip24hrDD, precip3hrDD, precip5minDD, precip6hrDD, 
     +     precipAccumDD, providerId, rawMessage, staticIds, 
     +     stationId, stationName, stationType, observationTime, 
     +     receivedTime, reportTime, riverReportChangeTime, 
     +     filterSetNum, firstInBin, firstOverflow, globalInventory, 
     +     invTime, inventory, isOverflow, lastInBin, lastRecord, 
     +     nStaticIds, numericWMOid, precip12hrICA, precip12hrICR, 
     +     precip12hrQCA, precip12hrQCR, precip1hrICA, precip1hrICR, 
     +     precip1hrQCA, precip1hrQCR, precip24hrICA, precip24hrICR, 
     +     precip24hrQCA, precip24hrQCR, precip3hrICA, precip3hrICR, 
     +     precip3hrQCA, precip3hrQCR, precip5minICA, precip5minICR, 
     +     precip5minQCA, precip5minQCR, precip6hrICA, precip6hrICR, 
     +     precip6hrQCA, precip6hrQCR, precipAccumICA, 
     +     precipAccumICR, precipAccumQCA, precipAccumQCR, 
     +     prevRecord, secondsStage1_2, secondsStage3,badflag)
C
      include 'netcdf.inc'
      integer ICcheckNum, QCcheckNum, maxStaticIds, nInventoryBins, 
     +     recNum,nf_fid, nf_vid, nf_status
      integer filterSetNum, firstInBin(nInventoryBins), firstOverflow,
     +     globalInventory, invTime(recNum), inventory(maxStaticIds),
     +     isOverflow(recNum), lastInBin(nInventoryBins),
     +     lastRecord(maxStaticIds), nStaticIds,
     +     numericWMOid(recNum), precip12hrICA(recNum),
     +     precip12hrICR(recNum), precip12hrQCA(recNum),
     +     precip12hrQCR(recNum), precip1hrICA(recNum),
     +     precip1hrICR(recNum), precip1hrQCA(recNum),
     +     precip1hrQCR(recNum), precip24hrICA(recNum),
     +     precip24hrICR(recNum), precip24hrQCA(recNum),
     +     precip24hrQCR(recNum), precip3hrICA(recNum),
     +     precip3hrICR(recNum), precip3hrQCA(recNum),
     +     precip3hrQCR(recNum), precip5minICA(recNum),
     +     precip5minICR(recNum), precip5minQCA(recNum),
     +     precip5minQCR(recNum), precip6hrICA(recNum),
     +     precip6hrICR(recNum), precip6hrQCA(recNum),
     +     precip6hrQCR(recNum), precipAccumICA(recNum),
     +     precipAccumICR(recNum), precipAccumQCA(recNum),
     +     precipAccumQCR(recNum), prevRecord(recNum),
     +     secondsStage1_2(recNum), secondsStage3(recNum)
      real elevation(recNum), latitude(recNum), longitude(recNum),
     +     precip12hr(recNum), precip12hrQCD( QCcheckNum, recNum),
     +     precip1hr(recNum), precip1hrQCD( QCcheckNum, recNum),
     +     precip24hr(recNum), precip24hrQCD( QCcheckNum, recNum),
     +     precip3hr(recNum), precip3hrQCD( QCcheckNum, recNum),
     +     precip5min(recNum), precip5minQCD( QCcheckNum, recNum),
     +     precip6hr(recNum), precip6hrQCD( QCcheckNum, recNum),
     +     precipAccum(recNum), precipAccumQCD( QCcheckNum, recNum),
     +     riverFlow(recNum), riverStage(recNum)
      double precision observationTime(recNum), receivedTime(recNum),
     +     reportTime(recNum), riverReportChangeTime(recNum)
      character precip5minDD(recNum)
      character precip12hrDD(recNum)
      character*72 ICT(ICcheckNum)
      character precip1hrDD(recNum)
      character*11 dataProvider(recNum)
      character*60 QCT(QCcheckNum)
      character*11 handbook5Id(recNum)
      character precip24hrDD(recNum)
      character precipAccumDD(recNum)
      character*51 stationName(recNum)
      character precip3hrDD(recNum)
      character*11 stationType(recNum)
      character precip6hrDD(recNum)
      character*256 rawMessage(recNum)
      character*24 staticIds(maxStaticIds)
      character*12 providerId(recNum)
      character*11 stationId(recNum)
      character*4 homeWFO(recNum)

      write(6,*)' subroutine read_madis_hydro_netcdf: recNum = ',recNum

C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     elevation     "elevation"
C
      nf_status=NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for elevation'
       print *,'Set elevation to badflag'
       elevation = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for elevation'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - elevation'
       else
        call ck_array_real(elevation,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - elevation'
       else
        call ck_array_real(elevation,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     latitude      "latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for latitude'
       print *,'Set latitude to badflag'
       latitude = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for latitude'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - latitude'
       else
        call ck_array_real(latitude,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - latitude'
       else
        call ck_array_real(latitude,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     longitude     "longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for longitude'
       print *,'Set longitude to badflag'
       longitude = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for longitude'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - longitude'
       else
        call ck_array_real(longitude,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - longitude'
       else
        call ck_array_real(longitude,recNum,valmis
     1                    ,badflag)
       endif
      endif

C
C     Variable        NETCDF Long Name
C     precip12hr    "12 hour precip accumulation "
C
      nf_status=NF_INQ_VARID(nf_fid,'precip12hr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip12hr'
       print *,'Set precip12hr to badflag'
       precip12hr = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip12hr)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip12hr'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip12hr'
       else
        call ck_array_real(precip12hr,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip12hr'
       else
        call ck_array_real(precip12hr,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip12hrQCD "12-hr precip amount QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip12hrQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip12hrQCD'
       print *,'Set precip12hrQCD to badflag'
       precip12hrQCD = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip12hrQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip12hrQCD'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip12hrQCD'
       else
        call ck_array_real(precip12hrQCD,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip12hrQCD'
       else
        call ck_array_real(precip12hrQCD,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip1hr     "1 hour precip accumulation "
C
      nf_status=NF_INQ_VARID(nf_fid,'precip1hr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip1hr'
       print *,'Set precip1hr to badflag'
       precip1hr = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip1hr)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip1hr'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip1hr'
       else
        call ck_array_real(precip1hr,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip1hr'
       else
        call ck_array_real(precip1hr,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip1hrQCD  "1-hr precip amount QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip1hrQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip1hrQCD'
       print *,'Set precip1hrQCD to badflag'
       precip1hrQCD = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip1hrQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip1hrQCD'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip1hrQCD'
       else
        call ck_array_real(precip1hrQCD,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip1hrQCD'
       else
        call ck_array_real(precip1hrQCD,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip24hr    "24 hour precip accumulation "
C
      nf_status=NF_INQ_VARID(nf_fid,'precip24hr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip24hr'
       print *,'Set precip24hr to badflag'
       precip24hr = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip24hr)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip24hr'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip24hr'
       else
        call ck_array_real(precip24hr,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip24hr'
       else
        call ck_array_real(precip24hr,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip24hrQCD "24-hr precip amount QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip24hrQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip24hrQCD'
       print *,'Set precip24hrQCD to badflag'
       precip24hrQCD = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip24hrQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip24hrQCD'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip24hrQCD'
       else
        call ck_array_real(precip24hrQCD,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip24hrQCD'
       else
        call ck_array_real(precip24hrQCD,recNum,valmis
     1                    ,badflag)
       endif
      endif

C
C     Variable        NETCDF Long Name
C     precip3hr     "3 hour precip accumulation "
C
      nf_status=NF_INQ_VARID(nf_fid,'precip3hr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip3hr'
       print *,'Set precip3hr to badflag'
       precip3hr = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip3hr)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip3hr'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip3hr'
       else
        call ck_array_real(precip3hr,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip3hr'
       else
        call ck_array_real(precip3hr,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip3hrQCD  "3-hr precip amount QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip3hrQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip3hrQCD'
       print *,'Set precip3hrQCD to badflag'
       precip3hrQCD = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip3hrQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip3hrQCD'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip3hrQCD'
       else
        call ck_array_real(precip3hrQCD,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip3hrQCD'
       else
        call ck_array_real(precip3hrQCD,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip5min    "5 minute precip accumulation "
C
      nf_status=NF_INQ_VARID(nf_fid,'precip5min',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip5min'
       print *,'Set precip5min to badflag'
       precip5min = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip5min)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip5min'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip5min'
       else
        call ck_array_real(precip5min,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip5min'
       else
        call ck_array_real(precip5min,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip5minQCD "5-min precip amount QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip5minQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip5minQCD'
       print *,'Set precip5minQCD to badflag'
       precip5minQCD = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip5minQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip5minQCD'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip5minQCD'
       else
        call ck_array_real(precip5minQCD,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip5minQCD'
       else
        call ck_array_real(precip5minQCD,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip6hr     "6 hour precip accumulation "
C
      nf_status=NF_INQ_VARID(nf_fid,'precip6hr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip6hr'
       print *,'Set precip6hr to badflag'
       precip6hr = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip6hr)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip6hr'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip6hr'
       else
        call ck_array_real(precip6hr,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip6hr'
       else
        call ck_array_real(precip6hr,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip6hrQCD  "6-hr precip amount QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip6hrQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip6hrQCD'
       print *,'Set precip6hrQCD to badflag'
       precip6hrQCD = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precip6hrQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip6hrQCD'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precip6hrQCD'
       else
        call ck_array_real(precip6hrQCD,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precip6hrQCD'
       else
        call ck_array_real(precip6hrQCD,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccum   "precip accumulation with an unknown time period"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipAccum'
       print *,'Set precipAccum to badflag'
       precipAccum = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipAccum)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipAccum'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precipAccum'
       else
        call ck_array_real(precipAccum,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precipAccum'
       else
        call ck_array_real(precipAccum,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccumQCD"precip amount QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipAccumQCD'
       print *,'Set precipAccumQCD to badflag'
       precipAccumQCD = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipAccumQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipAccumQCD'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precipAccumQCD'
       else
        call ck_array_real(precipAccumQCD,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precipAccumQCD'
       else
        call ck_array_real(precipAccumQCD,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     riverFlow     "River flow"
C
      nf_status=NF_INQ_VARID(nf_fid,'riverFlow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for riverFlow'
       print *,'Set riverFlow to badflag'
       riverFlow = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,riverFlow)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for riverFlow'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - riverFlow'
       else
        call ck_array_real(riverFlow,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - riverFlow'
       else
        call ck_array_real(riverFlow,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     riverStage    "River stage"
C
      nf_status=NF_INQ_VARID(nf_fid,'riverStage',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for riverStage'
       print *,'Set riverStage to badflag'
       riverStage = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,riverStage)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for riverStage'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - riverStage'
       else
        call ck_array_real(riverStage,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - riverStage'
       else
        call ck_array_real(riverStage,recNum,valmis
     1                    ,badflag)
       endif
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C     filterSetNum  
C
      nf_status=NF_INQ_VARID(nf_fid,'filterSetNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for filterSetNum'
       print *,'Set filterSetNum to -99'
       filterSetNum = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,filterSetNum)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for filterSetNum'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     firstInBin    
C
      nf_status=NF_INQ_VARID(nf_fid,'firstInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for firstInBin'
       print *,'Set firstInBin to -99'
       firstInBin = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,firstInBin)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for firstInBin'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     firstOverflow 
C
      nf_status=NF_INQ_VARID(nf_fid,'firstOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for firstOverflow'
       print *,'Set firstOverflow to -99'
       firstOverflow = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,firstOverflow)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for firstOverflow'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     globalInventory
C
      nf_status=NF_INQ_VARID(nf_fid,'globalInventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for globalInventory'
       print *,'Set globalInventory to -99'
       globalInventory = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,globalInventory)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for globalInventory'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     invTime       
C
      nf_status=NF_INQ_VARID(nf_fid,'invTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for invTime'
       print *,'Set invTime to -99'
       invTime = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,invTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for invTime'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     inventory     
C
      nf_status=NF_INQ_VARID(nf_fid,'inventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for inventory'
       print *,'Set inventory to -99'
       inventory = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,inventory)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for inventory'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     isOverflow    
C
      nf_status=NF_INQ_VARID(nf_fid,'isOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for isOverflow'
       print *,'Set isOverflow to -99'
       isOverflow = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,isOverflow)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for isOverflow'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     lastInBin     
C
      nf_status=NF_INQ_VARID(nf_fid,'lastInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for lastInBin'
       print *,'Set lastInBin to -99'
       lastInBin = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,lastInBin)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for lastInBin'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     lastRecord    
C
      nf_status=NF_INQ_VARID(nf_fid,'lastRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for lastRecord'
       print *,'Set lastRecord to -99'
       lastRecord = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,lastRecord)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for lastRecord'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     nStaticIds    
C
      nf_status=NF_INQ_VARID(nf_fid,'nStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for nStaticIds'
       print *,'Set nStaticIds to -99'
       nStaticIds = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,nStaticIds)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for nStaticIds'
       endif
      endif

C
C     Variable        NETCDF Long Name
C     numericWMOid  "numeric WMO identification"
C
      nf_status=NF_INQ_VARID(nf_fid,'numericWMOid',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for numericWMOid'
       print *,'Set numericWMOid to -99'
       numericWMOid = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numericWMOid)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for numericWMOid'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip12hrICA "12-hr precip amount IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip12hrICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip12hrICA'
       print *,'Set precip12hrICA to -99'
       precip12hrICA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip12hrICA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip12hrICA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip12hrICR "12-hr precip amount IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip12hrICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip12hrICR'
       print *,'Set precip12hrICR to -99'
       precip12hrICR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip12hrICR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip12hrICR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip12hrQCA "12-hr precip amount QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip12hrQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip12hrQCA'
       print *,'Set precip12hrQCA to -99'
       precip12hrQCA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip12hrQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip12hrQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip12hrQCR "12-hr precip amount QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip12hrQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip12hrQCR'
       print *,'Set precip12hrQCR to -99'
       precip12hrQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip12hrQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip12hrQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip1hrICA  "1-hr precip amount IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip1hrICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip1hrICA'
       print *,'Set precip1hrICA to -99'
       precip1hrICA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip1hrICA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip1hrICA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip1hrICR  "1-hr precip amount IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip1hrICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip1hrICR'
       print *,'Set precip1hrICR to -99'
       precip1hrICR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip1hrICR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip1hrICR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip1hrQCA  "1-hr precip amount QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip1hrQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip1hrQCA'
       print *,'Set precip1hrQCA to -99'
       precip1hrQCA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip1hrQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip1hrQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip1hrQCR  "1-hr precip amount QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip1hrQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip1hrQCR'
       print *,'Set precip1hrQCR to -99'
       precip1hrQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip1hrQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip1hrQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip24hrICA "24-hr precip amount IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip24hrICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip24hrICA'
       print *,'Set precip24hrICA to -99'
       precip24hrICA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip24hrICA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip24hrICA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip24hrICR "24-hr precip amount IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip24hrICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip24hrICR'
       print *,'Set precip24hrICR to -99'
       precip24hrICR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip24hrICR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip24hrICR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip24hrQCA "24-hr precip amount QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip24hrQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip24hrQCA'
       print *,'Set precip24hrQCA to -99'
       precip24hrQCA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip24hrQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip24hrQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip24hrQCR "24-hr precip amount QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip24hrQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip24hrQCR'
       print *,'Set precip24hrQCR to -99'
       precip24hrQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip24hrQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip24hrQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip3hrICA  "3-hr precip amount IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip3hrICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip3hrICA'
       print *,'Set precip3hrICA to -99'
       precip3hrICA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip3hrICA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip3hrICA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip3hrICR  "3-hr precip amount IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip3hrICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip3hrICR'
       print *,'Set precip3hrICR to -99'
       precip3hrICR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip3hrICR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip3hrICR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip3hrQCA  "3-hr precip amount QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip3hrQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip3hrQCA'
       print *,'Set precip3hrQCA to -99'
       precip3hrQCA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip3hrQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip3hrQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip3hrQCR  "3-hr precip amount QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip3hrQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip3hrQCR'
       print *,'Set precip3hrQCR to -99'
       precip3hrQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip3hrQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip3hrQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip5minICA "5-min precip amount IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip5minICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip5minICA'
       print *,'Set precip5minICA to -99'
       precip5minICA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip5minICA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip5minICA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip5minICR "5-min precip amount IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip5minICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip5minICR'
       print *,'Set precip5minICR to -99'
       precip5minICR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip5minICR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip5minICR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip5minQCA "5-min precip amount QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip5minQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip5minQCA'
       print *,'Set precip5minQCA to -99'
       precip5minQCA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip5minQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip5minQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip5minQCR "5-min precip amount QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip5minQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip5minQCR'
       print *,'Set precip5minQCR to -99'
       precip5minQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip5minQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip5minQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip6hrICA  "6-hr precip amount IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip6hrICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip6hrICA'
       print *,'Set precip6hrICA to -99'
       precip6hrICA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip6hrICA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip6hrICA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip6hrICR  "6-hr precip amount IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip6hrICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip6hrICR'
       print *,'Set precip6hrICR to -99'
       precip6hrICR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip6hrICR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip6hrICR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip6hrQCA  "6-hr precip amount QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip6hrQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip6hrQCA'
       print *,'Set precip6hrQCA to -99'
       precip6hrQCA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip6hrQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip6hrQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip6hrQCR  "6-hr precip amount QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip6hrQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip6hrQCR'
       print *,'Set precip6hrQCR to -99'
       precip6hrQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precip6hrQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip6hrQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccumICA"precip amount IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipAccumICA'
       print *,'Set precipAccumICA to -99'
       precipAccumICA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipAccumICA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipAccumICA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccumICR"precip amount IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipAccumICR'
       print *,'Set precipAccumICR to -99'
       precipAccumICR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipAccumICR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipAccumICR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccumQCA"precip amount QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipAccumQCA'
       print *,'Set precipAccumQCA to -99'
       precipAccumQCA = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipAccumQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipAccumQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccumQCR"precip amount QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipAccumQCR'
       print *,'Set precipAccumQCR to -99'
       precipAccumQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipAccumQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipAccumQCR'
       endif
      endif

C
C     Variable        NETCDF Long Name
C     prevRecord    
C
      nf_status=NF_INQ_VARID(nf_fid,'prevRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for prevRecord'
       print *,'Set prevRecord to -99'
       prevRecord = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,prevRecord)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for prevRecord'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     secondsStage1_2
C
      nf_status=NF_INQ_VARID(nf_fid,'secondsStage1_2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for secondsStage1_2'
       print *,'Set secondsStage1_2 to -99'
       secondsStage1_2 = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,secondsStage1_2)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for secondsStage1_2'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     secondsStage3 
C
      nf_status=NF_INQ_VARID(nf_fid,'secondsStage3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for secondsStage3'
       print *,'Set secondsStage3 to -99'
       secondsStage3 = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,secondsStage3)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for secondsStage3'
       endif
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C     observationTime"time of observation"
C
      nf_status=NF_INQ_VARID(nf_fid,'observationTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for observationTime'
       print *,'Set observationTime to badflag'
       observationTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,observationTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for observationTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - observationTime'
       else
        call ck_array_dble(observationTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - observationTime'
       else
        call ck_array_dble(observationTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     receivedTime  "time data was received"
C
      nf_status=NF_INQ_VARID(nf_fid,'receivedTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for receivedTime'
       print *,'Set receivedTime to badflag'
       receivedTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,receivedTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for receivedTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - receivedTime'
       else
        call ck_array_dble(receivedTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - receivedTime'
       else
        call ck_array_dble(receivedTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     reportTime    "time data was processed by the provider"
C
      nf_status=NF_INQ_VARID(nf_fid,'reportTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for reportTime'
       print *,'Set reportTime to badflag'
       reportTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reportTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for reportTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - reportTime'
       else
        call ck_array_dble(reportTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - reportTime'
       else
        call ck_array_dble(reportTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     riverReportChangeTime"time of last new river stage/flow rpt"
C
      nf_status=NF_INQ_VARID(nf_fid,'riverReportChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for riverReportChangeTime'
       print *,'Set riverReportChangeTime to badflag'
       riverReportChangeTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,riverReportChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for riverReportChangeTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - riverReportChangeTime'
       else
        call ck_array_dble(riverReportChangeTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - riverReportChangeTime'
       else
        call ck_array_dble(riverReportChangeTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C     ICT           "list of possible IC checks"
C
      nf_status=NF_INQ_VARID(nf_fid,'ICT',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for ICT'
       print *,'Set ICT to " "'
       ICT = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,ICT)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for ICT'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     QCT           "list of possible QC checks"
C
      nf_status=NF_INQ_VARID(nf_fid,'QCT',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for QCT'
       print *,'Set QCT to " "'
       QCT = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,QCT)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for QCT'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dataProvider  "LDAD data provider"
C
      nf_status=NF_INQ_VARID(nf_fid,'dataProvider',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dataProvider'
       print *,'Set dataProvider to " "'
       dataProvider = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,dataProvider)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dataProvider'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     handbook5Id   "handbook5 Id (AFOS or SHEF id)"
C
      nf_status=NF_INQ_VARID(nf_fid,'handbook5Id',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for handbook5Id'
       print *,'Set handbook5Id to " "'
       handbook5Id = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,handbook5Id)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for handbook5Id'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     homeWFO       "home WFO Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'homeWFO',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for homeWFO'
       print *,'Set homeWFO to " "'
       homeWFO = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,homeWFO)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for homeWFO'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip12hrDD  "12-hr precip amount QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip12hrDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip12hrDD'
       print *,'Set precip12hrDD to " "'
       precip12hrDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precip12hrDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip12hrDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip1hrDD   "1-hr precip amount QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip1hrDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip1hrDD'
       print *,'Set precip1hrDD to " "'
       precip1hrDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precip1hrDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip1hrDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip24hrDD  "24-hr precip amount QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip24hrDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip24hrDD'
       print *,'Set precip24hrDD to " "'
       precip24hrDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precip24hrDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip24hrDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip3hrDD   "3-hr precip amount QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip3hrDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip3hrDD'
       print *,'Set precip3hrDD to " "'
       precip3hrDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precip3hrDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip3hrDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip5minDD  "5-min precip amount QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip5minDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip5minDD'
       print *,'Set precip5minDD to " "'
       precip5minDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precip5minDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip5minDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precip6hrDD   "6-hr precip amount QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precip6hrDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precip6hrDD'
       print *,'Set precip6hrDD to " "'
       precip6hrDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precip6hrDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precip6hrDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccumDD "precip amount QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipAccumDD'
       print *,'Set precipAccumDD to " "'
       precipAccumDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precipAccumDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipAccumDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     providerId    "Data Provider station Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'providerId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for providerId'
       print *,'Set providerId to " "'
       providerId = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,providerId)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for providerId'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     rawMessage    "raw text LDAD hydro report"
C
      nf_status=NF_INQ_VARID(nf_fid,'rawMessage',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for rawMessage'
       print *,'Set rawMessage to " "'
       rawMessage = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,rawMessage)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for rawMessage'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     staticIds     
C
      nf_status=NF_INQ_VARID(nf_fid,'staticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for staticIds'
       print *,'Set staticIds to " "'
       staticIds = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,staticIds)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for staticIds'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationId     "alphanumeric station Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationId'
       print *,'Set stationId to " "'
       stationId = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationId)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationId'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationName   "alphanumeric station name"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationName'
       print *,'Set stationName to " "'
       stationName = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationName)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationName'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationType   "LDAD station type"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationType'
       print *,'Set stationType to " "'
       stationType = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationType)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationType'
       endif
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
