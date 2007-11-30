
      subroutine get_pirep_data_WFO(i4time_sys,ilaps_cycle_time,
     1                              filename,ext,NX_L,NY_L,istatus)

      character*(*) filename
      character*(*) ext

!.............................................................................

      include 'netcdf.inc'
      integer maxNumCloudLayers, recNum,nf_fid, nf_vid, nf_status
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',filename
      endif
C
C  Fill all dimension values
C
C
C Get size of maxNumCloudLayers
C
      nf_status = NF_INQ_DIMID(nf_fid,'maxNumCloudLayers',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxNumCloudLayers'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,maxNumCloudLayers)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxNumCloudLayers'
      endif
C
C Get size of recNum
C
      nf_status = NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif

      call pireps_WFO_sub(nf_fid, maxNumCloudLayers, recNum,
!.............................................................................
     1     ext,i4time_sys,ilaps_cycle_time,NX_L,NY_L,istatus)
      return
!.............................................................................
      end
C
C
      subroutine pireps_WFO_sub(nf_fid, maxNumCloudLayers, recNum,
!.............................................................................
     1     ext,i4time_sys,ilaps_cycle_time,NX_L,NY_L,istatus)
!.............................................................................

      include 'netcdf.inc'
      integer maxNumCloudLayers, recNum,nf_fid, nf_vid, 
     +        nf_status, cloudAmt(maxNumCloudLayers,recNum),
     +        written
      real lat(recNum), lon(recNum), 
     +     cloudBaseHt(maxNumCloudLayers,recNum),
     +     cloudTopHt(maxNumCloudLayers,recNum)
      double precision timeObs(recNum)

!.............................................................................

      integer i4time_sys, ilaps_cycle_time, NX_L, NY_L, istatus
      character*(*) ext
      character*9 a9_timeObs 
      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)

!.............................................................................
C
C      read latitude          "Latitude of report"
C
        nf_status = NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,lat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
C
C      read longitude          "Longitude of report"
C
        nf_status = NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,lon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
C
C      read timeObs      "Time of report"
C
        nf_status = NF_INQ_VARID(nf_fid,'timeObs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeObs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif
C
C      read cloudBaseHeight
C
        nf_status = NF_INQ_VARID(nf_fid,'cloudBaseHeight',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseHeight'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,cloudBaseHt)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseHeight'
      endif
C
C      read cloudTopHeight
C
        nf_status = NF_INQ_VARID(nf_fid,'cloudTopHeight',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudTopHeight'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,cloudTopHt)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudTopHeight'
      endif
C
C      read cloudAmount
C
        nf_status = NF_INQ_VARID(nf_fid,'cloudAmount',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudAmount'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,cloudAmt)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudAmount'
      endif
C
C The netcdf variables are filled 
C
!.............................................................................

      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a,
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_domain_perimeter'
          return
      endif

!     Write All Pireps to LAPS PIN file
      altitude = 0.

      r_nc_missing_data = 1e20

      write(6,*)' recNum = ',recNum

      num_pireps = recNum

      do i = 1,num_pireps

          write(6,*)
          write(6,*)' Pirep #',i

          if(lat(i) .ge. r_nc_missing_data)then
              write(6,*)' Missing first latitude',i
              goto 999
          endif
          if(lon(i) .ge. r_nc_missing_data)then
              write(6,*)' Missing first longitude',i
              goto 999
          endif
          if(timeObs(i) .ge. r_nc_missing_data)then
              write(6,*)' Missing ob time',i
              goto 999
          endif

!         Test to see how many lat/lons are present
          n_latitude_present = 0
          if(lat(i) .lt. r_nc_missing_data)then
            write(6,*)' lat/lon = ',lat(i),lon(i)
            n_latitude_present = n_latitude_present + 1
          endif

          write(6,*)' Num locations = ',n_latitude_present       

          if(n_latitude_present .gt. 1)then
              write(6,*)' Multiple locations, reject ob',i
              goto 999
          endif

          if(lat(i) .le. rnorth .and. lat(i) .ge. south .and.
     1       lon(i) .ge. west   .and. lon(i) .le. east      )then        
              write(6,*)' Pirep is inside domain lat/lon perimeter'
          else
              write(6,*)
     1            ' Outside domain lat/lon perimeter - reject'
              goto 999
          endif

!         Write Individual Pirep to LAPS PIN file
          call c_time2fname(nint(timeObs(i)),a9_timeObs)

          call cv_asc_i4time(a9_timeObs,i4time_ob)
          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. (ilaps_cycle_time / 2) )then
              write(6,*)' Outside time window - reject '
     1                              ,a9_timeObs,i4_resid
              goto 999        
          endif


!         Write out cloud base/top in feet and cloud amount in eighths
!         written = 0 if timeObs/lat/lon not written out, set to 1 once written to file
          written = 0
          do ilyr = 1,maxNumCloudLayers

!             Use -1000. for missing value of cloud base/top
              if ((cloudBaseHt(ilyr,i) .ge. 1e10) .or.
     +         (cloudBaseHt(ilyr,i) .eq. -9999.)) then
                  rbase = -1000.
              else
                  rbase = cloudBaseHt(ilyr,i)
              endif

              if((cloudTopHt(ilyr,i) .ge. 1e10) .or.
     +         (cloudTopHt(ilyr,i) .eq. -9999.)) then
                  rtop = -1000.
              else
                  rtop = cloudTopHt(ilyr,i)
              endif

!             Test for missing or unusable cloud cover
!             WFO bases cloudAmount on WFM BUFR Table 020011
              if((cloudAmt(ilyr,i) .gt. 8) .or.
     +           (cloudAmt(ilyr,i) .lt. 0)) then 
                  ieighths = -999
              else 
                  ieighths = cloudAmt(ilyr,i)
              endif                  

              write(6,3)rbase,rtop,ieighths
 3            format(' Cloud layer',2f8.0,i5)

              if ((rbase .ne. -1000.) .and. (rtop .ne. -1000.) .and.
     +            (ieighths .ne. -999)) then
                if (written .eq. 0) then   ! write timeObs, lat, lon, etc
                  call open_ext(31,i4time_sys,ext(1:3),istatus)

                  write(6,1)a9_timeObs 
                  write(31,1)a9_timeObs, a9_timeObs 
 1                format(' Time - prp/rcvd:'/1x,a9,2x,a9) 

                  write(6,2)lat(i),lon(i),altitude
                  write(31,2)lat(i),lon(i),altitude
 2                format(' Lat, lon, altitude'/f8.3,f10.3,f8.0)  

                  write(6,33)
                  write(31,33)
 33               format(' Cloud layer')
                  written = 1
                endif
!               write(6,*)' Above layer written to PIN file'
                write(31,3)rbase,rtop,ieighths
              endif

          enddo ! ilyr

          if (written .eq. 0) then ! none of layers in pirep contained info
            write(6,*)' Above pirep NOT written to PIN file - ',
     +'no usable data.'
          endif

 999      continue

      enddo ! i

!.............................................................................

      return
      end
