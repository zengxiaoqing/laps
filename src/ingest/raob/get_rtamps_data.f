      subroutine get_rtamps_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*170 filename

      integer numRawLevels,
     +     numRTamps, stdLevel, nfid, nf_vid,
     +     nf_status
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nfid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN 20031201_1200'
      endif
C
C  Fill all dimension values
C
C
C Get size of numRawLevels
C
      nf_status = NF_INQ_DIMID(nfid,'numRawLevels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim numRawLevels'
      endif
      nf_status = NF_INQ_DIMLEN(nfid,nf_vid,numRawLevels)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim numRawLevels'
      endif
C
C Get size of numRTamps
C
      nf_status = NF_INQ_DIMID(nfid,'numRTamps',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim numRTamps'
      endif
      nf_status = NF_INQ_DIMLEN(nfid,nf_vid,numRTamps)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim numRTamps'
      endif
C
C Get size of stdLevel
C
      nf_status = NF_INQ_DIMID(nfid,'stdLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim stdLevel'
      endif
      nf_status = NF_INQ_DIMLEN(nfid,nf_vid,stdLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim stdLevel'
      endif
C
C
      call read_rtamps_data(nfid, 
     +     numRawLevels, numRTamps, stdLevel, 
     +     i4time_sys, ilaps_cycle_time, NX_L, NY_L,
     +     i4time_earliest, i4time_latest, lun_out, istatus)

      return
      end
C
C
      subroutine read_rtamps_data(nfid, 
     +     numRawLevels, numRTamps, stdLevel, 
     +     i4time_sys, ilaps_cycle_time, NX_L, NY_L,
     +     i4time_earliest, i4time_latest, lun_out, istatus)


      include 'netcdf.inc'
      integer numRawLevels,
     +     numRTamps, stdLevel, nfid, nf_vid,
     +     nf_status

!     Linda's declarations
      real          geoHt(numRTamps,numRawLevels)
      real          windDir(numRTamps,numRawLevels)
      real          windSpeed(numRTamps,numRawLevels)
      real          temp(numRTamps,numRawLevels)
      real          dewPt(numRTamps,numRawLevels)
      real          baroPress(numRTamps,numRawLevels)
      integer       editFlag(numRTamps,numRawLevels)


!     Declarations for 'write_snd' call
      real stalat(numRawLevels),stalon(numRawLevels)
      character c5_staid*5,a9time_ob(numRawLevels)*9,c8_obstype*8
      real height_m(numRawLevels)
      real pressure_mb(numRawLevels)
      real temp_c(numRawLevels)
      real dewpoint_c(numRawLevels)
      real dir_deg(numRawLevels)
      real spd_mps(numRawLevels)

      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)

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

      call read_rtamps_netcdf(nfid, numRTamps, numRawLevels, 
     1                        missingData, geoHt, windDir,
     1                        windSpeed, temp, dewPt, baroPress, 
     1                        editFlag, istatus)
C
C The netcdf variables are filled - your snd write call may go here
C
      c8_obstype = 'RAOB    '

      do iob = 1,numRTamps
!         call 'write_snd' for a single ob
          call open_ext(lun_out,i4time_sys,'snd',istatus)

          call convert_array(geoHt(iob,:),height_m,numRawLevels
     1                      ,'none',r_missing_data,istatus)       

          call convert_array(baroPress(iob,:),pressure_mb,numRawLevels       
     1                      ,'pa_to_mb',r_missing_data,istatus)       

          call convert_array(temp(iob,:),temp_c,numRawLevels
     1                      ,'k_to_c',r_missing_data,istatus)       

          call convert_array(dewPt(iob,:),dewpoint_c,numRawLevels
     1                      ,'k_to_c',r_missing_data,istatus)       

          call convert_array(windDir(iob,:),dir_deg,numRawLevels
     1                      ,'none',r_missing_data,istatus)       

          call convert_array(windSpeed(iob,:),spd_mps,numRawLevels
     1                      ,'none',r_missing_data,istatus) 

          call write_snd(lun_out                         ! I
     +                  ,1,numRawLevels,1                ! I
     +                  ,iwmostanum                      ! I
     +                  ,stalat,stalon,staelev           ! I
     +                  ,c5_staid,a9time_ob,c8_obstype   ! I
     +                  ,nlvl                            ! I
     +                  ,height_m                        ! I
     +                  ,pressure_mb                     ! I
     +                  ,temp_c                          ! I
     +                  ,dewpoint_c                      ! I
     +                  ,dir_deg                         ! I
     +                  ,spd_mps                         ! I
     +                  ,istatus)                        ! O
      enddo ! iob
      return
      end


