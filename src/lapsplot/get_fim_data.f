      subroutine get_fim_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer latitude, longitude, time,nf_fid, nf_vid, nf_status
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
C Get size of latitude
C
      nf_status=NF_INQ_DIMID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim latitude'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim latitude'
      endif
C
C Get size of longitude
C
      nf_status=NF_INQ_DIMID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim longitude'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim longitude'
      endif
C
C Get size of time
C
      nf_status=NF_INQ_DIMID(nf_fid,'time',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim time'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,time)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim time'
      endif
      call read_fim_data(nf_fid, latitude, longitude, time,
     +     i4time_sys, ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, istatus)

      return
      end
C
C
      subroutine read_fim_data(nf_fid, latitude, longitude, time,
     +     i4time_sys, ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, istatus)


      include 'netcdf.inc'
      integer latitude, longitude, time,nf_fid, nf_vid, nf_status

      real CLWMR_10hybridlevel( longitude,  latitude, time),
     +     CLWMR_11hybridlevel( longitude,  latitude, time),
     +     CLWMR_12hybridlevel( longitude,  latitude, time),
     +     CLWMR_13hybridlevel( longitude,  latitude, time),
     +     CLWMR_14hybridlevel( longitude,  latitude, time),
     +     CLWMR_15hybridlevel( longitude,  latitude, time),
     +     CLWMR_16hybridlevel( longitude,  latitude, time),
     +     CLWMR_17hybridlevel( longitude,  latitude, time),
     +     CLWMR_18hybridlevel( longitude,  latitude, time),
     +     CLWMR_19hybridlevel( longitude,  latitude, time),
     +     CLWMR_1hybridlevel( longitude,  latitude, time),
     +     CLWMR_20hybridlevel( longitude,  latitude, time),
     +     CLWMR_21hybridlevel( longitude,  latitude, time),
     +     CLWMR_22hybridlevel( longitude,  latitude, time),
     +     CLWMR_23hybridlevel( longitude,  latitude, time),
     +     CLWMR_24hybridlevel( longitude,  latitude, time),
     +     CLWMR_25hybridlevel( longitude,  latitude, time),
     +     CLWMR_26hybridlevel( longitude,  latitude, time),
     +     CLWMR_27hybridlevel( longitude,  latitude, time),
     +     CLWMR_28hybridlevel( longitude,  latitude, time),
     +     CLWMR_29hybridlevel( longitude,  latitude, time),
     +     CLWMR_2hybridlevel( longitude,  latitude, time),
     +     CLWMR_30hybridlevel( longitude,  latitude, time),
     +     CLWMR_31hybridlevel( longitude,  latitude, time),
     +     CLWMR_32hybridlevel( longitude,  latitude, time),
     +     CLWMR_33hybridlevel( longitude,  latitude, time),
     +     CLWMR_34hybridlevel( longitude,  latitude, time),
     +     CLWMR_35hybridlevel( longitude,  latitude, time),
     +     CLWMR_36hybridlevel( longitude,  latitude, time),
     +     CLWMR_37hybridlevel( longitude,  latitude, time),
     +     CLWMR_38hybridlevel( longitude,  latitude, time),
     +     CLWMR_39hybridlevel( longitude,  latitude, time),
     +     CLWMR_3hybridlevel( longitude,  latitude, time),
     +     CLWMR_40hybridlevel( longitude,  latitude, time),
     +     CLWMR_41hybridlevel( longitude,  latitude, time),
     +     CLWMR_42hybridlevel( longitude,  latitude, time),
     +     CLWMR_43hybridlevel( longitude,  latitude, time),
     +     CLWMR_44hybridlevel( longitude,  latitude, time),
     +     CLWMR_45hybridlevel( longitude,  latitude, time),
     +     CLWMR_46hybridlevel( longitude,  latitude, time),
     +     CLWMR_47hybridlevel( longitude,  latitude, time),
     +     CLWMR_48hybridlevel( longitude,  latitude, time),
     +     CLWMR_49hybridlevel( longitude,  latitude, time),
     +     CLWMR_4hybridlevel( longitude,  latitude, time),
     +     CLWMR_50hybridlevel( longitude,  latitude, time),
     +     CLWMR_51hybridlevel( longitude,  latitude, time),
     +     CLWMR_52hybridlevel( longitude,  latitude, time),
     +     CLWMR_53hybridlevel( longitude,  latitude, time),
     +     CLWMR_54hybridlevel( longitude,  latitude, time),
     +     CLWMR_55hybridlevel( longitude,  latitude, time),
     +     CLWMR_56hybridlevel( longitude,  latitude, time),
     +     CLWMR_57hybridlevel( longitude,  latitude, time),
     +     CLWMR_58hybridlevel( longitude,  latitude, time),
     +     CLWMR_59hybridlevel( longitude,  latitude, time),
     +     CLWMR_5hybridlevel( longitude,  latitude, time),
     +     CLWMR_60hybridlevel( longitude,  latitude, time),
     +     CLWMR_61hybridlevel( longitude,  latitude, time),
     +     CLWMR_62hybridlevel( longitude,  latitude, time),
     +     CLWMR_63hybridlevel( longitude,  latitude, time),
     +     CLWMR_64hybridlevel( longitude,  latitude, time),
     +     CLWMR_6hybridlevel( longitude,  latitude, time),
     +     CLWMR_7hybridlevel( longitude,  latitude, time),
     +     CLWMR_8hybridlevel( longitude,  latitude, time),
     +     CLWMR_9hybridlevel( longitude,  latitude, time),
     +     PRES_10hybridlevel( longitude,  latitude, time),
     +     PRES_11hybridlevel( longitude,  latitude, time),
     +     PRES_12hybridlevel( longitude,  latitude, time),
     +     PRES_13hybridlevel( longitude,  latitude, time),
     +     PRES_14hybridlevel( longitude,  latitude, time),
     +     PRES_15hybridlevel( longitude,  latitude, time),
     +     PRES_16hybridlevel( longitude,  latitude, time),
     +     PRES_17hybridlevel( longitude,  latitude, time),
     +     PRES_18hybridlevel( longitude,  latitude, time),
     +     PRES_19hybridlevel( longitude,  latitude, time),
     +     PRES_1hybridlevel( longitude,  latitude, time),
     +     PRES_20hybridlevel( longitude,  latitude, time),
     +     PRES_21hybridlevel( longitude,  latitude, time),
     +     PRES_22hybridlevel( longitude,  latitude, time),
     +     PRES_23hybridlevel( longitude,  latitude, time),
     +     PRES_24hybridlevel( longitude,  latitude, time),
     +     PRES_25hybridlevel( longitude,  latitude, time),
     +     PRES_26hybridlevel( longitude,  latitude, time),
     +     PRES_27hybridlevel( longitude,  latitude, time),
     +     PRES_28hybridlevel( longitude,  latitude, time),
     +     PRES_29hybridlevel( longitude,  latitude, time),
     +     PRES_2hybridlevel( longitude,  latitude, time),
     +     PRES_30hybridlevel( longitude,  latitude, time),
     +     PRES_31hybridlevel( longitude,  latitude, time),
     +     PRES_32hybridlevel( longitude,  latitude, time),
     +     PRES_33hybridlevel( longitude,  latitude, time),
     +     PRES_34hybridlevel( longitude,  latitude, time),
     +     PRES_35hybridlevel( longitude,  latitude, time),
     +     PRES_36hybridlevel( longitude,  latitude, time),
     +     PRES_37hybridlevel( longitude,  latitude, time),
     +     PRES_38hybridlevel( longitude,  latitude, time),
     +     PRES_39hybridlevel( longitude,  latitude, time),
     +     PRES_3hybridlevel( longitude,  latitude, time),
     +     PRES_40hybridlevel( longitude,  latitude, time),
     +     PRES_41hybridlevel( longitude,  latitude, time),
     +     PRES_42hybridlevel( longitude,  latitude, time),
     +     PRES_43hybridlevel( longitude,  latitude, time),
     +     PRES_44hybridlevel( longitude,  latitude, time),
     +     PRES_45hybridlevel( longitude,  latitude, time),
     +     PRES_46hybridlevel( longitude,  latitude, time),
     +     PRES_47hybridlevel( longitude,  latitude, time),
     +     PRES_48hybridlevel( longitude,  latitude, time),
     +     PRES_49hybridlevel( longitude,  latitude, time),
     +     PRES_4hybridlevel( longitude,  latitude, time),
     +     PRES_50hybridlevel( longitude,  latitude, time),
     +     PRES_51hybridlevel( longitude,  latitude, time),
     +     PRES_52hybridlevel( longitude,  latitude, time),
     +     PRES_53hybridlevel( longitude,  latitude, time),
     +     PRES_54hybridlevel( longitude,  latitude, time),
     +     PRES_55hybridlevel( longitude,  latitude, time),
     +     PRES_56hybridlevel( longitude,  latitude, time),
     +     PRES_57hybridlevel( longitude,  latitude, time),
     +     PRES_58hybridlevel( longitude,  latitude, time),
     +     PRES_59hybridlevel( longitude,  latitude, time),
     +     PRES_5hybridlevel( longitude,  latitude, time),
     +     PRES_60hybridlevel( longitude,  latitude, time),
     +     PRES_61hybridlevel( longitude,  latitude, time),
     +     PRES_62hybridlevel( longitude,  latitude, time),
     +     PRES_63hybridlevel( longitude,  latitude, time),
     +     PRES_64hybridlevel( longitude,  latitude, time),
     +     PRES_65hybridlevel( longitude,  latitude, time),
     +     PRES_6hybridlevel( longitude,  latitude, time),
     +     PRES_7hybridlevel( longitude,  latitude, time),
     +     PRES_8hybridlevel( longitude,  latitude, time),
     +     PRES_9hybridlevel( longitude,  latitude, time)
!     double precision latitude(latitude), longitude(longitude),
!    +     time(time)

!     Declarations for 'write_zzz' call
      integer iwmostanum(recNum)
      character a9time_ob_r(recNum)*9
      logical l_closest_time, l_closest_time_i, l_in_domain
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

      call read_fim_netcdf(nf_fid, ! latitude, longitude, time, 
     +     CLWMR_10hybridlevel, CLWMR_11hybridlevel, 
     +     CLWMR_12hybridlevel, CLWMR_13hybridlevel, 
     +     CLWMR_14hybridlevel, CLWMR_15hybridlevel, 
     +     CLWMR_16hybridlevel, CLWMR_17hybridlevel, 
     +     CLWMR_18hybridlevel, CLWMR_19hybridlevel, 
     +     CLWMR_1hybridlevel, CLWMR_20hybridlevel, 
     +     CLWMR_21hybridlevel, CLWMR_22hybridlevel, 
     +     CLWMR_23hybridlevel, CLWMR_24hybridlevel, 
     +     CLWMR_25hybridlevel, CLWMR_26hybridlevel, 
     +     CLWMR_27hybridlevel, CLWMR_28hybridlevel, 
     +     CLWMR_29hybridlevel, CLWMR_2hybridlevel, 
     +     CLWMR_30hybridlevel, CLWMR_31hybridlevel, 
     +     CLWMR_32hybridlevel, CLWMR_33hybridlevel, 
     +     CLWMR_34hybridlevel, CLWMR_35hybridlevel, 
     +     CLWMR_36hybridlevel, CLWMR_37hybridlevel, 
     +     CLWMR_38hybridlevel, CLWMR_39hybridlevel, 
     +     CLWMR_3hybridlevel, CLWMR_40hybridlevel, 
     +     CLWMR_41hybridlevel, CLWMR_42hybridlevel, 
     +     CLWMR_43hybridlevel, CLWMR_44hybridlevel, 
     +     CLWMR_45hybridlevel, CLWMR_46hybridlevel, 
     +     CLWMR_47hybridlevel, CLWMR_48hybridlevel, 
     +     CLWMR_49hybridlevel, CLWMR_4hybridlevel, 
     +     CLWMR_50hybridlevel, CLWMR_51hybridlevel, 
     +     CLWMR_52hybridlevel, CLWMR_53hybridlevel, 
     +     CLWMR_54hybridlevel, CLWMR_55hybridlevel, 
     +     CLWMR_56hybridlevel, CLWMR_57hybridlevel, 
     +     CLWMR_58hybridlevel, CLWMR_59hybridlevel, 
     +     CLWMR_5hybridlevel, CLWMR_60hybridlevel, 
     +     CLWMR_61hybridlevel, CLWMR_62hybridlevel, 
     +     CLWMR_63hybridlevel, CLWMR_64hybridlevel, 
     +     CLWMR_6hybridlevel, CLWMR_7hybridlevel, 
     +     CLWMR_8hybridlevel, CLWMR_9hybridlevel, 
     +     PRES_10hybridlevel, PRES_11hybridlevel, 
     +     PRES_12hybridlevel, PRES_13hybridlevel, 
     +     PRES_14hybridlevel, PRES_15hybridlevel, 
     +     PRES_16hybridlevel, PRES_17hybridlevel, 
     +     PRES_18hybridlevel, PRES_19hybridlevel, PRES_1hybridlevel, 
     +     PRES_20hybridlevel, PRES_21hybridlevel, 
     +     PRES_22hybridlevel, PRES_23hybridlevel, 
     +     PRES_24hybridlevel, PRES_25hybridlevel, 
     +     PRES_26hybridlevel, PRES_27hybridlevel, 
     +     PRES_28hybridlevel, PRES_29hybridlevel, PRES_2hybridlevel, 
     +     PRES_30hybridlevel, PRES_31hybridlevel, 
     +     PRES_32hybridlevel, PRES_33hybridlevel, 
     +     PRES_34hybridlevel, PRES_35hybridlevel, 
     +     PRES_36hybridlevel, PRES_37hybridlevel, 
     +     PRES_38hybridlevel, PRES_39hybridlevel, PRES_3hybridlevel, 
     +     PRES_40hybridlevel, PRES_41hybridlevel, 
     +     PRES_42hybridlevel, PRES_43hybridlevel, 
     +     PRES_44hybridlevel, PRES_45hybridlevel, 
     +     PRES_46hybridlevel, PRES_47hybridlevel, 
     +     PRES_48hybridlevel, PRES_49hybridlevel, PRES_4hybridlevel, 
     +     PRES_50hybridlevel, PRES_51hybridlevel, 
     +     PRES_52hybridlevel, PRES_53hybridlevel, 
     +     PRES_54hybridlevel, PRES_55hybridlevel, 
     +     PRES_56hybridlevel, PRES_57hybridlevel, 
     +     PRES_58hybridlevel, PRES_59hybridlevel, PRES_5hybridlevel, 
     +     PRES_60hybridlevel, PRES_61hybridlevel, 
     +     PRES_62hybridlevel, PRES_63hybridlevel, 
     +     PRES_64hybridlevel, PRES_65hybridlevel, PRES_6hybridlevel, 
     +     PRES_7hybridlevel, PRES_8hybridlevel, PRES_9hybridlevel, 
     +     latitude, longitude, time)
C
C The netcdf variables are filled - your zzz write call may go here
C
!     Initial loop through obs to get times and stanums
      do iob = 1,recNum
          iwmostanum(iob) = 0
          if(abs(observationTime(iob)) .le. 1e10)then
              i4time_ob = idint(observationTime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

      do iob = 1,recNum
          height_m = r_missing_data
          l_closest_time = .true.

      enddo ! iob
      return
      end
C
C  Subroutine to read the file 
C
      subroutine read_fim_netcdf(nf_fid, ! latitude, longitude, time, 
     +     CLWMR_10hybridlevel, CLWMR_11hybridlevel, 
     +     CLWMR_12hybridlevel, CLWMR_13hybridlevel, 
     +     CLWMR_14hybridlevel, CLWMR_15hybridlevel, 
     +     CLWMR_16hybridlevel, CLWMR_17hybridlevel, 
     +     CLWMR_18hybridlevel, CLWMR_19hybridlevel, 
     +     CLWMR_1hybridlevel, CLWMR_20hybridlevel, 
     +     CLWMR_21hybridlevel, CLWMR_22hybridlevel, 
     +     CLWMR_23hybridlevel, CLWMR_24hybridlevel, 
     +     CLWMR_25hybridlevel, CLWMR_26hybridlevel, 
     +     CLWMR_27hybridlevel, CLWMR_28hybridlevel, 
     +     CLWMR_29hybridlevel, CLWMR_2hybridlevel, 
     +     CLWMR_30hybridlevel, CLWMR_31hybridlevel, 
     +     CLWMR_32hybridlevel, CLWMR_33hybridlevel, 
     +     CLWMR_34hybridlevel, CLWMR_35hybridlevel, 
     +     CLWMR_36hybridlevel, CLWMR_37hybridlevel, 
     +     CLWMR_38hybridlevel, CLWMR_39hybridlevel, 
     +     CLWMR_3hybridlevel, CLWMR_40hybridlevel, 
     +     CLWMR_41hybridlevel, CLWMR_42hybridlevel, 
     +     CLWMR_43hybridlevel, CLWMR_44hybridlevel, 
     +     CLWMR_45hybridlevel, CLWMR_46hybridlevel, 
     +     CLWMR_47hybridlevel, CLWMR_48hybridlevel, 
     +     CLWMR_49hybridlevel, CLWMR_4hybridlevel, 
     +     CLWMR_50hybridlevel, CLWMR_51hybridlevel, 
     +     CLWMR_52hybridlevel, CLWMR_53hybridlevel, 
     +     CLWMR_54hybridlevel, CLWMR_55hybridlevel, 
     +     CLWMR_56hybridlevel, CLWMR_57hybridlevel, 
     +     CLWMR_58hybridlevel, CLWMR_59hybridlevel, 
     +     CLWMR_5hybridlevel, CLWMR_60hybridlevel, 
     +     CLWMR_61hybridlevel, CLWMR_62hybridlevel, 
     +     CLWMR_63hybridlevel, CLWMR_64hybridlevel, 
     +     CLWMR_6hybridlevel, CLWMR_7hybridlevel, 
     +     CLWMR_8hybridlevel, CLWMR_9hybridlevel, 
     +     PRES_10hybridlevel, PRES_11hybridlevel, 
     +     PRES_12hybridlevel, PRES_13hybridlevel, 
     +     PRES_14hybridlevel, PRES_15hybridlevel, 
     +     PRES_16hybridlevel, PRES_17hybridlevel, 
     +     PRES_18hybridlevel, PRES_19hybridlevel, PRES_1hybridlevel, 
     +     PRES_20hybridlevel, PRES_21hybridlevel, 
     +     PRES_22hybridlevel, PRES_23hybridlevel, 
     +     PRES_24hybridlevel, PRES_25hybridlevel, 
     +     PRES_26hybridlevel, PRES_27hybridlevel, 
     +     PRES_28hybridlevel, PRES_29hybridlevel, PRES_2hybridlevel, 
     +     PRES_30hybridlevel, PRES_31hybridlevel, 
     +     PRES_32hybridlevel, PRES_33hybridlevel, 
     +     PRES_34hybridlevel, PRES_35hybridlevel, 
     +     PRES_36hybridlevel, PRES_37hybridlevel, 
     +     PRES_38hybridlevel, PRES_39hybridlevel, PRES_3hybridlevel, 
     +     PRES_40hybridlevel, PRES_41hybridlevel, 
     +     PRES_42hybridlevel, PRES_43hybridlevel, 
     +     PRES_44hybridlevel, PRES_45hybridlevel, 
     +     PRES_46hybridlevel, PRES_47hybridlevel, 
     +     PRES_48hybridlevel, PRES_49hybridlevel, PRES_4hybridlevel, 
     +     PRES_50hybridlevel, PRES_51hybridlevel, 
     +     PRES_52hybridlevel, PRES_53hybridlevel, 
     +     PRES_54hybridlevel, PRES_55hybridlevel, 
     +     PRES_56hybridlevel, PRES_57hybridlevel, 
     +     PRES_58hybridlevel, PRES_59hybridlevel, PRES_5hybridlevel, 
     +     PRES_60hybridlevel, PRES_61hybridlevel, 
     +     PRES_62hybridlevel, PRES_63hybridlevel, 
     +     PRES_64hybridlevel, PRES_65hybridlevel, PRES_6hybridlevel, 
     +     PRES_7hybridlevel, PRES_8hybridlevel, PRES_9hybridlevel, 
     +     latitude, longitude, time)
C
      include 'netcdf.inc'
      integer latitude, longitude, time,nf_fid, nf_vid, nf_status

      real CLWMR_10hybridlevel( longitude,  latitude, time),
     +     CLWMR_11hybridlevel( longitude,  latitude, time),
     +     CLWMR_12hybridlevel( longitude,  latitude, time),
     +     CLWMR_13hybridlevel( longitude,  latitude, time),
     +     CLWMR_14hybridlevel( longitude,  latitude, time),
     +     CLWMR_15hybridlevel( longitude,  latitude, time),
     +     CLWMR_16hybridlevel( longitude,  latitude, time),
     +     CLWMR_17hybridlevel( longitude,  latitude, time),
     +     CLWMR_18hybridlevel( longitude,  latitude, time),
     +     CLWMR_19hybridlevel( longitude,  latitude, time),
     +     CLWMR_1hybridlevel( longitude,  latitude, time),
     +     CLWMR_20hybridlevel( longitude,  latitude, time),
     +     CLWMR_21hybridlevel( longitude,  latitude, time),
     +     CLWMR_22hybridlevel( longitude,  latitude, time),
     +     CLWMR_23hybridlevel( longitude,  latitude, time),
     +     CLWMR_24hybridlevel( longitude,  latitude, time),
     +     CLWMR_25hybridlevel( longitude,  latitude, time),
     +     CLWMR_26hybridlevel( longitude,  latitude, time),
     +     CLWMR_27hybridlevel( longitude,  latitude, time),
     +     CLWMR_28hybridlevel( longitude,  latitude, time),
     +     CLWMR_29hybridlevel( longitude,  latitude, time),
     +     CLWMR_2hybridlevel( longitude,  latitude, time),
     +     CLWMR_30hybridlevel( longitude,  latitude, time),
     +     CLWMR_31hybridlevel( longitude,  latitude, time),
     +     CLWMR_32hybridlevel( longitude,  latitude, time),
     +     CLWMR_33hybridlevel( longitude,  latitude, time),
     +     CLWMR_34hybridlevel( longitude,  latitude, time),
     +     CLWMR_35hybridlevel( longitude,  latitude, time),
     +     CLWMR_36hybridlevel( longitude,  latitude, time),
     +     CLWMR_37hybridlevel( longitude,  latitude, time),
     +     CLWMR_38hybridlevel( longitude,  latitude, time),
     +     CLWMR_39hybridlevel( longitude,  latitude, time),
     +     CLWMR_3hybridlevel( longitude,  latitude, time),
     +     CLWMR_40hybridlevel( longitude,  latitude, time),
     +     CLWMR_41hybridlevel( longitude,  latitude, time),
     +     CLWMR_42hybridlevel( longitude,  latitude, time),
     +     CLWMR_43hybridlevel( longitude,  latitude, time),
     +     CLWMR_44hybridlevel( longitude,  latitude, time),
     +     CLWMR_45hybridlevel( longitude,  latitude, time),
     +     CLWMR_46hybridlevel( longitude,  latitude, time),
     +     CLWMR_47hybridlevel( longitude,  latitude, time),
     +     CLWMR_48hybridlevel( longitude,  latitude, time),
     +     CLWMR_49hybridlevel( longitude,  latitude, time),
     +     CLWMR_4hybridlevel( longitude,  latitude, time),
     +     CLWMR_50hybridlevel( longitude,  latitude, time),
     +     CLWMR_51hybridlevel( longitude,  latitude, time),
     +     CLWMR_52hybridlevel( longitude,  latitude, time),
     +     CLWMR_53hybridlevel( longitude,  latitude, time),
     +     CLWMR_54hybridlevel( longitude,  latitude, time),
     +     CLWMR_55hybridlevel( longitude,  latitude, time),
     +     CLWMR_56hybridlevel( longitude,  latitude, time),
     +     CLWMR_57hybridlevel( longitude,  latitude, time),
     +     CLWMR_58hybridlevel( longitude,  latitude, time),
     +     CLWMR_59hybridlevel( longitude,  latitude, time),
     +     CLWMR_5hybridlevel( longitude,  latitude, time),
     +     CLWMR_60hybridlevel( longitude,  latitude, time),
     +     CLWMR_61hybridlevel( longitude,  latitude, time),
     +     CLWMR_62hybridlevel( longitude,  latitude, time),
     +     CLWMR_63hybridlevel( longitude,  latitude, time),
     +     CLWMR_64hybridlevel( longitude,  latitude, time),
     +     CLWMR_6hybridlevel( longitude,  latitude, time),
     +     CLWMR_7hybridlevel( longitude,  latitude, time),
     +     CLWMR_8hybridlevel( longitude,  latitude, time),
     +     CLWMR_9hybridlevel( longitude,  latitude, time),
     +     PRES_10hybridlevel( longitude,  latitude, time),
     +     PRES_11hybridlevel( longitude,  latitude, time),
     +     PRES_12hybridlevel( longitude,  latitude, time),
     +     PRES_13hybridlevel( longitude,  latitude, time),
     +     PRES_14hybridlevel( longitude,  latitude, time),
     +     PRES_15hybridlevel( longitude,  latitude, time),
     +     PRES_16hybridlevel( longitude,  latitude, time),
     +     PRES_17hybridlevel( longitude,  latitude, time),
     +     PRES_18hybridlevel( longitude,  latitude, time),
     +     PRES_19hybridlevel( longitude,  latitude, time),
     +     PRES_1hybridlevel( longitude,  latitude, time),
     +     PRES_20hybridlevel( longitude,  latitude, time),
     +     PRES_21hybridlevel( longitude,  latitude, time),
     +     PRES_22hybridlevel( longitude,  latitude, time),
     +     PRES_23hybridlevel( longitude,  latitude, time),
     +     PRES_24hybridlevel( longitude,  latitude, time),
     +     PRES_25hybridlevel( longitude,  latitude, time),
     +     PRES_26hybridlevel( longitude,  latitude, time),
     +     PRES_27hybridlevel( longitude,  latitude, time),
     +     PRES_28hybridlevel( longitude,  latitude, time),
     +     PRES_29hybridlevel( longitude,  latitude, time),
     +     PRES_2hybridlevel( longitude,  latitude, time),
     +     PRES_30hybridlevel( longitude,  latitude, time),
     +     PRES_31hybridlevel( longitude,  latitude, time),
     +     PRES_32hybridlevel( longitude,  latitude, time),
     +     PRES_33hybridlevel( longitude,  latitude, time),
     +     PRES_34hybridlevel( longitude,  latitude, time),
     +     PRES_35hybridlevel( longitude,  latitude, time),
     +     PRES_36hybridlevel( longitude,  latitude, time),
     +     PRES_37hybridlevel( longitude,  latitude, time),
     +     PRES_38hybridlevel( longitude,  latitude, time),
     +     PRES_39hybridlevel( longitude,  latitude, time),
     +     PRES_3hybridlevel( longitude,  latitude, time),
     +     PRES_40hybridlevel( longitude,  latitude, time),
     +     PRES_41hybridlevel( longitude,  latitude, time),
     +     PRES_42hybridlevel( longitude,  latitude, time),
     +     PRES_43hybridlevel( longitude,  latitude, time),
     +     PRES_44hybridlevel( longitude,  latitude, time),
     +     PRES_45hybridlevel( longitude,  latitude, time),
     +     PRES_46hybridlevel( longitude,  latitude, time),
     +     PRES_47hybridlevel( longitude,  latitude, time),
     +     PRES_48hybridlevel( longitude,  latitude, time),
     +     PRES_49hybridlevel( longitude,  latitude, time),
     +     PRES_4hybridlevel( longitude,  latitude, time),
     +     PRES_50hybridlevel( longitude,  latitude, time),
     +     PRES_51hybridlevel( longitude,  latitude, time),
     +     PRES_52hybridlevel( longitude,  latitude, time),
     +     PRES_53hybridlevel( longitude,  latitude, time),
     +     PRES_54hybridlevel( longitude,  latitude, time),
     +     PRES_55hybridlevel( longitude,  latitude, time),
     +     PRES_56hybridlevel( longitude,  latitude, time),
     +     PRES_57hybridlevel( longitude,  latitude, time),
     +     PRES_58hybridlevel( longitude,  latitude, time),
     +     PRES_59hybridlevel( longitude,  latitude, time),
     +     PRES_5hybridlevel( longitude,  latitude, time),
     +     PRES_60hybridlevel( longitude,  latitude, time),
     +     PRES_61hybridlevel( longitude,  latitude, time),
     +     PRES_62hybridlevel( longitude,  latitude, time),
     +     PRES_63hybridlevel( longitude,  latitude, time),
     +     PRES_64hybridlevel( longitude,  latitude, time),
     +     PRES_65hybridlevel( longitude,  latitude, time),
     +     PRES_6hybridlevel( longitude,  latitude, time),
     +     PRES_7hybridlevel( longitude,  latitude, time),
     +     PRES_8hybridlevel( longitude,  latitude, time),
     +     PRES_9hybridlevel( longitude,  latitude, time)
!     double precision latitude(latitude), longitude(longitude),
!    +     time(time)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     CLWMR_10hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_10hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_10hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_10hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_10hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_11hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_11hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_11hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_11hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_11hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_12hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_12hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_12hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_12hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_12hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_13hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_13hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_13hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_13hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_13hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_14hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_14hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_14hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_14hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_14hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_15hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_15hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_15hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_15hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_15hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_16hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_16hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_16hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_16hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_16hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_17hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_17hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_17hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_17hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_17hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_18hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_18hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_18hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_18hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_18hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_19hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_19hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_19hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_19hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_19hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_1hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_1hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_1hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_1hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_1hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_20hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_20hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_20hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_20hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_20hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_21hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_21hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_21hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_21hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_21hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_22hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_22hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_22hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_22hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_22hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_23hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_23hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_23hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_23hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_23hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_24hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_24hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_24hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_24hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_24hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_25hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_25hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_25hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_25hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_25hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_26hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_26hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_26hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_26hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_26hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_27hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_27hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_27hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_27hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_27hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_28hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_28hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_28hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_28hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_28hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_29hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_29hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_29hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_29hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_29hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_2hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_2hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_2hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_2hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_2hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_30hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_30hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_30hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_30hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_30hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_31hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_31hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_31hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_31hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_31hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_32hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_32hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_32hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_32hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_32hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_33hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_33hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_33hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_33hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_33hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_34hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_34hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_34hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_34hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_34hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_35hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_35hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_35hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_35hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_35hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_36hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_36hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_36hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_36hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_36hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_37hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_37hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_37hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_37hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_37hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_38hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_38hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_38hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_38hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_38hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_39hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_39hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_39hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_39hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_39hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_3hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_3hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_3hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_3hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_3hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_40hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_40hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_40hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_40hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_40hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_41hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_41hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_41hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_41hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_41hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_42hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_42hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_42hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_42hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_42hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_43hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_43hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_43hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_43hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_43hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_44hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_44hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_44hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_44hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_44hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_45hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_45hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_45hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_45hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_45hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_46hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_46hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_46hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_46hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_46hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_47hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_47hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_47hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_47hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_47hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_48hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_48hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_48hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_48hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_48hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_49hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_49hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_49hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_49hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_49hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_4hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_4hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_4hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_4hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_4hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_50hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_50hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_50hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_50hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_50hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_51hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_51hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_51hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_51hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_51hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_52hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_52hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_52hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_52hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_52hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_53hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_53hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_53hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_53hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_53hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_54hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_54hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_54hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_54hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_54hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_55hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_55hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_55hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_55hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_55hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_56hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_56hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_56hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_56hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_56hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_57hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_57hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_57hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_57hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_57hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_58hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_58hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_58hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_58hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_58hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_59hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_59hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_59hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_59hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_59hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_5hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_5hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_5hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_5hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_5hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_60hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_60hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_60hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_60hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_60hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_61hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_61hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_61hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_61hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_61hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_62hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_62hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_62hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_62hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_62hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_63hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_63hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_63hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_63hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_63hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_64hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_64hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_64hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_64hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_64hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_6hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_6hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_6hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_6hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_6hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_7hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_7hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_7hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_7hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_7hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_8hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_8hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_8hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_8hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_8hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     CLWMR_9hybridlevel"Cloud Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'CLWMR_9hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for CLWMR_9hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,CLWMR_9hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for CLWMR_9hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_10hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_10hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_10hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_10hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_10hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_11hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_11hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_11hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_11hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_11hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_12hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_12hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_12hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_12hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_12hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_13hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_13hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_13hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_13hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_13hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_14hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_14hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_14hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_14hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_14hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_15hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_15hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_15hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_15hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_15hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_16hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_16hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_16hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_16hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_16hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_17hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_17hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_17hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_17hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_17hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_18hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_18hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_18hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_18hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_18hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_19hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_19hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_19hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_19hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_19hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_1hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_1hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_1hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_1hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_1hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_20hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_20hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_20hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_20hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_20hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_21hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_21hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_21hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_21hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_21hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_22hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_22hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_22hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_22hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_22hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_23hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_23hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_23hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_23hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_23hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_24hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_24hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_24hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_24hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_24hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_25hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_25hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_25hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_25hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_25hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_26hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_26hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_26hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_26hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_26hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_27hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_27hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_27hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_27hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_27hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_28hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_28hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_28hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_28hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_28hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_29hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_29hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_29hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_29hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_29hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_2hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_2hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_2hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_2hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_2hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_30hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_30hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_30hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_30hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_30hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_31hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_31hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_31hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_31hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_31hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_32hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_32hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_32hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_32hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_32hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_33hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_33hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_33hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_33hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_33hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_34hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_34hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_34hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_34hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_34hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_35hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_35hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_35hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_35hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_35hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_36hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_36hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_36hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_36hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_36hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_37hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_37hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_37hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_37hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_37hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_38hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_38hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_38hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_38hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_38hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_39hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_39hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_39hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_39hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_39hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_3hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_3hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_3hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_3hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_3hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_40hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_40hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_40hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_40hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_40hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_41hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_41hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_41hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_41hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_41hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_42hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_42hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_42hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_42hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_42hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_43hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_43hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_43hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_43hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_43hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_44hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_44hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_44hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_44hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_44hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_45hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_45hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_45hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_45hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_45hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_46hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_46hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_46hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_46hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_46hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_47hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_47hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_47hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_47hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_47hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_48hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_48hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_48hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_48hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_48hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_49hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_49hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_49hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_49hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_49hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_4hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_4hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_4hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_4hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_4hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_50hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_50hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_50hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_50hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_50hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_51hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_51hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_51hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_51hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_51hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_52hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_52hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_52hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_52hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_52hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_53hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_53hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_53hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_53hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_53hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_54hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_54hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_54hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_54hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_54hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_55hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_55hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_55hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_55hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_55hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_56hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_56hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_56hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_56hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_56hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_57hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_57hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_57hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_57hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_57hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_58hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_58hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_58hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_58hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_58hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_59hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_59hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_59hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_59hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_59hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_5hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_5hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_5hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_5hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_5hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_60hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_60hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_60hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_60hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_60hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_61hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_61hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_61hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_61hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_61hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_62hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_62hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_62hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_62hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_62hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_63hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_63hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_63hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_63hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_63hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_64hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_64hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_64hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_64hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_64hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_65hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_65hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_65hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_65hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_65hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_6hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_6hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_6hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_6hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_6hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_7hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_7hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_7hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_7hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_7hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_8hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_8hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_8hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_8hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_8hybridlevel'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PRES_9hybridlevel"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'PRES_9hybridlevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PRES_9hybridlevel'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PRES_9hybridlevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PRES_9hybridlevel'
       endif
      endif

C   Variables of type INT
C

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C     latitude      "latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for latitude'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,latitude)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for latitude'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     longitude     "longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for longitude'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,longitude)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for longitude'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     time          "verification time generated by wgrib2 function verftime()"
C
      nf_status=NF_INQ_VARID(nf_fid,'time',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for time'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,time)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for time'
       endif
      endif


C   Variables of type CHAR
C

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
