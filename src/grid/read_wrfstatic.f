
      subroutine read_wrfstatic(ni,nj,lat,lon,filename,topo,istatus)
      include 'netcdf.inc'
      integer Time, south_north, west_east, nf_fid, nf_vid, nf_status
      character*5 path
      character*255 filename
      data path/'data/'/

      real lat(ni,nj),lon(ni,nj),topo(ni,nj)
C
C  Open netcdf File for reading
C
      istatus = 1

      call s_len(filename,lenp)
      write(6,*)' read_wrfstatic - file is ',filename(1:lenp)

      nf_status = NF_OPEN(filename(1:lenp),NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN geo_em.d01.nc'
        istatus = 0
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of Time
C
      nf_status = NF_INQ_DIMID(nf_fid,'Time',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim Time'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,Time)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim Time'
      endif
C
C Get size of south_north
C
      nf_status = NF_INQ_DIMID(nf_fid,'south_north',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim south_north'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,south_north)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim south_north'
      endif
C
C Get size of west_east
C
      nf_status = NF_INQ_DIMID(nf_fid,'west_east',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim west_east'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,west_east)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim west_east'
      endif

      write(6,*)' WPS dims are ',west_east,south_north
      if(west_east .ne. ni .or. south_north .ne. nj)then
          write(6,*)' ERROR: WPS dims differ from LAPS'
          write(6,*)
     1    ' Only matching grids are currently supported: should be '
     1    ,ni,nj
          istatus = 0
          return
      endif

      call read_wrfstatic_sub(nf_fid , Time, south_north, west_east
     1                       ,topo,lat,lon,istatus)                     

      return
      end
C
C
      subroutine read_wrfstatic_sub(nf_fid, Time, south_north
     1                    ,west_east,HGT_M, XLAT_M, XLONG_M,istatus)
      include 'netcdf.inc'
      integer Time, south_north, west_east, nf_fid, nf_vid, nf_status
      real HGT_M( west_east,  south_north, Time), 
     +   XLAT_M( west_east,  south_north, Time), 
     +   XLONG_M( west_east,  south_north, Time)
      call read_netcdf(nf_fid , Time, south_north, west_east,
     +    HGT_M, XLAT_M, XLONG_M)
C
C The netcdf variables are filled - your code goes here
C
      return
      end

      subroutine read_netcdf(nf_fid , Time, south_north, west_east,
     +    HGT_M, XLAT_M, XLONG_M)
      include 'netcdf.inc'
      integer Time, south_north, west_east, nf_fid, nf_vid, nf_status

      real HGT_M( west_east,  south_north, Time), 
     +   XLAT_M( west_east,  south_north, Time), 
     +   XLONG_M( west_east,  south_north, Time)
C
C     Variable        NETCDF Long Name
C      HGT_M        
C
      nf_status = NF_INQ_VARID(nf_fid,'HGT_M',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var HGT_M'
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,HGT_M)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ HGT_M '
      endif
C
C     Variable        NETCDF Long Name
C      XLAT_M       
C
      nf_status = NF_INQ_VARID(nf_fid,'XLAT_M',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var XLAT_M'
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,XLAT_M)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ XLAT_M '
      endif
C
C     Variable        NETCDF Long Name
C      XLONG_M      
C
      nf_status = NF_INQ_VARID(nf_fid,'XLONG_M',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var XLONG_M'
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,XLONG_M)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ XLONG_M '
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
