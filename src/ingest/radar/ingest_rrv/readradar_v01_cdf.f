      subroutine read_radar_v01_cdf(filename,nxv01,nyv01,nzv01,
     &cref_com,cvel_com,cnyq_com,refd,veld,nyqd,istatus)
c
c note: nxv01= lon
c       nyv01= lat
c       nzv01= level
c
      real*4      refd(nxv01,nyv01,nzv01)
      real*4      veld(nxv01,nyv01,nzv01)
      real*4      nyqd(nxv01,nyv01,nzv01)

      character       cvel_com(nzv01)*126
      character       cref_com(nzv01)*126
      character       cnyq_com(nzv01)*126
      character*(*)   filename

      include 'netcdf.inc'
C
C  Open netcdf File for reading
C
      n=index(filename,' ')-1
      print*,'Open netCDF file',filename(1:n)
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',filename(1:n)
      endif
C
C  Fill all dimension values
C
C
C Get size of lat
C
      nf_status = NF_INQ_DIMID(nf_fid,'lat',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim lat'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,lat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim lat'
      endif
C
C Get size of level
C
      nf_status = NF_INQ_DIMID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim level'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,level)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim level'
      endif
C
C Get size of lon
C
      nf_status = NF_INQ_DIMID(nf_fid,'lon',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim lon'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,lon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim lon'
      endif
C
C Get size of n_2d_grids
C
      nf_status = NF_INQ_DIMID(nf_fid,'n_2d_grids',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim n_2d_grids'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,n_2d_grids)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim n_2d_grids'
      endif
      if(lon.ne.nxv01)then
         print *,'lon ne nxv01'
         stop
      endif
      if(lat.ne.nyv01)then
         print *,'lat ne nyv01'
         stop
      endif
      if(level.ne.nzv01)then
         print *,'level ne nzv01'
         stop
      endif

      call main_sub(nf_fid , nyv01, nzv01, nxv01, n_2d_grids,
     +cref_com,cvel_com,cnyq_com,refd,veld,nyqd)

      end
C
C
      subroutine main_sub(nf_fid, lat, level, lon, n_2d_grids,
     +refd_comment,veld_comment,nyqd_comment,refd,veld,nyqd)
      include 'netcdf.inc'
      integer lat, level, lon, n_2d_grids, nf_fid, nf_vid, nf_status
      character*18 asctime
      character*12 laps_domain_file
      character*132 model
      character*126 nyqd_comment(level)
      character*132 origin
      character*126 refd_comment(level)
      character*126 veld_comment(level)
      integer fctimes_var, imax, jmax, kdim, kmax, level_var(level),
     +     lvl(n_2d_grids), num_variables, nyqd_fcinv(level),
     +     refd_fcinv(level), veld_fcinv(level), version

      real nyqd( lon,  lat, level), refd( lon,  lat, level), veld(
     +     lon,  lat, level)

      call read_netcdf(nf_fid , lat, level, lon, n_2d_grids, asctime,
     +     fctimes_var, imax, jmax, kdim, kmax, laps_domain_file,
     +     level_var, lvl, model, num_variables, nyqd, nyqd_comment,
     +     nyqd_fcinv, origin, refd, refd_comment, refd_fcinv, veld,
     +     veld_comment, veld_fcinv, version)

C
C The netcdf variables are filled - your code goes here
C
      return
      end
      subroutine read_netcdf(nf_fid , lat, level, lon, n_2d_grids,
     +     asctime, fctimes_var, imax, jmax, kdim, kmax,
     +     laps_domain_file, level_var, lvl, model, num_variables,
     +     nyqd, nyqd_comment, nyqd_fcinv, origin, refd,
     +     refd_comment, refd_fcinv, veld, veld_comment, veld_fcinv,
     +     version)

      include 'netcdf.inc'
      integer lat, level, lon, n_2d_grids, nf_fid, nf_vid, nf_status

      character*18 asctime
      character*12 laps_domain_file
      character*132 model
      character*126 nyqd_comment(level)
      character*132 origin
      character*126 refd_comment(level)
      character*126 veld_comment(level)
      integer fctimes_var, imax, jmax, kdim, kmax, level_var(level),
     +     lvl(n_2d_grids), num_variables, nyqd_fcinv(level),
     +     refd_fcinv(level), veld_fcinv(level), version

      real nyqd( lon,  lat, level), refd( lon,  lat, level), veld(
     +     lon,  lat, level)

C
C     Variable        NETCDF Long Name
C      asctime
        nf_status = NF_INQ_VARID(nf_fid,'asctime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var asctime'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,asctime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ asctime '
      endif
C
C     Variable        NETCDF Long Name
C      laps_domain_file 
        nf_status = NF_INQ_VARID(nf_fid,'laps_domain_file',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var laps_domain_file'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,laps_domain_file)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ laps_domain_file '
      endif
C
C     Variable        NETCDF Long Name
C      model         
        nf_status = NF_INQ_VARID(nf_fid,'model',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var model'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,model)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ model '
      endif
C
C     Variable        NETCDF Long Name
C      nyqd_comment
       nf_status = NF_INQ_VARID(nf_fid,'nyqd_comment',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nyqd_comment'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,nyqd_comment)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ nyqd_comment '
      endif
C
C     Variable        NETCDF Long Name
C      origin    
        nf_status = NF_INQ_VARID(nf_fid,'origin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var origin'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,origin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ origin '
      endif
C
C     Variable        NETCDF Long Name
C      refd_comment 
        nf_status = NF_INQ_VARID(nf_fid,'refd_comment',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var refd_comment'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,refd_comment)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ refd_comment '
      endif
C
C     Variable        NETCDF Long Name
C      veld_comment    
        nf_status = NF_INQ_VARID(nf_fid,'veld_comment',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var veld_comment'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,veld_comment)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ veld_comment '
      endif
C
C     Variable        NETCDF Long Name
C      fctimes_var  "forecast times" 
C
        nf_status = NF_INQ_VARID(nf_fid,'fctimes',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fctimes'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,fctimes_var)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ fctimes_var '
      endif
C
C     Variable        NETCDF Long Name
C      imax
        nf_status = NF_INQ_VARID(nf_fid,'imax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var imax'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,imax)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ imax '
      endif
C
C     Variable        NETCDF Long Name
C      jmax
        nf_status = NF_INQ_VARID(nf_fid,'jmax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var jmax'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,jmax)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ jmax '
      endif
C
C     Variable        NETCDF Long Name
C      kdim
        nf_status = NF_INQ_VARID(nf_fid,'kdim',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var kdim'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,kdim)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ kdim '
      endif
C
C     Variable        NETCDF Long Name
C      kmax
        nf_status = NF_INQ_VARID(nf_fid,'kmax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var kmax'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,kmax)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ kmax '
      endif
C
C     Variable        NETCDF Long Name
C      level_var    "level" 
C
        nf_status = NF_INQ_VARID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var level'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,level_var)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ level_var '
      endif
C
C     Variable        NETCDF Long Name
C      lvl
        nf_status = NF_INQ_VARID(nf_fid,'lvl',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lvl'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,lvl)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ lvl '
      endif
C
C     Variable        NETCDF Long Name
C      num_variables
        nf_status = NF_INQ_VARID(nf_fid,'num_variables',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var num_variables'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,num_variables)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ num_variables '
      endif
C
C     Variable        NETCDF Long Name
C      nyqd_fcinv
        nf_status = NF_INQ_VARID(nf_fid,'nyqd_fcinv',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nyqd_fcinv'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,nyqd_fcinv)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ nyqd_fcinv '
      endif
C
C     Variable        NETCDF Long Name
C      refd_fcinv
        nf_status = NF_INQ_VARID(nf_fid,'refd_fcinv',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var refd_fcinv'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,refd_fcinv)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ refd_fcinv '
      endif
C
C     Variable        NETCDF Long Name
C      veld_fcinv
        nf_status = NF_INQ_VARID(nf_fid,'veld_fcinv',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var veld_fcinv'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,veld_fcinv)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ veld_fcinv '
      endif
C
C     Variable        NETCDF Long Name
C      version
        nf_status = NF_INQ_VARID(nf_fid,'version',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var version'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,version)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ version '
      endif
C
C     Variable        NETCDF Long Name
C      nyqd         "3D radar nyquist velocity" 
C
      print*,'Read nyqd'
      nf_status = NF_INQ_VARID(nf_fid,'nyqd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nyqd'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,nyqd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ nyqd '
      endif
C
C     Variable        NETCDF Long Name
C      refd         "3D radar" 
C
      print*,'Read refd'
      nf_status = NF_INQ_VARID(nf_fid,'refd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var refd'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,refd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ refd '
      endif
C
C     Variable        NETCDF Long Name
C      veld         "3D radar" 
C
      print*,'Read veld'
      nf_status = NF_INQ_VARID(nf_fid,'veld',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var veld'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,veld)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ veld '
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
