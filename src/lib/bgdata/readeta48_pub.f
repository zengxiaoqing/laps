      subroutine get_eta48_dims(filename,NX,NY,NZ,nf_status)
      implicit none
      include 'netcdf.inc'
      integer NZ, nrecs, NX, NY, nf_fid, nf_vid, nf_status
      character*(*) filename
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ', filename
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of isoLevel
C
      nf_status = NF_INQ_DIMID(nf_fid,'isoLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim isoLevel'
        return
      endif

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,NZ)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim isoLevel'
        return
      endif

C
C Get size of record
C
      nf_status = NF_INQ_DIMID(nf_fid,'record',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
        return
      endif

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,nrecs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
        return
      endif

      if(nrecs.ne.1) then
         PRINT *, 'ERROR IN ETA48 FORMAT'
         STOP
      endif
C
C Get size of x
C
      nf_status = NF_INQ_DIMID(nf_fid,'x',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
        return
      endif

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,NX)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
        return
      endif

C
C Get size of y
C
      nf_status = NF_INQ_DIMID(nf_fid,'y',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
        return
      endif

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,NY)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
        return
      endif


      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif


      return
      end
C




C
C
C  Subroutine to read the file "ETA 48 km AWIPS Regional CONUS " 
C
      subroutine read_eta_conusc(fname , NX,NY,NZ, ht, p, t, uw, vw,
     +     rh, ht_sfc, p_sfc, rh_sfc, t_sfc, uw_sfc, vw_sfc, mslp
     +     ,istatus)

      include 'netcdf.inc'
      integer NX,NY,NZ, nf_fid, nf_vid, nf_status, k
      character*(*) fname
      real MSLP( NX,  NY),    ht( NX,  NY,  NZ), 
     +     ht_sfc( NX,  NY),  p(NZ), p_sfc( NX,  NY), 
     +     rh( NX,  NY,  NZ), rh_sfc( NX,  NY), 
     +     t( NX,  NY,  NZ),  t_sfc( NX,  NY), 
     +     uw( NX,  NY,  NZ), uw_sfc( NX,  NY), 
     +     vw( NX,  NY,  NZ), vw_sfc( NX,  NY), tmp(nz)
      integer nxny,nxnynz
      logical reverse_fields
      data reverse_fields/.false./

      nf_status = NF_OPEN(fname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',fname
        return
      endif

      nxny = nx*ny
      nxnynz=nx*ny*nz
C
C     Variable        NETCDF Long Name
C      p            "isobaric levels" 
C
      call read_netcdf_real(nf_fid,'isoLevel',nz,p,0,0,nf_status)

c
c If p is going from top to bottom, reverse it and set a flag to 
c also reverse the other fields
c 
      if(p(1).lt.p(nz)) then
         do k=1,nz
            tmp(nz-k+1) = p(k)
         enddo
         do k=1,nz
            p(k) = tmp(k)
         enddo
         reverse_fields=.true.
      endif


C
C     Variable        NETCDF Long Name
C      ht           "geopotential height" 
C
      call read_netcdf_real(nf_fid,'gh',nxnynz,ht,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the height field is missing'
         print*, 'ABORTING file processing for eta file ',fname
         istatus=0
         return
      endif
         
      if(reverse_fields) call swap_array(nxny,nz,ht)
C
C     Variable        NETCDF Long Name
C      rh           "relative humidity" 
C
      call read_netcdf_real(nf_fid,'rh',nxnynz,rh,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the rh field is missing'
         print*, 'ABORTING file processing for eta file ',fname
         istatus=0
         return
      endif

      if(reverse_fields) call swap_array(nxny,nz,rh)

C
C     Variable        NETCDF Long Name
C      t            "temperature" 
C
      call read_netcdf_real(nf_fid,'t',nxnynz,t,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the temp field is missing'
         print*, 'ABORTING file processing for eta file ',fname
         istatus=0
         return
      endif


      if(reverse_fields) call swap_array(nxny,nz,t)
C
C     Variable        NETCDF Long Name
C      uw           "u-component of wind" 
C
      call read_netcdf_real(nf_fid,'uw',nxnynz,uw,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the u-wind field is missing'
         print*, 'ABORTING file processing for eta file ',fname
         istatus=0
         return
      endif


      if(reverse_fields) call swap_array(nxny,nz,uw)
C
C     Variable        NETCDF Long Name
C      vw           "v-component of wind" 
C
      call read_netcdf_real(nf_fid,'vw',nxnynz,vw,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the v-wind field is missing'
         print*, 'ABORTING file processing for eta file ',fname
         istatus=0
         return
      endif


      if(reverse_fields) call swap_array(nxny,nz,vw)
C
C     Variable        NETCDF Long Name
C      MSLP         "ETA mean sea level pressure" 
C

      call read_netcdf_real(nf_fid,'emspMSL',nxny,mslp,0,0,nf_status)

C
C     Variable        NETCDF Long Name
C      ht_sfc       "surface geopotential height" 
C
      call read_netcdf_real(nf_fid,'gh_sfc',nxny,ht_sfc,0,0,nf_status)

C
C     Variable        NETCDF Long Name
C      p_sfc        "surface pressure" 
C
      call read_netcdf_real(nf_fid,'p_sfc',nxny,p_sfc,0,0,nf_status)

C
C     Variable        NETCDF Long Name
C      rh_sfc       "relative humidity 2m fixed height abv ground" 
C
      call read_netcdf_real(nf_fid,'rh_2mFH',nxny,rh_sfc,0,0,nf_status)

C
C     Variable        NETCDF Long Name
C      t_sfc        "temperature 2m fixed height abv ground" 
C
      call read_netcdf_real(nf_fid,'t_2mFH',nxny,t_sfc,0,0,nf_status)

C
C     Variable        NETCDF Long Name
C      uw_sfc       "u wind component 10m fixed height abv ground" 
C
      call read_netcdf_real(nf_fid,'uw_10mFH',nxny,uw_sfc,0,0,nf_status)
C
C     Variable        NETCDF Long Name
C      vw_sfc       "v wind component 10m fixed height abv ground" 
C
      call read_netcdf_real(nf_fid,'vw_10mFH',nxny,vw_sfc,0,0,nf_status)

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        istatus = 0
        return
      endif
      istatus = 1

      return
      end


      subroutine swap_array(n1,n2,a1)
      integer i,j,n1,n2
      real a1(n1,n2),a2(n1,n2)
      do j=1,n2
         do i=1,n1
            a2(i,j)=a1(i,n2-j+1)
         enddo
      enddo
      do j=1,n2
         do i=1,n1
            a1(i,j)=a2(i,j)
         enddo
      enddo
      return
      end
