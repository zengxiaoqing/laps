      subroutine get_eta48_dims(filename,NX,NY,NZ
     &,StdLat1,StdLat2,Lon0,La1,Lo1,La2,Lo2,istatus)
c
c updated to get all nav info - JS 4-01
c
      implicit none
      include 'netcdf.inc'

      integer NZ, nrecs, NX, NY, nf_fid, nf_vid
      integer nf_status,istatus
      real    StdLat1,StdLat2
      real    Lon0
      real    La1,Lo1
      real    La2,Lo2

      character filename*200
C
C  Open netcdf File for reading
C
      print*,'filename =',filename

      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ', filename
      endif
C
C  Fill all dimension values
C
C
C  Get size of isoLevel (returned as NZ)
C
      nf_status = NF_INQ_DIMID(nf_fid,'isoLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim isoLevel'
      endif

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,NZ)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim isoLevel'
      endif

C
C Get size of record
C
      nf_status = NF_INQ_DIMID(nf_fid,'record',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
      endif

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,nrecs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
      endif

      if(nrecs.ne.1) then
         PRINT *, 'ERROR IN ETA48 FORMAT'
         STOP
      endif
C
C Get size of x (returned as NX)
C
      nf_status = NF_INQ_DIMID(nf_fid,'x',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,NX)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif

C
C Get size of y (returned as NY)
C
      nf_status = NF_INQ_DIMID(nf_fid,'y',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,NY)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif
C
C     Variable        NETCDF Long Name
C      IntLat1      "1st latitude of intersect earth-cone"
C                    returned as StdLat1
C
        nf_status = NF_INQ_VARID(nf_fid,'IntLat1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var IntLat1'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,StdLat1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var IntLat1'
      endif
C
C     Variable        NETCDF Long Name
C      IntLat2      "2nd latitude of intersect earth-cone"
C                    returned as StdLat2
C
        nf_status = NF_INQ_VARID(nf_fid,'IntLat2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var IntLat2'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,StdLat2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var IntLat2'
      endif
C
C     Variable        NETCDF Long Name
C      La1          "first latitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'La1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La1'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,La1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La1'
      endif
C
C     Variable        NETCDF Long Name
C      La2          "last latitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'La2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La2'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,La2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La2'
      endif
C
C     Variable        NETCDF Long Name
C      Lo1          "first longitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'Lo1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lo1'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Lo1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lo1'
      endif
C
C     Variable        NETCDF Long Name
C      Lo2          "last longitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'Lo2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lo2'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Lo2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lo2'
      endif
C
C     Variable        NETCDF Long Name
C      Lon0         "true longitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'Lon0',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lon0'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Lon0)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lon0'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      if(nf_status.ne.0)then
         print*,'Set istatus to Error condition'
         istatus = 0
         return
      endif

      istatus = 1
      return
      end
C
C
C  Subroutine to read the file "ETA 48 km AWIPS Regional CONUS " 
C
      subroutine read_eta_conusc(fname, NX,NY,NZ, ht,p,t,uw,vw,
     +     rh, pvv, ht_sfc, p_sfc, rh_sfc, td_sfc, t_sfc, uw_sfc, vw_sfc
     +     ,mslp,istatus)

c KML: CHANGES MADE APRIL 2004
c td_sfc (model 2m dew point) is now being read in
c td_sfc will ultimately be used in subroutine sfcbkgd
c KML: END

      include 'netcdf.inc'
      integer NX,NY,NZ, nf_fid, nf_status, k
      character*(*) fname
      real MSLP( NX,  NY),
     +    ht( NX, NY, NZ), ht_sfc( NX,  NY),
     +     p( NX, NY, NZ),  p_sfc( NX,  NY), 
     +    rh( NX,  NY,  NZ), rh_sfc( NX,  NY), 
     +                       td_sfc( NX,  NY),
     +     t( NX,  NY,  NZ),  t_sfc( NX,  NY), 
     +    uw( NX,  NY,  NZ), uw_sfc( NX,  NY), 
     +    vw( NX,  NY,  NZ), vw_sfc( NX,  NY),
     +    pvv(NX,  NY,  NZ),
     +    tmp(nz)
      real qsfc
      real make_ssh, make_td
      integer nxny,nxnynz
      character c8_project*8
      logical reverse_fields
      data reverse_fields/.false./

      istatus = 1

      nf_status = NF_OPEN(fname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',fname
        return
      endif

      nxny = nx*ny
      nxnynz=nx*ny*nz

      call get_c8_project(c8_project,istatus)
      if(istatus.ne. 1)then
         print*,'Error: returned from get_c8_project'
         print*,'Error: current routine: read_eta_conusc'
         return
      endif
C
C     Variable        NETCDF Long Name
C      p            "isobaric levels" 
C
      call read_netcdf_real(nf_fid,'isoLevel',nz,tmp,0,0,nf_status)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim isoLevel'
        return
      endif

c
c If p is going from top to bottom, reverse it and set a flag to 
c also reverse the other fields
c 
      if(tmp(1).lt.tmp(nz)) then

         do k=1,nz
         do j=1,ny
         do i=1,nx
            p(i,j,nz-k+1) = tmp(k)
         enddo
         enddo
         enddo
         reverse_fields=.true.
      else
         do k=1,nz
         do j=1,ny
         do i=1,nx
            p(i,j,k) = tmp(k)
         enddo
         enddo 
         enddo
      endif


C
C     Variable        NETCDF Long Name
C      ht           "geopotential height" 
C
      call read_netcdf_real(nf_fid,'gh',nxnynz,ht,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the height field is missing'
         print*, 'ABORTING file processing for eta file ',fname
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
         return
      endif


      if(reverse_fields) call swap_array(nxny,nz,vw)
C
C     Variable        NETCDF Long Name
C      pvv           "pressure vertical velocity"
C
      call read_netcdf_real(nf_fid,'pvv',nxnynz,pvv,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the w-wind field is missing'
         print*, 'ABORTING file processing for eta file ',fname
         return
      endif

      if(reverse_fields) call swap_array(nxny,nz,pvv)

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

c FSL netCDF file only has sfc rh
c ---------------------------------------
      if(c8_project.eq.'NIMBUS')then
C
C     Variable        NETCDF Long Name
C      rh_sfc       "relative humidity 2m fixed height abv ground" 
C
      print*,'FSL-NIMBUS: Read sfc RH'
      call read_netcdf_real(nf_fid,'rh_2mFH',nxny,rh_sfc,0,0,nf_status)
C
C     Variable        NETCDF Long Name
C      td_sfc        "dew point temp 2m fixed height abv ground"
C
      else
      print*,'Read td_2mFH: KML Upgrade'
      call read_netcdf_real(nf_fid,'td_2mFH',nxny,td_sfc,0,0,nf_status)
      endif
c ---------------------------------------
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
        return
      endif
C
C routines in file lib/bgdata/sfcbkgd.f require Td; thus, for FSL netcdf
C we must derive this from rh.
C
      if(c8_project.eq.'NIMBUS')then
       do j=1,ny
       do i=1,nx
         qsfc=make_ssh(p_sfc(i,j)/100.,t_sfc(i,j)-273.15,
     &rh_sfc(i,j)/100.,-132.)
         td_sfc(i,j)=make_td(p_sfc(i,j)/100.,t_sfc(i,j)-273.15,
     & qsfc,-132.)+273.15
       enddo
       enddo
      endif

      istatus = 0

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
