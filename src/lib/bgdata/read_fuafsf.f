      subroutine read_fuafsf_cdf(fullname, x, y, z, 
     +     ht, level, om, sh, t3, u3, v3,
     +     usfc, vsfc, tsfc, dsfc, psfc, mslp, zsfc,
     +     istatus)
C
      implicit none
      include 'netcdf.inc'
      character*(*) fullname
      character*200 cfname_int
      character*80 c8_proj
      integer x, y, z, nf_fid, nf_vid, nf_status
      integer lname,lname_in,lend,l,lenc8
      integer istatus

      real ht( x,  y,  z, 1)
     +,    om( x,  y,  z, 1)
     +,    sh( x,  y,  z, 1)
     +,    t3( x,  y,  z, 1)
     +,    u3( x,  y,  z, 1)
     +,    v3( x,  y,  z, 1)
     +,    w3( x,  y,  z, 1)

      real level(z)

      real tsfc( x,  y)
     +,    psfc( x,  y)
     +,    dsfc( x,  y)
     +,    mslp( x,  y)
     +,    usfc( x,  y)
     +,    vsfc( x,  y)
     +,    zsfc( x,  y)

      istatus = 0

      call get_c8_project(c8_proj,istatus)
      call downcase(c8_proj,c8_proj)
      call s_len(c8_proj,lenc8)
C
C  Open netcdf File for reading
C
      call s_len(fullname,lname_in)
      cfname_int=fullname(1:lname_in)//'.fua'
      call s_len(cfname_int,lname)
      print*,'Open/read ',cfname_int(1:lname)
      nf_status = NF_OPEN(cfname_int,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',cfname_int(1:lname)
        return
      endif
C
C     Variable        NETCDF Long Name
C      ht           "LAPS Fcst height"
C
      nf_status = NF_INQ_VARID(nf_fid,'ht',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ht'
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,ht)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ht'
        return
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(ht,x,y,z)
C
C     Variable        NETCDF Long Name
C      level        "level of data"
C
      nf_status = NF_INQ_VARID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var level'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,level)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var level'
        return
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(level,1,1,z)
C
C     Variable        NETCDF Long Name
C      om           "LAPS Fcst omega wind component"
C
      nf_status = NF_INQ_VARID(nf_fid,'om',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var om'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,om)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var om'
        return
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(om,x,y,z)
C
C     Variable        NETCDF Long Name
C      sh           "LAPS Fcst specific humidity"
C
      nf_status = NF_INQ_VARID(nf_fid,'sh',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var sh'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,sh)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var sh'
        return
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(sh,x,y,z)
C
C     Variable        NETCDF Long Name
C      t3           "LAPS Fcst temperature"
C
      nf_status = NF_INQ_VARID(nf_fid,'t3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var t3'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,t3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var t3'
        return
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(t3,x,y,z)
C
C     Variable        NETCDF Long Name
C      u3           "LAPS Fcst u wind component"
C
      nf_status = NF_INQ_VARID(nf_fid,'u3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var u3'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,u3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var u3'
        return
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(u3,x,y,z)
C
C     Variable        NETCDF Long Name
C      v3           "LAPS Fcst v wind component"
C
      nf_status = NF_INQ_VARID(nf_fid,'v3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var v3'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,v3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var v3'
        return
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(v3,x,y,z)

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif
c
c now fsf file:
c
      call get_directory_length(cfname_int,lend)
      search_dir: do l=lend,1,-1
         if(cfname_int(l:l).eq.'f')then
            if(cfname_int(l:l+2).eq.'fua')then
               exit search_dir
            endif
         endif
      enddo search_dir

      if(l.le.1)then
         print*,'didnt determine location of fua in'
         print*,'string for fsf. Return with no data'
         print*,'read_fuafsf.f: abort'
         return
      endif

      cfname_int(l:l+2)='fsf'

      cfname_int=cfname_int(1:lname_in)//'.fsf'

      call s_len(cfname_int,lname)
      print*,'Open/read ',cfname_int(1:lname)
      nf_status = NF_OPEN(cfname_int,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',cfname_int(1:lname)
        return
      endif

C
C     Variable        NETCDF Long Name
C      usf          "LAPS Fcst sfc u wind component"
C
      nf_status = NF_INQ_VARID(nf_fid,'usf',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var usf'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,usfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var usf'
        return
      endif
C
C     Variable        NETCDF Long Name
C      vsf          "LAPS Fcst sfc v wind component"
C
      nf_status = NF_INQ_VARID(nf_fid,'vsf',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vsf'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,vsfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vsf'
        return
      endif
C
C     Variable        NETCDF Long Name
C      vsf          "LAPS Fcst sfc temperature"
C
      nf_status = NF_INQ_VARID(nf_fid,'tsf',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tsf'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,tsfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tsf'
        return
      endif
C
C     Variable        NETCDF Long Name
C      vsf          "LAPS Fcst sfc dew point temperature"
C
      nf_status = NF_INQ_VARID(nf_fid,'dsf',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dsf'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,dsfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dsfc'
        return
      endif
C
C     Variable        NETCDF Long Name
C      vsf          "LAPS Fcst mean sea level pressure"
C
      nf_status = NF_INQ_VARID(nf_fid,'slp',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var slp'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,mslp)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var mslp'
        return
      endif
C
C     Variable        NETCDF Long Name
C      vsf          "LAPS Fcst sfc pressure"
C
      nf_status = NF_INQ_VARID(nf_fid,'psf',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var psf'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,psfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var psfc'
        return
      endif
C
C     Variable        NETCDF Long Name
C      ter          "LAPS sfc terrain"
C
      nf_status = NF_INQ_VARID(nf_fid,'ter',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ter'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,zsfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var zsfc'
        return
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 1
      return
      end
c
c =================================================
c
      subroutine swap_array_k(a1,nx,ny,nz)
c
      implicit none
      integer nx,ny,nz,k
      real  a1(nx,ny,nz)
      real, allocatable:: a2(:,:,:)

      allocate (a2(nx,ny,nz))

      do k=1,nz
         a2(:,:,k)=a1(:,:,nz-k+1)
      enddo
      a1=a2
      deallocate (a2)
      return
      end

