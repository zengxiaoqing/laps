C
C Get size of num_att
C
c     nf_status = NF_INQ_DIMID(nf_fid,'num_att',nf_vid)
c     if(nf_status.ne.NF_NOERR) then
c       print *, NF_STRERROR(nf_status)
c       print *,'dim num_att'
c     endif
c     nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,num_att)
c     if(nf_status.ne.NF_NOERR) then
c       print *, NF_STRERROR(nf_status)
c       print *,'dim num_att'
c     endif
c     call main_sub(nf_fid , num_att)

c     end
C
C
      subroutine read_orb_att(c_filespec,csatid, num_att, orb_att,
     &istatus)
      include 'netcdf.inc'
c     integer num_att, nf_fid, nf_vid, nf_status
      integer num_att, nf_fid, nf_status
      character*6 csatid
c     character*7 sat_name
      character*(*) c_filespec
      character*255 cfname
      double precision orb_att(num_att)

      istatus = 1
C
C  Open netcdf File for reading
C
      cfname=c_filespec//'/'//csatid//'_orbatt.dat'
      n=index(cfname,' ')-1
      nf_status = NF_OPEN(cfname,NF_NOWRITE,nf_fid)

      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',cfname(1:n)
        return
      endif

      call read_netcdf(nf_fid , num_att, orb_att,istatus)
C
C The netcdf variables are filled - your code goes here
C
      return
      end

      subroutine read_netcdf(nf_fid , num_att, orb_att,istatus)
      include 'netcdf.inc'
      integer num_att, nf_fid, nf_vid, nf_status

      character*7 sat_name
      double precision orb_att(num_att)

      istatus = 1
C
C     Variable        NETCDF Long Name
C      sat_name
      nf_status = NF_INQ_VARID(nf_fid,'sat_name',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var sat_name'
        return
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,sat_name)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ sat_name '
        return
      endif
      write(6,*)'sat name in netCDF file ',sat_name
C
C     Variable        NETCDF Long Name
C      orb_att      "orbit attitudes" 
C
        nf_status = NF_INQ_VARID(nf_fid,'orb_att',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var orb_att'
        return
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,orb_att)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ orb_att '
        return
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 0
      return
      end
