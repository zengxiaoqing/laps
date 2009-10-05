      subroutine read_netcdf_real(nf_fid,fname,n1,f,start,count,istatus)

      implicit none

      include 'netcdf.inc'
      include 'bgdata.inc'
      integer n1,i, nf_fid, nf_vid,istatus,nf_status
      integer start(10),count(10)
      real f(n1) , nfmissing
      character*(*) fname

      print*,'HERE: in read_netcdf_real A, n1=',n1

      istatus=0
      nf_status = NF_INQ_VARID(nf_fid,fname,nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ', fname
        istatus = 1
        return
      endif

      if(start(1).eq.0.and.count(1).eq.0) then
         print*,'HERE: in read_netcdf_real B1'
         nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,f)
      else
         print*,'HERE: in read_netcdf_real B2',start,count
         nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,f)
      endif
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ ', fname
        istatus = 1
        return
      endif

      print*,'HERE: in read_netcdf_real C'

      if(fname.ne.'isoLevel')then
         nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue'
     .,nfmissing)
         if(nf_status.ne.NF_NOERR) then
            print *, NF_STRERROR(nf_status)
         endif
      endif
      do i=1,n1
         if(f(i).eq.nfmissing) then
            f(i)=missingflag
            istatus=istatus-1
         endif
      enddo
      
      return
      end
