      subroutine get_CONUS3D_hdr(FILENAME
     &          ,nx_dim, ny_dim, nz_dim
     &          ,ISTAT)
c *
c * To read header info of CONUS3D netcdf format digital data
c * 4/20/04 Dongsoo Kim,  Original
c
c * input
      CHARACTER FILENAME*132
c * output
      integer nx_dim, ny_dim, nz_dim, ISTAT
C*************************************
c
      include '/usr/local/apps/netcdf-3.4/include/netcdf.inc'
c
      STATUS=NF_OPEN(FILENAME,NF_NOWRITE,NCID)
      if (STATUS .ne. NF_NOERR) then
           ISTAT = -1
           print *, FILENAME,' is not avail'
           return
      endif
C
C * statements to fill x, y and levels
C
      nf_status = NF_INQ_DIMID(NCID,'x',dimid)
      nf_status = NF_INQ_DIMLEN(NCID,dimid,nx_dim)

      nf_status = NF_INQ_DIMID(NCID,'y',dimid)
      nf_status = NF_INQ_DIMLEN(NCID,dimid,ny_dim)

      nf_status = NF_INQ_DIMID(NCID,'levels',dimid)
      nf_status = NF_INQ_DIMLEN(NCID,dimid,nz_dim)
C *
      nf_status = NF_CLOSE(NCID)
c *
      RETURN
      END
