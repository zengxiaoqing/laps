
c *
c * To read header info of wideband image data
c *
      subroutine get_wideband_hdr(FILENAME
     &          ,rad_dim, Z_dim, V_dim
     &          ,siteLat,siteLon,siteAlt,VCP
     &          ,ISTAT)
c
      parameter (rad_max = 999, ang_max = 999)
c * output
      integer rad_dim, Z_dim, V_dim, VCP
      real    siteLat,siteLon,siteAlt
C*************************************
      CHARACTER FILENAME*72
c
      include 'netcdf.inc'
c
      STATUS=NF_OPEN(FILENAME,NF_NOWRITE,NCID)
      if (STATUS .ne. NF_NOERR) then
           ISTAT = -1
c           print *, FILENAME,' is not avail'
           return
      endif
C
C * statements to fill Nx and Ny
C
      nf_status = NF_INQ_DIMID(NCID,'radial',nf_vid)
      nf_status = NF_INQ_DIMLEN(NCID,nf_vid,rad_dim)
      nf_status = NF_INQ_DIMID(NCID,'Z_bin',nf_vid)
      nf_status = NF_INQ_DIMLEN(NCID,nf_vid,Z_dim)
      nf_status = NF_INQ_DIMID(NCID,'V_bin',nf_vid)
      nf_status = NF_INQ_DIMLEN(NCID,nf_vid,V_dim)

      nf_status = NF_INQ_VARID(NCID,'VCP',nf_vid)
      nf_status = NF_GET_VAR_INT(NCID,nf_vid,VCP)
      nf_status = NF_INQ_VARID(NCID,'siteLat',nf_vid)
      nf_status = NF_GET_VAR_REAL(NCID,nf_vid,siteLat)
      nf_status = NF_INQ_VARID(NCID,'siteLon',nf_vid)
      nf_status = NF_GET_VAR_REAL(NCID,nf_vid,siteLon)
      nf_status = NF_INQ_VARID(NCID,'siteAlt',nf_vid)
      nf_status = NF_GET_VAR_REAL(NCID,nf_vid,siteAlt)
c
      nf_status = NF_CLOSE(NCID)
      ISTAT = 0
c *
      RETURN
      END
