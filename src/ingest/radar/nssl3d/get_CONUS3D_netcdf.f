      subroutine get_CONUS3D_netcdf(FILENAME
     &          ,nx_dim,ny_dim,nz_dim
     &          ,xxMin,xxMax,yyMin,yyMax, mrefl, ISTAT) 
c *
c * A subroutine to read 1 km resolution NSSL 3D in netcdf format
c * 4/19/2004 Dongsoo Kim, Original ingest routine
c *
c * input
      CHARACTER FILENAME*132
      integer*4 nx_dim, ny_dim, nz_dim
c * output
      real*4   mrefl(nx_dim, ny_dim, nz_dim)
      real*4   xxMin, xxMax, yyMin, yyMax
      integer  NCID
c
      include '/usr/local/apps/netcdf-3.4/include/netcdf.inc'
c
      STATUS=NF_OPEN(FILENAME,NF_NOWRITE,NCID)
      if (STATUS .ne. NF_NOERR) then
           ISTAT = -1
           print *, FILENAME,' is not avail'
           return
      endif
c
c  statements to fill image                          
c
      nf_status = NF_INQ_VARID(NCID,'mrefl_mosaic',nf_vid)
      nf_status = NF_GET_VAR_REAL(NCID,nf_vid,mrefl)
c
c  read global attributes
c
      nf_status = NF_GET_ATT_REAL(NCID,NF_GLOBAL,'xMin', xxMin)
      nf_status = NF_GET_ATT_REAL(NCID,NF_GLOBAL,'xMax', xxMax)
      nf_status = NF_GET_ATT_REAL(NCID,NF_GLOBAL,'yMin', yyMin)
      nf_status = NF_GET_ATT_REAL(NCID,NF_GLOBAL,'yMax', yyMax)
c
      nf_status = NF_CLOSE(NCID)
c
      RETURN
      END
