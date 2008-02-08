
      subroutine get_wideband_netcdf(FILENAME
     &          ,rad_dim,Z_dim,V_dim
     &          ,image_Z,image_V,image_W
     &          ,Azim_ang, Elev_ang, resolV
     &          ,ISTAT)
c
c * A subroutine to read wideband level II data
c   00 Nov 2001  Dongsoo Kim    Original
c 
      integer rad_max, ang_max
      parameter (rad_max = 999, ang_max = 999)
      integer rad_dim, Z_dim, V_dim
      byte  Z(Z_dim,rad_dim),W(V_dim,rad_dim)
     &     ,V(V_dim,rad_dim)
c
c * output
c
      real  Azim_ang(rad_dim), Elev_ang(rad_dim), resolV
      integer*2  image_Z(460,rad_max)
     &     ,image_V(920,rad_max)
     &     ,image_W(920,rad_max)
c
c *************************************
c
      CHARACTER FILENAME*72
c
      include 'netcdf.inc'
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
      nf_status = NF_INQ_VARID(NCID,'Z',nf_vid)
      nf_status = NF_GET_VAR_INT1(NCID,nf_vid,Z)
      nf_status = NF_INQ_VARID(NCID,'V',nf_vid)
      nf_status = NF_GET_VAR_INT1(NCID,nf_vid,V)
      nf_status = NF_INQ_VARID(NCID,'W',nf_vid)
      nf_status = NF_GET_VAR_INT1(NCID,nf_vid,W)
c
      nf_status = NF_INQ_VARID(NCID,'radialAzim',nf_vid)
      nf_status = NF_GET_VAR_REAL(NCID,nf_vid,Azim_ang)
      nf_status = NF_INQ_VARID(NCID,'radialElev',nf_vid)
      nf_status = NF_GET_VAR_REAL(NCID,nf_vid,Elev_ang)
      nf_status = NF_INQ_VARID(NCID,'resolutionV',nf_vid)
      nf_status = NF_GET_VAR_REAL(NCID,nf_vid,resolV)

      nf_status = NF_CLOSE(NCID)
c
      do j=1,rad_dim
      do i=1,Z_dim
         image_Z(i,j) = iand(Z(i,j),255)
      enddo
      enddo
c
      do j=1,rad_dim
      do i=1,V_dim
         image_V(i,j) = iand(V(i,j),255)
         image_W(i,j) = iand(W(i,j),255)
      enddo
      enddo
c
      RETURN
      END
