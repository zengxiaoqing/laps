      subroutine get_rams_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,pres_3d
     +                   ,ht_p
     +                   ,clwc_p
     +                   ,cice_p
     +                   ,rain_p
     +                   ,snow_p
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer x, y, z,nf_fid, nf_vid, nf_status
      real pres_3d(NX_L,NY_L,21)
      real ht_p(NX_L,NY_L,21)
      real clwc_p(NX_L,NY_L,21)
      real cice_p(NX_L,NY_L,21)
      real rain_p(NX_L,NY_L,21)
      real snow_p(NX_L,NY_L,21)
C
C  Open netcdf File for reading
C
      nf_status=NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),filename
        istatus=0
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of x
C
      nf_status=NF_INQ_DIMID(nf_fid,'x',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,x)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif
C
C Get size of y
C
      nf_status=NF_INQ_DIMID(nf_fid,'y',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,y)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif
C
C Get size of z
C
      nf_status=NF_INQ_DIMID(nf_fid,'z',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim z'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,z)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim z'
      endif
      call read_rams_data(nf_fid, x, y, z, i4time_sys,
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, ht_p, clwc_p,
     +     cice_p, rain_p, snow_p, istatus)

      return
      end
C
C
      subroutine read_rams_data(nf_fid, x, y, z, i4time_sys,
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, ht, lwc,
     +     ice, rai, sno, istatus)


      include 'netcdf.inc'
      integer x, y, z,nf_fid, nf_vid, nf_status
      integer ht_fcinv(z), ice_fcinv(z), lwc_fcinv(z), rai_fcinv(z),
     +     sno_fcinv(z)
      real glat( x, y), glon( x, y), ht( x,  y, z), ice( x,  y, z),
     +     level(z), lwc( x,  y, z), rai( x,  y, z), rtgt_ascii( x,
     +     y), sno( x,  y, z), topt( x, y), topt_ascii( x, y)


!     Declarations for 'write_zzz' call
      logical l_closest_time, l_closest_time_i, l_in_domain
      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting r_missing_data'
          return
      endif
      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a,
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_domain_perimeter'
          return
      endif

      call read_rams_netcdf(nf_fid, x, y, z, glat, glon, ht, ice, 
     +     level, lwc, rai, rtgt_ascii, sno, topt, topt_ascii, 
     +     ht_fcinv, ice_fcinv, lwc_fcinv, rai_fcinv, sno_fcinv)
C
C The netcdf variables are filled - your zzz write call may go here
C
!     Initial loop through obs to get times and stanums
      where (lwc .eq. r_missing_data)lwc = 0.
      where (ice .eq. r_missing_data)ice = 0.
      where (rai .eq. r_missing_data)rai = 0.
      where (sno .eq. r_missing_data)sno = 0.

!     Fill in missing heights
      do iz = 1,z-1
          where (ht(:,:,iz) .eq. r_missing_data)ht(:,:,iz) = iz * 100.
      enddo 

      where (ht(:,:,z) .eq. r_missing_data)ht(:,:,z) = 20000.

      return
      end
C
C  Subroutine to read the file "RAMS microphysics interpolated to constant pressure surfaces" 
C
      subroutine read_rams_netcdf(nf_fid, x, y, z, glat, glon, ht, 
     +     ice, level, lwc, rai, rtgt_ascii, sno, topt, topt_ascii, 
     +     ht_fcinv, ice_fcinv, lwc_fcinv, rai_fcinv, sno_fcinv)
C
      include 'netcdf.inc'
      integer x, y, z,nf_fid, nf_vid, nf_status
      integer ht_fcinv(z), ice_fcinv(z), lwc_fcinv(z), rai_fcinv(z),
     +     sno_fcinv(z)
      real glat( x, y), glon( x, y), ht( x,  y, z), ice( x,  y, z),
     +     level(z), lwc( x,  y, z), rai( x,  y, z), rtgt_ascii( x,
     +     y), sno( x,  y, z), topt( x, y), topt_ascii( x, y)



C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     glat          
C
      nf_status=NF_INQ_VARID(nf_fid,'glat',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for glat'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,glat)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for glat'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     glon          
C
      nf_status=NF_INQ_VARID(nf_fid,'glon',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for glon'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,glon)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for glon'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     ht            
C
      nf_status=NF_INQ_VARID(nf_fid,'ht',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for ht'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,ht)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for ht'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     ice           
C
      nf_status=NF_INQ_VARID(nf_fid,'ice',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for ice'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,ice)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for ice'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     level         
C
      nf_status=NF_INQ_VARID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for level'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,level)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for level'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     lwc           
C
      nf_status=NF_INQ_VARID(nf_fid,'lwc',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for lwc'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,lwc)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for lwc'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     rai           
C
      nf_status=NF_INQ_VARID(nf_fid,'rai',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for rai'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,rai)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for rai'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     rtgt_ascii    
C
      nf_status=NF_INQ_VARID(nf_fid,'rtgt_ascii',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for rtgt_ascii'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,rtgt_ascii)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for rtgt_ascii'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     sno           
C
      nf_status=NF_INQ_VARID(nf_fid,'sno',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for sno'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,sno)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for sno'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     topt          
C
      nf_status=NF_INQ_VARID(nf_fid,'topt',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for topt'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,topt)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for topt'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     topt_ascii    
C
      nf_status=NF_INQ_VARID(nf_fid,'topt_ascii',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for topt_ascii'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,topt_ascii)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for topt_ascii'
       endif
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C     ht_fcinv      
C
      nf_status=NF_INQ_VARID(nf_fid,'ht_fcinv',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for ht_fcinv'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,ht_fcinv)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for ht_fcinv'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     ice_fcinv     
C
      nf_status=NF_INQ_VARID(nf_fid,'ice_fcinv',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for ice_fcinv'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,ice_fcinv)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for ice_fcinv'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     lwc_fcinv     
C
      nf_status=NF_INQ_VARID(nf_fid,'lwc_fcinv',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for lwc_fcinv'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,lwc_fcinv)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for lwc_fcinv'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     rai_fcinv     
C
      nf_status=NF_INQ_VARID(nf_fid,'rai_fcinv',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for rai_fcinv'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,rai_fcinv)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for rai_fcinv'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     sno_fcinv     
C
      nf_status=NF_INQ_VARID(nf_fid,'sno_fcinv',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for sno_fcinv'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,sno_fcinv)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for sno_fcinv'
       endif
      endif

C   Variables of type DOUBLE
C


C   Variables of type CHAR
C

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
