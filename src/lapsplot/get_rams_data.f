      subroutine get_rams_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L,NZ_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lat,lon
     +                   ,pres_3d
     +                   ,ht_p
     +                   ,clwc_p
     +                   ,cice_p
     +                   ,rain_p
     +                   ,snow_p
     +                   ,aext_p
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer x, y, z,nf_fid, nf_vid, nf_status
      real pres_3d(NX_L,NY_L,NZ_L)
      real ht_p(NX_L,NY_L,NZ_L)
      real clwc_p(NX_L,NY_L,NZ_L)
      real cice_p(NX_L,NY_L,NZ_L)
      real rain_p(NX_L,NY_L,NZ_L)
      real snow_p(NX_L,NY_L,NZ_L)
      real aext_p(NX_L,NY_L,NZ_L)
      real glat(NX_L,NY_L)
      real glon(NX_L,NY_L)
      real lat(NX_L,NY_L)
      real lon(NX_L,NY_L)
      real gri(NX_L,NY_L)
      real grj(NX_L,NY_L)
      real buf_3d(NX_L,NY_L,NZ_L)
      real buf_2d(NX_L,NY_L)
      logical l_hinterp /.false./
      logical wrapped /.false./
C
C  Open netcdf File for reading
C
      write(6,*)' Opening ',trim(filename)

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

      if(x .eq. NX_L .and. y .eq. NY_L .and. z .eq. NZ_L)then
        write(6,*)' model grid dimensions match laps grid'
      else
        write(6,*)' model grid deviates from laps grid'
        write(6,*)x,y,z 
        write(6,*)NX_L,NY_L,NZ_L
        stop
      endif

      call read_rams_data(nf_fid, x, y, z, i4time_sys,
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, glat, glon, ht_p, clwc_p,
     +     cice_p, rain_p, snow_p, aext_p, istatus)

      if(istatus .eq. 1)then
        write(6,*)' model corners'
        write(6,*)' LL ',glat(1,1),glon(1,1)
        write(6,*)' UL ',glat(1,NY_L),glon(1,NY_L)
        write(6,*)' LR ',glat(NX_L,1),glon(NX_L,1)
        write(6,*)' UR ',glat(NX_L,NY_L),glon(NX_L,NY_L)

        write(6,*)' laps corners'
        write(6,*)' LL ',lat(1,1),lon(1,1)
        write(6,*)' UL ',lat(1,NY_L),lon(1,NY_L)
        write(6,*)' LR ',lat(NX_L,1),lon(NX_L,1)
        write(6,*)' UR ',lat(NX_L,NY_L),lon(NX_L,NY_L)

        ricen = (float(x)+1.)/2.
        rjcen = (float(y)+1.)/2.
        call bilinear_laps(ricen,rjcen,NX_L,NY_L,glat,glat_cen)
        call bilinear_laps(ricen,rjcen,NX_L,NY_L,glon,glon_cen)
        call bilinear_laps(ricen,rjcen,NX_L,NY_L,lat,lat_cen)
        call bilinear_laps(ricen,rjcen,NX_L,NY_L,lon,lon_cen)
        write(6,*)' RAMS center is ',glat_cen,glon_cen
        write(6,*)' LAPS center is ',lat_cen,lon_cen
      else
        write(6,*)' istatus from read_rams_data is ',istatus
        return
      endif

      if(l_hinterp)then
        call latlon_to_grij_new(lat,lon,nx_l,ny_l,
     1                          glat,glon,gri,grj,x,y,
     1                          istatus_grij)

        write(6,*)' rams gri,grj at laps corners'
        write(6,*)' LL ',gri(1,1),grj(1,1)
        write(6,*)' UL ',gri(1,NY_L),grj(1,NY_L)
        write(6,*)' LR ',gri(NX_L,1),grj(NX_L,1)
        write(6,*)' UR ',gri(NX_L,NY_L),grj(NX_L,NY_L)

        call stats_missing(pres_3d,x,y,nz_l,r_missing_data,nmiss)
        write(6,*)' nmiss for pres_3d is ',nmiss
        buf_3d = pres_3d
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,pres_3d
     1                       ,wrapped)
        call stats_missing(pres_3d,x,y,nz_l,r_missing_data,nmiss)
        write(6,*)' nmiss for pres_3d is ',nmiss

        stop

        buf_3d = ht_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,ht_p
     1                       ,wrapped)

        buf_3d = clwc_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,clwc_p
     1                       ,wrapped)

        buf_3d = cice_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,cice_p
     1                       ,wrapped)

        buf_3d = rain_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,rain_p
     1                       ,wrapped)

        buf_3d = rain_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,rain_p
     1                       ,wrapped)

        buf_3d = snow_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,snow_p
     1                       ,wrapped)

        buf_2d = glat
        call hinterp_field_2d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_2d,glat
     1                       ,wrapped)

        buf_2d = glon
        call hinterp_field_2d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_2d,glon
     1                       ,wrapped)

        write(6,*)' remapped glat/glon arrays at laps corners'
        write(6,*)' LL ',glat(1,1),glon(1,1)
        write(6,*)' UL ',glat(1,NY_L),glon(1,NY_L)
        write(6,*)' LR ',glat(NX_L,1),glon(NX_L,1)
        write(6,*)' UR ',glat(NX_L,NY_L),glon(NX_L,NY_L)

      endif

      return
      end
C
C
      subroutine read_rams_data(nf_fid, x, y, z, i4time_sys,
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, glat, glon, ht, lwc,
     +     ice, rai, sno, aext, istatus)


      include 'netcdf.inc'
      integer x, y, z,nf_fid, nf_vid, nf_status
      integer ht_fcinv(z), ice_fcinv(z), lwc_fcinv(z), rai_fcinv(z),
     +     sno_fcinv(z)
      real glat( x, y), glon( x, y), ht( x,  y, z), ice( x,  y, z),
     +     level(z), lwc( x,  y, z), rai( x,  y, z), rtgt_ascii( x,
     +     y), sno( x,  y, z), topt( x, y), topt_ascii( x, y)
      real aext(x,y,z)

!     Declarations for 'write_zzz' call
      logical l_closest_time, l_closest_time_i, l_in_domain
      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)

      write(6,*)' subroutine read_rams_data'

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting r_missing_data'
          return
      endif

!     call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a,
!    1            topo_a,1.0,rnorth,south,east,west,istatus)
!     if(istatus .ne. 1)then
!         write(6,*)' Error in get_domain_perimeter'
!         return
!     endif

      write(6,*)' calling read_rams_netcdf'

      call read_rams_netcdf(nf_fid, x, y, z, glat, glon, ht, ice, 
     +     level, lwc, rai, aext, rtgt_ascii, sno, topt, topt_ascii, 
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
     +     ice, level, lwc, rai, aext, rtgt_ascii, sno, topt, topt_ascii, 
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
C     aext
C
      nf_status=NF_INQ_VARID(nf_fid,'aext',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for aext'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,aext)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for aext'
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


C   Variables of type DOUBLE
C


C   Variables of type CHAR
C

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif
C
      return
      end

      subroutine stats_missing(array,ni,nj,nk,r_missing_data,nmiss)

      real array(ni,nj,nk)

      nmiss = 0

      do k = 1,nk
      do i = 1,ni
      do j = 1,nj
          if(array(i,j,k) .eq. r_missing_data)then
              nmiss = nmiss + 1
          endif
      enddo ! j
      enddo ! i
      enddo ! k

      return
      end
