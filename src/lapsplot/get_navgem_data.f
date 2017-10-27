      subroutine get_navgem_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L,NZ_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,pres_p
     +                   ,clwc_p
     +                   ,cice_p
     +                   ,rain_p
     +                   ,snow_p
!    +                   ,seaice,snow_depth
     +                   ,lun_out
     +                   ,istatus)

      use mem_allsky, ONLY: aod_3d

      include 'netcdf.inc'

      character*(*) filename

      integer latitude, level, longitude,nf_fid, nf_vid, nf_status

      parameter (LVLP=21)

      real ht_p(NX_L,NY_L,LVLP)
      real pres_p(NX_L,NY_L,LVLP)
      real clwc_p(NX_L,NY_L,LVLP)
      real cice_p(NX_L,NY_L,LVLP)
      real rain_p(NX_L,NY_L,LVLP)
      real snow_p(NX_L,NY_L,LVLP)
      real coarse_ext_p(NX_L,NY_L,LVLP)
      real fine_ext_p(NX_L,NY_L,LVLP)
      real total_ext_p(NX_L,NY_L,LVLP)
      real seaice(NX_L,NY_L)
      real snow_depth(NX_L,NY_L)

      write(6,*)' Entered get_navgem_data'

      if(NZ_L .ne. LVLP)then
        write(6,*)' incorrect vertical levels',NZ_L,LVLP
        istatus=0
        return
      endif

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
C Get size of latitude
C
      nf_status=NF_INQ_DIMID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim latitude'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim latitude'
      endif
C
C Get size of level
C
      nf_status=NF_INQ_DIMID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim level'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,level)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim level'
      endif
C
C Get size of longitude
C
      nf_status=NF_INQ_DIMID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim longitude'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim longitude'
      endif
      call read_navgem_data(nf_fid, latitude, level, longitude,
     +     i4time_sys, ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, clwc_p, cice_p, pres_p, 
     +     coarse_ext_p, fine_ext_p, total_ext_p, LVLP,
     +     seaice, snow_depth, istatus)

      if(.true.)then
         write(6,*)' Use NAVGEM aerosols'
         aod_3d(:,:,:) = total_ext_p(:,:,:)
      endif

      return
      end
C
C
      subroutine read_navgem_data(nf_fid, latitude, level, longitude,
     +     i4time_sys, ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, clw_p, ciw_p, pres_p,
     +     coarse_ext_p, fine_ext_p, total_ext_p, LVLP,
     +     seaice, snow_depth, istatus)


      include 'netcdf.inc'
      integer latitude, level, longitude,nf_fid, nf_vid, nf_status


      real ABF_aot( longitude, latitude), ABF_ext( longitude, 
     +     latitude, level), ABF_mass_concentration( longitude, 
     +     latitude, level), MASL( longitude,  latitude, level),
     +     SS_aot( longitude, latitude), SS_ext( longitude, 
     +     latitude, level), SS_mass_concentration( longitude, 
     +     latitude, level), cf( longitude,  latitude, level), ciw(
     +     longitude,  latitude, level), clw( longitude,  latitude,
     +     level), coarse_aot( longitude, latitude), coarse_ext(
     +     longitude,  latitude, level),
     +     coarse_mode_mass_concentration( longitude,  latitude,
     +     level), dust_aot( longitude, latitude), dust_ext(
     +     longitude,  latitude, level), dust_mass_concentration(
     +     longitude,  latitude, level), dz( longitude,  latitude,
     +     level), fine_aot( longitude, latitude), fine_ext(
     +     longitude,  latitude, level),
     +     fine_mode_mass_concentration( longitude,  latitude,
     +     level), height_AGL( longitude,  latitude, level),
     +     latitude_a(latitude), level_a(level), longitude_a(longitude),
     +     pressure( longitude,  latitude, level), ps( longitude,
     +     latitude), q( longitude,  latitude, level), rh( longitude,
     +     latitude, level), sigma_a(level), sigma_b(level),
     +     smoke_aot( longitude, latitude), smoke_ext( longitude, 
     +     latitude, level), smoke_mass_concentration( longitude, 
     +     latitude, level), temperature( longitude,  latitude,
     +     level), total_aot( longitude, latitude), total_ext(
     +     longitude,  latitude, level), tprec( longitude, latitude),
     +     u( longitude,  latitude, level), v( longitude,  latitude,
     +     level), w( longitude,  latitude, level), wind_speed(
     +     longitude,  latitude, level),
     +     seaice(longitude, latitude), snow_depth(longitude, latitude)
c
      real pres_p(longitude,latitude,LVLP)
      real ciw_p(longitude,latitude,LVLP)
      real clw_p(longitude,latitude,LVLP)
      real coarse_ext_p(longitude,latitude,LVLP)
      real fine_ext_p(longitude,latitude,LVLP)
      real total_ext_p(longitude,latitude,LVLP)

!     Declarations for 'write_nvg' call
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

      call read_navgem_netcdf(nf_fid, latitude, level, longitude, 
     +     ABF_aot, ABF_ext, ABF_mass_concentration, MASL, SS_aot, 
     +     SS_ext, SS_mass_concentration, cf, ciw, clw, coarse_aot, 
     +     coarse_ext, coarse_mode_mass_concentration, dust_aot, 
     +     dust_ext, dust_mass_concentration, dz, fine_aot, fine_ext, 
     +     fine_mode_mass_concentration, height_AGL, 
     +     latitude_a, level_a, longitude_a,
     +     pressure, ps, q, rh, sigma_a, sigma_b, 
     +     smoke_aot, smoke_ext, smoke_mass_concentration, 
     +     temperature, total_aot, total_ext, tprec, u, v, w, 
     +     wind_speed,seaice,snow_depth)
C
C The netcdf variables are filled - your nvg write call may go here
C

      write(6,*)' range of input ciw',minval(ciw)
     +                               ,maxval(ciw)

      write(6,*)' range of input clw',minval(clw)
     +                               ,maxval(clw)

      write(6,*)' range of input pressure',minval(pressure)
     +                                    ,maxval(pressure)

      write(6,*)' range of simgrid pres_p',minval(pres_p)
     +                                    ,maxval(pres_p)

!     Vertical Interpolation
      call vinterp_sub(r_missing_data,NX_L,NY_L,NX_L,NY_L
     .                ,1,NX_L,1,NY_L,LVLP,level
     .                ,pres_p,pressure,ciw,ciw_p)

!     Vertical Interpolation
      call vinterp_sub(r_missing_data,NX_L,NY_L,NX_L,NY_L
     .                ,1,NX_L,1,NY_L,LVLP,level
     .                ,pres_p,pressure,clw,clw_p)

!     Vertical Interpolation
      call vinterp_sub(r_missing_data,NX_L,NY_L,NX_L,NY_L
     .                ,1,NX_L,1,NY_L,LVLP,level
     .                ,pres_p,pressure,coarse_ext,coarse_ext_p)

!     Vertical Interpolation
      call vinterp_sub(r_missing_data,NX_L,NY_L,NX_L,NY_L
     .                ,1,NX_L,1,NY_L,LVLP,level
     .                ,pres_p,pressure,fine_ext,fine_ext_p)

!     Vertical Interpolation
      call vinterp_sub(r_missing_data,NX_L,NY_L,NX_L,NY_L
     .                ,1,NX_L,1,NY_L,LVLP,level
     .                ,pres_p,pressure,total_ext,total_ext_p)

      write(6,*)' range of interp ciw_p',minval(ciw_p)
     +                                  ,maxval(ciw_p)

      write(6,*)' range of interp clw_p',minval(clw_p)
     +                                  ,maxval(clw_p)

      write(6,*)' range of interp coarse_ext_p',minval(coarse_ext_p)
     +                                         ,maxval(coarse_ext_p)

      write(6,*)' range of interp fine_ext_p',minval(fine_ext_p)
     +                                       ,maxval(fine_ext_p)

      write(6,*)' range of interp total_ext_p',minval(total_ext_p)
     +                                        ,maxval(total_ext_p)

      return
      end
C
C  Subroutine to read the file 
C
      subroutine read_navgem_netcdf(nf_fid, latitude, level, 
     +     longitude, ABF_aot, ABF_ext, ABF_mass_concentration, MASL, 
     +     SS_aot, SS_ext, SS_mass_concentration, cf, ciw, clw, 
     +     coarse_aot, coarse_ext, coarse_mode_mass_concentration, 
     +     dust_aot, dust_ext, dust_mass_concentration, dz, fine_aot, 
     +     fine_ext, fine_mode_mass_concentration, height_AGL, 
     +     latitude_a, level_a, longitude_a,
     +     pressure, ps, q, rh, sigma_a, 
     +     sigma_b, smoke_aot, smoke_ext, smoke_mass_concentration, 
     +     temperature, total_aot, total_ext, tprec, u, v, w, 
     +     wind_speed,seaice,snow_depth)
C
      include 'netcdf.inc'
      integer latitude, level, longitude,nf_fid, nf_vid, nf_status

      real ABF_aot( longitude, latitude), ABF_ext( longitude, 
     +     latitude, level), ABF_mass_concentration( longitude, 
     +     latitude, level), MASL( longitude,  latitude, level),
     +     SS_aot( longitude, latitude), SS_ext( longitude, 
     +     latitude, level), SS_mass_concentration( longitude, 
     +     latitude, level), cf( longitude,  latitude, level), ciw(
     +     longitude,  latitude, level), clw( longitude,  latitude,
     +     level), coarse_aot( longitude, latitude), coarse_ext(
     +     longitude,  latitude, level),
     +     coarse_mode_mass_concentration( longitude,  latitude,
     +     level), dust_aot( longitude, latitude), dust_ext(
     +     longitude,  latitude, level), dust_mass_concentration(
     +     longitude,  latitude, level), dz( longitude,  latitude,
     +     level), fine_aot( longitude, latitude), fine_ext(
     +     longitude,  latitude, level),
     +     fine_mode_mass_concentration( longitude,  latitude,
     +     level), height_AGL( longitude,  latitude, level),
     +     latitude_a(latitude), level_a(level), longitude_a(longitude),
     +     pressure( longitude,  latitude, level), ps( longitude,
     +     latitude), q( longitude,  latitude, level), rh( longitude,
     +      latitude, level), sigma_a(level), sigma_b(level),
     +     smoke_aot( longitude, latitude), smoke_ext( longitude, 
     +     latitude, level), smoke_mass_concentration( longitude, 
     +     latitude, level), temperature( longitude,  latitude,
     +     level), total_aot( longitude, latitude), total_ext(
     +     longitude,  latitude, level), tprec( longitude, latitude),
     +     u( longitude,  latitude, level), v( longitude,  latitude,
     +     level), w( longitude,  latitude, level), wind_speed(
     +     longitude,  latitude, level),
     +     seaice(longitude, latitude), snow_depth(longitude, latitude)

C
C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     ABF_aot       "Anthropogenic and biogenic fine aerosol optical thickness at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'ABF_aot',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for ABF_aot'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,ABF_aot)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for ABF_aot'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     ABF_ext       "Anthropogenic and biogenic fine aerosol extinction at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'ABF_ext',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for ABF_ext'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,ABF_ext)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for ABF_ext'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     ABF_mass_concentration"Anthropogenic and biogenic fine_mass_concentration"
C
      nf_status=NF_INQ_VARID(nf_fid,'ABF_mass_concentration',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for ABF_mass_concentration'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,ABF_mass_concentration)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for ABF_mass_concentration'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     MASL          "meters above sea level"
C
      nf_status=NF_INQ_VARID(nf_fid,'MASL',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for MASL'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,MASL)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for MASL'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     SS_aot        "Sea Salt AOT"
C
      nf_status=NF_INQ_VARID(nf_fid,'SS_aot',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for SS_aot'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,SS_aot)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for SS_aot'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     SS_ext        "Sea Salt  aerosol extinction at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'SS_ext',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for SS_ext'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,SS_ext)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for SS_ext'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     SS_mass_concentration"Sea_Salt_mass_concentration"
C
      nf_status=NF_INQ_VARID(nf_fid,'SS_mass_concentration',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for SS_mass_concentration'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,SS_mass_concentration)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for SS_mass_concentration'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     cf            "cloud fraction"
C
      nf_status=NF_INQ_VARID(nf_fid,'cf',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for cf'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,cf)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for cf'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     ciw           "cloud ice water"
C
      nf_status=NF_INQ_VARID(nf_fid,'ciw',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for ciw'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,ciw)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for ciw'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     clw           "cloud liquid water"
C
      nf_status=NF_INQ_VARID(nf_fid,'clw',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for clw'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,clw)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for clw'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     coarse_aot    "coarse mode  aerosol optical thickness at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'coarse_aot',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for coarse_aot'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,coarse_aot)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for coarse_aot'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     coarse_ext    "coarse mode  aerosol extinction at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'coarse_ext',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for coarse_ext'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,coarse_ext)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for coarse_ext'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dust_aot      "dust aerosol optical thickness at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'dust_aot',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dust_aot'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dust_aot)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dust_aot'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dust_ext      "dust  aerosol extinction at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'dust_ext',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dust_ext'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dust_ext)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dust_ext'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dust_mass_concentration"dust_mass_concentration"
C
      nf_status=NF_INQ_VARID(nf_fid,'dust_mass_concentration',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dust_mass_concentration'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dust_mass_concentration)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dust_mass_concentration'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dz            "layer thickness"
C
      nf_status=NF_INQ_VARID(nf_fid,'dz',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dz'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dz)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dz'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     fine_aot      "fine mode  aerosol optical thickness at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'fine_aot',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for fine_aot'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,fine_aot)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for fine_aot'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     fine_ext      "fine mode  aerosol extinction at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'fine_ext',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for fine_ext'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,fine_ext)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for fine_ext'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     height_AGL    "height above ground level"
C
      nf_status=NF_INQ_VARID(nf_fid,'height_AGL',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for height_AGL'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,height_AGL)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for height_AGL'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     latitude      "latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for latitude'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,latitude_a)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for latitude'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     level         "model levels"
C
      nf_status=NF_INQ_VARID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for level'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,level_a)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for level'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     longitude     "longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for longitude'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,longitude_a)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for longitude'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressure      "Pressue"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for pressure'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,pressure)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for pressure'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     ps            "surface pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'ps',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for ps'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,ps)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for ps'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     q             "Specific Humidity"
C
      nf_status=NF_INQ_VARID(nf_fid,'q',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for q'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,q)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for q'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     rh            "relative humidity"
C
      nf_status=NF_INQ_VARID(nf_fid,'rh',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for rh'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,rh)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for rh'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     sigma_a       "sigma_a"
C
      nf_status=NF_INQ_VARID(nf_fid,'sigma_a',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for sigma_a'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,sigma_a)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for sigma_a'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     sigma_b       "sigma_b"
C
      nf_status=NF_INQ_VARID(nf_fid,'sigma_b',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for sigma_b'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,sigma_b)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for sigma_b'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     smoke_aot     "Smoke AOT"
C
      nf_status=NF_INQ_VARID(nf_fid,'smoke_aot',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for smoke_aot'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,smoke_aot)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for smoke_aot'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     smoke_ext     "Smoke  aerosol extinction at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'smoke_ext',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for smoke_ext'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,smoke_ext)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for smoke_ext'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     smoke_mass_concentration"smoke_mass_concentration"
C
      nf_status=NF_INQ_VARID(nf_fid,'smoke_mass_concentration',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for smoke_mass_concentration'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,smoke_mass_concentration)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for smoke_mass_concentration'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     temperature   "Temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for temperature'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for temperature'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     total_aot     "Total AOT"
C
      nf_status=NF_INQ_VARID(nf_fid,'total_aot',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for total_aot'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,total_aot)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for total_aot'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     total_ext     "Total  aerosol extinction at 550nm"
C
      nf_status=NF_INQ_VARID(nf_fid,'total_ext',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for total_ext'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,total_ext)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for total_ext'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     tprec         "total precipitation"
C
      nf_status=NF_INQ_VARID(nf_fid,'tprec',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for tprec'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,tprec)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for tprec'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     u             "eastward wind"
C
      nf_status=NF_INQ_VARID(nf_fid,'u',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for u'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,u)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for u'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     v             "northward wind"
C
      nf_status=NF_INQ_VARID(nf_fid,'v',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for v'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,v)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for v'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     w             "omega"
C
      nf_status=NF_INQ_VARID(nf_fid,'w',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for w'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,w)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for w'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     wind_speed    "wind speed (cup)"
C
      nf_status=NF_INQ_VARID(nf_fid,'wind_speed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for wind_speed'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,wind_speed)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for wind_speed'
       endif
      endif

C
C     Variable        NETCDF Long Name
C     seaice
C
      nf_status=NF_INQ_VARID(nf_fid,'seaice',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for seaice'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,seaice)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for seaice'
       endif
      endif

C
C     Variable        NETCDF Long Name
C     snow_depth    
C
      nf_status=NF_INQ_VARID(nf_fid,'snow_depth',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for snow_depth'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,snow_depth)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for snow_depth'
       endif
      endif


C   Variables of type INT
C

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
