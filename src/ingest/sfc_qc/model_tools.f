c
c
        subroutine model(bgmodel,grid_lat,grid_lon,lat_s,lon_s,elev_s,
     &                   i4time_cur,laps_cycle_time,n_obs,ni,nj,nk,
     &                   mxsta,dt,dtd,du,dv,dpmsl,dalt,badflag)
c
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          Completely new version.
c       06 Jan 2000  Peter Stamus, NOAA/FSL
c          Yet another new version.  This one uses actual model data.
c
c     Notes: Currently (Jan 2000), if any of the LGA fields is missing
c            we set all of the trends to zero.  Plans are to check for
c            each variable separately, and to use a simple model (like
c            diurnal cycle) to estimate the trend better.
c
c*********************************************************************
c
c.....  Arrays for the LAPS domain lat/lon
c
        real grid_lat(ni,nj), grid_lon(ni,nj)
c
c.....  Arrays for 3d fields.
c
        real heights_3d(ni,nj,nk), p_3d(ni,nj,nk)
        real u_3d(ni,nj,nk), v_3d(ni,nj,nk)
        real t_3d(ni,nj,nk), q_3d(ni,nj,nk)
        real p_1d(nk)
        character var_lga*3
c
c.....  Arrays for current run interpolation from 3d to station.
c
        real t_cur(mxsta), td_cur(mxsta)
        real u_cur(mxsta), v_cur(mxsta)
        real psfc_cur(mxsta), pmsl_cur(mxsta), alt_cur(mxsta)
c
c.....  Arrays for previous run interpolation from 3d to station.
c
        real t_1(mxsta), td_1(mxsta)
        real u_1(mxsta), v_1(mxsta)
        real psfc_1(mxsta), pmsl_1(mxsta), alt_1(mxsta)
c
c.....  Arrays for trends of each variable over the past cycle.
c.....  This is the output.
c
        real dt(mxsta), dtd(mxsta), du(mxsta), dv(mxsta)
        real dpmsl(mxsta), dalt(mxsta)
c
c.....  Arrays for station location info (lat/lon/elevation and i,j)
c
        real lat_s(mxsta), lon_s(mxsta), elev_s(mxsta)
        integer ii(mxsta), jj(mxsta)
c
c.....  Other stuff
c
        real    lapse_t, lapse_td
        integer bgmodel
c
c
c.....  Start.  Fill the trend arrays with something.
c
        do i=1,mxsta
           dt(i)    = badflag
           dtd(i)   = badflag
           du(i)    = badflag
           dv(i)    = badflag
           dpmsl(i) = badflag
           dalt(i)  = badflag
        enddo !i
c
        t_ref = -47.           !for specific humidity to Td calc
        lapse_t = -0.01167     !std atm lapse rate
        lapse_td = -0.007      !       "
c
c.....  Find and store the i,j location of each station in the LAPS
c.....  grid.
c
        do ista=1,n_obs
          call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista),grid_lat,
     &       grid_lon,ni,nj,rii,rjj,istatus)
          ii(ista) = rii
          jj(ista) = rjj
        enddo !ista
c
c.....  Start with the current hour.
c.....  Get the 3d fields for the current time.
c
        print *,' Getting LAPS pressure levels'
        call get_pres_3d(i4time_cur,ni,nj,nk,p_3d,istatus)
        do k=1,nk  !for now...grab a column of pressures
           p_1d(k) = p_3d(1,1,k) * 0.01  !conv from Pa to mb too...
        enddo !k
c
        print *,' Getting LGA 3-d Heights'
        var_lga = 'HT '
        call get_modelfg_3d(i4time_cur,var_lga,ni,nj,nk,heights_3d,
     &                      istatus)
c
        if(istatus .ne. 1)  then
           print *,' LGA HT not available. '
           go to 800
        else

        endif
c
        print *,' Getting LGA 3-d Temperatures'
        var_lga = 'T3 '
        call get_modelfg_3d(i4time_cur,var_lga,ni,nj,nk,t_3d,istatus)
c
        if(istatus .ne. 1)  then
           print *,' LGA T not available. '
           go to 800
        else

        endif
c
        print *,' Getting LGA 3-d Specific Humidity'
        var_lga = 'SH '  ! specific humidity 
        call get_modelfg_3d(i4time_cur,var_lga,ni,nj,nk,q_3d,istatus)
c
        if(istatus .ne. 1)  then
           print *,' LGA SH not available. '
           go to 800
        else

        endif
c
        print *,' Getting LGA 3-d U-wind component'
        var_lga = 'U3 '
        call get_modelfg_3d(i4time_cur,var_lga,ni,nj,nk,u_3d,istatus)
c
        if(istatus .ne. 1)  then
           print *,' LGA U not available. '
           go to 800
        else
           
        endif
c
        print *,' Getting LGA 3-d V-wind component'
        var_lga = 'V3 '
        call get_modelfg_3d(i4time_cur,var_lga,ni,nj,nk,v_3d,istatus)
c
        if(istatus .ne. 1)  then
           print *,' LGA V not available. '
           go to 800
        else

        endif
c
c.....  Now, for the current time, interpolate the 3-d variables to the
c.....  proper location of each station.
c
        do ista=1,n_obs
c
c.....  Is the station within the LAPS grid?  If so, continue.  Otherwise,
c.....  set the current variables to badflag.
c
           if(ii(ista).ge.1 .and. ii(ista).le.ni .and.
     &        jj(ista).ge.1 .and. jj(ista).le.nj) then
c
c.....  Find the vertical index for the station elevation.
c
              iii = ii(ista)
              jjj = jj(ista)
              zlow = height_to_zcoord2(elev_s(ista),heights_3d,ni,nj,nk,
     &                              iii,jjj,istatus)
              if(istatus .ne. 1)then
                 write(6,*) 
     &           ' Error in height_to_zcoord2 in interp_to_sfc', istatus       
                 write(6,*) iii,jjj,zlow,elev_s(ista),
     &                                 (heights_3d(iii,jjj,k),k=1,nk)     
                 return
              endif
c
              k = max(zlow, 1.)
c 
              t_in = t_3d(iii,jjj,k)   !K          -- temp
              q_in = q_3d(iii,jjj,k)   !kg/kg      -- specific hum in
              p_in = 1.                !dummy arg  -- nothing needed here
c
              call compute_sfc_bgfields(bgmodel,ni,nj,nk,iii,jjj,
     &               k,elev_s(ista),heights_3d,t_3d,p_1d,
     &               q_3d,t_ref,p_in,t_in,q_in)
c
              t_cur(ista)    = t_in !K   -- temp
              td_cur(ista)   = q_in !K   -- Td returned here
              psfc_cur(ista) = p_in !Pa  -- calc sfc p
c
              call interp_to_elevation(elev_s(ista),iii,jjj,
     &                u_3d,heights_3d,ni,nj,nk,badflag,u_cur(ista))
              call interp_to_elevation(elev_s(ista),iii,jjj,
     &                v_3d,heights_3d,ni,nj,nk,badflag,v_cur(ista))
c
c.....  Convert units to what LAPS is using.
c
              t_cur(ista) = ((t_cur(ista) - 273.17) * 1.8) + 32.   !K to F
              td_cur(ista) = ((td_cur(ista) - 273.17) * 1.8) + 32. !K to F
              u_cur(ista) = u_cur(ista) * 1.94254                  !m/s to kt
              v_cur(ista) = v_cur(ista) * 1.94254                  !m/s to kt
              psfc_cur(ista) = psfc_cur(ista) * 0.01               !Pa to mb
c
c.....  Use the interpolated values to calculate a MSL pressure and
c.....  altimeter.  Note that the routine used is the same as that
c.....  in the LAPS surface analysis, which may or may not be giving
c.....  good values in all areas (e.g., high terrain, very cold areas).
c     
              call reduce_p(t_cur(ista),td_cur(ista),psfc_cur(ista),
     &                   elev_s(ista),lapse_t,lapse_td,pmsl_cur(ista),
     &                   0.,badflag)
c
c.....  Altimeter is defined as being 10 feet (3.048 meters) above the 
c.....  station elevation.  So, we'll add that height and use the
c.....  same routine.  There are other differences in the way that
c.....  MSL and altimeter are calculated, but we'll ignore them here. 
c.....  The bias introduced should be small.
c
              elev_alt = elev_s(ista) + 3.048 !meters
              call reduce_p(t_cur(ista),td_cur(ista),psfc_cur(ista),
     &                   elev_alt,lapse_t,lapse_td,alt_cur(ista),
     &                   0.,badflag)
c
           else  !station outside of grid
c
              t_cur(ista)    = badflag
              td_cur(ista)   = badflag
              u_cur(ista)    = badflag
              v_cur(ista)    = badflag
              psfc_cur(ista) = badflag
              pmsl_cur(ista) = badflag
              alt_cur(ista)  = badflag
c
           endif
c
        enddo !ista
c
c.....  Now do all that again, but for the previous cycle.
c
        i4time_pre = i4time_cur - laps_cycle_time
c
c.....  Get the 3d fields for the previous cycle.
c
        print *,' Getting LAPS pressure levels'
        call get_pres_3d(i4time_pre,ni,nj,nk,p_3d,istatus)
        do k=1,nk  !for now...grab a column of pressures
           p_1d(k) = p_3d(1,1,k) * 0.01  !conv from Pa to mb too...
        enddo !k
c
        print *,' Getting LGA 3-d Heights'
        var_lga = 'HT '
        call get_modelfg_3d(i4time_pre,var_lga,ni,nj,nk,heights_3d,
     &                      istatus)
c
        if(istatus .ne. 1)  then
           print *,' LGA HT not available. '
           go to 800
        else

        endif
c
        print *,' Getting LGA 3-d Temperatures'
        var_lga = 'T3 '
        call get_modelfg_3d(i4time_pre,var_lga,ni,nj,nk,t_3d,istatus)
c
        if(istatus .ne. 1)  then
           print *,' LGA T not available. '
           go to 800
        else

        endif
c
        print *,' Getting LGA 3-d Specific Humidity'
        var_lga = 'SH '  ! specific humidity 
        call get_modelfg_3d(i4time_pre,var_lga,ni,nj,nk,q_3d,istatus)
c
        if(istatus .ne. 1)  then
           print *,' LGA SH not available. '
           go to 800
        else

        endif
c
        print *,' Getting LGA 3-d U-wind component'
        var_lga = 'U3 '
        call get_modelfg_3d(i4time_pre,var_lga,ni,nj,nk,u_3d,istatus)
c
        if(istatus .ne. 1)  then
           print *,' LGA U not available. '
           go to 800
        else

        endif
c
        print *,' Getting LGA 3-d V-wind component'
        var_lga = 'V3 '
        call get_modelfg_3d(i4time_pre,var_lga,ni,nj,nk,v_3d,istatus)
c
        if(istatus .ne. 1)  then
           print *,' LGA V not available. '
           go to 800
        else

        endif
c
c.....  Now, for the previous cycle, interpolate the 3-d variables to the
c.....  proper location of each station.
c
        do ista=1,n_obs
c
c.....  Is the station within the LAPS grid?  If so, continue.  Otherwise,
c.....  set the previous variables to badflag.
c
           if(ii(ista).ge.1 .and. ii(ista).le.ni .and.
     &        jj(ista).ge.1 .and. jj(ista).le.nj) then
c
c.....  Find the vertical index for the station elevation.
c
              iii = ii(ista)
              jjj = jj(ista)
c
              zlow = height_to_zcoord2(elev_s(ista),heights_3d,ni,nj,nk,
     &                              iii,jjj,istatus)
              if(istatus .ne. 1)then
                 write(6,*) 
     &           ' Error in height_to_zcoord2 in interp_to_sfc', istatus       
                 write(6,*) iii,jjj,zlow,elev_s(ista),
     &                               (heights_3d(iii,jjj,k),k=1,nk)     
                 return
              endif
c
              k = max(zlow, 1.)
c 
              t_in = t_3d(iii,jjj,k)    !K          -- temp
              q_in = q_3d(iii,jjj,k)    !kg/kg      -- specific hum in
              p_in = 1.                 !dummy arg  -- nothing needed here
c
              call compute_sfc_bgfields(bgmodel,ni,nj,nk,iii,jjj,
     &               k,elev_s(ista),heights_3d,t_3d,p_1d,
     &               q_3d,t_ref,p_in,t_in,q_in)
c
              t_1(ista)    = t_in      !K   -- temp
              td_1(ista)   = q_in      !K   -- Td returned here
              psfc_1(ista) = p_in      !Pa  -- calc sfc p
c
              call interp_to_elevation(elev_s(ista),iii,jjj,
     &                   u_3d,heights_3d,ni,nj,nk,badflag,u_1(ista))
              call interp_to_elevation(elev_s(ista),iii,jjj,
     &                   v_3d,heights_3d,ni,nj,nk,badflag,v_1(ista))
c
c.....  Convert units to what LAPS is using.
c
              t_1(ista) = ((t_1(ista) - 273.17) * 1.8) + 32.   !K to F
              td_1(ista) = ((td_1(ista) - 273.17) * 1.8) + 32. !K to F
              u_1(ista) = u_1(ista) * 1.94254                  !m/s to kt
              v_1(ista) = v_1(ista) * 1.94254                  !m/s to kt
              psfc_1(ista) = psfc_1(ista) * 0.01               !Pa to mb
c
c.....  Use the interpolated values to calculate a MSL pressure and
c.....  altimeter.  Note that the routine used is the same as that
c.....  in the LAPS surface analysis, which may or may not be giving
c.....  good values in all areas (e.g., high terrain, very cold areas).
c     
              call reduce_p(t_1(ista),td_1(ista),psfc_1(ista),
     &                   elev_s(ista),lapse_t,lapse_td,pmsl_1(ista),
     &                   0.,badflag)
c
c.....  Altimeter is defined as being 10 feet (3.048 meters) above the 
c.....  station elevation.  So, we'll add that height and use the
c.....  same routine.  There are other differences in the way that
c.....  MSL and altimeter are calculated, but we'll ignore them here. 
c.....  The bias introduced should be small.
c
              elev_alt = elev_s(ista) + 3.048 !meters
              call reduce_p(t_1(ista),td_1(ista),psfc_1(ista),
     &                   elev_alt,lapse_t,lapse_td,alt_1(ista),
     &                   0.,badflag)
c
           else  !station outside of grid
c
              t_1(ista)    = badflag
              td_1(ista)   = badflag
              u_1(ista)    = badflag
              v_1(ista)    = badflag
              psfc_1(ista) = badflag
              pmsl_1(ista) = badflag
              alt_1(ista)  = badflag
c
           endif
c

         enddo !ista

c
c.....  Now calculate the trends for each variable at each station.
c.....  Trend is defined as var(current) minus var(previous).
c.....  For now, if either estimate is bad, set the trend to zero
c.....  (which is the same as persistance).
c
         do i=1,n_obs
            if(t_cur(i).le.badflag .or. t_1(i).le.badflag) then
               dt(i) = 0.
            else
               dt(i) = t_cur(i) - t_1(i)                   !deg F
            endif
c
            if(td_cur(i).le.badflag .or. td_1(i).le.badflag) then
               dtd(i) = 0.
            else
               dtd(i) = td_cur(i) - td_1(i)                !deg F
            endif
c
            if(u_cur(i).le.badflag .or. u_1(i).le.badflag) then
               du(i) = 0.
            else
               du(i) = u_cur(i) - u_1(i)                   !kt
            endif
c
            if(v_cur(i).le.badflag .or. v_1(i).le.badflag) then
               dv(i) = 0.
            else
               dv(i) = v_cur(i) - v_1(i)                   !kt
            endif
c
            if(pmsl_cur(i).le.badflag .or. pmsl_1(i).le.badflag) then
               dpmsl(i) = 0.
            else
               dpmsl(i) = pmsl_cur(i) - pmsl_1(i)          !mb
            endif
c
            if(alt_cur(i).le.badflag .or. alt_1(i).le.badflag) then
               dalt(i) = 0.
            else
               dalt(i) = alt_cur(i) - alt_1(i)             !mb
            endif
         enddo !i
c
c..... That's all.
c
        return
c
c..... If missing the LGA from one or the other times, do something
c..... else to get the trends.  For now, just set them equal to
c..... zero (persistance forecast).
c
 800    continue
        print *,
     &    ' **WARNING. Problem getting LGA. Setting trends to zero.'
        do i=1,n_obs
           dt(i)    = 0.
           dtd(i)   = 0.
           du(i)    = 0.
           dv(i)    = 0.
           dpmsl(i) = 0.
           dalt(i)  = 0.
        enddo !i
        return
c
        end
c
c
        subroutine interp_to_elevation(pt_elev,pt_i,pt_j,field_3d,
     &                         heights_3d,ni,nj,nk,badflag,interp)
c
c=================================================================
c     
c     Routine to interpolate a 3-d field to a single point in 3-d
c     space.  The point can be the actual "surface" (the 
c     topography), or any other point within the grid.
c
c     Based on code in the LAPS wind analysis, Steve Albers, FSL,
c     and the "interp_to_sfc" routine.
c
c     Original: 01-06-00   Peter A. Stamus, NOAA/FSL
c     Changes:
c
c=================================================================
c
        real pt_elev                 !elevation of desired point
        real interp                  !output interpolated value
        real field_3d(ni,nj,nk),     !3d array of desired variable
     &       heights_3d(ni,nj,nk)    !3d height field
        integer pt_i, pt_j           !i,j of desired point
c
cc      i_sfc_bad = 0
c
c..... Interpolate from the 3-d grid to the elevation at the point.
c
        zlow = height_to_zcoord2(pt_elev,heights_3d,ni,nj,nk,
     &                           pt_i,pt_j,istatus)
        if(istatus .ne. 1)then
           write(6,*) 
     &       ' Error in height_to_zcoord2 in interp_to_sfc', istatus       
           write(6,*) pt_i,pt_j,zlow,pt_elev,(heights_3d(i,j,k),k=1,nk)     
           return
        endif
c
        klow = max(zlow, 1.)
        khigh = klow + 1
        fraclow = float(khigh) - zlow
        frachigh = 1.0 - fraclow
c
        if( field_3d(pt_i,pt_j,klow)  .eq. badflag .or.
     &      field_3d(pt_i,pt_j,khigh) .eq. badflag) then
           write(6,3333)
 3333      format(' Warning: cannot interpolate to the point') 
c          i_sfc_bad = 1
           interp = badflag
        else
           interp = field_3d(pt_i,pt_j,klow ) * fraclow  +  
     &              field_3d(pt_i,pt_j,khigh) * frachigh
        endif
c
c..... That's all.
c
        return
        end
c
c
	Subroutine perturb(ta,tb,dta,imax,m,offset,onoff)     
c       
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c     Notes:
c
c*********************************************************************
c
	parameter(badflag=-99.9)
	real ta(m),tb(m),dta(m)
	integer onoff
c
	if (onoff.eq.1) then
	   do i=1,imax
	      if(ta(i).ne.badflag.and.tb(i).ne.badflag) then
		 dta(i)=ta(i)-tb(i)+offset
	      else
		 dta(i)=badflag
	      endif
	   enddo !i
	else
	   do i=1,imax
	      ta(i)=tb(i)+(dta(i)-offset)
	   enddo !i
	endif
c
	return
	end
c
c        
	Subroutine weights(mwt,dwt,ua,va,theta,lat,lon,icycle,imax,m,
     &                      hcutoff,vcutoff)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  Peter Stamus, NOAA/FSL
c          New version (minor changes from old).
c
c     Notes:
c
c*********************************************************************
c
	common tab(10000),pi,re,rdpdg,reorpd
	real mwt(m,m),dwt(m,m),ua(m),va(m),theta(m),lat(m),lon(m)
        integer icycle
c       
        cycle = float( icycle )
	expsclr=1000.
c       hcutoff= specify  !dist in m for wt e**-1
c       vcutoff=specify	 !deg C for vertical wt to go to e**-1
	C3=2.		 !wt on wind factor
	c1=sqrt(expsclr)/hcutoff**2
	c1m=sqrt(expsclr)/hcutoff**2
	c2=hcutoff**2/vcutoff**2
	c2m=hcutoff**2/vcutoff**2
c
c corrects for wind component along station/station vector
c
	do i=1,imax
	do j=1,imax
	   ang=(lat(i)+lat(j))*.5*rdpdg
	   dy=(lat(j)-lat(i))*reorpd
	   dx=(lon(j)-lon(i))*reorpd*cos(ang)
	   rr=(dx*dx+dy*dy)+1.e-20
	   r=sqrt(rr)
	   uav=dx*(ua(i)+ua(j))*.5/r
	   vav=dy*(va(i)+va(i))*.5/r
	   dxx=dx+uav*cycle/c3
	   dyy=dy+vav*cycle/c3
	   dz=theta(j)-theta(i)
	   iii=int(c1*((dx*dx+dy*dy) +c2*dz*dz))+1
	   iiii=int(c1m*((dxx*dxx+dyy*dyy) +c2m*dz*dz))+1
	   if(iii.gt.10000) iii=10000
	   if(iiii.gt.10000) iiii=10000 
	   dwt(i,j)=tab(iii)
	   mwt(i,j)=tab(iiii)
	enddo !j
	enddo !i
c
c normalize dwt,mwt
c
	do i=1,imax
	   sum=0.
	   summ=0.
	   do j=1,imax
	      sum=dwt(i,j)+sum
	      summ=mwt(i,j)+summ
	   enddo !on j
	   do j=1,imax
	      dwt(i,j)=dwt(i,j)/sum
	      mwt(i,j)=mwt(i,j)/summ
	   enddo !on j
	enddo !on i
c
	return
	end
c
c
        subroutine project(ord,y,monster,nv,nvar,maxsta,m,ncycles,
     &              nn,atime,it,icyc,oberr,badthr,bmonster,ihr)     
c
c*********************************************************************
c
c     This routine uses Taylor series of order 'ord' to forecast 
c     the variable 'nv' one cycle.  
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          Completely new version.
c       Bias correction (bmonster) added to self trend. Includes
c      a cycle dependent bias correction based on z time (ihr)
c      and model performance
c
c     Notes:
c          For the time being (14 Dec 1999) 'ord' is limited to 2.
c
c*********************************************************************
c
        parameter(badflag=-99.9)
        common tab(10000),pi,re,rdpdg,reorpd
        integer ord
        real y(m), monster(m,ncycles,nvar),bmonster(m,24,nvar)
        character*24 atime
c
c.....  Apply a filter to each derivative
c
        em1=exp(-1.)
        em1s=em1*em1
        do i=1,maxsta
           a=1.
           b=1.
           if(monster(i,2,nv).eq.badflag) a=0.
           if(monster(i,3,nv).eq.badflag) b=0.
           sum0 = monster(i,1,nv)*(1.+a*em1+b*em1s*.5) + 
     &            a*monster(i,2,nv)*(-a*em1-b*em1s) + 
     &            monster(i,3,nv)*(b*em1s*.5)    
           y(i)=sum0-bmonster(i,ihr,nv)
        enddo !on i
c
        return 
        end

