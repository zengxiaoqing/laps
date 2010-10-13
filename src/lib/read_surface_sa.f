cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
c
c
        subroutine read_surface_sa(i4time,maxsta,                     ! I
     &   n_obs_b,stn,reptype,atype,                                   ! O
     &   lat,lon,elev,wx,t,td,                                        ! O
     &   kloud,store_amt,store_hgt,                                   ! O
     &   solar,solar_ea,obstime,istatus)                              ! O
c
cdoc    This routine calls 'read_surface_data' and is used primarily to read
cdoc    in cloud info along the lines of the arrays in the "old" LSO format.
cdoc    This is called only from the cloud analysis at present.
c
c
        real          badflag
c
c.....  Input arrays (for new format LSO)
c
	real lat(maxsta), lon(maxsta), elev(maxsta)
	real t(maxsta), t_ea(maxsta), max24t(maxsta), min24t(maxsta)
	real td(maxsta), td_ea(maxsta), rh(maxsta), rh_ea(maxsta)
	real dd(maxsta), ddg(maxsta), dd_ea(maxsta)
	real ff(maxsta), ffg(maxsta), ff_ea(maxsta)
	real alt(maxsta), alt_ea(maxsta), delp(maxsta)
	real pstn(maxsta), pmsl(maxsta), p_ea(maxsta)
	real vis(maxsta), vis_ea(maxsta)
	real solar(maxsta), solar_ea(maxsta)
	real sfct(maxsta), sfct_ea(maxsta)
	real sfcm(maxsta), sfcm_ea(maxsta)
	real pcp1(maxsta), pcp3(maxsta), pcp6(maxsta), pcp24(maxsta)
	real snow(maxsta), snow_ea(maxsta), pcp_ea(maxsta)
	real store_hgt(maxsta,5)
c
	integer i4time, wmoid(maxsta), jstatus
	integer time(maxsta), delpch(maxsta), kloud(maxsta)
c
	character filetime*9, infile*256 
	character stations(maxsta)*20, provider(maxsta)*11
        character reptype(maxsta)*6, atype(maxsta)*6
	character wx_in(maxsta)*25, store_amt(maxsta,5)*4
c
c.....  Output arrays (as old format LSO if different)
c
        Integer   obstime(maxsta)
c
        Character   atime*24,stn(maxsta)*3     
        character   store_emv(maxsta,5)*1, wx(maxsta)*8
c
c
c.....  Start here.  Set the status to nothing, zero out the cloud storage
c.....  and character arrays.
c
        call get_sfc_badflag(badflag,istatus)

        istatus = 0
        ibadflag = int(badflag)
c
        do j=1,5
        do i=1,maxsta
          store_emv(i,j) = ' '
          store_amt(i,j) = '    '
          store_hgt(i,j) = badflag
        enddo !i
        enddo !j
c
        do i=1,maxsta
           stn(i)(1:3) = '   '
           wx(i)(1:8) = '        '
        enddo !i


        call make_fnam_lp(i4time,filetime,istatus)

        if(.false.)then
	   call read_surface_data(i4time,atime,n_obs_g,n_obs_b,time,
     &    wmoid,stations,provider,wx_in,reptype,atype,lat,lon,
     &    elev,t,td,rh,dd,ff,ddg,ffg,alt,pstn,pmsl,delpch,delp,vis,
     &    solar,
     &    sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kloud,max24t,min24t,t_ea,
     &    td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,sfct_ea,
     &    sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,
     &    maxsta,jstatus)
c

        else
           call read_cloud_obs(i4time,maxsta,                      ! I
     &      n_obs_b,stations,reptype,                              ! O
     &      atype,                                                 ! O
     &      lat,lon,elev,wx_in,t,td,vis,solar,                     ! O
     &      kloud,store_amt,store_hgt,obstime,jstatus)             ! O

        endif

        if(jstatus .ne. 1 .and. jstatus .ne. -1) then
           print *,' ERROR: No valid LSO file found for ', filetime       
           istatus = -1
           return
        endif
c
c.....  Shuffle data for the differences between old and new formats.
c.....  Now the station data.
c
        do i=1,n_obs_b
           wx(i)(1:8) = wx_in(i)(1:8)
c
           if(reptype(i)(1:4) .eq. 'LDAD') then
              if(provider(i)(1:4) .eq. 'CDOT') then
                 stn(i)(1:1) = stations(i)(1:1)
                 stn(i)(2:2) = stations(i)(3:3) ! skip the '-'
                 if(ichar(stations(i)(4:4)) .lt. 32) then
                    stn(i)(3:3) = ' '
                 else
                    stn(i)(3:3) = stations(i)(4:4)
                 endif
              else
                 stn(i)(1:3) = stations(i)(1:3)
              endif

           else
              call right_justify(stations(i))
              stn(i)(1:3) = stations(i)(18:20)

           endif
c
        enddo !i
c
c..... End of data gathering. Let's go home...
c
        istatus = 1             ! everything's ok...
        print *, ' Normal completion of new READ_SURFACE_SA'
c
        return
        end

c
c
        subroutine read_cloud_obs(i4time,maxsta,                      ! I
     &   n_obs_b_out,stations_out,reptype_out,                        ! O
     &   autostntype_out,                                             ! O
     &   lat_s_out,lon_s_out,elev_s_out,wx_s_out,t_s_out,td_s_out,    ! O
     &   vis_s_out,solar_s_out,                                       ! O
     &   kloud_s_out,store_amt_out,store_hgt_out,obstime_out,istatus) ! O
c
c       The argument list is or should be consistent with 'read_sfc.inc' except
c       that a duplicate subset of '_out' arrays are used
c
cdoc    This routine calls 'read_surface_data' and is used primarily to read
cdoc    in cloud info from the current LSO file and potentially SYNOP obs from
cdoc    a wider time window. The initial design is to call this from the
cdoc    'read_surface_sa' routine.

        include 'read_sfc.inc'
c
c.....  Output arrays (duplicate declarations with _out suffix)
c
	real lat_s_out(maxsta), lon_s_out(maxsta), elev_s_out(maxsta)
	real t_s_out(maxsta), td_s_out(maxsta), vis_s_out(maxsta)
        real solar_s_out(maxsta)

	real store_hgt_out(maxsta,5) 

	character store_amt_out(maxsta,5)*4
        character stations_out(maxsta)*20
        character reptype_out(maxsta)*6, autostntype_out(maxsta)*6
        character wx_s_out(maxsta)*25 

	integer kloud_s_out(maxsta), obstime_out(maxsta)

c       End of output arrays

        real          badflag
c
	character filetime*9
c
c.....  Start here.  Set the status to nothing, zero out the cloud storage
c.....  and character arrays.
c
        call get_sfc_badflag(badflag,istatus)

        istatus = 0
        ibadflag = int(badflag)
c
c.....  Figure out the i4time and call the read routine.
c
        call read_surface_data(i4time,atime_s,n_obs_g,n_obs_b, !regular LSO
     &         obstime,wmoid,stations,provider,wx_s,reptype,autostntype,       
     &         lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &         alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &         pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &         td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &         sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,
     &         maxsta,jstatus)
c
        call make_fnam_lp(i4time,filetime,istatus)
        if(jstatus .ne. 1 .and. jstatus .ne. -1) then
           print *,' ERROR: No valid LSO file found for ', filetime       
           istatus = -1
           return
        endif
c
c.....  Place main cloud obs into output arrays

        n_obs_b_out = n_obs_b
        stations_out = stations
        reptype_out = reptype
        autostntype_out = autostntype
        lat_s_out = lat_s
        lon_s_out = lon_s
        elev_s_out = elev_s
        wx_s_out = wx_s
        t_s_out = t_s
        td_s_out = td_s
        vis_s_out = vis_s
        solar_s_out = solar_s
        kloud_s_out = kloud_s
        store_amt_out = store_amt
        store_hgt_out = store_hgt
        obstime_out = obstime

        write(6,*)' n_obs_b_out = ',n_obs_b_out       

        n_obs_b_ontime = n_obs_b_out
c
c.....  Figure out the i4time and call the read routine for SYNOPs.
c
        i4time_synop = (i4time / 10800) * 10800

        if(i4time_synop .ne. i4time .and. .false.)then

           write(6,*)' Reading 3 hourly sfc data for SYNOPs...'

           call read_surface_data(i4time_synop,atime_s,n_obs_g,n_obs_b, !regular LSO
     &         obstime,wmoid,stations,provider,wx_s,reptype,autostntype,       
     &         lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &         alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &         pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &         td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &         sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,      
     &         maxsta,jstatus)
c
           if(jstatus .ne. 1 .and. jstatus .ne. -1) then
              print *,' ERROR: No valid LSO file found for ', filetime       
              istatus = -1
              return
           endif

c
c.....     Place synop cloud obs into output arrays

           do i=1,n_obs_b
              if(reptype(i)(1:5) .eq. 'SYNOP') then

c                Check that station is not near other "on time" stations

                 n_obs_b_out = n_obs_b_out + 1
                 stations_out(n_obs_b_out) = stations(i)
                 reptype_out(n_obs_b_out) = reptype(i)
                 autostntype_out(n_obs_b_out) = autostntype(i)
                 lat_s_out(n_obs_b_out) = lat_s(i)
                 lon_s_out(n_obs_b_out) = lon_s(i)
                 elev_s_out(n_obs_b_out) = elev_s(i)
                 wx_s_out(n_obs_b_out) = wx_s(i)
                 t_s_out(n_obs_b_out) = t_s(i)
                 td_s_out(n_obs_b_out) = td_s(i)
                 vis_s_out(n_obs_b_out) = vis_s(i)
                 solar_s_out(n_obs_b_out) = solar_s(i)
                 kloud_s_out(n_obs_b_out) = kloud_s(i)
                 store_amt_out(n_obs_b_out,:) = store_amt(i,:)
                 store_hgt_out(n_obs_b_out,:) = store_hgt(i,:)
                 obstime_out(n_obs_b_out) = obstime(i)
              endif
           enddo !i

           write(6,*)' n_obs_b_out = ',n_obs_b_out       

        endif ! SYNOP time is different from current systime
c
c.....End of data gathering. Let's go home...
c
        istatus = 1             ! everything's ok...
        print *, ' Normal completion of new READ_CLOUD_OBS'
c
        return
        end

