!       Declarations for surface obs

	real lat_s(maxsta), lon_s(maxsta), elev_s(maxsta)
	real t_s(maxsta), t_ea(maxsta), max24t(maxsta), min24t(maxsta)
        real td_s(maxsta), td_ea(maxsta), rh_s(maxsta), rh_ea(maxsta)

        real dd_s(maxsta), ddg_s(maxsta), dd_ea(maxsta)
        real ff_s(maxsta), ffg_s(maxsta), ff_ea(maxsta)

        real alt_s(maxsta), alt_ea(maxsta), delp(maxsta)
	real pstn_s(maxsta), pmsl_s(maxsta), p_ea(maxsta)
        real pred_s(maxsta)      

	real store_hgt(maxsta,5) 

        real vis_s(maxsta), vis_ea(maxsta)
        real solar_s(maxsta), solar_ea(maxsta)

        real sfct(maxsta), sfct_ea(maxsta)
        real sfcm(maxsta), sfcm_ea(maxsta)
        real pcp1(maxsta), pcp3(maxsta), pcp6(maxsta), pcp24(maxsta)
        real snow(maxsta), snow_ea(maxsta), pcp_ea(maxsta)

	character atime_s*24
	character store_amt(maxsta,5)*4
        character stations(maxsta)*20, provider(maxsta)*11,stn20*20      
        character reptype(maxsta)*6, autostntype(maxsta)*6
        character wx_s(maxsta)*25 

	integer kloud_s(maxsta), obstime(maxsta)
        integer wmoid(maxsta), delpch(maxsta)
	integer ii(maxsta), jj(maxsta)
