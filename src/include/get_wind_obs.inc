cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
!       LAPS Grid Dimensions
                                                   
        real*4 lat(NX_L,NY_L)                      
        real*4 lon(NX_L,NY_L)                      

!       Sectorized Uniform Wind

!       Profiler Stuff
        real lat_pr(MAX_PR)                        
        real lon_pr(MAX_PR)                        


!       Profiler Observations

        integer nlevels_obs_pr(MAX_PR)             
        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)    
        real ob_pr_di_obs(MAX_PR,MAX_PR_LEVELS)    
        real ob_pr_sp_obs(MAX_PR,MAX_PR_LEVELS)    
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS)     
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS)     
        real ob_pr_ht(MAX_PR,NZ_L)                 
        real ob_pr_di(MAX_PR,NZ_L)                 
        real ob_pr_sp(MAX_PR,NZ_L)                 
        real ob_pr_u (MAX_PR,NZ_L) ! Vertically interpolated Profiler wind
        real ob_pr_v (MAX_PR,NZ_L) ! Vertically interpolated Profiler wind
        real ob_pr_r (MAX_PR,NZ_L) ! Vertically interpolated Profiler wind
        real ob_pr_t (MAX_PR,NZ_L) ! Vertically interpolated Profiler wind


!       Barnes Profiler analysis

        real weights_pr(NX_L,NY_L)                                      
        real max_weights_pr(NX_L,NY_L,NZ_L)                             

        real sum_pr_u (NX_L,NY_L,NZ_L)                                  
        real sum_pr_v (NX_L,NY_L,NZ_L)                                  
        real sum_wt_pr (NX_L,NY_L,NZ_L)                                 

!       Laps Analysis Grids
        real grid_laps_wt(NX_L,NY_L,NZ_L)                               
        real grid_laps_u(NX_L,NY_L,NZ_L)                                
        real grid_laps_v(NX_L,NY_L,NZ_L)                                

!       Data flags for LAPS analysis
        logical L_profiler
        parameter (L_profiler = .true.)
