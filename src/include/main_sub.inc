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
        real*4 topo(NX_L,NY_L)                

!       Profiler Stuff
        real lat_pr(MAX_PR)                   
        real lon_pr(MAX_PR)                   


!       Profiler Observations

        integer nlevels_obs_pr(MAX_PR)        
        real ob_pr_ht(MAX_PR,NZ_L)            
        real ob_pr_di(MAX_PR,NZ_L)            
        real ob_pr_sp(MAX_PR,NZ_L)            
        real ob_pr_u (MAX_PR,NZ_L) ! Vertically interpolated Profiler wind
        real ob_pr_v (MAX_PR,NZ_L) ! Vertically interpolated Profiler wind
        real ob_pr_r (MAX_PR,NZ_L) ! Vertically interpolated Profiler wind
        real ob_pr_t (MAX_PR,NZ_L) ! Vertically interpolated Profiler wind


!       Mesonet

        integer meso_i(N_MESO) ! X mesonet coordinates                
        integer meso_j(N_MESO) ! Y mesonet coordinates                
        integer meso_k(N_MESO) ! Z mesonet coordinates                
        real    meso_u(N_MESO) ! u mesonet component                  
        real    meso_v(N_MESO) ! v mesonet component                  

!       SAO

        integer sao_i(N_SAO) ! X sao coordinates                      
        integer sao_j(N_SAO) ! Y sao coordinates                      
        integer sao_k(N_SAO) ! Z sao coordinates                      
        real    sao_u(N_SAO) ! u sao component                        
        real    sao_v(N_SAO) ! v sao component                        

!       Pireps

        integer pirep_i(N_PIREP) ! i pirep coordinates                
        integer pirep_j(N_PIREP) ! j pirep coordinates                
        integer pirep_k(N_PIREP) ! k pirep coordinates                
        real    pirep_u(N_PIREP) ! u pirep component                  
        real    pirep_v(N_PIREP) ! v pirep component                  


!       Laps Analysis Grids
        real grid_laps_wt(NX_L,NY_L,NZ_L)                             
        real grid_laps_u(NX_L,NY_L,NZ_L)                              
        real grid_laps_v(NX_L,NY_L,NZ_L)                              

!       Data flags for LAPS analysis
        logical L_profiler
     1  ,laps_rotate_winds,l_derived_output,l_3pass
        parameter (L_profiler = .true.)
        parameter (laps_rotate_winds = .false.)! True North Winds
        parameter (l_derived_output = .true.)
        parameter (l_3pass = .true.)

        logical L_correct_unfolding
        parameter (L_correct_unfolding = .true.)
