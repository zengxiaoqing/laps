!dis    Forecast Systems Laboratory
!dis    NOAA/OAR/ERL/FSL
!dis    325 Broadway
!dis    Boulder, CO     80303
!dis 
!dis    Forecast Research Division
!dis    Local Analysis and Prediction Branch
!dis    LAPS 
!dis 
!dis    This software and its documentation are in the public domain and 
!dis    are furnished "as is."  The United States government, its 
!dis    instrumentalities, officers, employees, and agents make no 
!dis    warranty, express or implied, as to the usefulness of the software 
!dis    and documentation for any purpose.  They assume no responsibility 
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis     technical support to users.
!dis    
!dis    Permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  All modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making 
!dis    the modifications.  If significant modifications or enhancements 
!dis    are made to this software, the FSL Software Policy Manager  
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis 
!dis 
!dis 
!dis 
!dis 
!dis 

!       Data flags for LAPS analysis
        logical L_profiler,l_derived_output,l_3pass
        parameter (L_profiler = .true.)
        parameter (l_derived_output = .true.)
        parameter (l_3pass = .true.)

        logical L_correct_unfolding
        parameter (L_correct_unfolding = .true.)

!       If this is true, then we rotate the first guess and obs arrays to
!       Grid North U+V just prior to the call to 'laps_anl'. 
!       The winds are rotated back to True North after the 'laps_anl' call 
!       and before the LW3 file is written out (depending on l_grid_north_out).

        logical l_grid_north
        logical l_grid_north_bkg, l_grid_north_anal, l_grid_north_out       
        parameter (l_grid_north      = .true.)
        parameter (l_grid_north_bkg  = .true.)
        parameter (l_grid_north_anal = .true.)
        parameter (l_grid_north_out  = .true.)

