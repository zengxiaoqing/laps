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

