 &surface_analysis
 use_lso_qc = 0,
 skip_internal_qc = 0,
 l_require_lso = .false.,
 itheta=-1,
 redp_lvl=1500.,
 del=3.24e6,
 gam=.0001,
 ak=1.e-6, 
 /
c..... This is the namelist for the LAPS surface analysis
c..... process (LSX).  Switches and similar things can go
c..... here, and are read at runtime (rather than requiring
c..... a recompile.
c
c..... Current switches and their default values:
c
c..... use_lso_qc = 0, (a "1" tells LSX to use the quality-
c.....                  controlled version of LSO (lso_qc),
c.....                  a "0" uses the normal LSO file. Note
c.....                  that setting this to "1" - using the
c.....                  QC'd LSO file - turns off the internal
c.....                  LSX QC). 
c.....                  
c
c..... skip_internal_qc = 0, (a "1" tells LSX to skip it's
c.....                        internal QC routine; a "0" uses
c.....                        it.  Note that the internal QC is potentially 
c.....                        used only when "use_lso_qc" is set to zero.)
c.....
c
c      l_require_lso: boolean flag to indicate whether an LSO file (with obs)
c                     is required in order to generate an LSX surface analysis
c
c.......... itheta,  
c
c.......... Surface Theta check:  Check to see that the surface potential
c..........     temperatures are not greater than the potential temperature
c..........     at an upper level.  Set this variable equal to the desired
c..........     upper level:
c
c..........     -1 = auto-set itheta (5/7) based on terrain (centr gridpt>1000m)
c..........      0 = No sfc theta check done
c..........      7 = Use 700 mb level
c..........      5 = Use 500 mb level
c
c..........     Recommended:  Use 700 mb most places, 500 mb over higher
c..........                   terrain areas (like Colorado).
c
c.......... redp_lvl (Pressure reduction):  The main pressure analysis that 
c               LAPS produces is a reduction to this elevation (m).  For 
c               example, the Colorado LAPS uses 1500 m, about the elevation of
c               Denver, since it is representative of the elevations in the 
c               region of interest.
c
c..... comments on del, gam, ak are in the surface code
c
c          wse is assumed to be low .50 so that winds are retained and 
c          most adjustment is in p ...
c          del is sqd error of wind/sqd error in eqn of motion
c          gam is sqd error of wind/sqd error of press
c          with mslp error at 50pa,eqn of motion residual (1m/s/hr)^2
c
c..... if del=0., then variational section would be skipped for (u,v,p)
c
c..... ANALYSIS QC THREHOLDS
c
c      The following list represents the default QC thresholds. These 
c      parameters can be added to this namelist to override the default values
c      set in the code...
c
c       QC parms: # of standard deviations 
c       bad_p  = 3.0 	        ! for reduced pressure
c       bad_mp = 4.0 	        ! for MSL pressure
c       bad_t  = 2.5 	        ! for temperature
c       bad_td = 2.0 	        ! for dewpoint
c       bad_u  = 4.0 	        ! for u-wind
c       bad_v  = 4.0 	        ! for v-wind
c       bad_th = 3.5 	        ! for theta
c       bad_the = 2.5           ! for theta-e
c       bad_vis = 500. 	        ! for visibility
c       bad_tb8 = 5.0 	        ! for tb8 Brightness temps.
c
