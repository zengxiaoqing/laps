 &surface_analysis
 use_lso_qc = 0,
 skip_internal_qc = 0,
 l_dev_ck = .true.,
 l_require_lso = .false.,
 itheta=-1,
 del=3.24e6,
 gam=.0001,
 ak=1.e-6, 
 bad_t=2.5,
 bad_td=1.7,
 bad_u=3.0,
 bad_v=3.0,
 bad_p=2.5,
 bad_mp=5.0,
 bad_th=3.5,
 bad_the=2.5,
 bad_vis=500,
 bad_tb8=5.0,
 bad_tgd_land=3.0,
 bad_tgd_water=3.0,
 thresh_t=30.,
 thresh_td=40.,
 thresh_mslp=10.,
 rms_wind=1.0,
 rms_temp=1.0,
 rms_dewpoint=1.2,
 rms_pres=1.0,
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
c      l_dev_ck: boolean flag that turns on or off the "Hartsough" QC checks
c                that are a subset of the "internal" QC
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
c..... comments on del, gam, ak are in the surface code
c
c          wse is assumed to be low .50 so that winds are retained and 
c          most adjustment is in p ...
c          del is sqd error of wind/sqd error in eqn of motion
c          gam is sqd error of wind/sqd error of press
c              high values favor the pressure, lower values favor the wind
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
c       bad_t   	        ! for temperature
c       bad_td  	        ! for dewpoint
c       bad_u   	        ! for u-wind
c       bad_v   	        ! for v-wind
c       bad_p   	        ! for reduced pressure
c       bad_mp  	        ! for MSL pressure
c       bad_th  	        ! for theta
c       bad_the                 ! for theta-e
c       bad_vis    	        ! for visibility
c       bad_tb8    	        ! for tb8 Brightness temps.
c
c       These parameters should be defined in this namelist as they aren't
c       initialized in the code
c       
c       QC parms: # of standard deviations 
c       bad_tgd_land            ! for ground temperature
c       bad_tgd_water           ! for water/sea surface temperature
c
c       QC parms: threshold checks
c       thresh_t                ! for temperature (deg F)
c       thresh_td               ! for dewpoint (deg F)
c       thresh_mslp             ! for MSL pressure (millibars)
c
c       The following parametsrs will adjust the RMS analysis fit to the obs,
c       given the idealized case of flat terrain. They are multiplied by the 
c       default instrument error (set at 1.5 deg F for T and Td, and 1.5kt for
c       the U and V wind components).
c
c       rms_wind                ! scaling factor for wind rms threshold
c       rms_temp                ! scaling factor for temperature rms threshold
c       rms_dewpoint            ! scaling factor for dewpoint rms threshold
c       rms_pres                ! scaling factor for MSLP rms threshold
