 &surface_analysis
 use_lso_qc = 0,
 skip_internal_qc = 0,
 
 /
c
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
c.....                  that setting this to one--using the
c.....                  QC'd LSO file--turns off the internal
c.....                  LSX QC). 
c.....                  
c
c..... skip_internal_qc = 0, (a "1" tells LSX to skip it's
c.....                        internal QC routine; a "0" uses
c.....                        it.  Note that this is only used
c.....                        if "use_lso_qc" is set to zero.)
c.....

