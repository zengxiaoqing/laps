
!       RADAR
        integer*4  MAX_RADAR_FILES
        parameter (MAX_RADAR_FILES = 20000)     

!       Background
        integer*4  MAX_BACKGROUND_FILES
        parameter (MAX_BACKGROUND_FILES = 6000)     

!       Ingest
        integer*4  MAX_INGEST_FILES
        parameter (MAX_INGEST_FILES = 10000)     

        integer      MAXBGMODELS
        parameter    (MAXBGMODELS=50) 

!       Maximum number of LAPS grid levels
        integer*4 MAX_LVLS
        parameter (MAX_LVLS=150)

!       Century time cutoff. LAPS will be set to run with filenames having
!       two digits for the year between the actual years of 'iyear_earliest'
!       and 'iyear_earliest+99'. Valid values are from 1901 to 1999. Note that
!       other time constraints may limit the usable time span of LAPS to less
!       than a century.
        integer*4 iyear_earliest
        parameter (iyear_earliest=1950)
