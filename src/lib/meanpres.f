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
cdis
c
c
        subroutine mean_pres(num_sfc,p_s,pbar)
c
c*******************************************************************************
c
c       Routine to calculate the mean surface or reduced pressure.
c
c       Changes:
c               P. A. Stamus    12-21-88        Original
c
c       Inputs/Outputs:
c
c          Variable     Var Type     I/O     Description
c         ----------   ----------   -----   -------------
c          num_sfc         I          I      Number of surface stations.
c          p_s             RA         I      Surface (or reduced) pressure.
c          pbar            R          O      Mean pressure of all the stations.
c
c       User Notes:
c
c       1.  Units are not changed in this routine.
c
c*******************************************************************************
c
        real p_s(num_sfc)
c
        sump = 0.
        cntp = 0.
c
        do 1 i=1,num_sfc
          if(p_s(i) .le. 0.) go to 1
          sump = sump + p_s(i)
          cntp = cntp + 1.
1       continue
c
        pbar = sump / cntp
c
        return
        end
