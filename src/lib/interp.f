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
        subroutine interp (li,l1,l2,v1,v2,output)
c       simple linear interpolator
c       birkenheuer 0ct 1990

c               l1...........l2
c               v1...........v2
c                   li
c                    .
c                    .
c                   vi   <-----computed value (output)

        implicit none

        real li,l1,l2,v1,v2,output, range
        
        range = l2-l1
        if (abs(range/li) .le.  1.0e-4) then
           output = (v1+v2)*0.5
        else
           output = (li-l1)/range * v2 + (l2-li)/range*v1
        endif
        return
        end
