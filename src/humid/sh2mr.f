cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis



        subroutine sh2mr (sh,mr,istatus)


c       $log: sh2mr.for,v $
c revision 1.1  1996/08/30  20:57:36  birk
c initial revision
c

c       given the specific humidity (g/g)
c       this routine computes the mixing ratio (g/g)

        real sh,mr
        real epsilon
        integer istatus

        epsilon = 0.622 !(ratio of molecular weight of water ~ 18,
c                       over dry air~29)

        istatus  = 1

        if (sh/epsilon .eq. 1.) then
                write(6,*) 'error in input to routine sh2mr'
                write(6,*) 'divide by zero has occurred'
                istatus = 0
        endif


        if (sh/epsilon .gt. 1.) then
                write(6,*) 'error in input to routine sh2mr'
                write(6,*) 'impossible condition has occurred'
                write(6,*) 'sh has value of ', sh
                istatus = 0
        endif

        mr = sh / (1.-sh/epsilon)

        return
        end





