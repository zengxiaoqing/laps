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


        subroutine plotfield (data_anal,igrid,jgrid)


c       $log: plotfield.for,v $
c revision 1.1  1996/08/30  20:48:28  birk
c initial revision
c

c ported to the unix environment 1 oct 1993

        implicit none

        integer igrid,jgrid


        real data_anal(igrid,jgrid)

        real finc




        integer ioffp,ioffm
        real spval,epsval,cntmin,cntmax,cntint

        common /conre1/ioffp,spval,epsval,cntmin,cntmax,cntint,ioffm

        ioffp = 1   ! says that the missing value exists.
        spval = 1.0e+37  ! sets the missing value flag for plots
        ioffm = 0 !omit the conrec message on the plot -- doesn't work





c        call opngks


        ioffm = 1




        call conrec (data_anal,igrid,igrid,jgrid,0.,0.,finc,-1,0,0)
        call frame

c        call clsgks


        end
