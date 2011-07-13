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
        function make_td (p,t,q,t_ref)

c   This function is designed to compute (td)
c   from basic variablesp (mb), t(c) and q (g/kg)
c   to give td in (c).  The reference temperature
c   t_ref (c) is used to describe the temperature
c   at which the liquid and ice phase change
c   occurs. Basically this routine inverts
c   routine ssh2 that discriminates between ice
c   and liquid phases.
c   
c   Author:  Dan Birkenheuer
c   
c   Date:    28 September 1993
c
c   QC checks added by S. Albers in 2010
c   
c   

        implicit none

        real p  !ambient pressure (mb)
        real t   !ambient temp (c)
        real q  !specific humidity in g/kg
        real t_ref  !phase reference temp (c)
        real pc,tc,qc,t_refc
        common /root/ pc,tc,qc,t_refc
        real internal_q

c       type functions
        real make_td
        real func
        real ssh2
        real rootfind


        external ssh2
        external func
        external rootfind


c   ---------------------------------------------
c   
c   
c   basic algorithm:
c   
c   compare the input q to the saturation ssh2
c   computed value and then iterate until 
c   ssh2 - q = 0.0 
c   this is a simple function zero solver.
c   
c   if the ssh2 value at the ambient temperature t
c   is less than the input q then we must generate
c   an error since we cannot super saturate in the
c   module.
c   
c   ---------------------------------------------
c   
c   
c   This module conforms to the other modules that
c   refer to phase of water in computing the vapor
c   pressure functions.  Related modules are:
c   
c       Make_rh
c       Make_ssh
c       ssh2
c   
c   These are the only functions that have been
c   assembled with the phase of the condensate in
c   mind.  These draw upon the "best" modules in
c   the mthermo library written by Tom Schlatter
c   and Don Baker.  The modules referenced for
c   the vapor computations take 2 things into
c   consideration.  1) accuracy (compared to
c   Smithsonian Tables.  and 2) computational
c   speed.
c   
c   first part of code is to check input q
c   against saturation q calling ssh2 in this

        internal_q = ssh2 (p,t,t,t_ref)
        if (q .gt. internal_q) then ! saturate only
                make_td = t
                q = internal_q
        return
        endif

c second test for too little a q for routine to work.
        if(q .lt. 7.214e-26) then
        print *, 'input value of q to module MAKE_TD was lower than 
     1 7.214e-26 ',q
        print *, ' assigning returned td the values of -199. c'
        make_td = -199.
        return
        endif

c now find the appropriate td that gives us internal_q = q

c assign common variables for func call
        if(p .gt. 1500. .or. p .le. 0.)then
            write(6,*)' p out of bounds in make_td ',p
        endif
        if(t .gt. 1000. .or. t .lt. -274.)then
            write(6,*)' t out of bounds in make_td ',t
        endif
        
        pc = p
        tc = t
        qc = q
        t_refc = t_ref
        make_td = rootfind (func,t,-199.,.0000001)

        return
        end


        function func (x)

c       This is a coding construct to call the
c       rootfinder func is a routine with only one
c       input and interfaces with the modules that
c       are needed in this call.

c       In this case, the input is the dewpoint
c       temperature that is varied by routine
c       rootfind until there is a zero found.
c       then that value is the dewpoint
c       temperature.  In order to interface with
c       the routines ssh2, we must pass the other
c       necessary variables through a common block.

c       coded by Dan Birkenheuer, September 28 1993

        implicit none

        real x
        real func
        real ssh2 !function type
        real pc,tc,qc,t_refc
        common /root/ pc,tc,qc,t_refc

        func = qc-ssh2(pc,tc,x,t_refc)  ! function will be zero when ssh2=q
        return
        end
