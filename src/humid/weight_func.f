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
       subroutine weight_func (tau,p_model,nm,w_dm,p_dm,nw)

       real tau (nm)
       real p_model(nm)
       integer nm
       real w_dm(nw), p_dm(nw) !wieghts and pressures (derivative model)
       integer nw

       integer i
       real delta_log_p(100)  ! 
       real sum
       integer first_time
       data first_time /1/


c  check array dimensions

       if(nm-nw .ne. 1) then
           write(6,*) 'check array dimensions on tau and w_model'
           write(6,*) 'they should differ by 1, but they dont'
           write(6,*) 'Program halted'
           stop
       endif

c check first time

        if(first_time.eq.1) then !first time
           first_time = 0
           do i = 1,nw
        p_dm(i) = sqrt( p_model(i) * p_model(i+1) )  ! central pressure
        delta_log_p(i) = log(p_model(i)) - log(p_model(i+1))
           enddo
        endif

       
c compute central derivative pressures and delta log pressure

        sum = 0.0

        do i =1,nw

        w_dm(i) =  - tau(i) + tau(i+1)
        w_dm(i) = w_dm(i) / delta_log_p(i)
        sum = sum + w_dm(i) 

        enddo  

c normalize

        do i = 1,nw
        w_dm(i) = w_dm(i) / sum
        enddo


        return

        end

  
