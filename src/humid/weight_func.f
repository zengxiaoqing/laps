cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
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

  
