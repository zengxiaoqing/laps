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

  
