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
C
C==========================================
C
        function ishow_timer()

        integer sec_elapsed

        character*24 atime

        common /timer/ i4time_start

!       save i4time_start

        i4time_now = i4time_now_gg()

        i4time_elapsed = i4time_now - i4time_start

        min_elapsed = i4time_elapsed / 60

        sec_elapsed = i4time_elapsed - min_elapsed * 60

        call cv_i4tim_asc_lp(i4time_now,atime,istatus)

        if(sec_elapsed .ge. 10)then
            write(6,1,err=99)min_elapsed,sec_elapsed,atime(1:20)
 1          format(1x,'Elapsed time -',i6,':',i2,5x,a20)
        elseif(sec_elapsed .ge. 0)then
            write(6,2,err=99)min_elapsed,sec_elapsed,atime(1:20)
 2          format(1x,'Elapsed time -',i6,':0',i1,5x,a20)
        else
            write(6,*)'Elapsed time (sec) = ',sec_elapsed
        endif

        ishow_timer = i4time_elapsed

        return

 99     write(6,*)' Error in ishow_timer: i4time_now, i4time_start = '
     1                  ,i4time_now,i4time_start
        write(6,*)' This should not happen, stopping program'
        stop

        return
        end
C
C==========================================
C
        function init_timer()

        common /timer/ i4time_start

!       save i4time_start

        i4time_start = i4time_now_gg()

        write(6,*)' Initializing elapsed timer at ',i4time_start

        init_timer = 1

        return
        end
C
C==========================================
C
        subroutine LIB$SHOW_TIMER(my_show_timer)
!	Linda Wharton 26 Oct 1998 commented declaration of istatus never used
C       integer*4 istatus
C
        write(6,10)
10      format('Called Lib$show_Timer')
C

C       istatus = system('date')
        i4time = i4time_now_gg()
        write(6,*)' TIMER: ',i4time

        return
        end
C
C==========================================
C
        FUNCTION LIB$INIT_TIMER()
C
c       WRITE(*,20)
20      format('Called LIB$INIT_TIMER')
C
        LIB$INIT_TIMER = 1

        return
        end
C
C==========================================
C
        FUNCTION lib$set_logical()
C
c       WRITE(*,20)
20      format('Called LIB$INIT_TIMER')
C
        lib$set_logical = 1

        return
        end
C
C==========================================
C
        FUNCTION ltest_log_gg()
C
c       WRITE(*,20)
20      format('Called LIB$INIT_TIMER')
C
        ltest_log_gg = 1

        return
        end
C
C==========================================
C
        FUNCTION lib$set_symbol()
C
c       WRITE(*,20)
20      format('Called LIB$INIT_TIMER')
C
        lib$set_symbol = 1

        return
        end
C
C==========================================
C
C       Subroutine lib$find_file(C_filespec,C_fname,I_context)
C
Cd      write(*,30)
C30     format('Called Lib$find_file')
C
C       return
C       end
C
C==========================================
C
C       Subroutine lib$signal(istatus)
C       integer*4 istatus
C
Cd      write(*,40)
C40     format('Called lib$signal')
C
C       return
C       end
C
C==========================================
C
C       subroutine lib$find_file_end(i_context)
C
Cd      write(*,50)
C50     format('Called lib$find_file_end')
C
C       return
C       end
C
C==========================================
C
        subroutine notify_exec_gs(izero,             ! INPUT
     1                  prod_array,          ! INPUT
     1                  i4time_array,        ! INPUT
     1                  istatflag,           ! OUTPUT?
     1                  j_status)            ! INPUT

        integer*4 prod_array(10)
C
Cd      write(*,60)
C60     format('Called lib$find_file_end')
C
        return
        end
