
        function a10_to_a9(a10_time,istatus)

!       Convert a10_time (yyMMddhhmm) to a9_time (yydddhhmm)
!       Steve Albers 1998

        character*10 a10_time
        character*9 a10_to_a9, a8_to_a9, a9_time
        character*8 a8_time

        a8_time = a10_time(1:8)
        a9_time = a8_to_a9(a8_time) 
        a9_time = a9_time(1:7)//a10_time(9:10)

        a10_to_a9 = a9_time

        istatus = 1

        return
        end


        function a8_to_a9(a8_time)

!       Convert a8_time (yyMMddhh) to a9_time (yydddhhmm)
!       Steve Albers 1998

        character*9 a8_to_a9
        character*9 a9
        character*8 a8_time

        integer*4 imon_a(12)

        data imon_a/0,31,59,90,120,151,181,212,243,273,304,334/

        read(a8_time,1)iyr,imn,idy,ih
1       format(i2,i2,i2,i2)

        id = imon_a(imn)

        idays = id + idy

!       Decide whether to add a day for leap year.
        if(iyr .eq. (iyr / 4) * 4 )then
            if(imn .ge. 3)then
                idays = idays + 1
            endif
        endif

        write(a9,2)iyr,idays,ih,im
2       format(i2,i3,i2,i2)

        if(a9(3:3) .eq. ' ')a9(3:3) = '0'
        if(a9(4:4) .eq. ' ')a9(4:4) = '0'
        if(a9(6:6) .eq. ' ')a9(6:6) = '0'
        if(a9(8:8) .eq. ' ')a9(8:8) = '0'

        a8_to_a9 = a9

        return
        end

