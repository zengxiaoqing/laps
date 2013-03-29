
        subroutine get_interval_precip(c_pcp,ext,i4time_start,i4time_end       
     1                                ,laps_cycle_time,r_missing_data
     1                                ,ni,nj,nf,pcp_2d,istatus)

!       This is written for just one field at present, though it can be
!       readily expanded to two fields. 

        real pcp_2d(ni,nj,nf)
        real pcp_buf_2d(ni,nj,nf)
        character*3 var_2d(nf)
        character*(*)ext
        character*1 c_pcp ! Used only if a single field is chosen
        character*10  units_2d(nf)
        character*125 comment_2d(nf)

        pcp_2d = 0.

        if(trim(ext) .eq. 'st4')then
            var_2d(1) = 'ppt'       
            ipcp_cycle_time = 3600         
        elseif(nf .eq. 2)then
            var_2d(1) = 'R01'
            var_2d(2) = 'S01'
            ipcp_cycle_time = laps_cycle_time
        else
            var_2d(1) = c_pcp(1:1)//'01'
            ipcp_cycle_time = laps_cycle_time
        endif

        i4time_interval = i4time_end - i4time_start
C    abdel       
        if (ipcp_cycle_time.eq.0) i4time_interval = 0
	    
        if(i4time_interval .ne. 
     1       (i4time_interval/ipcp_cycle_time) *ipcp_cycle_time)then
            write(6,*)' ERROR: Interval is not an integral # of cycles'
            pcp_2d = r_missing_data
            istatus = 0
            return
        endif

        do i4time = i4time_start+ipcp_cycle_time,i4time_end
     1             ,ipcp_cycle_time

!           if(nf .le. 2)then
            if(.false.)then
                call get_laps_multi_2d(i4time,ext,var_2d
     1                                ,units_2d,comment_2d,ni,nj,nf
     1                                ,pcp_buf_2d,istatus_file)
            else ! this can be turned on for testing if needed
                call get_laps_2d(i4time,ext,var_2d
     1                          ,units_2d,comment_2d,ni,nj
     1                          ,pcp_buf_2d,istatus_file)

            endif
            if(istatus_file .ne. 1 .and. istatus_file .ne. -1)then
                write(6,*)' ERROR: Missing precip file: ',istatus_file 
                pcp_2d = r_missing_data
                istatus = 0
                return
            endif                

            do i = 1,nf
                call array_range(pcp_buf_2d,ni,nj,rmin,rmax
     1                          ,r_missing_data)
                write(6,*)' get_interval_precip range: ',nf,rmin,rmax 
            enddo ! i

            do i = 1,nf
                call add_miss(pcp_2d(1,1,nf),pcp_buf_2d(1,1,nf)
     1                       ,pcp_2d(1,1,nf),ni,nj)
!               call add_miss(pcp_2d,pcp_buf_2d
!    1                       ,pcp_2d,ni,nj)
            enddo ! i

        enddo ! i4time

        do i = 1,nf
            call array_range(pcp_2d,ni,nj,rmin,rmax,r_missing_data)
            write(6,*)' get_interval_precip total range: ',nf,rmin,rmax
        enddo ! i

        return
        end
