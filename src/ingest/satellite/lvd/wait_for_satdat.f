
       subroutine wait_for_satdat(c_channel_type,i4time_data,ntm,
     &data_dir,max_channels,nchannels,c_mtype,lvis_flag)
c
c
c
       implicit none

       integer max_channels
       integer nchannels
       integer ntm
       integer i4time_data
       integer i4_check_interval
       integer i4_thresh_age
       integer i4_total_wait
       integer i4time_nearest
       integer istatus
       integer i,j,n,nn
       integer lend

       character*3   c_mtype(max_channels)
       character*3   c_channel_type(max_channels)
c
       character*5   c_ext(max_channels)
       character*200 data_dir
       character*150 static_dir
       character*255 path

       logical   foundit(max_channels)
       logical   lvis_flag
c
c get wait parameters in case we need to wait for the data
c
       call get_directory('static',static_dir,lend)
       path=static_dir(1:lend)//'lvd/'
       lend=index(path,' ')-1
       open(15,file=path(1:lend)//'lvd_wait.parms',
     &form='formatted',err=995)
       read(15,87)i4_check_interval
       read(15,88)i4_total_wait
       read(15,89)i4_thresh_age
       close(15)
87     format(i2)
88     format(i3)
89     format(i4)
       write(6,*)'check-interval, total_wait, thresh age'
       write(6,*)i4_check_interval,i4_total_wait,i4_thresh_age
c
c set "found" flag. set extensions (because the filenames are not in the same order
c                   every time.
c
       do i=1,nchannels
          foundit(i)=.false.
       enddo

       do i=1,nchannels
          do j=1,ntm
             if(c_mtype(j).eq.c_channel_type(i))then
                foundit(i)=.true.
             endif
          enddo

          if(.not.foundit(i))then
             c_ext(i)='_'//c_channel_type(i)
          endif

       enddo
c
c this disables waiting for visible data when solar_altitude is low
c
       if(lvis_flag)then
          do i=1,nchannels
             if(c_ext(i)(1:4).eq.'_vis')then
                foundit(i)=.true.
             endif
          enddo
       endif
c
       n=index(data_dir,' ')-1
c
       do i=1,nchannels

          if(.not.foundit(i))then
             nn=index(c_ext(i),' ')-1
             path=data_dir(1:n)//'*'//c_ext(i)(1:nn)
             nn=index(path,' ')
             write(6,*)'Waiting for ',path(1:nn)
             write(6,*)'---------------------'
             call wait_for_data(path,i4time_data
     1               ,i4_check_interval,i4_total_wait
     1               ,i4_thresh_age
     1               ,istatus)
             if(istatus .ne. 1)then
                write(6,*)'No data found: ',c_ext(i)
                foundit(i)=.true.
             else
                call get_latest_file_time(path,i4time_nearest)
                if(i4time_nearest.eq.i4time_data)then
                   ntm=ntm+1
                   c_mtype(ntm)=c_channel_type(i)
                   foundit(i)=.true.
                else
                   write(6,*)'data newer than that waited for'
                   foundit(i)=.true.
                endif
             end if
          endif

       enddo

       goto 999

995    write(6,*)'Error getting wait parms'

999    return
       end
