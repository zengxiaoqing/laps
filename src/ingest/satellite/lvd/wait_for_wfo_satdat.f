
       subroutine wait_for_wfo_satdat(chtype,i4time_data,ntm,data_dir,
     &max_channels,nchannels,c_mtype,lvis_flag)
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
       integer lstatus
       integer i,j,n,nn
       integer lend
       integer ispec

       character*3   c_mtype(max_channels)
       character*3   chtype(max_channels)
c
       character*200 data_dir(max_channels)
       character*255 path

       logical   foundit(max_channels)
       logical   lvis_flag
c
c get wait parameters in case we need to wait for the data
c
       call get_directory('static',path,lend)
       path=path(1:lend)//'lvd/'
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
       write(6,*)'In wait-for-wfo-data '
c
c set "found" flag. set extensions (because the filenames are not in the same order
c
       do i=1,nchannels
          foundit(i)=.false.
       enddo

       do i=1,nchannels
          do j=1,ntm
             if(c_mtype(j).eq.chtype(i))then
                foundit(i)=.true.
             endif
          enddo
       enddo
c
c for wfo data, there is no filename extension indicative of sat data
c type. This is contained in the database directory structure.
c
c this disables waiting for visible data when solar_altitude is low
c
       if(lvis_flag)then
          do i=1,nchannels
             if(chtype(i).eq.'vis'.and.
     &(.not.foundit(i)))then
                foundit(i)=.true.
             endif
          enddo
       endif
c
       do i=1,nchannels

          call lvd_file_specifier(chtype(i),ispec,lstatus)

          if(.not.foundit(ispec))then
             nn=index(data_dir(ispec),' ')-1
             path=data_dir(ispec)(1:nn)//'*'
             nn=index(path,' ')
             write(6,*)'Waiting for ',path(1:nn)
             write(6,*)'---------------------'
             call wait_for_data(path,i4time_data
     1               ,i4_check_interval,i4_total_wait
     1               ,i4_thresh_age
     1               ,istatus)
             if(istatus .ne. 1)then
                write(6,*)'No data found: ',path(1:nn)
             else
                call get_latest_file_time(path,i4time_nearest)
                if(i4time_nearest.eq.i4time_data)then
                   ntm=ntm+1
                   c_mtype(ntm)=chtype(i)
                else
                   write(6,*)'data newer than that waited for'
                endif
             end if
          endif

       enddo

       goto 999

995    write(6,*)'Error getting wait parms'

999    return
       end
