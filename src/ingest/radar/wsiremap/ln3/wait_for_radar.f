
       subroutine wait_for_wsi_3d_radar(nfiles,n_data_types,
     &c_data_types,c_type_found,i4time_data)
c
c
       implicit none

       integer msng
       integer nfiles
       integer n_data_types
       integer i4time_data(n_data_types)
       integer i4time_cur
       integer i4_check_interval
       integer i4_thresh_age
       integer i4_total_wait
       integer i4time_nearest
       integer istatus
       integer i,j,n
       integer id
       integer lend

       character*2   c_data_types(n_data_types)
       character*2   c_type_found(n_data_types)
       character*150 cdir_static
       character*200 c_data_path
       character*255 path

       logical   foundit(4)
c
c get wait parameters in case we need to wait for the data
c
       call get_ln3_parameters(c_data_path,msng,i4_check_interval,
     +i4_total_wait,i4_thresh_age,id,id,id,id,istatus)
       if(istatus.ne.1)then
          print*,'Error getting wait parms from ln3.nl'
          goto 999
       endif
       write(6,*)'check-interval, total_wait, thresh age'
       write(6,*)i4_check_interval,i4_total_wait,i4_thresh_age
c
       do i=1,n_data_types
          foundit(i)=.false.
          do j=1,nfiles
             if(c_type_found(j).eq.c_data_types(i))then
                foundit(i)=.true.
                i4time_cur=i4time_data(j)
             endif
          enddo
       enddo
c
       n=index(c_data_path,' ')

       do i = 1,n_data_types

       if(.not.foundit(i))then
          path=c_data_path(1:n-1)//'*_'//c_data_types(i)
          write(6,*)
          write(6,*)'Waiting for ',c_data_types(i),' data '
          write(6,*)'-------------------------------------'
          call wait_for_data(path,i4time_cur
     1               ,i4_check_interval,i4_total_wait
     1               ,i4_thresh_age
     1               ,istatus)
          if(istatus .ne. 1)then
             write(6,*)'No data found in Wait-for-data',
     &c_data_types(i)
          else
             call get_latest_file_time(path,i4time_nearest)
             if(i4time_nearest.eq.i4time_cur)then
                nfiles=nfiles+1
                i4time_data(nfiles)=i4time_nearest
                c_type_found(nfiles)=c_data_types(i)
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
