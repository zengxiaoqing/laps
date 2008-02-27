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
       Program ln3_driver
c
c Program transforms WSI/NOWRAD nexrad layer composite/echo tops/VIL 
c to the laps domain (subroutine NOWRAD_to_LAPS). The data
c have been mosaiced by WSI corporation and FSL's Facility Division stores them
c FD decodes these WSI files and stores them as netCDF files
c with names yyjjjhhmm_ll _lm and _lh, for layer low (ll), layer middle (lm) and
c layer high (lh). Echo Tops (et) and VIL (vi).
c

       call get_grid_dim_xy(nx_l,ny_l,istatus)
       if(istatus.eq.1)then
         write(*,*)'nx_l and ny_l Parameters obtained'
       else
          write(*,*)'IStatus = ',IStatus,' Error - Get_laps_config'
          write(*,*)'Terminating WSI 3d radar ingest (ln3).'
          stop
       end if

       call ln3_driver_sub(NX_L,NY_L)
       stop
       end
c
c------------------------------------------------------------------------------
c
       subroutine ln3_driver_sub(nx_l,ny_l)

       implicit none
       integer   nx_l,ny_l
       real  extreme_thrsh_70
       real  extreme_thrsh_47
       integer n_data_types
       integer max_files
       integer nlevs

       parameter (max_files=100,
     &            extreme_thrsh_47=0.30,
     &            extreme_thrsh_70=0.10,
     &            n_data_types=6)

       character*150 dir_ln3
       character*31 ext_ln3
       character*150 dir_static

       character*125 comment_ll(2)
       character*125 comment_ln3(n_data_types)
       character*10 units_ll(2)
       character*10 units_ln3(n_data_types)
       character*3 var_ll(2)
       character*3 var_ln3(n_data_types)
       character*11 laps_dom_file
       character*4 lvl_coord_2d(n_data_types)

       real lat(nx_l,ny_l)
       real lon(nx_l,ny_l)
       real remapped_prod(nx_l,ny_l)
       real laps_data(nx_l,ny_l,n_data_types)
       real grid_spacing
       real data(nx_l,ny_l,2)

       real percent_extreme_47
       real percent_extreme_70
       real baddata
       real r_missing_data
       real ref_base
       real rdum

       integer   n,k,j,i,nn,id
       integer   i4time_cur
       integer   i4time_diff
       integer   min_i4time
       integer   min_i4time_diff
       integer   i4time_nearest
       integer   i4time_nearest_ln3
       integer   i4time_data(n_data_types)
       integer   i4time_data_p(n_data_types)
       integer   i4time_now_gg
       integer   i4_validTime
       integer   istatus
       integer   kstatus
       integer   lstatus
       integer   istatus_io
       integer   lvl_2d(n_data_types)
       integer   nfiles
       integer   nfiles_p
       integer   n_extreme
       integer   lend
       integer   istart,jstart
       integer   iend,jend
       integer   msng_radar
       integer   lines,elements

       character c_gridfname*50
       character c_generic_dataroot*255
       character c_filetime*9
       character ctype_data*2
       character c_data_types(n_data_types)*2
       character c_type_found(n_data_types)*2
       character c_type_found_p(n_data_types)*2
       character ctype_cur*2
       character c_filespec*255
       character cpathwsi3d*200
       character c_filenames_proc(max_files)*200
c
c these could be runtime parameters placed in
c static/ln3/ln3_parms.dat. 
c
       data c_data_types /'ll','lm','lh','et','cr','vi'/

       call get_ln3_parameters(cpathwsi3d,msng_radar
     +,id,id,id,istatus)

       call get_r_missing_data(r_missing_data,istatus)
       call get_ref_base(ref_base,istatus)
c
       i4time_cur = i4time_now_gg()
c
c get the latest ln3 product (ln3 is laps nexrad 3d incoming from WSI)
c
       call get_directory('ln3',dir_ln3,lend)
       call get_file_time(dir_ln3,i4time_cur,i4time_nearest_ln3) 
c
c get file closest in time to the current time.
c
       n=index(cpathwsi3d,' ')
       write(*,*)'Data pathname: ',cpathwsi3d(1:n-1)

       nfiles=0
       do k=1,n_data_types

          write(6,*)'WSI 3d data type: ',c_data_types(k)
          c_filespec = cpathwsi3d(1:n-1)//'*'//c_data_types(k)
          call get_file_time(c_filespec
     1          ,i4time_cur,i4time_nearest)
          if(i4time_nearest .gt. i4time_nearest_ln3)then
c     &.and.(i4time_nearest-i4time_nearest_ln3).lt.1500)then
             nfiles=nfiles+1
             c_type_found(nfiles)=c_data_types(k)
             i4time_data(nfiles)=i4time_nearest
          endif

       enddo

       if(nfiles.lt.n_data_types.and.nfiles.gt.0)then !we know that some data_types have not updated

          call wait_for_wsi_3d_radar(nfiles,n_data_types,
     &c_data_types,c_type_found,i4time_data)
 
          nfiles_p=nfiles
          do i=1,nfiles_p
             c_type_found_p(i)=c_type_found(i)
             i4time_data_p(i)=i4time_data(i)
          enddo

       elseif(nfiles.eq.n_data_types)then    !still need to determine if all are the same time.

          min_i4time = 2000000000
          do i=1,nfiles
             min_i4time = min(i4time_data(i),min_i4time)
          enddo
 
          do i=1,n_data_types
             ctype_cur=c_data_types(i)
             do j=1,nfiles
                if(c_type_found(j).eq.ctype_cur)then
                   c_filespec=cpathwsi3d(1:n-1)//'*_'//ctype_cur
                   call get_file_time(c_filespec,min_i4time,
     &i4time_nearest)
                   if(i4time_nearest.eq.min_i4time)then
                      i4time_data(i)=i4time_nearest
                   else
                      c_type_found(j)='  '
                   endif
                endif
             enddo
          enddo 
                   
          nfiles_p=0
          do i=1,n_data_types
             if(c_type_found(i).ne.'  ')then
                nfiles_p=nfiles_p+1
                i4time_data_p(nfiles_p)=i4time_data(i) 
                c_type_found_p(nfiles_p)=c_type_found(i)
              endif
          enddo

       else     !nfiles_p .eq. 0

          write(6,*)'No Current WSI 3d data - Terminating'
          return

       endif

       do k=1,nfiles_p

          call make_fnam_lp(i4time_data_p(k),c_filetime,ISTATUS)

          c_filenames_proc(k)=cpathwsi3d(1:n-1)//c_filetime//'_'//
     &c_type_found_p(k)
          nn=index(c_filenames_proc(k),' ')
          write(6,*)'Filename Data: ',c_filenames_proc(k)(1:nn-1)
 
       enddo
c
c This for output.  LAPS ln3 files as indicated.
c
       ext_ln3 = 'ln3'

       var_ln3(1)='r04'
       var_ln3(2)='r48'
       var_ln3(3)='r8c'
       var_ln3(4)='et'
       var_ln3(5)='rco'
       var_ln3(6)='vil'

       comment_ln3(1) = 'WSI NEXRAD low-level layer composite
     &reflectivity (dBZ)'
       comment_ln3(2) = 'WSI NEXRAD mid-level layer composite
     &reflectivity (dBZ)'
       comment_ln3(3) = 'WSI NEXRAD high-level layer composite
     &reflectivity (dBZ)'
       comment_ln3(4) = 'WSI NEXRAD echo tops height (m)'
       comment_ln3(5) = 'WSI NEXRAD composite reflectivity (dBZ)'
       comment_ln3(6) = 'WSI NEXRAD vertically integrated liquid
     &(g/m3)'
c
c Definitions needed for acquiring LAPS latitude and longitude arrays.
c
       call get_directory('static',dir_static,lend)
       var_ll(1)='LAT'
       var_ll(2)='LON'

       call find_domain_name(c_generic_dataroot,c_gridfname,istatus)

       call rd_laps_static(dir_static, c_gridfname, nx_l, ny_l, 2,
     &     var_ll, units_ll, comment_ll, data, grid_spacing,
     &     istatus)

       if(istatus.eq.1)then
          write(*,*)'LAPS lat/lon grid obtained'
          write(*,*)
          do j=1,ny_l
             do i=1,nx_l
                lat(i,j)=data(i,j,1)
                lon(i,j)=data(i,j,2)
             end do
          end do
       else
          write(6,*)'Unable to get lat/lon data'
          write(6,*)'NOWRAD-ln3 process terminating'
          stop
       end if
c
c Initialize the laps output data array
c
       do k=1,n_data_types
       do j=1,ny_l
       do i=1,nx_l
          laps_data(i,j,k)=r_missing_data
       enddo
       enddo
       enddo

       call get_wsi_3drad_dims(c_filenames_proc(1),nlevs,
     .                       elements,lines,istatus)
c
c
c *****************************************************************************
c process the WSI NEXRAD data. Remap to LAPS domain.  Generate output.
c *****************************************************************************
c
       do k=1,nfiles_p

          ctype_data=c_filenames_proc(k)(nn-2:nn-1)
          write(6,*)'ctype data: ',ctype_data

          call NEXRADWSI_to_LAPS(c_filenames_proc(k),
     &                 msng_radar,
     &                 lines,elements,nlevs,
     &                 nx_l,ny_l,
     &                 lat,
     &                 lon,
     &                 i4_validTime,
     &                 remapped_prod,
     &                 istatus)
          if(istatus .eq. 1)then
             write(6,*)'WSI data properly remapped'
          else
             write(6,*)'Problem remapping the data'
          endif
c
          percent_extreme_70=0.0
          percent_extreme_47=0.0
c
c c_data_types(4)='et': c_data_types(6)='vi' {echo tops and VIL}
c
          if( (ctype_data.ne.c_data_types(4)) .and.
     &         ctype_data.ne.c_data_types(6))then

             n_extreme=0
             do j=1,ny_l
             do i=1,nx_l
                if(remapped_prod(i,j) .ge. 70.)then
                   n_extreme=n_extreme+1
                end if
             end do
             end do
             percent_extreme_70 = n_extreme/(nx_l*ny_l)
             n_extreme=0
             do j=1,ny_l
             do i=1,nx_l
                if(remapped_prod(i,j) .ge. 47.)then
                   n_extreme=n_extreme+1
                 end if
             end do
             end do
             percent_extreme_47 = n_extreme/(nx_l*ny_l)

             if(percent_extreme_70 .le. extreme_thrsh_70 .and.
     &percent_extreme_47 .le. extreme_thrsh_47)then

                baddata=80.   !this is reflectivity

                call check_radar_data(nx_l,ny_l,baddata,remapped_prod,
     &kstatus)
                if(kstatus.lt.0)then
                   write(6,*)'Found Bad data points: ',kstatus
                endif
c
c load the appropriate array for output
c
                if(ctype_data.eq.c_data_types(1))then
                   call move(remapped_prod,laps_data(1,1,1),nx_l,ny_l)
                elseif(ctype_data.eq.c_data_types(2))then
                   call move(remapped_prod,laps_data(1,1,2),nx_l,ny_l)
                elseif(ctype_data.eq.c_data_types(3))then
                   call move(remapped_prod,laps_data(1,1,3),nx_l,ny_l)
                elseif(ctype_data.eq.c_data_types(5))then
                   call move(remapped_prod,laps_data(1,1,5),nx_l,ny_l)
                endif

             else

                write(6,*)'% extreme 70 = ',percent_extreme_70
                write(6,*)'% extreme 47 = ',percent_extreme_47
                write(6,*)'Not using this layer data'
                write(*,*)'i4Time: ',i4time_data(k),
     &              ' fTime: ',c_filenames_proc(k)(1:nn-1)

             endif

          elseif(ctype_data.eq.c_data_types(4))then
c
c echo tops
c
             baddata=75.   !this is kft. Max value in WSI data = 70,000 ft.
             call check_radar_data(nx_l,ny_l,baddata,
     &remapped_prod,kstatus)
             if(kstatus.lt.0)then
                write(6,*)'Found Bad data points: ',kstatus
             endif
             do j=1,ny_l
             do i=1,nx_l
                if(remapped_prod(i,j).ne.r_missing_data.and.
     &remapped_prod(i,j).ge.0.0)then
                   laps_data(i,j,4)=remapped_prod(i,j)*304.8  !1000.0*.3048 convert from kft to m.
                elseif(remapped_prod(i,j).lt.0.0)then
                   laps_data(i,j,4)=0.0
                endif
             enddo
             enddo

          elseif(ctype_data.eq.c_data_types(6))then
c
c vertically integrated liquid
c
             baddata=1000.   !this is kg/m**2
             call check_radar_data(nx_l,ny_l,baddata,
     &remapped_prod,kstatus)
             if(kstatus.lt.0)then
                write(6,*)'Found Bad data points: ',kstatus
             endif
             do j=1,ny_l
             do i=1,nx_l
                if(remapped_prod(i,j).ne.r_missing_data.and.
     &remapped_prod(i,j).ge.0.0)then
                   laps_data(i,j,6)=remapped_prod(i,j)*1000.0   !convert kg/m**2 to g/m**2
                elseif(remapped_prod(i,j).lt.0.0)then
                   laps_data(i,j,6)=0.0
                endif
             enddo
             enddo

          endif

       enddo

       do i=1,n_data_types
          lvl_2d(i)=0
       enddo

       call write_laps_data(i4time_data(1),
     &                      dir_ln3,
     &                      ext_ln3,
     &                      nx_l,ny_l,
     &                      n_data_types,
     &                      n_data_types,
     &                      var_ln3,
     &                      lvl_2d,
     &                      lvl_coord_2d,
     &                      units_ln3,
     &                      comment_ln3,
     &                      laps_data,
     &                      istatus)

       if(istatus.eq.1)then
          write(*,*)'LN3 file successfully written'
          do i=1,nfiles
             write(*,*)'for: ',c_filenames_proc(i)(1:nn-1)
          enddo

          write(*,*)'i4 time: ',i4time_data(1)
       else
          write(*,*)'Error writing LN3 file - Terminating'
          write(*,*)'i4Time: ',i4time_data(1),
     &             ' fTime: ',c_filenames_proc(k)(1:nn-1)

       end if

1000   write(6,*)'finished in WSI-3dprod-driver'

16     stop
       end
