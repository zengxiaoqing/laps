       Subroutine getlapsvxx(imax,jmax,kmax,maxradar,c_radar_id,
     &      n_radars,c_extension_proc,i4timefile_proc,i4_tol,rheight_3d,      
     &      lat,lon,topo,
     &      rlat_radar,rlon_radar,rheight_radar,n_valid_radars,
     &      grid_ra_ref,grid_ra_vel,istatus)
c
       Integer       imax,jmax,kmax  !same as imax,jmax,kmax in lapsparms.for
       Integer       maxradar
       Integer       maxfiles
       Integer       i,j,k
       Integer       lvl_3d(kmax)
       Integer       len_dir
       Integer       i4time
       Integer       i4_tol
       Integer       n_radars
       Integer       n_ref_grids
       Integer       istatus
       Integer       istatus_2dref
       Integer       istatus_3dref
       Integer       level

       Integer       i4timefile_proc

       Real*4          grid_ra_ref(imax,jmax,kmax,maxradar)
       Real*4          grid_ra_vel(imax,jmax,kmax,maxradar)
       Real*4          rheight_3d (imax,jmax,kmax)
       Real*4          lat(imax,jmax)
       Real*4          lon(imax,jmax)
       Real*4          topo(imax,jmax)
       Real*4          rlat_radar(maxradar)
       Real*4          rlon_radar(maxradar)
       Real*4          rheight_radar(maxradar)
       Real*4          r_missing_data

       Real*4          zcoord_of_level
c
c readlaps stuff
c
       Character ext*31, directory*50, var_3d(kmax)*3,
     &lvl_coord_3d(kmax)*4,units_3d(kmax)*10,
     &comment_3d(kmax)*125

       Character       c_extension_proc(maxradar)*3
       Character*31    ext_a(maxradar)
       Character*4     c_radar_id(maxradar)
       Character*1     c_again

       Logical l_apply_map
       Logical l_low_fill
       Logical l_high_fill

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)return

c ========================================================
c Read intermediate radar file
c
      l_apply_map=.true.
      l_low_fill = .true.
      l_high_fill= .true.

      write(6,*)
      write(6,*)'get_laps_vxx: Reading v-file Reflectivity, ',
     1          '# of potential radars = ',n_radars

      n_valid_radars = 0

      do k=1,n_radars

         call get_directory(c_extension_proc(k),directory,len_dir)

         EXT = c_extension_proc(k)

         call get_file_time(directory,i4timefile_proc,i4time_nearest)       
         i4_diff = i4timefile_proc - i4time_nearest
         write(6,*)' radar #, i4_diff = ',k,i4_diff

         call read_radar_3dref(i4timefile_proc,
!    1   i4_tol,i4_ret,
     1   l_apply_map,r_missing_data,
     1   imax,jmax,kmax,ext,
     1   lat,lon,topo,l_low_fill,l_high_fill,
     1   rheight_3d,
     1   grid_ra_ref(1,1,1,k),
     1   rlat_radar(k),rlon_radar(k),rheight_radar(k),c_radar_id(k),
     1   n_ref_grids,istatus_2dref,istatus_3dref)

c check laps analysis values
c100       write(6,*)'Which output level do you want? [1-21]'
c          read(5,*)level
         if(istatus_3dref.ne.1 .and. istatus_3dref.ne.-1)then
            write(6,*)'Error reading radar ',ext
            goto 44
         else
            write(6,*)'Successful reading radar ',ext
            write(6,*)'radar lat/lon/elev: ',rlat_radar(k),
     &                         rlon_radar(k),rheight_radar(k)
            n_valid_radars = n_valid_radars + 1
         endif

         level=9
         write(6,*)
         write(6,*)'Level ',level,' Analysis output'
         write(6,*)'------------------------------'
29       format(1x,'  i  j    lat   lon     topo    ref ')
         write(6,29)

         do j=1,jmax,10
         do i=1,imax,10
            write(6,30)i,j,lat(i,j),lon(i,j),topo(i,j)
     &                ,grid_ra_ref(i,j,level,k)
         end do
         end do

c        write(6,*)'Another level [Y/y or N/n]?'
c        read(5,38)c_again
c38       format(a1)
c        if(c_again.eq.'Y'.or.c_again.eq.'y')goto 100
30       format(1x,2i3,1x,3f7.1,1x,2(f8.1,1x))

44       if(k.lt.n_radars)write(6,*)'Ok, next radar '

      enddo

      return
      end
