       Subroutine getlapsvxx(imax,jmax,kmax,maxradar,c_radar_id,         ! I
     &      n_radars,c_extension_proc,i4timefile_proc,i4_tol,rheight_3d, ! I   
     &      lat,lon,topo,i4_file_closest,                                ! I
     &      rlat_radar,rlon_radar,rheight_radar,n_valid_radars,          ! O
     &      grid_ra_ref,grid_ra_vel,istatus)                             ! O
c
       Integer       imax,jmax,kmax  !same as imax,jmax,kmax in lapsparms.for
       Integer       maxradar
       Integer       maxfiles
       Integer       i,j,k
       Integer       lvl_3d(kmax)
       Integer       i4time
       Integer       i4_tol
       Integer       n_radars
       Integer       n_ref_grids
       Integer       istatus
       Integer       istatus_2dref
       Integer       istatus_3dref
       Integer       level

       Integer       i4timefile_proc
       Integer       i4_file_closest(n_radars)

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
       Character ext*31, var_3d(kmax)*3,
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
c Read intermediate radar files. This will return the array only with those
c valid radars within the internal 20min time window of 'read_radar_3dref'.
c
      l_apply_map=.true.
      l_low_fill = .true.
      l_high_fill= .true.

      write(6,*)
      write(6,*)'get_laps_vxx: Reading v-file Reflectivity, ',
     1          '# of potential radars = ',n_radars

      k = 0

      do kcount=1,n_radars

         write(6,*)

         I4_elapsed = ishow_timer()

         EXT = c_extension_proc(kcount)

         i4_diff = i4timefile_proc - i4_file_closest(kcount)

         k = k+1      ! Used for output arrays
         write(6,*)' radar #, i4_diff = ',kcount,k,i4_diff

         call read_radar_3dref(i4_file_closest(kcount),               ! I
!    1   i4_tol,i4_ret,
     1   l_apply_map,r_missing_data,
     1   imax,jmax,kmax,ext,                                          ! I
     1   lat,lon,topo,l_low_fill,l_high_fill,
     1   rheight_3d,
     1   grid_ra_ref(1,1,1,k),                                        ! O
     1   rlat_radar(k),rlon_radar(k),rheight_radar(k),c_radar_id(k),  ! O
     1   n_ref_grids,istatus_2dref,istatus_3dref)                     ! O

c check laps analysis values
         if(istatus_3dref.ne.1 .and. istatus_3dref.ne.-1)then
            write(6,*)'ERROR: Unsuccessful reading radar ',kcount,k,ext       
            k=k-1
         else
            write(6,*)'Successful reading radar ',kcount,k,ext
            write(6,*)'radar lat/lon/elev: ',rlat_radar(k),
     &                         rlon_radar(k),rheight_radar(k)
            level=9
            write(6,*)
            write(6,*)'Level ',level,' Analysis output'
            write(6,*)'------------------------------'
29          format(1x,'  i  j    lat   lon     topo    ref ')
            write(6,29)

            do j=1,jmax,10
            do i=1,imax,10
               write(6,30)i,j,lat(i,j),lon(i,j),topo(i,j)
     &                   ,grid_ra_ref(i,j,level,k)
            end do
            end do

30          format(1x,2i3,1x,3f7.1,1x,2(f8.1,1x))

            I4_elapsed = ishow_timer()

         endif

      enddo

      n_valid_radars = k

      return
      end
