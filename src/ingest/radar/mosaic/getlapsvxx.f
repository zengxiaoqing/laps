       Subroutine getlapsvxx(imax,jmax,kmax,maxradar,c_radar_id,         ! I
     &      n_radars,c_extension_proc,i4timefile_proc,i4_tol,rheight_3d, ! I   
     &      lat,lon,topo,i4_file_closest,                                ! I
     &      nx_r,ny_r,igrid_r,                                           ! I
     &      rlat_radar,rlon_radar,rheight_radar,n_valid_radars,          ! O
     &      grid_ra_ref,maxradarg,                                       ! O
     &      grid_ra_ref_offset,ioffset,joffset,maxradaro,                ! O
     &      l_offset,                                                    ! I
     &      istatus)                                                     ! O
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
       Integer       ioffset(maxradar)
       Integer       joffset(maxradar)

!      Note that only one of maxradarg, maxradaro is non-zero
       Real          grid_ra_ref(imax,jmax,kmax,maxradarg)
       Real          grid_ra_ref_offset(nx_r,ny_r,kmax,maxradaro)
       Real          grid_ra_ref_3d(imax,jmax,kmax)             
       Real          rheight_3d (imax,jmax,kmax)
       Real          lat(imax,jmax)
       Real          lon(imax,jmax)
       Real          topo(imax,jmax)
       Real          rlat_radar(maxradar)
       Real          rlon_radar(maxradar)
       Real          rheight_radar(maxradar)
       Real          r_missing_data

       Real          zcoord_of_level
c
c readlaps stuff
c
       Character ext*31, var_3d(kmax)*3,
     &lvl_coord_3d(kmax)*4,units_3d(kmax)*10,
     &comment_3d(kmax)*125

       Character       c_extension_proc(maxradar)*4
       Character*31    ext_a(maxradar)
       Character*4     c_radar_id(maxradar)
       Character*1     c_again

       Logical l_apply_map
       Logical l_low_fill
       Logical l_high_fill
       Logical l_offset

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

!     Number of valid radars - change this index to 'l'?
      k = 0

      do kcount=1,n_radars

         write(6,*)

         I4_elapsed = ishow_timer()

         EXT = trim(c_extension_proc(kcount))

         i4_diff = i4timefile_proc - i4_file_closest(kcount)

         k = k+1      ! Used for output arrays
         write(6,*)' radar #, i4_diff = ',kcount,k,i4_diff

         i4_tol_radar = 0 

         call read_radar_3dref(i4_file_closest(kcount),               ! I
     1   i4_tol_radar,i4_ret,                                         ! I/O
     1   l_apply_map,r_missing_data,
     1   imax,jmax,kmax,ext,                                          ! I
     1   lat,lon,topo,l_low_fill,l_high_fill,
     1   rheight_3d,
     1   grid_ra_ref_3d,                                              ! O
     1   rlat_radar(k),rlon_radar(k),rheight_radar(k),c_radar_id(k),  ! O
     1   n_ref_grids,istatus_2dref,istatus_3dref)                     ! O

c check laps analysis values
         if(istatus_3dref.ne.1 .and. istatus_3dref.ne.-1)then
            write(6,*)'ERROR: Unsuccessful reading radar ',kcount,k,ext       
            k=k-1
         else
            write(6,*)'Successful reading radar ',kcount,k,ext
            I4_elapsed = ishow_timer()
            write(6,*)'radar lat/lon/elev: ',rlat_radar(k),
     &                         rlon_radar(k),rheight_radar(k)
            level=9
            write(6,*)
            write(6,*)'Level ',level,' Analysis output'
            write(6,*)'------------------------------'
29          format(1x,'  i  j    lat   lon     topo    ref ')
            write(6,29)

!           Move radar data to second array via offset
            if(l_offset)then
              if(rlat_radar(k) .eq. r_missing_data .or.
     1           rlon_radar(k) .eq. r_missing_data      )then
                write(6,*)' No valid or single lat/lon for radar ',k
!    1                   ,' ',radar_name(k)
!               l_valid_latlon(k) = .false.

                ioffset(k) = 0
                joffset(k) = 0

              else
                call latlon_to_rlapsgrid(rlat_radar(k),
     &                                   rlon_radar(k),
     &                                   lat,lon,
     &                                   imax,jmax,
     &                                   ri,rj,
     &                                   jstatus)
                if(jstatus.ne.1)then
                    write(6,*)
     1               'computing ri/rj for radar (outside domain)'    
                endif
!               write(6,*)'Name: ',radar_name(k),ri(k),rj(k),k
!               l_valid_latlon(k) = .true.

!               Offset is location of lower left corner of small array in the large array
                ioffset(k) = (nint(ri) - igrid_r) - 1
                joffset(k) = (nint(rj) - igrid_r) - 1

              endif

              I4_elapsed = ishow_timer()

              write(6,*)' getlapsvxx - offset info '
     1               ,'ri,rj,ioffset(k),joffset(k),igrid_r : '
     1                ,ri,rj,ioffset(k),joffset(k),igrid_r

              nfill = 0  

              do jo = 1,ny_r
                j = jo + joffset(k)
                if(j .ge. 1 .and. j .le. jmax)then
                  do io = 1,nx_r
                    i = io + ioffset(k)
                    if(i .ge. 1 .and. i .le. imax)then
                      grid_ra_ref_offset(io,jo,:,k) = 
     1                grid_ra_ref_3d(i,j,:)
                      nfill = nfill + 1
                    endif ! in i bounds
                  enddo ! io
                endif ! in j bounds
              enddo ! j

              write(6,*)' nfill = ',nfill

            else ! fill 4D array with 3D array contents
              grid_ra_ref(:,:,:,k) = grid_ra_ref_3d(:,:,:)

            endif ! l_offset

!           Write sample of radar data
            do j=1,jmax,20
            do i=1,imax,20
               write(6,30)i,j,lat(i,j),lon(i,j),topo(i,j)
     &                   ,grid_ra_ref_3d(i,j,level)
            end do
            end do

30          format(1x,2i5,1x,3f7.1,1x,2(f8.1,1x))

            write(6,*)'midpoint column ',grid_ra_ref_3d(imax/2,jmax/2,:)      

            I4_elapsed = ishow_timer()

         endif

      enddo

      n_valid_radars = k

      return
      end
