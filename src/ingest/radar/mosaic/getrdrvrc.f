      subroutine get_laps_rdr(nx_l,ny_l,nz_l,z,record,i_ra_count,
     & cfname,c_radar_id,rlat_radar,rlon_radar,rheight_radar,
     & n_valid_radars,grid_ra_ref,istatus)

      implicit none

      integer nx_l,ny_l,nz_l
      integer record, x, y, z
      integer nf_fid, nf_vid, nf_status
      integer Nx, Ny, nc
      integer imax, jmax, kdim, kmax
      integer ref_fcinv( z, record)
      integer i,j, i_ra_count, n_valid_radars
      integer n_ref_grids(i_ra_count)
      integer istatus

      real Dx, Dy, La1, Latin1, Latin2, Lo1, LoV, level, r_missing_data       
      real ref(nx_l,ny_l,z,record)
      real grid_ra_ref(nx_l, ny_l, nz_l, i_ra_count) 
      real rlat_radar(i_ra_count)
      real rlon_radar(i_ra_count)
      real rheight_radar(i_ra_count)

      double precision reftime(record), valtime(record)

      character*132 origin_name
      character*132 x_dim
      character*132 y_dim
      character*132 earth_shape
      character*132 asctime(record)
      character*132 grid_name
      character*132 process_name
      character*132 ref_comment( z, record)
      character*132 grid_type
      character*(*) cfname(i_ra_count)
      character*(*) c_radar_id(i_ra_count)
c
c ------------------ start ----------------------
c

      write(6,*)' Subroutine get_laps_rdr:'

      istatus = 0
 
      n_valid_radars = 0

      call get_r_missing_data(r_missing_data,istatus)

      do i=1,i_ra_count

         write(6,*)' Reading radar # ',i

         call s_len(cfname(i),nc)

         call read_rdr_vrc(cfname(i),nx_l,ny_l,z,record, Nx, Ny,
     + imax, jmax, kdim, kmax, ref_fcinv, Dx, Dy, La1,
     + Latin1, Latin2, Lo1, LoV, level, ref, reftime, valtime,
     + asctime, earth_shape, grid_name, grid_type, origin_name,
     + process_name, ref_comment, x_dim, y_dim, istatus)

         if(istatus .ne. 1)then
            print*,'Error reading ',cfname(i)(1:nc)
            return
         else
            print*,'Success reading ',cfname(i)(1:nc)
            n_valid_radars = n_valid_radars + 1
         endif

         grid_ra_ref(:,:,:,i) = r_missing_data  ! Initialize this 3D radar ref

         call move(ref(1,1,1,1),grid_ra_ref(1,1,1,i),nx_l,ny_l)

c extract parameters from ref_comment ... format(2f9.3,f8.0,i7,a4,3x)
c
         read(ref_comment(1,1)(1:9),'(f9.3)')rlat_radar(i)
         read(ref_comment(1,1)(10:18),'(f9.3)')rlon_radar(i)
         read(ref_comment(1,1)(19:26),'(f8.0)')rheight_radar(i)
         read(ref_comment(1,1)(27:33),'(i7)')n_ref_grids(i)
         c_radar_id(i)(1:4)=ref_comment(1,1)(34:37)

      enddo

      return
      end
