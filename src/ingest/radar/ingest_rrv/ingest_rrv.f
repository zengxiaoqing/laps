      program ingest_rrv

      include 'lapsparms.cmn'
      include 'ingest_rrv_dims.inc'
      include 'ingest_rrv_common.inc'
c
c ==================================================================
c
      Call get_laps_config('nest7grid',IStatus)
      if(istatus.eq.1)then
         write(*,*)'LAPS Parameters obtained'
      else
         write(*,*)'IStatus = ',IStatus,'Error - Get_LAPS_Config'
         write(*,*)'Terminating LAPS-VRC. WSI remapping'
         stop
      end if

      call config_ingest_rrv_parms(nxv01,nyv01,nzv01,istatus)

      call ingest_rrvv01_sub(NX_L_CMN,NY_L_CMN,NK_LAPS,
     &nxv01,nyv01,nzv01,istatus)

      if(istatus.ne.1)then
         write(6,*)'Error returned from ingest_rrvv01_sub'
      else
         write(6,*)'Finished processing in ingest_rrvv01_sub'
      endif

      stop
      end
c
c =====================================================================
c
      subroutine ingest_rrvv01_sub(nx_l,ny_l,nz_l,
     &nxv01,nyv01,nzv01,istatus)

      implicit none
c
      include 'ingest_rrv_dims.inc'
      include 'ingest_rrv_common.inc'

      integer  max_fields
      parameter (max_fields = 3)
      integer  max_files
      parameter (max_files  = 1000)

c     integer  nxv01,nyv01,nzv01
      integer  nx_l,ny_l,nz_l
      integer  istatus
      integer  i,j,k,l,m,n,nf,nr,nt
      integer  jj
      integer  lend
      integer  i4time_cur
      integer  i4time_data_rrv(mxf,max_radars_rrv)
      integer  i4time_now_gg
      integer  i4times_vxx(max_files,max_radars_rrv)
      integer  i4times_rrv(max_files,max_radars_rrv)
      integer  irrv(max_radars_rrv)
      integer  irrv_index(mxf,max_radars_rrv)
      integer  i_nbr_files_out_rrv(max_radars_rrv)
      integer  i_nbr_files_out_vxx(max_radars_rrv)

      integer  ishow_timer
      integer  init_timer
      integer  itstatus

      real*4     ri(nx_l,ny_l,max_radars_rrv)
      real*4     rj(nx_l,ny_l,max_radars_rrv)
      real*4     nyqd_v01(nxv01,nyv01,nzv01,mxf,max_radars_rrv)
      real*4     refd_v01(nxv01,nyv01,nzv01,mxf,max_radars_rrv)
      real*4     veld_v01(nxv01,nyv01,nzv01,mxf,max_radars_rrv)
      real*4     out_array_4d(nx_l,ny_l,nz_l,max_fields)
      real*4     rdummy(nx_l,ny_l)
      real*4     r_grid_ratio
      real*4     rlaps_grid_spacing
      real*4     r_missing_data

      logical    found_match
      logical    founddata

      character  cradar_fname(mxf,max_radars_rrv)*200
      character  cdir_vxx(max_radars_rrv)*255
      character  cfnames_vxx(max_files,max_radars_rrv)*200
      character  cfnames_rrv(max_files,max_radars_rrv)*200
      character  cfname*255

      character  ext*31
      character  var_a(max_fields)*3
      character  comment_a(max_fields)*125
      character  units_a(max_fields)*10

      character  cfname9*9
      character  c_ext(max_radars_rrv)*3
      character  c_num*2

      character  cref_comment(nzv01,mxf,max_radars_rrv)*126
      character  cvel_comment(nzv01,mxf,max_radars_rrv)*126
      character  cnyq_comment(nzv01,mxf,max_radars_rrv)*126

c ======================== START ==============================
c -------------------------------------------------------------
c search for fresh v01 files. Assumming that the order of radars in
c the data statements above matches the files in v01, v02, and v03.
c That is, alphabetical (KCYS, KFTG, KGLD) -> v01, v02, v03.
c
c 2/15/98 (JRS) -> For testing we now have KFTG as the 4th radar and
c the second radar does not exist! (KCYS, , KGLD, KFTG) -> v01, v02, v03, V04).
c
      i4time_cur = i4time_now_gg()

      write(*,*)'Start rrv_ingest'
      itstatus=init_timer()
      itstatus=ishow_timer()

      write(*,*)'Get lapsprd/vxx filesnames'

      do i=1,max_radars_rrv
         write(c_num,100)i 
100      format(i2)
         if(c_num(1:1).eq.' ')then
            c_num(1:1)='0'
         endif
         c_ext(i)='v'//c_num
         call get_directory(c_ext(i),cdir_vxx(i),lend)
c        call get_file_time(cdir_vxx(i),i4time_cur,i4time_vxx(i))
         call get_file_times(cdir_vxx(i),max_files,cfnames_vxx(1,i)
     1                      ,i4times_vxx(1,i),i_nbr_files_out_vxx(i)
     1                      ,istatus)
      enddo

      write(*,*)
      write(*,*)'Find new rrv filenames'
      n=index(path_to_raw_rrv,' ')-1
      do i=1,max_radars_rrv
         cfname=path_to_raw_rrv(1:n)//c_radar_name(i)//'/*'
         call get_directory_length(cfname,nr)
         write(6,*)'Search for data ',cfname(1:nr)

         call get_file_times(cfname,max_files,cfnames_rrv(1,i)
     1                      ,i4times_rrv(1,i),i_nbr_files_out_rrv(i)
     1                      ,istatus)

c        call get_file_time(cfname,i4time_cur,i4timedata)

      enddo
c
c compare the vxx and rrv file results and identify rrv data that has not
c yet been processed into vxx.
c
      do m=1,max_radars_rrv
         do j=1,i_nbr_files_out_rrv(m)

            found_match=.false.
            do i=1,i_nbr_files_out_vxx(m)
               if(i4times_rrv(j,m).eq.i4times_vxx(i,m))then
                  found_match=.true.
               endif
            enddo

            if(.not.found_match)then
               irrv(m)=irrv(m)+1
               irrv_index(irrv(m),m)=j
            endif

         enddo
      enddo
      do m=1,max_radars_rrv
       do i=1,irrv(m)
         founddata=.true.
       enddo
      enddo
      if(.not.founddata)goto 999
c 
c this is new data

      do m=1,max_radars_rrv
       do j=1,irrv(m)
        jj=irrv_index(j,m)
        cradar_fname(j,m)=cfnames_rrv(jj,m)
        i4time_data_rrv(j,m)=i4times_rrv(jj,m)
       enddo
      enddo
c
c read the new data files
c -----------------------
      do i=1,max_radars_rrv
       do j=1,irrv(i)
         call read_radar_v01_cdf(cradar_fname(j,i),nxv01,nyv01,nzv01,
     &cref_comment(1,j,i),cvel_comment(1,j,i),cnyq_comment(1,j,i),
     &refd_v01(1,1,1,j,i),veld_v01(1,1,1,j,i),nyqd_v01(1,1,1,j,i),
     &istatus)
         if(istatus.ne.1)then
            write(*,*)'Error returned from read_radar_v01_cdf'
         else
            write(*,*)'Success'
         endif

c        call check_v01_cdf(veld_v01(1,1,1,i),refd_v01,results,istatus)
       enddo
      enddo
c
c compute mapping look-up-table
c -----------------------------
      call genv01lut_sub(nx_l,ny_l,max_radars_rrv,radlat,radlon,
     &c_radar_name,radar_la1,radar_lo1,nxv01,nyv01,dxv01,dyv01,ri,rj,
     &istatus)
c
c remap the rrv-v01 files to laps-vxx where xx is 01, 02, 03, ..., v20.
c ---------------------------------------------------------------------
      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         write(6,*)'Error getting r_missing_data'
         goto 1000
      endif
      call get_grid_spacing(rlaps_grid_spacing,istatus)
      if(istatus.ne.1)then
         write(6,*)'Error getting grid_spacing'
         goto 1000
      endif

      nf=3
      var_a(1) = 'REF'
      var_a(2) = 'VEL'
      var_a(3) = 'NYQ'
      units_a(1) = 'dBZ'
      units_a(2) = 'M/S'
      units_a(3) = 'M/S'

      r_grid_ratio=((dxv01+dyv01)/2.)/rlaps_grid_spacing
      write(6,*)'r_grid_ratio = ',r_grid_ratio

      do m=1,max_radars_rrv
       do n=1,irrv(m)

         do l=1,max_fields
         do k=1,nz_l
         do j=1,ny_l
         do i=1,nx_l
            out_array_4d(i,j,k,l)=r_missing_data
         enddo
         enddo
         enddo
         enddo

         call get_time_length(cradar_fname(n,m),nt)
         cfname9=cradar_fname(n,m)(nr+1:nt)
         print*,'========================================='
         print*,'Radar Name = ',c_radar_name(m), ' Time = ',cfname9
         print*,'========================================='

         do k=1,nzv01      !nzv01 better = nz_l

            write(*,*)'Processing Level ',k
            print*,' Type = VEL'
            call  rrvdat2laps_v01(nx_l,ny_l,
     &                  r_grid_ratio,
     &                  r_missing_data,
     &                  veld_v01(1,1,k,n,m),
     &                  ri(1,1,m),
     &                  rj(1,1,m),
     &                  nyv01,nxv01, ! input array dimensions
     &                  out_array_4d(1,1,k,1),
     &                  rdummy,rdummy,
     &                  istatus)

            if(istatus .ne. 1)then
               write(6,*)'Error returned from rrvdat2laps_v01'
               goto 900
            endif 

            print*,' Type = REF'
            call  rrvdat2laps_v01(nx_l,ny_l,
     &                  r_grid_ratio,
     &                  r_missing_data,
     &                  refd_v01(1,1,k,n,m),
     &                  ri(1,1,m),
     &                  rj(1,1,m),
     &                  nyv01,nxv01, ! input array dimensions
     &                  out_array_4d(1,1,k,2),
     &                  rdummy,rdummy,
     &                  istatus)

            if(istatus .ne. 1)then
               write(6,*)'Error returned from rrvdat2laps_v01'
               goto 900
            endif

            print*,' Type = NYQ'
            call  rrvdat2laps_v01(nx_l,ny_l,
     &                  r_grid_ratio,
     &                  r_missing_data,
     &                  nyqd_v01(1,1,k,n,m),
     &                  ri(1,1,m),
     &                  rj(1,1,m),
     &                  nyv01,nxv01, ! input array dimensions
     &                  out_array_4d(1,1,k,3),
     &                  rdummy,rdummy,
     &                  istatus)

            if(istatus .ne. 1)then
               write(6,*)'Error returned from rrvdat2laps_v01'
               goto 900
            endif
            print*,'*****************************************'
            print*

900      enddo
c
c output
c
         if(istatus.eq.1)then
            ext=c_ext(m)

            call build_comment(cref_comment(1,n,m),cvel_comment(1,n,m),
     &cnyq_comment(1,n,m),nxv01,nyv01,nzv01,max_fields,comment_a,
     &out_array_4d,istatus)

            call put_laps_multi_3d(i4time_data_rrv(n,m),ext,var_a,
     1              units_a,comment_a,out_array_4d,NX_L,NY_L,NZ_L,nf,
     1              istatus)

            if(istatus.ne.1)then
               write(6,*)'Error status returned from put_laps_multi_3d'
            endif
            print*
         else
            write(6,*)'Not writing output due to error in rrvdat2laps'
         endif

       enddo   !irrv = number of new rrv files for this radar
      enddo    !max_radars_rrv

      itstatus=ishow_timer()
      write(*,*)'Elapsed time (sec): ',itstatus
      goto 1000

999   print*,'No new rrv data'

1000  return
      end
