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

      call ingest_rrvv01_sub(NX_L_CMN,NY_L_CMN,NK_LAPS,max_rrv_radars,
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
      subroutine ingest_rrvv01_sub(nx_l,ny_l,nz_l,max_radars,nxv01,
     &nyv01,nzv01,istatus)

      implicit none
c
c acquired using an i/o routine in lapsgrid.f
c
c     integer*4  nradars
c     parameter (nradars=4)
c     real*4     dxv01,dyv01
c     parameter (dxv01=5079.3,dyv01=5079.3)

      integer*4  max_fields
      parameter (max_fields = 3)

c     integer*4  nxv01,nyv01,nzv01
      integer*4  nx_l,ny_l,nz_l
      integer*4  max_radars
      integer*4  istatus
      integer*4  i,j,k,l,m,n,nf,nr
      integer*4  nradars_proc
      integer*4  lend
      integer*4  i4time_cur
      integer*4  i4timedata
      integer*4  i4time_now_gg
      integer*4  i4time_vxx(max_radars)
      integer*4  i4time_data_v01(max_radars)
      integer*4  iradar_index(max_radars)

      integer*4  ishow_timer
      integer*4  init_timer
      integer*4  itstatus

c     real*4     radar_la1(max_radars)
c     real*4     radar_lo1(max_radars)
c     real*4     radar_la1_info(nradars)
c     real*4     radar_lo1_info(nradars)

      real*4     ri(nx_l,ny_l,max_radars)
      real*4     rj(nx_l,ny_l,max_radars)
      real*4     nyqd_v01(nxv01,nyv01,nzv01,max_radars)
      real*4     refd_v01(nxv01,nyv01,nzv01,max_radars)
      real*4     veld_v01(nxv01,nyv01,nzv01,max_radars)
      real*4     out_array_4d(nx_l,ny_l,nz_l,max_fields)
      real*4     rdummy(nx_l,ny_l)
      real*4     r_grid_ratio
      real*4     rlaps_grid_spacing
      real*4     r_missing_data

      character  cradar_fname(max_radars)*255
c     character  cpath_to_raw_v01*255
      character  cdir_vxx(max_radars)*255
      character  cfname*255

      character  ext*31
      character  var_a(max_fields)*3
      character  comment_a(max_fields)*125
      character  units_a(max_fields)*10

      character  cfname9*9
c     character  c_radar_name(max_radars)*4
c     character  c_radar_name_info(nradars)*4
      character  c_ext(max_radars)*3
      character  c_num*2

      character  cref_comment(nzv01,max_radars)*126
      character  cvel_comment(nzv01,max_radars)*126
      character  cnyq_comment(nzv01,max_radars)*126

      include   'ingest_rrv_dims.inc'
      include   'ingest_rrv_common.inc'
c
c this now in ingest_rrv_constants.inc
c
c     data radar_la1_info/38.2892,0.0,36.55,36.9252/
c     data radar_lo1_info/-107.948,0.0,-104.8627,-107.948/
c     data c_radar_name_info/'KCYS','NULL','KGLD','KFTG'/
c     data cpath_to_raw_v01/'/data/lapb/import/lapsdat/radar/rrv/'/
c
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

      do i=1,max_radars
         write(c_num,100)i 
100      format(i2)
         if(c_num(1:1).eq.' ')then
            c_num(1:1)='0'
         endif
         c_ext(i)='v'//c_num
         call get_directory(c_ext(i),cdir_vxx(i),lend)
         call get_file_time(cdir_vxx(i),i4time_cur,i4time_vxx(i))
      enddo

c src/include/ingest_rrv_common.inc + ingest_rrv_constants.dat
c takes of this now (via block_data.f and ingest_rrv.nl).
c     do i=1,max_radars
c        radar_la1(i)=radar_la1_info(i)
c        radar_lo1(i)=radar_lo1_info(i)
c        c_radar_name(i)=c_radar_name_info(i)
c     enddo

      write(*,*)
      write(*,*)'Find new rrv filenames'
      n=index(path_to_raw_rrv,' ')-1
      do i=1,max_radars
         cfname=path_to_raw_rrv(1:n)//c_radar_name(i)
         nr=index(cfname,' ')-1
         write(6,*)'Search for data ',cfname(1:nr)
         call get_file_time(cfname,i4time_cur,i4timedata)
         if(i4timedata.gt.i4time_vxx(i))then
c this is new data
            nradars_proc=nradars_proc+1
            call make_fnam_lp(i4timedata,cfname9,istatus)
            iradar_index(nradars_proc)=i
          cradar_fname(i)=cfname(1:nr)//'/'//cfname9//'.v01'
            i4time_data_v01(i)=i4timedata
         endif
      enddo

      if(nradars_proc.eq.0)then
         write(6,*)'*********************************'
         write(6,*)'No new rrv/v01 files: Terminating'
         write(6,*)'*********************************'
         goto 1000
      endif
c
c read the new data files
c -----------------------
      do i=1,nradars_proc

         j=iradar_index(i)
         call read_radar_v01_cdf(cradar_fname(j),nxv01,nyv01,nzv01,
     &cref_comment(1,i),cvel_comment(1,i),cnyq_comment(1,i),
     &refd_v01(1,1,1,i),veld_v01(1,1,1,i),nyqd_v01(1,1,1,i),istatus)
         if(istatus.ne.1)then
            write(*,*)'Error returned from read_radar_v01_cdf'
         else
            write(*,*)'Success'
         endif

c        call check_v01_cdf(veld_v01(1,1,1,i),refd_v01,results,istatus)

      enddo
c
c compute mapping look-up-table
c -----------------------------
      call genv01lut_sub(nx_l,ny_l,max_radars,nradars_proc,iradar_index,
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

      do l=1,nradars_proc    !nradars_proc represents the # of new radars found.

c
c initialize
c
         do m=1,max_fields
         do k=1,nz_l
         do j=1,ny_l
         do i=1,nx_l
            out_array_4d(i,j,k,m)=r_missing_data
         enddo
         enddo
         enddo
         enddo

         j=iradar_index(l)

         do k=1,nzv01      !nzv01 better = nz_l

            write(*,*)'Processing Level ',k
            call  rrvdat2laps_v01(nx_l,ny_l,
     &                  r_grid_ratio,
     &                  r_missing_data,
     &                  veld_v01(1,1,k,l),
     &                  ri(1,1,j),
     &                  rj(1,1,j),
     &                  nyv01,nxv01, ! input array dimensions
     &                  out_array_4d(1,1,k,1),
     &                  rdummy,rdummy,
     &                  istatus)

            if(istatus .ne. 1)then
               write(6,*)'Error returned from rrvdat2laps_v01'
               goto 900
            endif 

            call  rrvdat2laps_v01(nx_l,ny_l,
     &                  r_grid_ratio,
     &                  r_missing_data,
     &                  refd_v01(1,1,k,l),
     &                  ri(1,1,j),
     &                  rj(1,1,j),
     &                  nyv01,nxv01, ! input array dimensions
     &                  out_array_4d(1,1,k,2),
     &                  rdummy,rdummy,
     &                  istatus)

            if(istatus .ne. 1)then
               write(6,*)'Error returned from rrvdat2laps_v01'
               goto 900
            endif

            call  rrvdat2laps_v01(nx_l,ny_l,
     &                  r_grid_ratio,
     &                  r_missing_data,
     &                  nyqd_v01(1,1,k,l),
     &                  ri(1,1,j),
     &                  rj(1,1,j),
     &                  nyv01,nxv01, ! input array dimensions
     &                  out_array_4d(1,1,k,3),
     &                  rdummy,rdummy,
     &                  istatus)

            if(istatus .ne. 1)then
               write(6,*)'Error returned from rrvdat2laps_v01'
               goto 900
            endif

900      enddo
c
c output
c
         if(istatus.eq.1)then
            ext=c_ext(j)

            call build_comment(cref_comment(1,l),cvel_comment(1,l),
     &cnyq_comment(1,l),nxv01,nyv01,nzv01,max_fields,comment_a,
     &out_array_4d,istatus)

            call put_laps_multi_3d(i4time_data_v01(j),ext,var_a,
     1              units_a,comment_a,out_array_4d,NX_L,NY_L,NZ_L,nf,
     1              istatus)

            if(istatus.ne.1)then
               write(6,*)'Error status returned from put_laps_multi_3d'
            endif
         else
            write(6,*)'Not writing output due to error in rrvdat2laps'
         endif

      enddo    !nradars_proc

      itstatus=ishow_timer()
      write(*,*)'Elapsed time (sec): ',itstatus

1000  return
      end
