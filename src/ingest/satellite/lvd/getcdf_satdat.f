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
       subroutine getcdf_satdat(c_sat_id,c_sat_type,
     &                      nchannels,chtype,path_to_raw_sat,
     &                      c_fname_in,lvis_flag,i_delta_t,
     &                      n_ir_lines, n_ir_elem,
     &                      n_vis_lines,n_vis_elem,
     &                      n_wv_lines,n_wv_elem,
     &                      max_channel,n_images,
     &                      nft,ntm,c_mtype,max_files,
     &                      image_ir,image_vis,
     &                      image_12,image_39,image_67,
     &                      image_lat_ir,image_lon_ir,          ! O
     &                      scale_img,                          ! O
     &                      i4time_data,
     &                      istatus)
c
c Subroutine determines files available for input time and returns
c those image files. Input time (c_fname_in) is from i4time_now_gg
c
c J Smart               1/95        
c J Smart              11/95    removed routine that set r_missing_data. This
c                               now done in main program.
c J Smart              12/95    nearly total rewrite to accomodate processing using the
c                               "at" command.
c J Smart               2/96    incorporate lvis_flag indicating visible availability
c
c J Smart               2/96    added 6.7 (channel 3) to data gathering process. This is
c                               the /public file with _wv
c J Smart               4/96    added max_channel to discriminate the number of satellite channels
c 				to process for a given satellite time period.
c J Smart               6/96    add get_static_info routine for runtime parameters
c
c J Smart               7/96    Another major rewrite allowing for the process to acquire ground
c                               station sat data. This has a different set of channel type nharacter
c                               identifiers. Also consolidate to have only image_vis and image_ir arrays.
c J Smart               9/97    Dynamic array mods
c   "                   2/98    modifications associated with satellite_dims.inc and satellite_common.inc.
      implicit none
c
c                                IMPORTANT: adjust i_delta_t to conform to the amount of time
c                                between the process run time (specified by c_fname_in) and the
c                                expected time of the data to process.
c				 For example, if n_images = 1 (only one set of satellite data to
c				 process per run), and the data become available 15 minutes after the
c				 data valid time (ie, the filename time), then i_delta_t should be just
c 				 greater than 900 seconds. If n_images = 2, then i_delta_t should be
c				 just greater than 1800 seconds.
c				 Conversely, if the data frequency is, say 5 minutes and they become
c				 available 4 minutes after their valid time n_images should = 1,
c				 i_delta_t should be greater than 240 seconds but less than 300 seconds.
c
      integer max_files
      integer i_delta_t
      integer nchannels
      integer n_wv_lines,n_wv_elem
      integer n_ir_lines,n_ir_elem
      integer n_vis_lines,n_vis_elem
      integer max_channel
      integer n_images

      real image_vis (n_vis_elem,n_vis_lines,n_images)
      real image_ir  (n_ir_elem,n_ir_lines,n_images)
      real image_12  (n_ir_elem,n_ir_lines,n_images)
      real image_39  (n_ir_elem,n_ir_lines,n_images)
      real image_67  (n_wv_elem,n_wv_lines,n_images)

      real image_lat_ir  (n_ir_elem,n_ir_lines)
      real image_lon_ir  (n_ir_elem,n_ir_lines)

      real scale_img

      integer i,j,k,n,jj,il
      integer in(max_channel)
      integer ispec
      integer nft,ntm(max_files)
      integer intm
      integer ntmnew(max_files)
      integer ith_ch
      integer ith_pre_ch
      integer ith_cur_file
      integer ith_pre_file
      integer icnt
      integer newcnt
      integer nft_tmp
      integer nfiles_lvd
      integer ifiles_sat
      integer ifiles_sat_raw
      integer nfiles_sat(max_channel)
      integer i4time_lvd(max_files)
      integer i4time_sat(max_files,max_channel)
      integer i4time_sat_raw(max_files)
      integer i4time_data(max_files)
      integer i4time_array(max_channel,max_files)
      integer i4time_data_new(max_files)
      integer cvt_wfo_fname13_i4time
      integer i4time_now
      integer i4time_now_gg
      integer i4time_in
      integer i4time
      integer itm
      integer ivar
      integer lend
      integer laps_cycle_time

      logical first_time
      logical found_lvd_match
      logical found_it
      logical lvis_flag
      logical l_parse

      integer istatus
      integer jstatus
      integer lstatus
      integer i4status
      integer gfn_status

      character     c_fname_in*9
      character     cfname9*9
      character     wfo_fname13_to_fname9*9
      character     wfo_fname13_in*13
      character     wfo_fname13*13
      character     fname9_to_wfo_fname13*13
      character     cvt_i4time_wfo_fname13*13
      character     lvd_dir*150
      character     chtype(max_channel)*3
      character     c_mtype(max_channel,max_files)*3
      character     c_mtype_new(max_channel,max_files)*3
      character     cmt(max_channel)*3
      character     c_type_sat(max_files,max_channel)*3
      character     c_sat_type*3
      character     c_sat_id*6
      character     c_fname_data(max_files)*9
      character     c_fname_data_new(max_files)*9
      character     cfd*9
      character     path_to_raw_sat(max_channel)*200
      character     pathname*255
      character     c_filename(max_files)*255
      character     c_filename_sat(max_files)*255
      character     c_filename_lvd(max_files)*255

c
c ************************************************************
c
      istatus = 1

      call get_directory('lvd',lvd_dir,lend)
      lvd_dir=lvd_dir(1:lend)//trim(c_sat_id)//'/'

      do j=1,nchannels
         call lvd_file_specifier(chtype(j),ispec,lstatus)
         call s_len(path_to_raw_sat(ispec),in(j))
      enddo
c
c Adjust time when just past top of hour to ensure processing files that are
c just before top of hour.
 
      if(c_sat_type.eq.'wfo'.or.c_sat_type.eq.'ncp')then
         print*,'Convert wfo or ncp to 9 char time'
         wfo_fname13_in=fname9_to_wfo_fname13(c_fname_in)
         if(wfo_fname13_in(12:12).eq.'0')then
            write(6,*)'Adjusting i4time'
            i4time_in=cvt_wfo_fname13_i4time(wfo_fname13_in)
C           i4time_in = i4time_in-480
C            480 seconds is not enough to subtract if running at
C            more than 7 minutes past the hour...make it
C            10 minutes to be safe...BLS 19 Jul 2002
C            i4time_in = i4time_in - 600
C Use namelist parameter for this time offset
            i4time_in = i4time_in-i_delta_t
            wfo_fname13_in = cvt_i4time_wfo_fname13(i4time_in)
            write(6,*)'New time: ',wfo_fname13_in
         endif
      else
        i4time_now = i4time_now_gg()
        call cv_asc_i4time(c_fname_in,i4time_in)
        call get_laps_cycle_time(laps_cycle_time,istatus)
        if(i4time_now-i4time_in .lt. 2*laps_cycle_time)then ! real-time
c          if(c_fname_in(8:8).eq.'0')then
            print*
!           if(c_sat_type.eq.'rll')then ! global satellites can have latency
!               write(6,*)'Adjusting I4time_In (for rll data) by ',3600
!               i4time_in = i4time_in-3600         
!           else
                write(6,*)'Adjusting I4time_In by i_delta_t ',i_delta_t       
                i4time_in = i4time_in-i_delta_t
!           endif
            call make_fnam_lp (i4time_in, c_fname_in, i4status)
            if(i4status.ne.1)then
               write(6,*)'Error converting i4time_in to new c_fname_in'
               goto 992
            endif
            write(6,*)'New time: ',c_fname_in
c          endif
        else
           write(6,*)' No time adjustment for i4time_in (non-realtime)'
        endif
      endif

      ifiles_sat_raw = 0
c
c Find raw satellite data files in the path_to_raw_sat.
c
      if(c_sat_type.eq.'wfo'.or.c_sat_type.eq.'ncp')then  

         do j=1,nchannels

            call lvd_file_specifier(chtype(j),ispec,lstatus)
      pathname=path_to_raw_sat(ispec)(1:in(j))//wfo_fname13_in(1:9)//'*'
         call s_len(pathname,n)
c           n=index(pathname,' ')
            print*,'Data pathname: ',TRIM(pathname)
 
            call get_file_names(pathname,
     &                     ifiles_sat,
     &                     c_filename,
     &                     max_files,
     &                     gfn_status)

            if(gfn_status.eq.1)then
               write(*,*)'Success in GFN (Satellite): ifiles_sat = '
     1                                               ,ifiles_sat
            endif
            do i=1,ifiles_sat
               ifiles_sat_raw=ifiles_sat_raw+1
               wfo_fname13=c_filename(i)(in(j)+1:in(j)+13)
               cfname9=wfo_fname13_to_fname9(wfo_fname13)
         c_filename_sat(ifiles_sat_raw)=c_filename(i)(1:in(j))//cfnam
     &e9//'_'//chtype(j)
            i4time_sat_raw(ifiles_sat_raw)=cvt_wfo_fname13_i4time(wfo_fn
     &ame13)
            enddo

         enddo

         if(ifiles_sat_raw .le. 0)then
            write(*,*)'+++ No Data Available +++ '
            goto 998
         else
            call get_file_names(lvd_dir,
     &                     nfiles_lvd,
     &                     c_filename_lvd,
     &                     max_files,
     &                     gfn_status)
            if(gfn_status.eq.1)then
               write(*,*)'Success in GFN (lvd)'
            else
               write(6,*)'Error GFN (lvd) ',lvd_dir
               istatus=-1
               goto 996
            endif
         endif
c
c ------   End of WFO switch  ------
c
      else
c
c all satellite channels are in files within same directory.
c assume j=1 represents this for the minimum # of channels to process.
c

         pathname=path_to_raw_sat(1)(1:in(1))//c_fname_in(1:7)//'*'
         call s_len(pathname,n)
c        n=index(pathname,' ')
         print*,'Data pathname: ',TRIM(pathname)
c
c get filenames for the indicated in pathname
c
         call get_file_names(pathname,
     &                     ifiles_sat_raw,
     &                     c_filename_sat,
     &                     max_files,
     &                     gfn_status)

         if(gfn_status.eq.1)then
            write(*,*)'Success in GFN (Satellite): ifiles_sat_raw = '
     1                                            ,ifiles_sat_raw

            if(ifiles_sat_raw .le. 0)then
               write(*,*)'+++ No Data Available +++ '
               goto 998
            else
               call get_file_names(lvd_dir,
     &                     nfiles_lvd,
     &                     c_filename_lvd,
     &                     max_files,
     &                     gfn_status)
               if(gfn_status.eq.1)then
                  write(*,*)'Success in GFN (lvd)'
               else
                  write(6,*)'Error GFN (lvd) ',lvd_dir
                  istatus=-1
                  goto 996
               endif
            endif
         endif
c
         il=9
         if(c_sat_type.eq.'ncp')il=13
         do i=1,ifiles_sat_raw
          write(6,*)' Raw filename: ',c_filename_sat(i)
          write(6,*)' Raw a9time: ',c_filename_sat(i)(in(1)+1:in(1)+il)
          call i4time_fname_lp(c_filename_sat(i)(in(1)+1:in(1)+il),
     &i4time_sat_raw(i),jstatus)
         enddo

      endif

      do i=1,nfiles_lvd
         n=index(c_filename_lvd(i),' ')-1
         call i4time_fname_lp(c_filename_lvd(i)(n-12:n-4),
     &i4time_lvd(i),jstatus)
      enddo

      do i=1,nchannels
         nfiles_sat(i)=0
      enddo
c
c determine which raw satellite data files have already been
c processed by lvd.
c
      do i=1,ifiles_sat_raw
         found_lvd_match=.false.
         do j=1,nfiles_lvd
            if(i4time_sat_raw(i).eq.i4time_lvd(j))then
               found_lvd_match=.true.
            endif
         enddo
         if(.not.found_lvd_match)then

!           Determine j as the beginning of the channel part
            do k=255,1,-1
               if(c_filename_sat(i)(k:k).eq.'_')then
                  j=k+1
                  goto 121
               endif
            enddo

!           Determine jj as the end of the channel part
121         jj=index(c_filename_sat(i),' ')-1

!           So far jj (end of file index) is assumed to be the end of the
!           channel part. We now strip off the '.nc' if present.
            if(c_filename_sat(i)(jj-2:jj) .eq. '.nc')then
               if(l_parse(c_filename_sat(i),'_vis.nc') 
     1                                        .eqv. .true.)then     
                  jj = jj-3 ! (e.g. read 'vis' part of filename)
               else
                  jj = jj-4 ! (e.g. read '10p' part of filename)
               endif
            endif

            do k=1,nchannels
               if(c_filename_sat(i)(j:jj).eq.chtype(k))then
                  nfiles_sat(k) = nfiles_sat(k) + 1
                  i4time_sat(nfiles_sat(k),k)=i4time_sat_raw(i)
                  c_type_sat(nfiles_sat(k),k)=c_filename_sat(i)(j:jj)
               endif
            enddo

         endif
      enddo
c
c categorize the files by time and by the max number of channels selected.
c
c initialize
c
      ntm=0
      first_time=.true.
      nft=0

      do ith_ch=1,nchannels

         if(first_time)then

            do ith_cur_file=1,nfiles_sat(ith_ch)

               nft=nft+1                        !# of files total (# of distinct time periods)
               ntm(nft)=ntm(nft)+1              !# of time matches (for this nft)
               c_mtype(ntm(nft),nft)=c_type_sat(ith_cur_file,ith_ch)
               i4time_array(ntm(nft),nft)=i4time_sat(ith_cur_file,ith_ch
     &)

            enddo

            first_time=.false.

         else

            ith_pre_ch = ith_ch-1
            do ith_cur_file = 1,nfiles_sat(ith_ch)

               found_it = .false.
               do while(.not.found_it .and. ith_pre_ch.ge.1)

                  do ith_pre_file = 1,nfiles_sat(ith_pre_ch)

                     if(i4time_sat(ith_pre_file,ith_pre_ch).eq.i4time_sa
     &t(ith_cur_file,ith_ch))then

                        i=1
                        do while(.not.found_it.and.i.le.nft)

                           if(i4time_array(1,i).eq.i4time_sat(ith_cur_fi
     &le,ith_ch) )then
                              ntm(i)=ntm(i)+1
                              i4time_array(ntm(i),i)=i4time_sat(ith_cur_
     &file,ith_ch)
                              c_mtype(ntm(i),i)=c_type_sat(ith_cur_file,
     &ith_ch)
                              found_it=.true.
                           else
                              i=i+1
                           endif

                        enddo
                     endif

                  enddo

                  ith_pre_ch = ith_pre_ch-1

               enddo

               if(.not.found_it)then
                  nft=nft+1
                  ntm(nft)=ntm(nft)+1
                  i4time_array(ntm(nft),nft)=i4time_sat(ith_cur_file,
     &ith_ch)
                  c_mtype(ntm(nft),nft)=c_type_sat(ith_cur_file,ith_c
     &h)
               endif 

               ith_pre_ch = ith_ch-1

            enddo

         endif

      enddo
c
c accumulate i4time into output array. 
c
       do i =1,nft
       do j =1,ntm(i)

          i4time_data(i)=i4time_array(j,i)

       enddo
       enddo
c
c list the files found that need lvd processing
c
      write(6,*)'Number of new files found ',nft
      write(6,*)'--------------------------------'
      do i=1,nft
         call make_fnam_lp(i4time_data(i),c_fname_data(i),
     &istatus)
      enddo
c
c sort times in ascending order
c
      do j=2,nft
         itm=i4time_data(j)
         cfd=c_fname_data(j)
         intm=ntm(j)
         do k=1,ntm(j)
            cmt(k)=c_mtype(k,j)
         enddo
         do i=j-1,1,-1
            if(i4time_data(i).le.itm)goto 10
            i4time_data(i+1)=i4time_data(i)
            c_fname_data(i+1)=c_fname_data(i)
            ntm(i+1)=ntm(i)
            do k=1,ntm(i)
               c_mtype(k,i+1)=c_mtype(k,i)
            enddo
         enddo
         i=0
10       i4time_data(i+1)=itm
         c_fname_data(i+1)=cfd
         ntm(i+1)=intm
         do k=1,intm
            c_mtype(k,i+1)=cmt(k)
         enddo
      enddo

      do i=1,nft
         write(6,*)'Time: ',c_fname_data(i)
         write(6,*)'Number of matches: ',ntm(i)
         ivar=ntm(i)
         write(6,101)(c_mtype(j,i),j=1,ivar)
101      format(5x,5(:,a3,1x))
      enddo
c
c This write up must now consider that we are looking for 6.7 micron data too.
c This discussion must consider parameter max_channel. (4-96)
c
c This part tests for a full complement of files for lvd processing. We are
c expecting at most 5 files (ir, 12, 39, 67, and vis). If, for the latest time,
c there are less than 5 then subroutine wait_for_data is used to wait a preset
c amount of time for this data to come in.  If the data do not come in then
c only those files that have been found are processed by lvd. Otherwise, wait_for
c data will allow processing of all 5 files.
c
c The following discussion must consider parameter max_image.
c Before this, however, it is necessary to check for a large number of files to
c be processed.  This will happen for a cold lvd start or if laps goes down and files
c continue to build up in /public
c for a long time. Limit the number of new files processed to those within parameter
c i_delta_t (sec)
c of the current processsing time.  This should limit the number of files processed.
c Note: these may not be cronological.
c
      if(nft.gt.n_images)then      !nft is the number of file times and is related to
c                                     parameter n_images
         icnt=0
         call i4time_fname_lp(c_fname_in,i4time,jstatus)
c
c since the times and other arrays have been sorted in ascending time
c order
c
         do j=nft,1,-1
            icnt=icnt+1
            if(icnt.le.n_images)then
               newcnt=icnt
               ntmnew(newcnt)=ntm(j)
               i4time_data_new(newcnt)=i4time_data(j)
               c_fname_data_new(newcnt)=c_fname_data(j)
               do i=1,ntmnew(newcnt)
                  c_mtype_new(i,newcnt)=c_mtype(i,j)
               enddo
            endif
         enddo
         nft=newcnt

      elseif(nft.gt.0)then
         icnt=0
         call i4time_fname_lp(c_fname_in,i4time,jstatus)
         do i = 1,nft
c           if( abs(i4time-i4time_data(i)).lt.i_delta_t)then
               icnt=icnt+1
               ntmnew(icnt)=ntm(i)
               i4time_data_new(icnt)=i4time_data(i)
               c_fname_data_new(icnt)=c_fname_data(i)
               do j=1,ntmnew(icnt)
                  c_mtype_new(j,icnt)=c_mtype(j,i)
               enddo
c           endif
         enddo
         nft=icnt

      elseif(nft.eq.0)then

         write(6,*)'No new satellite data files found'
         write(6,*)'Terminating lvd'
         istatus=-1
         goto 16

      endif
c
c move stuff back into the original arrays
c but first reinitialize the channel type (c_mtype) array!
c     Need to have a second check to see if nft is > 0 since if
c     data are too old then icnt = 0 as does ntf

      if(nft.gt.0)then

         do i=1,nft
            ntm(i)=ntmnew(i)
            i4time_data(i)=i4time_data_new(i)
            c_fname_data(i)=c_fname_data_new(i)
            i4time_data(i)=i4time_data_new(i)
            do j=1,ntm(i)
               c_mtype(j,i)=c_mtype_new(j,i)
            enddo
         enddo

      else

         write(6,*)'Data too old to process'
         write(6,*)'Terminating lvd'
         istatus=-1
         goto 16

      endif
c
c this part determines if the last file time (nft) has the max number of channels.
c
      nft_tmp=nft
      if(ntm(nft).lt.nchannels) nft_tmp=nft_tmp-1

      if(nft_tmp.eq.0)then

         write(6,*)'Not enough files found for this time'
         if(c_sat_type.eq.'wfo')then
            call wait_for_wfo_satdat(chtype,i4time_data(nft),ntm(nft),
     &path_to_raw_sat,max_channel,nchannels,c_mtype(1,nft),lvis_flag)

         else
            call wait_for_satdat(chtype,i4time_data(nft),ntm(nft),
     &path_to_raw_sat(1),max_channel,nchannels,c_mtype(1,nft),lvis_flag)
         endif

      elseif(nft_tmp.lt.nft)then

         write(6,*)'Not enough files found for latest time'
         if(c_sat_type.eq.'wfo')then
            call wait_for_wfo_satdat(chtype,i4time_data(nft),ntm(nft),
     &path_to_raw_sat,max_channel,nchannels,c_mtype(1,nft),lvis_flag)
         else
            call wait_for_satdat(chtype,i4time_data(nft),ntm(nft),
     &path_to_raw_sat(1),max_channel,nchannels,c_mtype(1,nft),lvis_flag)
         endif

      else

         write(6,*)'Full set of files found for all times!'

      endif
c
c get the data.
c
      call readsatdat(c_sat_id,
     &                c_sat_type,
     &                path_to_raw_sat,
     &                c_fname_data,
     &                c_mtype,
     &                nft,ntm,
     &                max_channel,
     &                max_files,
     &                n_images,
     &                n_wv_lines,n_wv_elem,
     &                n_ir_lines,n_ir_elem,
     &                n_vis_lines,n_vis_elem,
     &                image_67,
     &                image_ir,image_12,
     &                image_39,image_vis,
     &                image_lat_ir,image_lon_ir,         ! O
     &                scale_img,
     &                istatus)
 
      goto 16

992   write(6,*)'Error adjusting time'
      goto 16

996   write(6,*)'Error getting lvd filenames'
      goto 16

997   write(6,*)'Error getting static info'
      istatus = -1
      goto 16

 998  write(*,*)' +++ neither IR or VIS data available'
      istatus = -1

 16   write(*,*)' Normal completion in GETCDF_SATDAT'
      return
      end
