      subroutine get_acceptable_files(i4time_now,bgpath,bgmodel,names
     +     ,max_files,oldest_forecast,max_forecast_delta,use_analysis
     +     ,bg_files,accepted_files,forecast_length,cmodel,NX,NY,NZ
     +     ,rejected_files,rejected_cnt)

      implicit none

      include 'netcdf.inc'

      integer bgmodel
      integer max_files, NX,NY,NZ,rejected_cnt
      character*(*) names(max_files), cmodel, bgpath
     +     ,rejected_files(max_files)
      integer oldest_forecast, bg_files,forecast_length
      integer i, j, k, kk
      integer max_forecast_delta
      integer ntbg
      parameter (ntbg=100)
      integer ivaltimes(ntbg)
      integer nvaltimes
      character*4   af,c4valtime,c4_FA_valtime
      character*3   c_fa_ext
      character*200 bg_names(max_files), fullname
      character*256 cfilespec
      integer istatus
      integer i4time_fa
      integer sbnvaltimes
      logical use_analysis
      character*9   fname,wfo_fname13_to_fname9,fname9
      integer itimes(max_files)
      integer i4time_now,bgtime,bgtime2,previous_time,next_time,len
      integer lenf,bg_len,ifile
      integer bigint, ihour, n, accepted_files, final_time, lend
      parameter(bigint=2000000000)
      logical print_message
      data print_message/.true./
      integer nlapsprds,lentodot,nc
      integer record
      character*3 lapsprds(4)
      character*2 cwb_model_type
      parameter (nlapsprds=4)
      data lapsprds/'lq3','lt1','lsx','lw3'/

C
C if forecast_length < 0 return only the file which matches i4time_now 
C if forecast_length = 0 return files for i4time_now and the preceeding
C forecast.  if forecast_length > 0 return all files for i4time_now
C to >= i4time_now+forecast_length
C      
      call s_len(bgpath,len)
      if(bgpath(len:len).ne.'/')then
         len=len+1
         bgpath(len:len)='/'
      endif

      cfilespec=bgpath(1:len)//'*'

      do ifile=1,max_files
         do i=1,200
            bg_names(ifile)(i:i) = char(0)
         enddo
      enddo

      if(bgmodel.eq.0) then
         do i=max(0,forecast_length),0,-1
            bg_files=bg_files+1
            call make_fnam_lp(i4time_now+3600*i
     +           ,bg_names(bg_files),istatus)
            bg_names(bg_files)(10:13) = '0000'
            bg_names(bg_files)(14:14) = char(0)
         enddo
         call s_len(bg_names(bg_files),j)

      elseif(bgmodel.eq.3)then

         call s_len(cmodel,nc)
         cwb_model_type = cmodel(nc-1:nc)
         call downcase(cwb_model_type,cwb_model_type)
         next_time=bigint
         final_time=i4time_now+3600*max(0,forecast_length)
         call get_file_times(cfilespec,max_files,names,itimes
     +,bg_files,istatus)
         if(istatus.ne.1)then
            print*,'error status returned: get_file_names'
     +,' in get_acceptable_files'
            return
         endif

         print*,'bg_files = ',bg_files
         do i=1,bg_files
            call get_directory_length(names(i),lend)
            call s_len(names(i),j)
            if(names(i)(lend+1:lend+2).eq.cwb_model_type)then
c              print*
c              print*,'i/names(i) = ',i,names(i)(1:j)
               call i4time_fname_lp(names(i),i4time_fa,istatus)
               call make_fnam_lp(i4time_fa,fname9,istatus)
               lentodot=index(names(i),'.')
               c_fa_ext=names(i)(lentodot+1:j)
               c4valtime=c4_FA_valtime(c_fa_ext)
               bg_names(i)=fname9//c4valtime
c              print*,'i/bg_names(i) ',i,bg_names(i)(1:14)
            endif
         enddo

      else

         next_time=bigint
         final_time = i4time_now+3600*max(0,forecast_length)
         call get_file_times(cfilespec,max_files,names,itimes
     +,bg_files,istatus)
         if(istatus.ne.1)then
            print*,'error status returned: get_file_names'
     +,' in get_acceptable_files'
            return
         endif

         nvaltimes = 1
         do i=1,bg_files
            call s_len(names(i),j)
            j=j-13
            if (j .ge. 0) then
               if(index(names(i)(j+1:j+13),'/').eq.0 .and.
     +              names(i)(j:j).eq.'/') then
                  if (bgmodel .eq. 4) then

c     print *, 'SBN file:',names(i)(j+1:j+13)

                     fname=wfo_fname13_to_fname9(names(i)(j+1:j+13))
                     
                     call get_nvaltimes(names(i),sbnvaltimes,ivaltimes
     +                                 ,istatus)
                     if(istatus.ne.1) then
                        print*,'error returned from get_sbn_model_id '
     +                         ,fname
                        return
                     else
                        kk=0
                        do k=nvaltimes,nvaltimes+sbnvaltimes
                         kk=kk+1
                         if(kk.le.sbnvaltimes)then
                           write(af,'(i4.4)') ivaltimes(kk)/3600
                           bg_names(k)=fname//af
c     print*,'SBN: ',bg_names(k),k
                         endif
                        enddo

                        nvaltimes = k-1
                     endif
                  else 
c     if(names(i)(j:j) .eq. '/') then
                     bg_names(i)=names(i)(j+1:j+13)
c     print*,'NOTSBN: ',bg_names(i),bg_files
                  endif
               endif
            endif
         enddo

         if(bgmodel.eq.4)then
            bg_files = nvaltimes-1
         endif
      
      endif

c      print *,bg_names(bg_files),bg_files

      bgtime2=0
      accepted_files = 0
      n=bg_files+1
      previous_time = 0
      do while((previous_time.eq.0
     +     .or.(forecast_length.lt.0.and.accepted_files.lt.1))
     +     .and.bg_files.gt.0.and.n.gt.1)
         n=n-1

         do i=1,rejected_cnt
            if(bg_names(n).eq.rejected_files(i)) goto 40
         enddo
            
         call s_len(bg_names(n),bg_len)
	 if(bg_len.gt.0.and.bg_len.le.13)then
	    if(bg_names(n)(1:1).eq.'0'.or.
     +         bg_names(n)(1:1).eq.'1'.or.
     +         bg_names(n)(1:1).eq.'2')then
                  fname=bg_names(n)(1:9)
                  af=bg_names(n)(10:13)
                  read(af,'(i4)',err=888) ihour
	    else
	       print*,'Ignoring weird filename ',bg_names(n)(1:bg_len)
	    endif 
         endif

         call i4time_fname_lp(fname,bgtime,istatus)

         if(bgtime.lt.bgtime2.and. .not.use_analysis) then
            print*, 'Trying earlier forecast ',bgtime,bgtime2,fname,af
            bgtime2=bgtime
            accepted_files=0
            previous_time=0
            next_time=bigint
         endif
c     
c *** File names are returned sorted from newest to oldest ***
c     process only the newest which is at least as old as the laps time
c
      
         if(bgtime.gt.i4time_now .and. .not. use_analysis) then
            print*,
     +           'Background model newer than requested - skipping'
            print*,'This behavior can be changed in the namelist '
            print*,'parameter use_analysis '
            goto 40
         endif  
         
c
c ****** Do NOT process model fcst if fcst is greater than oldest_forecast hours.
c
         if (ihour .gt. oldest_forecast)then
            if(print_message)then
               print *,'IHOUR > ',oldest_forecast,' no file created.'
               print*,'oldest_forecast is a namelist parameter'
               print_message=.false.
            endif
            goto 40
         endif
c     
         if(forecast_length.lt.0) then
            if(accepted_files.le.0.and.
     +           bgtime+ihour*3600.eq.i4time_now) then
               accepted_files=1
               names(accepted_files) = bg_names(n)
               previous_time=bigint
               next_time=0
            endif
         else
            if(bgtime+ihour*3600.ge.final_time.and.
     +           bgtime+ihour*3600.lt.next_time) then
               next_time = bgtime+ihour*3600
               accepted_files=1
               names(accepted_files) = bg_names(n)
               bgtime2=bgtime
cc               print *, '1: ',next_time,accepted_files,bg_names(n)
            else if(bgtime+ihour*3600.ge.i4time_now.and.
     +              bgtime+ihour*3600.lt.next_time) then
               next_time = bgtime+ihour*3600
               if(bgtime+ihour*3600.lt.final_time) then
                  accepted_files=accepted_files+1
               endif
               names(accepted_files) = bg_names(n)
cc               print *, '2: ',next_time,accepted_files,bg_names(n)
            endif
c            print*,bg_names(n),bgtime+ihour*3600,i4time_now,n
            if(bgtime+ihour*3600.lt.i4time_now.and.
     +           bgtime+ihour*3600.gt.previous_time) then
               if(forecast_length.gt.0) then
c
c If this is a forecast request we want the first forecast > i4time_now
c to be the last file on the list 
c
                  if(next_time.ge.i4time_now) then  !was .gt.
                     previous_time=next_time
                  else
                     names(accepted_files) = ' '
                     accepted_files=accepted_files-1
                     previous_time=next_time
                  endif
                     
c              else if(next_time.lt.bigint.and. 

                  if(next_time.lt.bigint.and. 
     +              (next_time-(bgtime+ihour*3600))
     +              .le.max_forecast_delta*3600) then
                     previous_time = bgtime+ihour*3600
                     accepted_files=accepted_files+1
                     names(accepted_files) = bg_names(n)
                  endif
                  
               else
c     
c     otherwise the latest model is incomplete and doesn't have a current forecast 
c     look for an older one
c
               endif
cc       print *, '3: ',next_time,accepted_files,bg_names(n),previous_time
            endif
         endif
      
 40      continue

      enddo

c      print*,accepted_files,previous_time,
c     +     forecast_length,bg_files,n


      if (previous_time.eq.0 .or. next_time.eq.bigint) then
        print*, 'Did not find acceptable background files:'
     +          ,previous_time,next_time
      endif



      do i=1,accepted_files
         print*,'Accepted name ',i, '= ',names(i)(1:j)
      enddo

      if(accepted_files.gt.0)then

       call get_fname_length(names(1),lenf)
       fullname=bgpath(1:len)//names(1)(1:lenf)
       if(bgmodel.eq.0) then
          call get_laps_dimensions(nz,istatus)
          call get_grid_dim_xy(nx,ny,istatus)
       else if(bgmodel.eq.1) then
          NX = 81
          NY = 62
          NZ = 25        
       else if(bgmodel.eq.7) then
          NX = 185
          NY = 129
          NZ = 42
       else if(bgmodel.eq.9) then
          call get_conus_dims(fullname,NX,NY,NZ)
          print *,'nws:',NX,NY,NZ
       endif

      else

       print*
       rejected_cnt=bg_files
       do j=1,bg_files
          rejected_files(j)=bg_names(j)
       enddo
       print*,'return to main'
       print*

      endif
      return

888   print*,'Error decoding fname ',fname
      print*,'n/bg_names(n) = ',n,bg_names(n)
      print*,'af = ',af

      return 
      end
