      subroutine get_acceptable_files(i4time_now,bgpath,bgmodel,names
     +     ,max_files,oldest_forecast,max_forecast_delta,use_analysis
     +     ,bg_files,forecast_length,cmodel,ntbg,NX,NY,NZ
     +     ,rejected_files,rejected_cnt)

      implicit none

      include 'netcdf.inc'

      integer bgmodel
      integer max_files, NX,NY,NZ,rejected_cnt
      character*(*) names(max_files), cmodel, bgpath
     +     ,rejected_files(rejected_cnt)
      integer oldest_forecast, bg_files,forecast_length, i, j, k
     + , max_forecast_delta
      integer ntbg
      integer ivaltimes(ntbg)
      character*4   af,c4valtime,c4_FA_valtime
      character*3   c_fa_ext
      character*100 bg_names(max_files), fullname
      integer istatus
      integer i4time_fa
      logical use_analysis
      character*9   fname,wfo_fname13_to_fname9,fname9
      integer i4time_now,bgtime,bgtime2,previous_time,next_time,len
      integer lenf
      integer bigint, ihour, n, accepted_files, final_time, lend
      parameter(bigint=2000000000)
      logical print_message
      data print_message/.true./
      integer nlapsprds,lens,lentodot,nc
      integer record
      character*3 lapsprds(4)
      character*2 cwb_model_type
      parameter (nlapsprds=4)
      data lapsprds/'lq3','lt1','lsx','lw3'/

C
C if forecast_length < 0 return only the file which matchs i4time_now 
C if forecast_length = 0 return files for i4time_now and the preceeding
C forecast.  if forecast_length > 0 return all files for i4time_now
C to >= i4time_now+forecast_length
C      
      previous_time = 0
      call s_len(bgpath,len)
      if(bgpath(len:len).ne.'/')then
         len=len+1
         bgpath(len:len)='/'
      endif

      do bg_files=1,max_files
         do i=1,100
            bg_names(i:i) = char(0)
         enddo
      enddo
      bg_files=0


      if(bgmodel.eq.0) then
         do i=max(0,forecast_length),0,-1
            bg_files=bg_files+1
            call make_fnam_lp(i4time_now+3600*i
     +           ,bg_names(bg_files),istatus)
            bg_names(bg_files)(10:13) = '0000'
            bg_names(bg_files)(14:14) = char(0)
         enddo

      elseif(bgmodel.eq.3)then

         call s_len(cmodel,nc)
         cwb_model_type = cmodel(nc-1:nc)
         call downcase(cwb_model_type,cwb_model_type)
         next_time=bigint
         final_time=i4time_now+3600*max(0,forecast_length)
         call get_file_names(bgpath,bg_files,names,max_files,istatus)
         if(istatus.ne.1)then
            print*,'error status returned: get_file_names'
     +,' in get_acceptable_files'
            return
         endif

         do i=1,bg_files
            call get_directory_length(names(i),lend)
            call s_len(names(i),lens)
            if(names(i)(lend+1:lend+2).eq.cwb_model_type)then
               call i4time_fname_lp(names(i),i4time_fa,istatus)
               call make_fnam_lp(i4time_fa,fname9,istatus)
               lentodot=index(names(i),'.')
               c_fa_ext=names(i)(lentodot+1:lens)
               c4valtime=c4_FA_valtime(c_fa_ext)
               bg_names(i)=fname9//c4valtime
            endif
         enddo

      else

         next_time=bigint
         final_time = i4time_now+3600*max(0,forecast_length)
         call get_file_names(bgpath,bg_files,names,max_files,istatus)
         if(istatus.ne.1)then
            print*,'error status returned: get_file_names'
     +,' in get_acceptable_files'
            return
         endif
         bg_files=0

         do i=1,max_files
            call s_len(names(i),j)
            j=j-13
            if (j .ge. 0) then
               if(index(names(i)(j+1:j+13),'/').eq.0 .and.
     +              names(i)(j:j).eq.'/') then
                  if (bgmodel .eq. 4) then
                     print *, 'SBN file:',names(i)(j+1:j+13)
                     fname=wfo_fname13_to_fname9(names(i)(j+1:j+13))
                     
                     call get_sbn_model_id(names(i),cmodel,ivaltimes,
     +                    ntbg,istatus)
                     if(istatus.eq.0) then
                        print*,'error returned from get_sbn_model_id '
     +,fname
                     else
                        do k=1,ntbg
                           write(af,'(i4.4)') ivaltimes(k)/3600
                           bg_files=bg_files+1
                           bg_names(bg_files)=fname//af
c     print*,'SBN: ',bg_names(bg_files),bg_files
                        enddo
                     endif
                  else 
c     if(names(i)(j:j) .eq. '/') then
                     bg_files=bg_files+1
                     bg_names(bg_files)=names(i)(j+1:j+13)
c     print*,'NOTSBN: ',bg_names(bg_files),bg_files
                  endif
               endif
            endif
         enddo

      
      endif
c      print *,bg_names(bg_files),bg_files
      bgtime2=0
      accepted_files = 0
      n=bg_files+1
      do while((previous_time.eq.0
     +     .or.(forecast_length.lt.0.and.accepted_files.lt.1))
     +     .and.bg_files.gt.0.and.n.gt.1)
         n=n-1

         do i=1,rejected_cnt
            if(bg_names(n).eq.rejected_files(i)) goto 40
         enddo
            
         fname=bg_names(n)(1:9)
         af=bg_names(n)(10:13)
         read(af,'(i4)') ihour
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
         if (ihour .gt. oldest_forecast.and.print_message) then
            print *,'IHOUR > ',oldest_forecast,' no file created.'
            print*,'oldest_forecast is a namelist parameter'
            print_message=.false.
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
                  if(next_time.gt.i4time_now) then
                     previous_time=next_time
                  else
                     names(accepted_files) = ' '
                     accepted_files=accepted_files-1
                     previous_time=next_time
                  endif
                     
               else if(next_time.lt.bigint.and. 
     +              (next_time-(bgtime+ihour*3600))
     +              .le.max_forecast_delta*3600) then
                  previous_time = bgtime+ihour*3600
                  accepted_files=accepted_files+1
                  names(accepted_files) = bg_names(n)
                  
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
         
         bg_names(n)=' '

      enddo

c      print*,accepted_files,previous_time,
c     +     forecast_length,bg_files,n


      if (previous_time.eq.0 .or. next_time.eq.bigint) then
        print*, 'Did not find two files:',previous_time,next_time
      endif



      do i=1,accepted_files
         print*,  names(i)
      enddo

      bg_files = accepted_files

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
       else if(bgmodel.eq.2) then
          call get_eta48_dims(fullname,NX,NY,NZ,istatus)
       else if(bgmodel.eq.3) then
          if(cmodel.eq.'CWB_20FA_LAMBERT_RE')then
             NX = 91
             NY = 91
             NZ = 16
          elseif(cmodel.eq.'CWB_20FA_LAMBERT_NF')then
             NX=191
             NY=127
             NZ=16
          else
             print*,'Unknow type for bgmodel = ',bgmodel,
     +'and cmodel = ',cmodel(1:nc)
          endif
       else if(bgmodel.eq.4) then
          NX = 93
          NY = 65
          NZ = 19        
       else if(bgmodel.eq.5) then
          call get_ruc2_dims(fullname,NX,NY,NZ,istatus)

c     NX = 151
c     NY = 113
c     NZ = 40   
          print*, NX,NY,NZ    
       else if(bgmodel.eq.6.or.bgmodel.eq.8) then
          if(cmodel .eq. 'AVN_FSL_NETCDF')then
             call readavnpublicdims(fullname,NX,NY,NZ,record,istatus)
             if(istatus .ne. 0)then
                print*,'Error reading AVN public dims'
                return
             else
                print *,'AVN_FSL_NETCDF Dimensions from file',
     &' ',fullname(1:len+lenf+1),NX,NY,NZ
             endif
          elseif(cmodel.eq.'AVN_AFWA_DEGRIB')then         !this switch for AVN AFWA DEGRIB
             NX = 360
             NY = 181
             NZ = 26 
          else
             NX = 360
             NY = 181
             NZ = 16
          endif
       else if(bgmodel.eq.7) then
          NX = 185
          NY = 129
          NZ = 42
       else if(bgmodel.eq.9) then
          call get_conus_dims(fullname,NX,NY,NZ)
          print *,'nws:',NX,NY,NZ
       endif

      else
       print*,'Names array is empty. Cannot set nx/ny/nz'
      endif
      
      return 
      end
