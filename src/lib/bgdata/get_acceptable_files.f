
      subroutine get_acceptable_files(i4time_anal,bgpath,bgmodel,names
     +     ,max_files,oldest_forecast,max_forecast_delta,use_analysis
     +     ,bg_files,accepted_files,forecast_length,cmodel,NX,NY,NZ
     +     ,rejected_files,rejected_cnt)

      implicit none

      include 'netcdf.inc'

      integer bgmodel
      integer max_files, NX,NY,NZ,rejected_cnt
      character*256 names(max_files)
      character*256 rejected_files(max_files)
      character*132 cmodel
      integer oldest_forecast, bg_files,forecast_length
      integer i, j, k, kk, ij, jj
      integer max_forecast_delta
      integer ntbg,nvt
      parameter (ntbg=100)
      integer ivaltimes(ntbg)
      integer nvaltimes
      character*4   af,c4valtime,c4_FA_valtime
      character*3   c_fa_ext
      character*200 fullname
      character*256 cfilespec
      integer istatus
      integer i4time_fa
      integer sbnvaltimes
      logical use_analysis
      character*9   fname,wfo_fname13_to_fname9,fname9
      integer itimes(max_files)
      integer i4time_anal
      integer bgtime_init
      integer len,lenf,bg_len,ifile
      integer ihour, n, accepted_files, final_time, lend
      integer lentodot,nc
      character*9 bkgd(3500)
      character*4 fcst(3000,120),finit1,finit2
      integer  ifcst,ibkgd
      integer  ifcst_bkgd(3000)
      integer  i4timeinit(3000)
      integer  valid_time_1,valid_time_2
      integer  i4time_min_diff
      integer  indx_for_best_init
      integer  indx_for_best_fcst
      character*2 cwb_model_type

      character(len=256), allocatable :: bg_names(:)
      character(len=256), allocatable :: bgnames_tmp(:)

      character*256 bgpath

C
C if forecast_length < 0 return only the file which matches i4time_anal 
C if forecast_length = 0 return files for i4time_anal and the preceeding
C forecast.  if forecast_length > 0 return all files for i4time_anal
C to >= i4time_anal+forecast_length
C      
      print*, '-----------------------------'
      print*, 'get_acceptable_files: 2-27-04'
      print*, '-----------------------------'

      if(.not.allocated(bg_names))allocate(bg_names(max_files))

      bg_files=0
       
      call s_len(bgpath,len)
      if(bgpath(len:len).ne.'/')then
         len=len+1
         bgpath(len:len)='/'
      endif

      cfilespec=bgpath(1:len)//'*'

      do ifile=1,max_files
         do i=1,256
            bg_names(ifile)(i:i) = char(0)
         enddo
      enddo

c     if(bgmodel.eq.0 .and. cmodel.ne.'LAPS_FUA'.and.
c    +   cmodel.ne.'MODEL_FUA')then
c        do i=max(0,forecast_length),0,-1
c           bg_files=bg_files+1
c           call make_fnam_lp(i4time_anal+3600*i
c    +           ,bg_names(bg_files),istatus)
c           bg_names(bg_files)(10:13) = '0000'
c           bg_names(bg_files)(14:14) = char(0)
c        enddo
c        call s_len(bg_names(bg_files),j)

c     elseif(bgmodel.eq.3)then

      if(bgmodel.eq.3)then

         nvt=0
         call s_len(cmodel,nc)
         cwb_model_type = cmodel(nc-1:nc)
         call downcase(cwb_model_type,cwb_model_type)
         final_time=i4time_anal+3600*max(0,forecast_length)
c        final_time=i4time_anal+3600*max(0,oldest_forecast)
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
               nvt=nvt+1
c              print*
c              print*,'i/names(i) = ',i,names(i)(1:j)
               call i4time_fname_lp(names(i),i4time_fa,istatus)
               call make_fnam_lp(i4time_fa,fname9,istatus)
               lentodot=index(names(i),'.')
               c_fa_ext=names(i)(lentodot+1:j)
               c4valtime=c4_FA_valtime(c_fa_ext)
               bg_names(nvt)=fname9//c4valtime
c              print*,'nvt/bg_names(nvt) ',i,bg_names(nvt)(1:14)
            endif
         enddo

         bg_files=nvt

      else

         final_time = i4time_anal+3600*max(0,forecast_length)
         call get_file_times(cfilespec,max_files,names,itimes
     +,bg_files,istatus)
         if(istatus.ne.1)then
            print*,'error status returned: get_file_names'
     +,' in get_acceptable_files'
            return
         endif

         nvaltimes = 1
         do i=1,bg_files

            if(bgmodel .eq. 0)then
               call get_time_length(names(i),j)
            else    
               call s_len(names(i),j)
            endif

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

c ok, if we do not want 0-hr fcst files (analysis files) then filter them.
c ---------------------------------------------------------------------------
      ij=0
      if(bg_files .gt.0)then
         if(.not.use_analysis)then
          allocate(bgnames_tmp(bg_files))
          call s_len(bg_names(1),j)
          do i=1,bg_files
             if(bg_names(i)(j-3:j).ne.'0000')then
                ij=ij+1
                bgnames_tmp(ij)=bg_names(i)
             endif
          enddo
          print*,'Removed ',bg_files-ij,' initial cond files '
          do i=1,bg_files
             bg_names(i)=' '
          enddo

          if(ij.eq.0)then
             print*,'All files must be initial files. Return to main'
             deallocate(bgnames_tmp)
             return
          endif

          bg_files=ij
          do i=1,bg_files
             bg_names(i)=bgnames_tmp(i)
          enddo
          deallocate(bgnames_tmp)
         endif
      else
          print*,'Variable bg_files = 0, return to main'
          print*
          return
      endif

      print*,TRIM(bg_names(bg_files)),bg_files
c
c categorize file by init and fcst times. convert init time to i4time
c -------------------------------------------------------------------
      ij=0
      ifcst=1
      ibkgd=1
      do n=1,bg_files-1
         finit1=bg_names(n)(6:9)
         finit2=bg_names(n+1)(6:9)
         if(finit1 .eq. finit2)then
            fcst(ibkgd,ifcst)=bg_names(n)(10:13) !4 character time (ffff) of forecast time assoc with initial.
            ifcst=ifcst+1
         else
            if(bgmodel.eq.0 .and. cmodel.eq.'LAPS_FUA'.or.
     +cmodel.eq.'MODEL_FUA')then
               fcst(ibkgd,ifcst)=bg_names(n)(10:11)
            else
               fcst(ibkgd,ifcst)=bg_names(n)(10:13)
            endif
            bkgd(ibkgd)=bg_names(n)(1:9) !9 character name of the background initial time
            ifcst_bkgd(ibkgd)=ifcst        !number of fcsts for this initial background time
            call i4time_fname_lp(bkgd(ibkgd),i4timeinit(ibkgd),istatus)
            ibkgd=ibkgd+1
            ifcst=1
         endif
      enddo

      ifcst_bkgd(ibkgd)=ifcst
      fcst(ibkgd,ifcst)=bg_names(bg_files)(10:13)
      bkgd(ibkgd)=bg_names(bg_files)(1:9)
      call i4time_fname_lp(bkgd(ibkgd),i4timeinit(ibkgd),istatus)

      print*,'Found ',ibkgd,' initial backgrounds '

c     do n=1,ibkgd,10
c        print*,'bkgd ',n, bkgd(n),' has ',ifcst_bkgd(n),' fcsts'
c     enddo

c
c this section determines only the two fcsts bounding the analysis (i4time_anal) time.
c -----------------------------------------------------------------------------------
      n=1
      i4time_min_diff=100000
      do while(n.le.ibkgd)
         if(i4timeinit(n).le.i4time_anal .and.
     +      i4timeinit(n)+forecast_length*3600 .gt. i4time_anal)then   !let "forecast_length" window on init time qualify

            print*,'Found bkgd init that corresp to anal: ',bkgd(n)
            print*,'Num of fcst = ',ifcst_bkgd(n)

c this separate (near identical section) is due to local model having ffff = hhmm format
c whereas the second section below is ffff = hhhh.

            if(bgmodel.eq.0.and.
     +(cmodel.eq.'LAPS_FUA'.or.cmodel.eq.'MODEL_FUA'))then
               do jj=2,ifcst_bkgd(n)
                  af=fcst(n,jj-1)(1:2)
                  read(af,'(i2)',err=888) ihour
                  valid_time_1=i4timeinit(n)+ihour*3600
                  af=fcst(n,jj)(1:2)
                  read(af,'(i2)',err=888) ihour
                  valid_time_2=i4timeinit(n)+ihour*3600
                  if(valid_time_1.le.i4time_anal.and.
     +valid_time_2.ge.i4time_anal)then
                   if(abs(i4timeinit(n)-i4time_anal).lt.i4time_min_diff)
     +             then
                      i4time_min_diff=abs(i4timeinit(n)-i4time_anal)
                      indx_for_best_init=n
                      indx_for_best_fcst=jj-1
                   endif
                   print*,'Found fcsts bounding anal'
c                  print*,'Full name 1: ', bkgd(n),fcst(n,jj-1)
c                  print*,'Full name 2: ', bkgd(n),fcst(n,jj)
                   print*
                  endif
               enddo

            else

             do jj=2,ifcst_bkgd(n)
                af=fcst(n,jj-1)
                read(af,'(i4)',err=888) ihour
                valid_time_1=i4timeinit(n)+ihour*3600
                af=fcst(n,jj)
                read(af,'(i4)',err=888) ihour
                valid_time_2=i4timeinit(n)+ihour*3600
                if(valid_time_1.le.i4time_anal.and.
     +             valid_time_2.ge.i4time_anal)then
                 if(abs(i4timeinit(n)-i4time_anal).lt.i4time_min_diff)
     +then
                    i4time_min_diff=abs(i4timeinit(n)-i4time_anal)
                    indx_for_best_init=n
                    indx_for_best_fcst=jj-1
                 endif
                 print*,'Found fcsts bounding anal'
c                print*,'Full name 1: ', bkgd(n),fcst(n,jj-1)
c                print*,'Full name 2: ', bkgd(n),fcst(n,jj)
                 print*
              endif
             enddo
            endif

         endif
         n=n+1
      enddo

c
c now restore the filename corresponding to the actual original "raw" file name 
c -----------------------------------------------------------------------------
      if(indx_for_best_init.ne.0.and.indx_for_best_fcst.ne.0)then
         accepted_files = 2
         names(2) = bkgd(indx_for_best_init)//fcst(indx_for_best_init
     +,indx_for_best_fcst)
         names(1) = bkgd(indx_for_best_init)//fcst(indx_for_best_init
     +,indx_for_best_fcst+1)
         print*,'Accepted file 1: ',TRIM(names(1))
         print*,'Accepted file 2: ',TRIM(names(2))
      else
         print*,'!******************************************'
         print*,'!*** WARNING: Did not find acceptable files'
         print*,'!*** -------> Returning to main'
         print*,'!******************************************'
         accepted_files = 0
         return
      endif

c     if(bg_len.gt.0.and.bg_len.le.13)then
c        if(bg_names(n)(1:1).eq.'0'.or.
c    +         bg_names(n)(1:1).eq.'1'.or.
c    +         bg_names(n)(1:1).eq.'2'.or.
c    +         bg_names(n)(1:1).eq.'9')then
c                 fname=bg_names(n)(1:9)
c                 af=bg_names(n)(10:13)
c                 read(af,'(i4)',err=888) ihour
c                 if(bgmodel.eq.0.and.cmodel.eq.'LAPS_FUA'.or.
c    +               cmodel.eq.'MODEL_FUA')ihour=ihour/100
c        else
c           print*,'Ignoring unexpected filename ',TRIM(bg_names(n))
c        endif 
c     endif

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

      return

888   print*,'Error decoding fname ',fname
      print*,'n/bg_names(n) = ',n,bg_names(n)
      print*,'af = ',af

      return 
      end
