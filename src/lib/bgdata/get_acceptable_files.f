
      subroutine get_acceptable_files(i4time_anal,bgpath,bgmodel,names
     +     ,max_files,use_analysis,use_forecast
     +     ,bg_files,accepted_files,forecast_length,cmodel,NX,NY,NZ
     +     ,rejected_files,rejected_cnt)

      implicit none

      include 'netcdf.inc'

      integer bgmodel
      integer max_files, NX,NY,NZ,rejected_cnt
      character*256 names(max_files)           ! should be returned as basenames
      character*256 names_tmp(max_files)
      character*256 rejected_files(max_files)
      character*132 cmodel
      integer bg_files,forecast_length
      integer i, j, k, kk, ij, jj,l,nclen
      integer ntbg,nvt
      parameter (ntbg=100)
      integer ivaltimes(100)
      integer nvaltimes
      character*4   af,c4valtime,c4_FA_valtime
      character*3   c_fa_ext
      character*200 fullname
      character*256 cfilespec
      integer istatus, idebug
      integer i4time_fa
      integer sbnvaltimes
      logical use_analysis, use_forecast, l_parse
      character*9   fname,wfo_fname13_to_fname9,fname9
      integer itimes(max_files)
      integer i4time_anal
      integer bgtime_init
      integer len,lenf,bg_len,ifile
      integer ihour, imin, n, accepted_files, final_time, lend
      integer lentodot,nc
      character*9 bkgd(3500)
      character*4 fcst(3000,120)
      character*9 finit1,finit2
      integer  ibkgd
      integer  ifcst_bkgd(3000)
      integer  i4timeinit(3000)
      integer  valid_time_1,valid_time_2
      integer  i4time_min_diff
      integer  indx_for_best_init1
      integer  indx_for_best_init2
      integer  indx_for_best_init
      integer  indx_for_best_fcst
      integer  laps_cycle_time
      character*2 cwb_model_type
      character*3 clapsdirs(4)
      data        clapsdirs/'lsx','lt1','lw3','lq3'/

      character(len=256), allocatable :: bg_names(:)
      character(len=256), allocatable :: bgnames_tmp(:)

      character*256 bgpath, null_string

C
C if forecast_length < 0 return only the file which matches i4time_anal 
C if forecast_length = 0 return files for i4time_anal and the preceeding
C forecast.  if forecast_length > 0 return all files for i4time_anal
C to >= i4time_anal+forecast_length
C      
      print*, '----------------------------------------------'
      print*, 'get_acceptable_files: use_analysis / use_forecast = '
     1                              ,use_analysis,use_forecast
      print*, '----------------------------------------------'

      idebug = 0

!     If later on we want to set use_forecast to .false. then we can use 
!     analyses only similar to what is presently being done with the 'LAPS' 
!     option.

      if((use_forecast .eqv. .true.) .OR. 
     1   (use_analysis .eqv. .true.)      )then       
          continue
      else
          write(6,*)' ERROR: use_forecast AND use_analysis are false'       
          return
      endif

      call get_laps_cycle_time(laps_cycle_time,istatus)

      if(.not.allocated(bg_names))allocate(bg_names(max_files))
   
      call s_len(cmodel,nclen)

      bg_files=0
       
      call s_len(bgpath,len)
      if(len .gt. 0)then
         if(bgpath(len:len).ne.'/')then
            len=len+1
            bgpath(len:len)='/'
         endif

         cfilespec=bgpath(1:len)//'*'
    
      else
         cfilespec=' '

      endif

      do i=1,256
         null_string(i:i) = char(0)
      enddo

      do ifile=1,max_files
         bg_names(ifile) = null_string
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

         search_mdltype: do l=nclen,1,-1
            if(cmodel(l:l).eq.'_')then
               exit search_mdltype
            endif
         enddo search_mdltype
         if(cmodel(l+1:nclen).eq.'NF'.or.
     &      cmodel(l+1:nclen).eq.'NF15'.or.
     &      cmodel(l+1:nclen).eq.'NF45')then
            cwb_model_type='nf'
         elseif(cmodel(l+1:nclen).eq.'GFS')then
            cwb_model_type='gb'
         elseif(cmodel(l+1:nclen).eq.'TFS')then
            cwb_model_type='sb'
         else
            print*,'Not able to determine CWB model type'
            print*,'cmodel = ',cmodel(1:nclen)
            stop
         endif

         nvt=0
         final_time=i4time_anal+3600*max(0,forecast_length)
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
!              if(i .eq. 1)then
!                  print*
!                  print*,'i/names(i) = ',i,names(i)(1:j)
!              endif
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

      elseif(cmodel.ne.'LAPS')then

         final_time = i4time_anal+3600*max(0,forecast_length)
         call get_file_times(cfilespec,max_files,names,itimes
     +,bg_files,istatus)
         if(istatus.ne.1)then
            print*,'error status returned: get_file_times'
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

            if(i .eq. 1)then
                write(6,*)' get_acceptable_files: i,names(i) ',
     1                                            i,names(i)(1:j)
            endif

            j=j-13
            if (j .ge. 0) then
!              if(index(names(i)(j+1:j+13),'/').eq.0 .and.
!    +              names(i)(j:j).eq.'/') then

!              Test for "grib" extension in the name
               if(l_parse(names(i)(j+1:j+13),'grib'))then
                 write(6,*)' ERROR: filename has .grib extension: '
     1                    ,trim(names(i))
                 write(6,*)
     1             ' remove .grib extension if present and/or use links'
               endif

               if(l_parse(names(i)(j+1:j+13),'/') .eqv. .false. .AND.
     +              names(i)(j:j).eq.'/') then
                
                  if(i .eq. 1)write(6,*)
     1            ' index no slash in last 13 chars test was true:'
     1            ,j,names(i)(j:j),' ',names(i)(j+1:j+13)

CWNI-BLS ... Support for bgmodel = 10, assuming one uses the YYYYMMDD_HHMM
CWNI-BLS     convention
                  if ((bgmodel .eq. 4).OR.(bgmodel.eq.10)) then
c     print *, 'SBN file:',names(i)(j+1:j+13)
                     fname=wfo_fname13_to_fname9(names(i)(j+1:j+13))
CWNI-BLS ... Call to new version of get_nvaltimes for Unidata
                     IF (bgmodel .EQ. 4) THEN                     
                       call get_nvaltimes(names(i),sbnvaltimes,ivaltimes
     +                                 ,istatus)
                     ELSEIF (bgmodel.EQ.10) THEN
                       print *,"Reading Unidata: ",TRIM(names(i))
                       call get_nvaltimes_unidata(names(i),sbnvaltimes,
     +                                 ivaltimes,istatus)
                     ENDIF
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
                           bg_names(k)=TRIM(fname)//af
                           names_tmp(k) = names(i)
                           nvaltimes = nvaltimes + 1
                         endif
                        enddo
                        
                     endif
                  else 
                     bg_names(i)=names(i)(j+1:j+13)
                  endif
               else
                  if(i .eq. 1)write(6,*)
     1            ' index no slash in last 13 chars test was false:'
     1            ,j,names(i)(j:j),' ',names(i)(j+1:j+13)
               endif
            endif
         enddo

         if((bgmodel.eq.4) .or. (bgmodel .eq. 10))then
            bg_files = nvaltimes-1
            DO i=1,bg_files
              names(i) = names_tmp(i)
            ENDDO
         endif
      else ! c_model = 'LAPS' (exclusively using analyses)
         print*,'LAPS analysis selected as background'
         print*,'Checking lapsprd/lt1 subdirectory'
         cfilespec=bgpath(1:len)//'/lt1/*'
         call get_file_times(cfilespec,max_files,names,itimes
     +,bg_files,istatus)
         if(istatus.ne.1)then
            print*,'error status returned: get_file_times'
     +,' in get_acceptable_files',cfilespec(1:len+6)
            return
         endif

         do i=1,bg_files
            call s_len(names(i),j)
            call get_directory_length(names(i),lend)
            bg_names(i)=names(i)(lend+1:lend+9)//'0000'
         enddo         
       
      endif

      print *, "Total (nvaltimes) = ",nvaltimes
c ok, if we do not want 0-hr fcst files (analysis files) then filter them.
c ---------------------------------------------------------------------------
      if(cmodel.eq.'LAPS' .and. (.not.use_analysis))then
         print*,'Error in namelist variables'
         print*,'use_analysis must = true when cmodel = LAPS'
         return
      endif
      
      ij=0
      if(bg_files .gt.0)then
         if(.not.use_analysis)then
          allocate(bgnames_tmp(bg_files))
          call s_len(bg_names(1),j)
          do i=1,bg_files
             if(bg_names(i)(j-3:j).ne.'0000')then
                ij=ij+1
                bgnames_tmp(ij)=bg_names(i)
                IF ((bgmodel .EQ. 4).OR.(bgmodel.EQ.10))THEN
                  names_tmp(ij) = names(i)
                endif 
             endif
          enddo
          print*,'Removed ',bg_files-ij,' initial cond files '
          do i=1,bg_files
             bg_names(i)=' '
          enddo

          if(ij.eq.0)then
             print*,
     1       'WARNING: all files must be initial files. Return to main'
             print*,'Consider setting use_analysis parameter to true'
             print*,'Forecasts should also be added to the analyses'
             deallocate(bgnames_tmp)
             accepted_files = 0
             return
          endif

          bg_files=ij
          do i=1,bg_files
             bg_names(i)=bgnames_tmp(i)
             IF ((bgmodel .EQ. 4).OR.(bgmodel.EQ.10))THEN
               names(i) = names_tmp(i)
             ENDIF
             print *, TRIM(bg_names(i)),":", TRIM(names(i))
          enddo
          deallocate(bgnames_tmp)
         endif
      else
          print*,'Variable bg_files = 0, return to main'
          print*
          accepted_files=0
          return
      endif
c
c categorize file by init and fcst times. convert init time to i4time
c -------------------------------------------------------------------
      ij=0
      ibkgd=1
cwni      do n=1,bg_files-1
      do n=1,bg_files
C -- WNI-BLS ... Reworked this whole section
         if (n .EQ. 1) then
            bkgd(ibkgd)=bg_names(n)(1:9)
            ifcst_bkgd(ibkgd)=1
            call i4time_fname_lp(bkgd(ibkgd),i4timeinit(ibkgd),istatus)
            print *, "First: ",ibkgd, ifcst_bkgd(ibkgd)
         else
           finit1=bg_names(n-1)(1:9)
           finit2=bg_names(n)(1:9)
           if(finit1 .eq. finit2)then
C            Same as previous init tim
             ifcst_bkgd(ibkgd) = ifcst_bkgd(ibkgd) + 1
           else
C            Start a new init time
             ibkgd = ibkgd + 1
             bkgd(ibkgd)=bg_names(n)(1:9)
             ifcst_bkgd(ibkgd)=1 
             call i4time_fname_lp(bkgd(ibkgd),i4timeinit(ibkgd),istatus)
           endif

         endif

Curiously had to comment these lines as we now need 4 char to retain full filename
c        if(bgmodel.eq.0 .and. cmodel.eq.'LAPS_FUA'.or.
c    +       cmodel.eq.'MODEL_FUA'.or.cmodel.eq.'LAPS')then
c            fcst(ibkgd,ifcst_bkgd(ibkgd))=bg_names(n)(10:11)
c        else
             fcst(ibkgd,ifcst_bkgd(ibkgd))=bg_names(n)(10:13)
             IF (bgmodel .EQ. 4 .OR. bgmodel .EQ. 10) THEN
               IF (ifcst_bkgd(ibkgd).EQ.1)THEN
                 names_tmp(ibkgd) = names(n)
                 print *,"**** ",TRIM(names_tmp(ibkgd))
               ENDIF
             endif
c        endif

      enddo ! n

      print *,'Found ',ibkgd,' initial backgrounds '
      do n=1,ibkgd
        print *,'bkgd ',n, bkgd(n),' has ',ifcst_bkgd(n),' fcsts'
      enddo

c
c this section determines only the two fcsts bounding the analysis (i4time_anal) time.
c -----------------------------------------------------------------------------------
      n=1
      indx_for_best_init=0
      indx_for_best_init1=0
      indx_for_best_init2=0
      indx_for_best_fcst=0
      i4time_min_diff=1000000
      do while(n.le.ibkgd)
         if(i4timeinit(n).le.i4time_anal .and.
     +i4timeinit(n)+forecast_length*3600 .gt. i4time_anal)then   !let "forecast_length" window on init time qualify

            if(idebug .ge. 1)print*  
            print*,'Found bkgd init that corresp to anal: ',bkgd(n)
     1            ,' i4timeinit = ',i4timeinit(n)
            print*,'Num of fcst = ',ifcst_bkgd(n)

c this separate (near identical section) is due to local model having ffff = hhmm format
c whereas the second section below is ffff = hhhh.

            if(cmodel.eq.'LAPS_FUA'.or.cmodel.eq.'MODEL_FUA'
     +                             .or.cmodel.eq.'HRRR'      
     +                             .or.cmodel.eq.'RRs')then

c    +       (cmodel.eq.'LAPS') )then ! all cmodel types with bgmodel = 0 need this switch.

               do jj=2,ifcst_bkgd(n)
                  af=fcst(n,jj-1)
                  read(af,'(2i2)',err=888) ihour,imin
                  valid_time_1=i4timeinit(n)+ihour*3600+imin*60

                  af=fcst(n,jj)
                  read(af,'(2i2)',err=888) ihour,imin
                  valid_time_2=i4timeinit(n)+ihour*3600+imin*60

                  if(valid_time_1.le.i4time_anal.and.
     +               valid_time_2.ge.i4time_anal)then

                   if(abs(i4timeinit(n)-i4time_anal).lt.i4time_min_diff)
     +             then
                      i4time_min_diff=abs(i4timeinit(n)-i4time_anal)
                      indx_for_best_init=n
                      indx_for_best_fcst=jj-1
                   endif
                   print*,'Found fcsts bounding anal init/fcst='
     1                   , indx_for_best_init                
     1                   , indx_for_best_fcst                
                   print*,'Full name 1: ', bkgd(n),fcst(n,jj-1)
                   print*,'Full name 2: ', bkgd(n),fcst(n,jj)
                   print*,'i4times: fcst1/anal/fcst2',valid_time_1
     1                                   ,i4time_anal,valid_time_2
                   print*
                  endif
               enddo

            elseif(cmodel.ne.'LAPS')then

!!! add Huiling Yuan 20120904  AAA001, tested with RUC, use_analysis=true
             if ((use_analysis .eqv. .true.) .AND. 
     1           (use_forecast .eqv. .false.)      )then
              jj=1

              if(n.lt.ibkgd)then
                valid_time_1=i4timeinit(n)
                valid_time_2=i4timeinit(n+1)

                if(idebug.ge.1)write(6,*)bkgd(n),fcst(n,jj),' '
     1                                  ,bkgd(n),fcst(n+1,jj)  ,' '
     1                        ,valid_time_1,i4time_anal,valid_time_2

                if(valid_time_1.le.i4time_anal.and.
     +             valid_time_2.ge.i4time_anal)then

                 if(abs(i4timeinit(n)-i4time_anal).lt.i4time_min_diff)
     +then
                    i4time_min_diff=abs(i4timeinit(n)-i4time_anal)
                    indx_for_best_init=n+1
                    indx_for_best_fcst=jj
                 endif

                 print*,'Found fcsts bounding anal'
                 print*, "Index/init/fcst = ",n+1
     1                   , indx_for_best_init
     1                   , indx_for_best_fcst
                 print*,'Full name 1: ', bkgd(n),fcst(n,jj)
                 print*,'Full name 2: ', bkgd(n+1),fcst(n+1,jj)
                 write(6,*)valid_time_1,i4time_anal,valid_time_2
                 print*

                endif ! valid time

              endif   ! (n.lt.ibkgd)

!!!  partial if block, use_analysis=true,  by Huiling Yuan,  AAA001_A
             else   !! add by Huiling Yuan, AAA001_B

!             use_analysis = .false. .OR. use_forecast = .true.

              do jj=2,ifcst_bkgd(n)
                af=fcst(n,jj-1)
                read(af,'(i4)',err=888) ihour
                valid_time_1=i4timeinit(n)+ihour*3600
                af=fcst(n,jj)
                read(af,'(i4)',err=888) ihour
                valid_time_2=i4timeinit(n)+ihour*3600
                if(idebug.ge.1)write(6,*)bkgd(n),fcst(n,jj-1),' '
     1                                  ,bkgd(n),fcst(n,jj)  ,' '
     1                        ,valid_time_1,i4time_anal,valid_time_2

                if(valid_time_1.le.i4time_anal.and.
     +             valid_time_2.ge.i4time_anal)then

                 if(abs(i4timeinit(n)-i4time_anal).lt.i4time_min_diff)
     +then
                    i4time_min_diff=abs(i4timeinit(n)-i4time_anal)
                    indx_for_best_init=n
                    indx_for_best_fcst=jj-1
                 endif

                 print*,'Found fcsts bounding anal'
                 print*, "Index/init/fcst = ",n
     1                   , indx_for_best_init                
     1                   , indx_for_best_fcst                
                 print*,'Full name 1: ', bkgd(n),fcst(n,jj-1)
                 print*,'Full name 2: ', bkgd(n),fcst(n,jj)
                 write(6,*)valid_time_1,i4time_anal,valid_time_2
                 print*

                endif
              enddo

             endif   !use_analysis.eqv..true., end Huiling Yuan,  AAA001

            else  !must be LAPS

             if(n.gt.1)then
              valid_time_1=i4timeinit(n)
              valid_time_2=i4timeinit(n-1)
              if(valid_time_1.eq.i4time_anal)then
               indx_for_best_init1=n
               print*,'Found prev analysis corresp to current need'
               print*,'Full name 1: ', bkgd(indx_for_best_init1)
              endif
              if(valid_time_2.eq.i4time_anal)then
               indx_for_best_init2=n-1
               print*,'Found prev analysis corresp to current need'
               print*,'Full name 2: ', bkgd(indx_for_best_init2)
              endif

c             endif

              if(indx_for_best_init2.eq.0.and.
     +         indx_for_best_init1.eq.0)then
               if(valid_time_1+laps_cycle_time.eq.i4time_anal)then
                indx_for_best_init2=n
                indx_for_best_init1=ibkgd+1
                call i4time_fname_lp(names(i),i4time_fa,istatus)
                call make_fnam_lp(valid_time_1+laps_cycle_time
     +,bkgd(indx_for_best_init1),istatus)
 
                print*,'Found analyses corresp to anal+cycle_time'
                print*,'Full name 1: ', bkgd(indx_for_best_init1)
                print*,'Full name 2: ', bkgd(indx_for_best_init2)
               endif
              endif ! indx_for_best_init2

             endif ! n

            endif ! cmodel
         endif
         n=n+1
      enddo

c
c now restore the filename corresponding to the actual original "raw" file name 
c -----------------------------------------------------------------------------
      if(cmodel.ne.'LAPS')then
         if(indx_for_best_init.ne.0.and.indx_for_best_fcst.ne.0)then
            accepted_files = 2
           if ((use_analysis .eqv. .true.) .AND.
     1         (use_forecast .eqv. .false.))then   !! add  if block Huiling Yuan 20120905, AAA002
!!!!  names(1) is the earlier file, and names(2) is the latter file, Huiling Yuan, 20120905
            names(2)=bkgd(indx_for_best_init-1)//fcst(indx_for_best_init
     +        ,indx_for_best_fcst)
            names(1)=bkgd(indx_for_best_init)//fcst(indx_for_best_init
     +        ,indx_for_best_fcst)
           else
            names(2) = bkgd(indx_for_best_init)//fcst(indx_for_best_init
     +        ,indx_for_best_fcst)
            names(1) = bkgd(indx_for_best_init)//fcst(indx_for_best_init
     +        ,indx_for_best_fcst+1)
           endif     !! endif Huiling Yuan AAA002

            print*,'Accepted file 1: ',TRIM(names(1))
            print*,'Accepted file 2: ',TRIM(names(2))
         else
            print*,'!******************************************'
            print*,'!*** WARNING: Did not find acceptable files'
            print*,'!*** -------> Returning to main'
            print*,'!******************************************'
            print*,'indx_for_best_init/indx_for_best_fcst = '
     1            , indx_for_best_init,indx_for_best_fcst 
            accepted_files = 0
            bg_files = 0
            return
         endif

      else  ! cmodel = LAPS 
c At t = cycle time or t+1 cycle time we need to use previous analyses to
c advance fields in lga.
         if(indx_for_best_init1.ne.0.or.indx_for_best_init2.ne.0)then
            accepted_files=2
            IF (bgmodel.EQ.4 .or. bgmodel .eq. 10) THEN
              names(2) = names_tmp(indx_for_best_init2)
              names(1) = names_tmp(indx_for_best_init1)
            ELSE
              names(2)=bkgd(indx_for_best_init2)
              names(1)=bkgd(indx_for_best_init1)
            ENDIF
            print*,'Accepted file 1: ',TRIM(names(1))
            print*,'Accepted file 2: ',TRIM(names(2))
         else
            print*,'!******************************************'
            print*,'!*** WARNING: Did not find acceptable files'
            print*,'!*** -------> Returning to main'
            print*,'!******************************************'
            accepted_files = 0
            bg_files = 0
            return
         endif
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

      if(bgmodel .eq. 13 .and. lenf .ne. 13 .and. lenf .ne. 14)then
         write(6,*)' WARNING: filename length appears incorrect ',lenf
         write(6,*)' for example ',fullname(1:len+lenf)
         write(6,*)' remove .grib extension if present and/or use links'       
      endif

      write(6,*)
     1    ' Normal finish of get_acceptable_files, lenf/names(1) = '  
     1    ,lenf,TRIM(names(1))
      return

888   print*,'Error decoding fname ',fname
      print*,'n/bg_names(n) = ',n,bg_names(n)
      print*,'af = ',af

      return 
      end
