
      subroutine check_for_new_ctp(itime_ctp_window,istatus_ctp)

      implicit none
      character*200  cdir_ctp
      character*9    a9_systime
      integer        itime_ctp_window
      integer        istatus
      integer        istatus_ctp
      integer        i4time_latest
      integer        i4timedif
      integer        i4time_sys
      integer        ldctp


      istatus_ctp=0
      call get_directory('ctp',cdir_ctp,ldctp)
      cdir_ctp=cdir_ctp(1:ldctp)//'*'
      call get_latest_file_time(cdir_ctp,i4time_latest)
      call get_systime(i4time_sys,a9_systime,istatus)
      if(i4time_latest.eq.0)i4time_latest=i4time_sys  !if = 0 then no files in directory
      i4timedif=abs(i4time_sys-i4time_latest)         !if negative then maybe a case rerun
      print*,'Age of ctp data: ',i4timedif/60,' minutes' 
      if(i4timedif.le.itime_ctp_window)then
         print*,'Found new data'
         istatus_ctp=1
      endif
      return
      end
c
c--------------------------------------------------------
c
      subroutine read_cld_top_p(nxl,nyl,path_to_ctp
     &,rlctp,rlca,rlct,ri4time_ob,iwindow_ctp_s,i4time_data
     &,a9time_data,istatus)

      implicit none

      integer   max_files
      parameter (max_files=20000)
c
c the max_ctp parameter is set based upon the size of the
c cld-top-p files we read from FSL's /public data base.
      integer   max_ctp
      parameter (max_ctp=110000)

      character path_to_ctp*256
      character c_filenames(max_files)*256
      character cheader*256
      character c_dataroot*200
      character a24time_data*24
      character c_domain_name*10
      character c9timestring*9
      character a9_systime*9
      character a9time_data*9
      character cfname9*9
      character cyrjday*7
      character chrminsec*6

      integer   nxl,nyl
      integer   i,j,ii,jj
      integer   isave
      integer   isec
      integer   iobs
      integer   iwindow_ctp_s
      integer   i4timemin
      integer   i4timedif
      integer   i4time_data
      integer   i4time_sys
      integer   i4time_obs
      integer   i4time_file
      integer   i4timefile_min
      integer   i4timemx,i4timemn
      integer   inobs
      integer   lenp,lenf,lend
      integer   numoffiles
      integer   istatus
      integer   jstatus

      integer   lday(max_ctp)
      integer   ltim(max_ctp)
      integer   ca(max_ctp)

      real      lat(max_ctp)
      real      lon(max_ctp)
      real      pct(max_ctp)
      real      tc(max_ctp)

      real      rlat(nxl,nyl)
      real      rlon(nxl,nyl)
      real      rlctp(nxl,nyl)
      real      rlct(nxl,nyl)
      real      rlca(nxl,nyl)
      real      topo(nxl,nyl)
      real      ri,rj
      real*8    rtime
      real*8    rtimemx
      real*8    rtimemn
      real      r_missing_data
      real      pctobs
      real      ri4time_ob(nxl,nyl)
 
      call get_systime(i4time_sys,a9_systime,istatus)
      if(istatus .ne. 1)then
         print*,'error returned from get_systime'
         stop
      endif

      call get_file_names(path_to_ctp,numoffiles,c_filenames,
     +     max_files,istatus)
 
      call s_len(path_to_ctp,lenp)

      if(istatus.ne.1)then
         print*,'error returned from get_file_names:: ',
     +'path to cloud top pressure = ',path_to_ctp(1:lenp)
         return
      endif

      print*,'numoffiles returned from get_file_names = '
     +,numoffiles 

      call s_len(c_filenames(1),lenf)
c     print*,'first file name: ',c_filenames(1)(1:lenf)
c     print*,'last file name:  ',c_filenames(numoffiles)(1:lenf)
      print*
      print*,'find the file closest in time to systime'

      call get_directory_length(c_filenames(1),lend)

      isave = 0

      i4timemin=2000000000
      do i=1,numoffiles
         cfname9=c_filenames(i)(lend+1:lend+9)
         call i4time_fname_lp (cfname9,i4time_file,istatus)
         if(istatus.ne.1)then
            print*,'unable to convert to i4time'
            stop
         endif
         i4timedif=abs(i4time_sys-i4time_file)
         if(i4timedif.lt.i4timemin)then
            isave=i
            i4timefile_min=i4time_file
            i4timemin=i4timedif
         endif
      enddo

      print*,'index for minimum time: i= ',isave
      if(isave .gt. 0)then
          print*,'current cld-top-p file: ',c_filenames(isave)(1:lenf)
      else
          print*,'no current cld-top-p data - returning '
          istatus = 0
          return
      endif

c     if(i4timefile_min.ge.i4time_sys-iwindow_ctp_s.and.
c    &   i4timefile_min.le.i4time_sys+iwindow_ctp_s)then

        print*
        print*,'read current cld top p file'

        open(15,file=c_filenames(isave),form='formatted'
     &,status='old')

        iobs=0
        read(15,*,err=9)cheader
        do while (iobs.lt.max_ctp)
           iobs=iobs+1
           read(15,*,end=20,err=99)lday(iobs),ltim(iobs),lat(iobs)
     +,lon(iobs),ca(iobs),pct(iobs),tc(iobs)
           lon(iobs)=-1*lon(iobs)
        enddo
        if(iobs.eq.max_ctp)then
           print*,'WARNING: Entire file may not have been read'
        endif
20      iobs=iobs-1
        print*,'Done reading:',iobs, ' records' 
c
c determine which observations are within the analysis domain
c
        close(15)

        call find_domain_name(c_dataroot,c_domain_name,istatus)
        if(istatus.ne.1)then
           print*,'Error returned: find_domain_name'
           return
        endif

        call get_laps_domain(nxl,nyl,c_domain_name
     +,rlat,rlon,topo,istatus)
        if (istatus.lt.1)then
            print *,'Error reading lat, lon, topo data.'
            stop
        endif

c determine cld top p data spacing

        call get_r_missing_data(r_missing_data,istatus)

        rlctp=r_missing_data
        rlct=r_missing_data
        rlca=r_missing_data
        i4timemx=0
        i4timemn=2000000000

        do ii=1,iobs
      
         Call latlon_to_rlapsgrid(lat(ii),
     &                            lon(ii),
     &                            rlat,rlon,         !LAPS lat/lon arrays
     &                            nxl,nyl,           !LAPS horiz domain
     &                            ri,rj,             !Output: real i,j
     &                            jstatus)

         if(jstatus .eq. 1)then

c           call gdtost(pct(ii),1,1,ri,rj,rlctp(i,j),0)
c           call gdtost(float(ca(ii))/100.(ii),1,1,ri,rj
c    &,rlca(i,j),0)
c           call gdtost(tc(ii),1,1,ri,rj,rlct(i,j),0)

            i = nint(ri)
            j = nint(rj)

            rlctp(i,j)=pct(ii)*100.      !pa
            rlca(i,j)=float(ca(ii))/100. !unitless
            rlct(i,j)=tc(ii)

c make and save i4time
            write(cyrjday,'(i7.7)')lday(ii)
            do jj=1,7
               if(cyrjday(jj:jj).eq.' ')cyrjday(jj:jj)='0'
            enddo
            cyrjday=cyrjday(3:7)

            write(chrminsec,'(i6.6)')ltim(ii)
            do jj=1,6
               if(chrminsec(jj:jj).eq.' ')chrminsec(jj:jj)='0'
            enddo
            c9timestring=cyrjday(1:5)//chrminsec(1:4)
            read(chrminsec(5:6),'(i2)')isec

            call cv_asc_i4time(c9timestring,i4time_obs)
            ri4time_ob(i,j)=float(i4time_obs+isec)
            i4timemx=max(int(ri4time_ob(i,j)),i4timemx)
            i4timemn=min(int(ri4time_ob(i,j)),i4timemn)

         endif

        enddo

        inobs=0
        do j=1,nyl
        do i=1,nxl
           if(rlctp(i,j).ne.r_missing_data)then
              inobs=inobs+1
           endif
        enddo
        enddo
        pctobs=float(inobs)/float(nxl*nyl)
        if(pctobs.gt.0.0)then
           rtimemx=float(i4timemx)/1.0e6
           rtimemn=float(i4timemn)/1.0e6
           i4time_data=int(((rtimemx+rtimemn)*1.0e6)/2.0)
           call cv_i4tim_asc_lp(i4time_data,a24time_data,istatus)
           call make_fnam_lp(i4time_data,a9time_data,istatus)
           print*,'num/pct obs for domain ',inobs,pctobs
           print*,'time of data: ',a24time_data,' ',i4time_data
           print*,'returning data to main: fname data = ',a9time_data
        else
           istatus=0
           print*,'no cloud top p data in domain'
           print*,'no data processed for ctp: return to main'
        endif

c     else
c       istatus=0
c       print*,'no current cloud top p files: return to main'
c       print*,'no data processed'
c     endif

      goto 1000

9     print*,'Error reading cloud top p header'
      istatus = 0
      goto 1000
99    print*,'Error reading cloud top p file'
      istatus = 0

1000  return
      end
