      subroutine get_sounding_data_cdf(csat_id,
     &                         c_sounding_path,
     &                         i4time_latest,
     &                         c_filename_sat,
     &                         ires_x,ires_y,
     &                         nelems,nlines,nch,
     &                         i4snddata,
     &                         wavelength,
     &                         scalingBias,
     &                         scalingGain,
     &                         nw_pix,nw_line,
     &                         se_pix,se_line,
     &                         ewCycles,
     &                         ewIncs,
     &                         nsCycles,
     &                         nsIncs,
     &                         f_time,
     &                         lineTimeBeg,lineTimeEnd,
     &                         imcI4,
     &                         orbitAttitude,
     &                         ndimx,ndimy,ndimch,
     &                         istatus)
c
c subroutine determines if goes8/9 satellite sounder data are available
c in FSL's /public data base. These are raw GVAR satellite sounder
c files, 16 bit, 8 km resolution.
c It is then determined which sat files have not yet been processed
c for laps and returns these files for further mapping/navigation.
c It turns out that only the goes9 eastern sector covers the ROC-laps
c domain.
c 
c Any port of this software to another location, or perhaps a larger
c size than ROC-laps will require adjustments (software) to determine
c which sector properly maps into the domain.
c
c There is no preliminary installation step ATM to make this sector
c determination.
c
c J. Smart   2/96     Original Version
c J. Smart   10/97    Modify for newlaps structure;
c                     Modify for dynamic memory;
c                     Modify to accept both goes8 and goes9.
      implicit none

      integer     max_files
      parameter    (max_files=100)
      integer     i_snddata_age
      parameter    (i_snddata_age=20000)
      integer     nch
      integer     nelems,nlines

      character*255 c_filename_sat
      character*255 c_filename_lsr
      character*255 c_filename(max_files)
      character*150 c_lsr_dir
      character*255 pathname
      character*255 c_sounding_path
      character*9   c_filetime_sat
      character*9   c_filetime_lsr
      character*9   c_fname
      character*6   csat_id

      integer     i,j,k,n,jj
      integer     istatus
      integer     mstatus
      integer     gfn_status
      integer     nfiles_sat
      integer     numoffiles
c     integer*4     nfiles_lsr
c     integer*4     icnt
      integer     lend
      integer     ndimy,ndimch
      integer     ndimx(nlines)
      integer     imcI4
c     integer*4     nlsr
      integer     itstatus
      integer     init_timer
      integer     ishow_timer

      integer     i4time
      integer     i4time_proc
      integer     i4time_latest
      integer     i4time_nearest_lsr
      integer     i4time_now_gg
c     integer*4     i4time_lsr(max_files)
c     integer*4     i4time_sat(max_files)

      Real*8        orbitAttitude(336)
      Real*8        f_time
      Real*8        lineTimeBeg(nlines,nch)
      Real*8        lineTimeEnd(nlines,nch)
      Real*4        scalingBias(nlines,nch)
      Real*4        scalingGain(nlines,nch)

      Integer     ewCycles,ewIncs
      Integer     nsCycles,nsIncs
      Integer     nw_pix,nw_line
      Integer     se_pix,se_line
      Integer     ires_x,ires_y
c
      Integer     i4snddata(nelems,nlines,nch)
      REAL*8        wavelength(nch)
      Character*1   imc(4)
c
c ----------------------------- START --------------------------------------
c
      istatus = 1
c
c Find files for GOES-8 sounding data sector over Colorado. This is
c the _sdr file at ~ 45 or 46 past the hour..
c
      i4time_proc=i4time_now_gg()

      call get_directory('lsr',c_lsr_dir,lend)
      c_lsr_dir=c_lsr_dir(1:lend)//'/'//csat_id//'/'

      n=index(c_sounding_path,' ')-1
      c_filename_sat=c_sounding_path(1:n)//'*_sdr'
      call get_file_names(c_filename_sat,numoffiles,c_filename
     1        ,max_files,istatus)
      i4time_latest=0
      do i = 1,numoffiles
         j=index(c_filename(i),' ')-5
         jj=j-8
         if(csat_id.eq.'goes08')then
            if(c_filename(i)(j-1:j-1).eq.'4')then
               call cv_asc_i4time(c_filename(i)(jj:j),i4time)
               if(i4time.gt.i4time_latest)then
                  i4time_latest = i4time
               endif
            endif
         elseif(csat_id.eq.'goes09')then
            if(c_filename(i)(j-1:j-1).eq.'0')then
               call cv_asc_i4time(c_filename(i)(jj:j),i4time)
               if(i4time.gt.i4time_latest)then
                  i4time_latest = i4time
               endif
            endif
         else
            write(6,*)'Unknown Satellite Type - returning to main'
            goto 1000
         endif
      enddo
      call make_fnam_lp(i4time_latest,c_filetime_sat
     +,mstatus)
      c_filename_sat=c_sounding_path(1:n)//c_filetime_sat//'_sdr'

      if(i4time_proc-i4time_latest.gt.i_snddata_age)then
c
c put in wait for data here
c
         write(6,*)'Either data is too old or the correct sector'
     &,' is not available.' 
         write(6,*)'Returning to main without new sounder data'
         istatus = -1
         goto 1000
      endif
c
c We have "current" data. Determine if it has already been processed.
c Get latest lsr file (ie., in lapsprd/lsr)
c
      call get_file_time(c_lsr_dir,i4time_proc,
     +i4time_nearest_lsr)

      if(i4time_latest-i4time_nearest_lsr.le.0)then
         write(6,*)'lsr has the same time as the sat data.'
         write(6,*)'Returning to Main without new sounding data.'
         istatus = -1
         goto 1000
      endif

      call make_fnam_lp (i4time_nearest_lsr, c_filetime_lsr, mstatus)

c     call get_file_names(c_lsr_dir,
c    &                    nfiles_lsr,
c    &                    c_filename_lsr,
c    &                    max_files,
c    &                    gfn_status)
c     if(gfn_status.eq.1)then
c        write(*,*)'Success in GFN (laps lsr)'
c     else
c        write(6,*)'Error GFN (laps lsr)'
c        istatus=-1
c        goto 996
c     endif
c
c Need i4time.
c
c     nlsr = index(c_filename_lsr,' ')
c     do i = 1,nfiles_lsr
c        c_fname = c_filename_lsr(i)(nlsr-12:nlsr-4)
c        call cv_asc_i4time(c_fname,i4time_lsr(i))
c     enddo
c
c Get the times for the current sounding data
c
      n=index(c_filename_lsr,' ')
      write(6,*)'Filetime lsr: ',c_filetime_lsr
      n=index(c_filename_sat,' ')
      write(6,*)'Filetime Sat: ',c_filetime_sat

      n=index(c_sounding_path,' ')
      pathname=c_sounding_path(1:n-1)//c_filetime_sat//'_sdr'
      n=index(pathname,' ')
      write(6,*)'Get full satellite filename'
      write(*,*)'Satellite data filename: ',pathname(1:n-1)
c
c----------------------------------------------
c read the data
c
      write(6,*)'Read the sounding database '

      itstatus=init_timer()
      itstatus=ishow_timer()
      do k=1,nch
      do j=1,nlines
      do i=1,nelems
         i4snddata(i,j,k)=0.0
      enddo
      enddo
      enddo
         
      Call Read_sounder_db_cdf(c_filename_sat,
     &                         nelems,nlines,nch,
     &                         i4snddata,
     &                         wavelength,
     &                         scalingBias,
     &                         scalingGain,
     &                         nw_pix,nw_line,
     &                         se_pix,se_line,
     &                         ewCycles,
     &                         ewIncs,
     &                         nsCycles,
     &                         nsIncs,
     &                         f_time,
     &                         lineTimeBeg,lineTimeEnd,
     &                         imc,ires_x,ires_y,
     &                         orbitAttitude,
     &                         ndimx,ndimy,ndimch,
     &                         istatus)
      itstatus=ishow_timer()
      write(6,*)'Elapsed time (sec): ',itstatus
c
c !for sounder data the image motion compensation is off (ie.,=1)
c
      imcI4 = 1

      if(istatus.ne.1)then
         write(6,*)'Error occurred reading sounding db'
         write(6,*)'Getlsrdata.f'
         return
      endif
      write(6,*)'Sounding data obtained'

      goto 1000

996   write(6,*)'Error'

1000  return
      end
