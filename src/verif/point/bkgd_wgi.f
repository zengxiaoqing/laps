      subroutine bkgd_wgi(a9_time,i4time_init,i4time_fcst
     +,ext,cmdltype,balance,istatus)

      implicit none

      include 'grid_fname.cmn'

      character  clogdir*200
      character  cfilespec*200
      character  cline_wgi*200
      character  cmdltype*20
      character  a9_time*9
      character  a9time_init*9
      character  a9time_fcst*9
      character  ext*3
      character  cpads_type*3
      
      integer    istatus
      integer    balance
      integer    i4time
      integer    i4time_sys
      integer    i4time_init
      integer    i4time_fcst
      integer    lend,lens,lenw
      integer    hh,mm,hhmm
      integer    i,ii,lun
      integer    adv_anal_by_t_min
      logical    lexist,lrunbal,lfndch
      
c     call get_laps_config(grid_fnam_common,istatus)
c     if(istatus .ne. 1)then
c        print*,' error in get_laps_config'
c        goto 999
c     endif

      call get_directory('log',clogdir,lend)

      if(balance.eq.1)then
         call get_balance_nl(lrunbal,adv_anal_by_t_min,cpads_type
     1,istatus)
         call i4time_fname_lp(a9_time,i4time,istatus)
         i4time=i4time-adv_anal_by_t_min*60
         call make_fnam_lp(i4time,a9_time,istatus)
      endif

c this only works temporarily for Airdrop since there
c is only one cron job which is delayed purposefully.
c In the future, a9_time will be input from the verif
c software (ie., it computes the delay). Thus, we'll
c know which "wgi" file to open.
c     call get_systime(i4time_sys,a9_time,istatus)
c     if(istatus.ne.1)then
c        print*,'Error returned: get_systime'
c        return
c     endif

      cfilespec=clogdir(1:lend)//'wind.wgi.'//a9_time

      call s_len(cfilespec,lens)
      inquire(file=cfilespec,exist=lexist)
      if(lexist)then
         lun=20
         open(lun,file=cfilespec,status='old',form='formatted',
     +err=990)
         do i=1,7
            read(lun,100)cline_wgi(1:100)
            call downcase(cline_wgi,cline_wgi)
            if(cline_wgi(7:16).eq."determine")
     .      then
               istatus = 0
               print*,'Warning: bkgd not determined in ',cfilespec
               return
            endif
         enddo
100      format(100a)
         lfndch=.false.
         i=80
         do while(i.ge.1.and.(.not.lfndch))
            if(cline_wgi(i:i).ne.' ')then
               lenw=i
               lfndch=.true.
            endif
            i=i-1
         enddo
         print*,cline_wgi(lenw-12:lenw)
         close(lun)
      endif

c first, check to be sure a model background was determined in wgi
      if(cline_wgi(1:lenw).eq.'BACKGROUND FIELDS:')then
         print*,'Bkgd field not avail in ',cfilespec(1:lens)
         istatus=0
         return
      endif

      a9time_init=cline_wgi(lenw-12:lenw-4)  !a9_time(1:5)//cline_wgi(69:72)
      read(cline_wgi(lenw-3:lenw-2),'(i2)')hh
      read(cline_wgi(lenw-1:lenw),'(i2)')mm
      hhmm=hh*3600+mm*60

      call i4time_fname_lp(a9time_init,i4time_init,istatus)
      if(istatus.ne.1)goto 995
      i4time_fcst=i4time_init+hhmm
      call make_fnam_lp(i4time_fcst,a9time_fcst,istatus)
      if(istatus.ne.1)goto 995

      print*,'init  time ',a9time_init, i4time_init
      print*,'valid time ',a9time_fcst, i4time_fcst

      ext = cline_wgi(1:3)
      cmdltype=cline_wgi(5:lenw-14)
      istatus=1
      return  
995   print*,'Error converting a9 time string'
      print*,'a9time = ',a9time_init
      return
990   print*,'Error opening file ',cfilespec(1:lens)
      return
      end
