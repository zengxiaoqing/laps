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
      subroutine vinterp(nx,ny,nz_bg,nz_laps,prlaps,
     .                   prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .                   htvi,tpvi,shvi,uwvi,vwvi)
c
      implicit none
      include 'bgdata.inc'
c
      integer nx,ny,nz_bg,nz_laps
c
c *** Input background variables.
c
      real*4 prbg(nx,ny,nz_bg),   !pressure (mb)
     .       tpbg(nx,ny,nz_bg),   !temperature (K)
     .       htbg(nx,ny,nz_bg),   !height (m)
     .       shbg(nx,ny,nz_bg),   !specific humidity (kg/kg)
     .       uwbg(nx,ny,nz_bg),   !u-wind (m/s)
     .       vwbg(nx,ny,nz_bg)    !v-wind (m/s)
c
c *** Output vertically interpolated variables.
c
      real*4 tpvi(nx,ny,nz_laps), !temperature (K)
     .       htvi(nx,ny,nz_laps), !height (m)
     .       shvi(nx,ny,nz_laps), !specific humidity (kg/kg)
     .       uwvi(nx,ny,nz_laps), !u-wind (m/s)
     .       vwvi(nx,ny,nz_laps)  !v-wind (m/s)
c
      real*4 prlaps(nz_laps),prilaps,fact1,fact2,fact3
      real*4 datmsg
      integer i,j,k,kk
c_______________________________________________________________________________
c

      do k=1,nz_laps
         prilaps=1./prlaps(k)
         do j=1,ny
            do i=1,nx
               do kk=1,nz_bg


                  if (prlaps(k) .gt. prbg(i,j,1)) then
                     datmsg = max(htbg(i,j,1),tpbg(i,j,1),shbg(i,j,1),
     +                    uwbg(i,j,1),vwbg(i,j,1))
                     if (datmsg .lt. missingflag) then
                        fact2=14.642857*alog(prbg(i,j,1)*prilaps)
                        tpvi(i,j,k)=tpbg(i,j,1)
     +                       +(prlaps(k)-prbg(i,j,1))*0.056
                        htvi(i,j,k)=htbg(i,j,1)
     +                       +(tpvi(i,j,k)+tpbg(i,j,1))*fact2
                        shvi(i,j,k)=shbg(i,j,1)
                        uwvi(i,j,k)=uwbg(i,j,1)
                        vwvi(i,j,k)=vwbg(i,j,1)
                     else
                        htvi(i,j,k)=missingflag
                        tpvi(i,j,k)=missingflag
                        shvi(i,j,k)=missingflag
                        uwvi(i,j,k)=missingflag
                        vwvi(i,j,k)=missingflag
                     endif
                     goto 10
                  elseif (prlaps(k) .lt. prbg(i,j,nz_bg)) then
                     datmsg = max(htbg(i,j,nz_bg),tpbg(i,j,nz_bg)
     +                ,shbg(i,j,nz_bg),uwbg(i,j,nz_bg),vwbg(i,j,nz_bg))
                     if (datmsg .lt. missingflag) then
                        fact2=29.285714*alog(prbg(i,j,nz_bg)*prilaps)
                        tpvi(i,j,k)=tpbg(i,j,nz_bg)
                        htvi(i,j,k)=htbg(i,j,nz_bg)
     +                       +tpbg(i,j,nz_bg)*fact2
                        shvi(i,j,k)=shbg(i,j,nz_bg)
                        uwvi(i,j,k)=uwbg(i,j,nz_bg)
                        vwvi(i,j,k)=vwbg(i,j,nz_bg)
                     else
                        htvi(i,j,k)=missingflag
                        tpvi(i,j,k)=missingflag
                        shvi(i,j,k)=missingflag
                        uwvi(i,j,k)=missingflag
                        vwvi(i,j,k)=missingflag
                     endif
                     goto 10
                  elseif (prlaps(k) .eq. prbg(i,j,kk)) then
                     datmsg = max(htbg(i,j,kk),tpbg(i,j,kk),
     +                    shbg(i,j,kk),uwbg(i,j,kk),vwbg(i,j,kk))
                     if (datmsg .lt. missingflag) then
                        htvi(i,j,k)=htbg(i,j,kk)
                        tpvi(i,j,k)=tpbg(i,j,kk)
                        shvi(i,j,k)=shbg(i,j,kk)
                        uwvi(i,j,k)=uwbg(i,j,kk)
                        vwvi(i,j,k)=vwbg(i,j,kk)
                     else
                        htvi(i,j,k)=missingflag
                        tpvi(i,j,k)=missingflag
                        shvi(i,j,k)=missingflag
                        uwvi(i,j,k)=missingflag
                        vwvi(i,j,k)=missingflag
                     endif
                     goto 10
                  elseif (prlaps(k) .lt. prbg(i,j,kk) .and. 
     .                    prlaps(k) .gt. prbg(i,j,kk+1)) then
                     if (datmsg .lt. missingflag) then
                        fact1=alog(prlaps(k)/prbg(i,j,kk))/
     .                       alog(prbg(i,j,kk+1)/prbg(i,j,kk))
                        fact2=14.642857*alog(prbg(i,j,kk)*prilaps)
                        fact3=(prlaps(k)-prbg(i,j,kk))/
     .                       (prbg(i,j,kk+1)-prbg(i,j,kk))
                        tpvi(i,j,k)=tpbg(i,j,kk)
     .                       +(tpbg(i,j,kk+1)-tpbg(i,j,kk))*fact1
                        htvi(i,j,k)=htbg(i,j,kk)
     .                       +(tpvi(i,j,k)+tpbg(i,j,kk))*fact2
                        shvi(i,j,k)=shbg(i,j,kk)
     .                       +(shbg(i,j,kk+1)-shbg(i,j,kk))*fact3
                        uwvi(i,j,k)=uwbg(i,j,kk)
     .                       +(uwbg(i,j,kk+1)-uwbg(i,j,kk))*fact3
                        vwvi(i,j,k)=vwbg(i,j,kk)
     .                       +(vwbg(i,j,kk+1)-vwbg(i,j,kk))*fact3
                     else
                        htvi(i,j,k)=missingflag
                        tpvi(i,j,k)=missingflag
                        shvi(i,j,k)=missingflag
                        uwvi(i,j,k)=missingflag
                        vwvi(i,j,k)=missingflag
                     endif
                     goto 10
                  endif
               enddo
 10            continue
            enddo
         enddo
      enddo
c
      return
c
      end
c
c===============================================================================
c

      subroutine time_interp(dir,ext,nx,ny,nz,kdim,pr,cycle_time,
     .                       time1,fcst1,time2,fcst2)
c
      implicit none
c
      integer nx,ny,nz,kdim
c
      integer time1,time2,
     .          fcst1,fcst2,
     .          cycle_time,
     .          ip(kdim),
     .          newfcst,
     .          imin,ihour,
     .          i,j,k,kk,istatus
c
      real*4 pr(nz),weight,
     .       grid1(nx,ny,kdim),
     .       grid2(nx,ny,kdim),
     .       gridn(nx,ny,kdim)
c
      character*(*)  dir
      character*(*)  ext
      character*3   var(kdim)
      character*4   lvl_coord(kdim)
      character*10  units(kdim)
      character*125 comment(kdim)
      character*9   fname
      character*4   af

c_______________________________________________________________________________
c
c *** Read in the two existing .lga files.
c
      if(ext.eq.'lga') then
         do k=1,nz
            ip(k)=int(pr(k))
            var(k)='HT '
            kk=k+nz
            ip(kk)=int(pr(k))
            var(kk)='T3 '
            kk=k+2*nz
            ip(kk)=int(pr(k))
            var(kk)='SH '
            kk=k+3*nz
            ip(kk)=int(pr(k))
            var(kk)='U3 '
            kk=k+4*nz
            ip(kk)=int(pr(k))
            var(kk)='V3 '
         enddo
      else
         do k=1,6
            ip(k)=0
         enddo
         var(1)='USF'
         var(2)='VSF'
         var(3)='TSF'
         var(4)='PSF'
         var(5)='SLP'
         var(6)='RSF'
      endif
         
c
      call read_laps(time1,time1+fcst1,dir,ext,
     .               nx,ny,kdim,kdim,var,
     .               ip,lvl_coord,units,comment,grid1,istatus)
c
      if(istatus.ne.1) then
         print *, 'ERROR returned from read_laps'
         stop 'lga_interp'
      endif
      call read_laps(time2,time2+fcst2,dir,ext,
     .               nx,ny,kdim,kdim,var,
     .               ip,lvl_coord,units,comment,grid2,istatus)
      if(istatus.ne.1) then
         print *, 'ERROR returned from read_laps'
         stop 'lga_interp'
      endif
c
c *** Do interpolation with time for each new file.
c
      newfcst=fcst2-cycle_time
10    continue
      weight=float(newfcst-fcst1)/float(fcst2-fcst1)
c
      do k=1,kdim
         do j=1,ny
            do i=1,nx
               gridn(i,j,k)=(1.-weight)*grid1(i,j,k)+weight*grid2(i,j,k)

               if(gridn(i,j,k) .gt. 1.e35) then
                  print *,'SEVERE Error in time_interp',
     +                 i,j,k,grid1(i,j,k),grid2(i,j,k),fcst1,fcst2
                  call erase_file(time1,fcst1,dir,ext)
                  call erase_file(time2,fcst2,dir,ext)

                  stop
               endif

            enddo
         enddo
      enddo
c
c *** Write out file.
c
      call make_fnam_lp(time1,fname,istatus)
      imin=mod(newfcst,3600)/60
      ihour=newfcst/3600
      write(af,'(2i2.2)') ihour,imin
      print *,'Writing - ',fname//af,'.'//ext//' (Backfill)'
      comment(1)='LI: '//comment(1)

      call write_laps(time1,time1+newfcst,dir,ext,
     .                nx,ny,nz,kdim,var,
     .                ip,lvl_coord,units,comment,gridn,istatus)
c
      newfcst=newfcst-cycle_time
      if (newfcst .gt. fcst1) goto 10
c
      return
      end
      subroutine erase_file(inittime,validtime,dir,ext)
      integer inittime,validtime, istatus, rename
      character*(*) dir, ext
      character*200 filename
      character*13 fname

      call make_fnam13_lp(inittime,validtime, fname, istatus)

      print*,inittime,validtime, fname

      call s_len(dir,len_dir)
      if(dir(len_dir:len_dir) .ne. '/') then
         dir(len_dir:len_dir)='/'
         len_dir=len_dir+1
      endif
      write(filename,*) dir(1:len_dir)//fname//'.'
      
      istatus = rename(filename(1:len_dir+15)//ext(1:3)
     +     ,filename(1:len_dir+15)//'bad')
      
      print*,'renamed ',filename(1:len_dir+15)//ext(1:3),
     +     ' with istatus= ',istatus

      return
      end
