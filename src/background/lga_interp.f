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
      subroutine vinterp(nz_laps,nx,ny,
     .	nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww,
     .  prlaps, prbght,prbgsh,prbguv,prbgww,
     .  htbg,tpbg,shbg,uwbg,vwbg,wwbg,
     .  htvi,tpvi,shvi,uwvi,vwvi,wwvi)
c
      implicit none
      include 'bgdata.inc'
c
      integer nx,ny
      integer nzbg_ht
      integer nzbg_sh
      integer nzbg_uv
      integer nzbg_ww

      integer nzbght
      integer nzbgsh
      integer nzbguv
      integer nzbgww

      integer nz_laps
c
c *** Input background 3D variable.
c

      real*4 prbght(nx,ny,nzbg_ht),   !pressure (mb) of levels for hgt and T
     .       prbgsh(nx,ny,nzbg_sh),   !pressure (mb) of levels for rh
     .       prbguv(nx,ny,nzbg_uv),   !pressure (mb) of levels for u/v comps
     .       prbgww(nx,ny,nzbg_ww),   !pressure (mb) of levels for vertical wind comp
     .       tpbg(nx,ny,nzbg_ht),     !temperature (K)
     .       htbg(nx,ny,nzbg_ht),     !height (m)
     .       shbg(nx,ny,nzbg_sh),     !specific humidity (kg/kg)
     .       uwbg(nx,ny,nzbg_uv),     !u-wind (m/s)
     .       vwbg(nx,ny,nzbg_uv),     !v-wind (m/s)
     .       wwbg(nx,ny,nzbg_ww)      !w-wind (omega [pa/s])

c
c *** Output vertically interpolated variables.
c
      real*4 tpvi(nx,ny,nz_laps), !temperature (K)
     .       htvi(nx,ny,nz_laps), !height (m)
     .       shvi(nx,ny,nz_laps), !specific humidity (kg/kg)
     .       uwvi(nx,ny,nz_laps), !u-wind (m/s)
     .       vwvi(nx,ny,nz_laps), !v-wind (m/s)
     .       wwvi(nx,ny,nz_laps)  !w-wind (omega [pa/s])
c
      real*4 prlaps(nz_laps),prilaps,fact1,fact2
      real*4 datmsg,datmsg1,datmsg2
      integer i,j,k,kk,lencm

      interface
         subroutine vinterp_sub(msngflag
     .,nx,ny,nz,nzbg,pr,prbg,bgdata,bgdatavi)
         integer    nx,ny,nz
         integer    nzbg
         real  ::   pr(nz)
         real  ::   prbg(nx,ny,nzbg)
         real  ::   bgdata(nx,ny,nzbg)
         real  ::   bgdatavi(nx,ny,nz)
         real       msngflag
         end subroutine
      end interface  

c_______________________________________________________________________________
c
c first loop is required for getting the heights and temps.
c currently only SBN grids have variable pressure levels for
c individual fields (like sh, u/v and ww).


      nzbght=nzbg_ht
      nzbgsh=nzbg_sh
      nzbguv=nzbg_uv
      nzbgww=nzbg_ww

      datmsg = 0.
      do k=1,nz_laps
         prilaps=1./prlaps(k)
         do j=1,ny
            do i=1,nx
               do kk=1,nzbght

c lowest bg pressure level is above analysis lowest pressure levels
                  if (prlaps(k) .gt. prbght(i,j,1)) then
                     datmsg = max(htbg(i,j,1),tpbg(i,j,1))

                     if (datmsg .lt. missingflag) then
                        fact2=14.642857*alog(prbght(i,j,1)*prilaps)

                        tpvi(i,j,k)=tpbg(i,j,1)
     +                       +(prlaps(k)-prbght(i,j,1))*0.056

                        htvi(i,j,k)=htbg(i,j,1)
     +                       +(tpvi(i,j,k)+tpbg(i,j,1))*fact2

                     else
                        htvi(i,j,k)=missingflag
                        tpvi(i,j,k)=missingflag
                     endif
                     goto 10

c highest bg pressure level is below analysis highest pressure level
                  elseif (prlaps(k) .lt. prbght(i,j,nzbght)) then

                     datmsg = max(htbg(i,j,nzbght),tpbg(i,j,nzbght))

                     if (datmsg .lt. missingflag) then
                      fact2=29.285714*alog(prbght(i,j,nzbght)*prilaps)
                        tpvi(i,j,k)=tpbg(i,j,nzbght)
                        htvi(i,j,k)=htbg(i,j,nzbght)
     +                       +tpbg(i,j,nzbg_ht)*fact2
                     else
                        htvi(i,j,k)=missingflag
                        tpvi(i,j,k)=missingflag
                     endif
                     goto 10

c analysis pressure of level equals bg pressure of level
                  elseif (prlaps(k) .eq. prbght(i,j,kk)) then
                     datmsg = max(htbg(i,j,kk),tpbg(i,j,kk))
                     if (datmsg .lt. missingflag) then
                        htvi(i,j,k)=htbg(i,j,kk)
                        tpvi(i,j,k)=tpbg(i,j,kk)
                     else
                        htvi(i,j,k)=missingflag
                        tpvi(i,j,k)=missingflag
                     endif
                     goto 10

c analysis pressure of level is inbetween bg pressures of levels kk and kk+1
                  elseif (prlaps(k) .lt. prbght(i,j,kk) .and. 
     +                    prlaps(k) .gt. prbght(i,j,kk+1)) then

                     datmsg1 = max(htbg(i,j,kk),tpbg(i,j,kk))
                     datmsg2 = max(htbg(i,j,kk+1),tpbg(i,j,kk+1))

                     if (datmsg1 .lt. missingflag.and.
     .                   datmsg2 .lt. missingflag)then

                        fact1=alog(prlaps(k)/prbght(i,j,kk))/
     .                       alog(prbght(i,j,kk+1)/prbght(i,j,kk))
                        fact2=14.642857*alog(prbght(i,j,kk)*prilaps)

                        tpvi(i,j,k)=tpbg(i,j,kk)
     .                       +(tpbg(i,j,kk+1)-tpbg(i,j,kk))*fact1
                        htvi(i,j,k)=htbg(i,j,kk)
     .                       +(tpvi(i,j,k)+tpbg(i,j,kk))*fact2
                     else
                        htvi(i,j,k)=missingflag
                        tpvi(i,j,k)=missingflag
                     endif
                     goto 10
                  endif
               enddo
 10            continue
            enddo
         enddo
      enddo
c
c second loops for remaining variables
c
      call vinterp_sub(missingflag,nx,ny,nz_laps,nzbgsh
     .                     ,prlaps,prbgsh,shbg,shvi)
      call vinterp_sub(missingflag,nx,ny,nz_laps,nzbguv
     .                     ,prlaps,prbguv,uwbg,uwvi)
      call vinterp_sub(missingflag,nx,ny,nz_laps,nzbguv
     .                     ,prlaps,prbguv,vwbg,vwvi)
      call vinterp_sub(missingflag,nx,ny,nz_laps,nzbgww
     .                     ,prlaps,prbgww,wwbg,wwvi)

      return
      end

      subroutine vinterp_sub(msngflag,nx,ny,nz,nzbg
     .          ,pr,prbg,bgdata,bgdatavi)

      implicit none

      integer  nx,ny,nz
      integer  nzbg

      real, intent(in)  ::   pr(nz)
      real, intent(in)  ::   prbg(nx,ny,nzbg)
      real, intent(in)  ::   bgdata(nx,ny,nzbg)
      real, intent(out) ::   bgdatavi(nx,ny,nz)

      real     msngflag
      real     fact

      integer  i,j,k,kk

      do k=1,nz
         do j=1,ny
            do i=1,nx
               do kk=1,nzbg


c lowest bg pressure level is above analysis lowest pressure levels
                  if (pr(k) .gt. prbg(i,j,1)) then

                     if (bgdata(i,j,1) .lt. msngflag) then
                        bgdatavi(i,j,k)=bgdata(i,j,1)
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     goto 20


c highest bg pressure level is below analysis highest pressure level
                  elseif (pr(k) .lt. prbg(i,j,nzbg)) then

                     if (bgdata(i,j,nzbg) .lt. msngflag) then
                        bgdatavi(i,j,k)=bgdata(i,j,nzbg)
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     goto 20

c analysis pressure of level equals bg pressure of level
                  elseif (pr(k) .eq. prbg(i,j,kk)) then
                     if (bgdata(i,j,kk) .lt. msngflag) then
                        bgdatavi(i,j,k)=bgdata(i,j,kk)
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     goto 20

c analysis pressure of level is inbetween bg pressures of levels kk and kk+1
                  elseif (pr(k) .lt. prbg(i,j,kk) .and.
     +                    pr(k) .gt. prbg(i,j,kk+1)) then

                     if (bgdata(i,j,kk)   .lt. msngflag.and.
     .                   bgdata(i,j,kk+1) .lt. msngflag)then

                        fact=(pr(k)-prbg(i,j,kk))/
     .                       (prbg(i,j,kk+1)-prbg(i,j,kk))

                        bgdatavi(i,j,k)=bgdata(i,j,kk)
     .                       +(bgdata(i,j,kk+1)-bgdata(i,j,kk))*fact
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     goto 20
                  endif
               enddo
 20            continue
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
      subroutine time_interp(dir,ext,nx,ny,nz,ngrids,pr,
     .                  i4time_valid1,i4time_valid2,i4time_now,
     .                  time1,fcst1,time2,fcst2)
c
      implicit none
      include 'bgdata.inc'
c
      integer nx,ny,nz
      integer n,ngrids
      integer warncnt
c
      integer   i4time_valid1,
     .          i4time_valid2,
     .          i4time_now,
     .          time1,time2,
     .          fcst1,fcst2,
     .          ip(nz),
     .          newfcst,
     .          imin,ihour,
     .          i,j,k,kk,
     .          istatus,nstatus
c
      real*4 pr(nz),weight,
     .       grid1(nx,ny,nz),
     .       grid2(nx,ny,nz),
     .       gridn(nx,ny,nz)
c
      integer nan
      character*(*)  dir
      character*(*)  ext
      character*3    var(nz,ngrids)
      character*4    lvl_coord(nz)
      character*10   units(nz)
      character*125  comment(nz)
      character*9    fname9
      character*4    af

c_______________________________________________________________________________
c
      if(ext.eq.'lga') then
         do k=1,nz
            ip(k)=int(pr(k))
            var(k,1)='HT '
            var(k,2)='T3 '
            var(k,3)='SH '
            var(k,4)='U3 '
            var(k,5)='V3 '
            var(k,6)='OM '
         enddo
      else
         do k=1,nz
            ip(k)=0
         enddo
         var(1,1)='USF'
         var(1,2)='VSF'
         var(1,3)='TSF'
         var(1,4)='PSF'
         var(1,5)='SLP'
         var(1,6)='RSF'
         var(1,7)='DSF'
      endif

      newfcst=fcst1-(i4time_valid2-i4time_now)
      weight=float(i4time_valid2-i4time_now)/
     .       float(i4time_valid2-i4time_valid1)
      print*,'Time interp weight = ',weight

      write(af,'(i4.4)') fcst1/3600
      call  make_fnam_lp(time1, fname9, nstatus)
      print*,'Reading: ',fname9,af,'.'//ext
      write(af,'(i4.4)') fcst2/3600
      call  make_fnam_lp(time2, fname9, nstatus)
      print*,'Reading: ',fname9,af,'.'//ext

      call make_fnam_lp(time1,fname9,istatus)
      imin=mod(newfcst,3600)/60
      ihour=newfcst/3600
      write(af,'(2i2.2)') ihour,imin
      print *,'Writing - ',fname9//af,'.'//ext//' (Backfill)'

      do n=1,ngrids

         call read_laps(time1,time1+fcst1,dir,ext,
     .        nx,ny,nz,nz,var(1,n),ip,lvl_coord,units,comment,
     .        grid1,istatus)
c
         if(istatus.ne.1) then
            print *, 'ERROR returned from read_laps'
            stop 'lga_interp'
         endif

         call read_laps(time2,time2+fcst2,dir,ext,
     .        nx,ny,nz,nz,var(1,n),ip,lvl_coord,units,comment,
     .        grid2,istatus)
         if(istatus.ne.1) then
            print *, 'ERROR returned from read_laps'
            stop 'lga_interp'
         endif
c
c *** Do interpolation with time for each new file.
c
         warncnt = 0
         do k=1,nz
         do j=1,ny
            do i=1,nx
               if(nan(grid1(i,j,k))+nan(grid2(i,j,k)).gt.0 .or.
     +              grid1(i,j,k).ge.missingflag .or.
     +              grid2(i,j,k).ge.missingflag) then

                    if(warncnt.eq. 0)then
                       print*,'Missingflag at ',i,j,k,
     +                 grid1(i,j,k),grid2(i,j,k)
                       warncnt = 1
                    endif
                    gridn(i,j,k) = missingflag
               else
               
                  gridn(i,j,k)= (1.- weight)*grid1(i,j,k) +
     +                               weight *grid2(i,j,k)

               endif

            enddo
         enddo
         enddo
c
         comment(1) = 'Time Interpolated: '//comment(1)
         if(ext.eq.'lga')then
            call write_laps(time1,time1+newfcst,dir,ext,
     .           nx,ny,nz,nz,var(1,n),ip,lvl_coord,units,comment,
     .           gridn,istatus)
         else
            call write_laps(time1,time1+newfcst,dir,ext,
     .           nx,ny,1,1,var(1,n),ip,lvl_coord,units,comment,
     .           gridn,istatus)
         endif
          
      enddo
c
      return
      end

      subroutine erase_file(inittime,validtime,dir,ext)
      integer inittime,validtime, istatus, rename
      character*(*) dir, ext
      character*256 filename
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
