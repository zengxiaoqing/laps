
      subroutine vinterp(nz_laps,nx,ny,nx_pr,ny_pr,
     .  ixmin,ixmax,iymin,iymax,
     .	nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,
     .  prlaps, prbght,prbgsh,prbguv,prbgww,
     .  htbg,tpbg,shbg,uwbg,vwbg,wwbg,
     .  htvi,tpvi,shvi,uwvi,vwvi,wwvi)
c
      implicit none
      include 'bgdata.inc'
c
      integer nx,ny
      integer nx_pr,ny_pr ! either model background dims, or 1,1 based on 
                          ! vertical grid
      integer ip,jp
      integer ixmin, ixmax, iymin, iymax ! bounding box of needed gridpoints
      integer ixtest, iytest
      integer nzbg_ht
      integer nzbg_tp
      integer nzbg_sh
      integer nzbg_uv
      integer nzbg_ww

      integer nzbght
      integer nzbgtp
      integer nzbgsh
      integer nzbguv
      integer nzbgww

      integer nz_laps
c
c *** Input background 3D variable.
c

      real   prbght(nx,ny,nzbg_ht),   !pressure (mb) of levels for hgt and T
     .       prbgsh(nx,ny,nzbg_sh),   !pressure (mb) of levels for rh
     .       prbguv(nx,ny,nzbg_uv),   !pressure (mb) of levels for u/v comps
     .       prbgww(nx,ny,nzbg_ww),   !pressure (mb) of levels for vertical wind comp
     .       tpbg(nx,ny,nzbg_tp),     !temperature (K)
     .       htbg(nx,ny,nzbg_ht),     !height (m)
     .       shbg(nx,ny,nzbg_sh),     !specific humidity (kg/kg)
     .       uwbg(nx,ny,nzbg_uv),     !u-wind (m/s)
     .       vwbg(nx,ny,nzbg_uv),     !v-wind (m/s)
     .       wwbg(nx,ny,nzbg_ww)      !w-wind (omega [pa/s])

c
c *** Output vertically interpolated variables.
c
      real   tpvi(nx,ny,nz_laps), !temperature (K)
     .       htvi(nx,ny,nz_laps), !height (m)
     .       shvi(nx,ny,nz_laps), !specific humidity (kg/kg)
     .       uwvi(nx,ny,nz_laps), !u-wind (m/s)
     .       vwvi(nx,ny,nz_laps), !v-wind (m/s)
     .       wwvi(nx,ny,nz_laps)  !w-wind (omega [pa/s])
c
!     3D pressures on the model grid (or input as a 1D constant pressure array)
      real   prlaps(nx_pr,ny_pr,nz_laps) ! pressure (mb)
      real   prilaps,fact1,fact2
      real   datmsg,datmsg1,datmsg2
      integer i,j,k,kk,lencm,istatus,ishow_timer
      integer kkguess,nguess_eq,nguess_int,noguess,nmiss_write,mode

      interface
         subroutine vinterp_sub(msngflag
     .,nx,ny,nx_pr,ny_pr,ixmin,ixmax,iymin,iymax
     .,nz,nzbg,pr,prbg,bgdata,bgdatavi)
         integer    nx,nx_pr,ny_pr,ny,nz
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

      write(6,*)' Start vinterp, # lvls (bg/laps) = ',nzbg_ht,nz_laps       

      if(nzbg_ht.ne.nzbg_tp)then
         print*,'vinterp requires nzbg_ht=nzbg_tp'
         print*,'values are: ',nzbg_ht,nzbg_tp       
         print*,'no interp performed. terminating'
         stop
      endif

      nzbght=nzbg_ht
      nzbgtp=nzbg_tp
      nzbgsh=nzbg_sh
      nzbguv=nzbg_uv
      nzbgww=nzbg_ww

      datmsg = 0.

      noguess = 0
      nguess_eq = 0
      nguess_int = 0

      nmiss_write = 0
      htvi = missingflag
      tpvi = missingflag

      do k=1,nz_laps
         kkguess = 1 ! default value
         do j=iymin,iymax ! 1,ny
            do i=ixmin,ixmax ! 1,nx
               ip = min(i,nx_pr) ! Collapse indices to 1,1 for 1D 'prlaps' input
               jp = min(j,ny_pr) 
               prilaps=1./prlaps(ip,jp,k)

c guessed pressure level
               kk=kkguess
               mode = 0 

c analysis pressure of level equals bg pressure of guessed level
               if (prlaps(ip,jp,k) .eq. prbght(i,j,kk)) then
                  datmsg = max(htbg(i,j,kk),tpbg(i,j,kk))
                  if (datmsg .lt. missingflag) then
                     htvi(i,j,k)=htbg(i,j,kk)
                     tpvi(i,j,k)=tpbg(i,j,kk)
                  else
                     htvi(i,j,k)=missingflag
                     tpvi(i,j,k)=missingflag
                  endif
                  nguess_eq = nguess_eq + 1
                  mode = 1
                  goto 10

c analysis pressure of level is inbetween bg pressures of guessed levels kk and kk+1
               elseif (prlaps(ip,jp,k) .lt. prbght(i,j,kk) .and. 
     +                 prlaps(ip,jp,k) .gt. prbght(i,j,kk+1)) then

                  datmsg1 = max(htbg(i,j,kk),tpbg(i,j,kk))
                  datmsg2 = max(htbg(i,j,kk+1),tpbg(i,j,kk+1))

                  if (datmsg1 .lt. missingflag.and.
     .                datmsg2 .lt. missingflag)then

                     fact1=alog(prlaps(ip,jp,k)/prbght(i,j,kk))/
     .                     alog(prbght(i,j,kk+1)/prbght(i,j,kk))
                     fact2=14.642857*alog(prbght(i,j,kk)*prilaps)

                     tpvi(i,j,k)=tpbg(i,j,kk)
     .                    +(tpbg(i,j,kk+1)-tpbg(i,j,kk))*fact1
                     htvi(i,j,k)=htbg(i,j,kk)
     .                    +(tpvi(i,j,k)+tpbg(i,j,kk))*fact2
                  else
                     htvi(i,j,k)=missingflag
                     tpvi(i,j,k)=missingflag
                  endif
                  nguess_int = nguess_int + 1
                  mode = 2
                  goto 10

               endif

               noguess = noguess + 1

c analysis pressure level is below lowest bg pressure level
               if (prlaps(ip,jp,k) .gt. prbght(i,j,1)) then
                  datmsg = max(htbg(i,j,1),tpbg(i,j,1))

                  if (datmsg .lt. missingflag) then ! extrapolate for T, Ht
                     fact2=14.642857*alog(prbght(i,j,1)*prilaps)

                     tpvi(i,j,k)=tpbg(i,j,1)
     +                    +(prlaps(ip,jp,k)-prbght(i,j,1))*0.056

                     htvi(i,j,k)=htbg(i,j,1)
     +                    +(tpvi(i,j,k)+tpbg(i,j,1))*fact2

                  else
                     htvi(i,j,k)=missingflag
                     tpvi(i,j,k)=missingflag
                  endif
                  mode = 3
                  goto 10

c analysis pressure level is above highest bg pressure level 
               elseif (prlaps(ip,jp,k) .lt. prbght(i,j,nzbght)) then       

                  datmsg = max(htbg(i,j,nzbght),tpbg(i,j,nzbght))

                  if (datmsg .lt. missingflag) then
                   fact2=29.285714*alog(prbght(i,j,nzbght)*prilaps)
                     tpvi(i,j,k)=tpbg(i,j,nzbght)
                     htvi(i,j,k)=htbg(i,j,nzbght)
     +                    +tpbg(i,j,nzbg_ht)*fact2
                  else
                     htvi(i,j,k)=missingflag
                     tpvi(i,j,k)=missingflag
                  endif
                  mode = 4
                  goto 10

               endif

               do kk=1,nzbght
c analysis pressure of level equals bg pressure of level
                  if (prlaps(ip,jp,k) .eq. prbght(i,j,kk)) then
                     datmsg = max(htbg(i,j,kk),tpbg(i,j,kk))
                     if (datmsg .lt. missingflag) then
                        htvi(i,j,k)=htbg(i,j,kk)
                        tpvi(i,j,k)=tpbg(i,j,kk)
                     else
                        htvi(i,j,k)=missingflag
                        tpvi(i,j,k)=missingflag
                     endif
                     kkguess=kk
                     mode = 5
                     goto 10

c analysis pressure of level is inbetween bg pressures of levels kk and kk+1
                  elseif (prlaps(ip,jp,k) .lt. prbght(i,j,kk) .and. 
     +                    prlaps(ip,jp,k) .gt. prbght(i,j,kk+1)) then

                     datmsg1 = max(htbg(i,j,kk),tpbg(i,j,kk))
                     datmsg2 = max(htbg(i,j,kk+1),tpbg(i,j,kk+1))

                     if (datmsg1 .lt. missingflag.and.
     .                   datmsg2 .lt. missingflag)then

                        fact1=alog(prlaps(ip,jp,k)/prbght(i,j,kk))/
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
                     kkguess=kk
                     mode = 6
                     goto 10
                  endif
               enddo ! kk
 10            continue

               if(htvi(i,j,k) .eq. missingflag)then
                  if(nmiss_write .le. 50)then
                     write(6,*)' WARNING: missing htvi ',i,j,k,mode
                     nmiss_write = nmiss_write + 1
                  endif
               endif

               if(tpvi(i,j,k) .lt. 150.)then
                  if(nmiss_write .le. 50)then
                     write(6,*)' WARNING: tpvi < 150. ',i,j,k,mode
                     nmiss_write = nmiss_write + 1
                  endif
               endif

            enddo ! i
         enddo ! j

         ixtest = (ixmin+ixmax)/2
         iytest = (iymin+iymax)/2
         if(ixtest .lt.  1 .OR. iytest .lt.  1 .OR. 
     1      ixtest .gt. nx .OR. iytest .gt. ny      )then
             write(6,*)' WARNING: ixtest/iytest out of bounds',ixtest
     1                                                        ,iytest
         else
             write(6,*)k,htvi(ixtest,iytest,k),tpvi(ixtest,iytest,k)
         endif

      enddo ! k

      write(6,*)' nguess_eq/nguess_int/noguess = ',nguess_eq,nguess_int
     1                                            ,noguess

      istatus=ishow_timer()

      write(6,*)' calls to vinterp_sub'
c
c second loops for remaining variables
c
      call vinterp_sub(missingflag,nx,ny,nx_pr,ny_pr
     .                     ,ixmin,ixmax,iymin,iymax,nz_laps,nzbgsh
     .                     ,prlaps,prbgsh,shbg,shvi)
      call vinterp_sub(missingflag,nx,ny,nx_pr,ny_pr
     .                     ,ixmin,ixmax,iymin,iymax,nz_laps,nzbguv
     .                     ,prlaps,prbguv,uwbg,uwvi)
      call vinterp_sub(missingflag,nx,ny,nx_pr,ny_pr
     .                     ,ixmin,ixmax,iymin,iymax,nz_laps,nzbguv
     .                     ,prlaps,prbguv,vwbg,vwvi)
      call vinterp_sub(missingflag,nx,ny,nx_pr,ny_pr
     .                     ,ixmin,ixmax,iymin,iymax,nz_laps,nzbgww
     .                     ,prlaps,prbgww,wwbg,wwvi)

      return
      end

      subroutine vinterp_sub(msngflag,nx,ny,nx_pr,ny_pr
     .          ,ixmin,ixmax,iymin,iymax,nz,nzbg
     .          ,pr,prbg,bgdata,bgdatavi)

      implicit none

      integer  nx,ny,nx_pr,ny_pr,nz
      integer  ixmin,ixmax,iymin,iymax ! bounding box of needed gridpoints
      integer  ip,jp
      integer  nzbg

      real, intent(in)  ::   pr(nx_pr,ny_pr,nz)
      real, intent(in)  ::   prbg(nx,ny,nzbg)
      real, intent(in)  ::   bgdata(nx,ny,nzbg)
      real, intent(out) ::   bgdatavi(nx,ny,nz)

      real     msngflag
      real     fact

      integer  i,j,k,kk,kkguess,nguess_eq,nguess_int,noguess

      noguess = 0
      nguess_eq = 0
      nguess_int = 0

      do k=1,nz
         kkguess = 1 ! default value
         do j=iymin,iymax ! 1,ny
            do i=ixmin,ixmax ! 1,nx
               ip = min(i,nx_pr) ! Collapse indices to 1,1 for 1D 'prlaps' input
               jp = min(j,ny_pr)

c guessed pressure level
               kk=kkguess

c analysis pressure of level equals bg pressure of guessed level
               if (pr(ip,jp,k) .eq. prbg(i,j,kk)) then
                     if (bgdata(i,j,kk) .lt. msngflag) then
                         bgdatavi(i,j,k)=bgdata(i,j,kk)
                     else
                         bgdatavi(i,j,k)=msngflag
                     endif
                     nguess_eq = nguess_eq + 1
                     goto 20

c analysis pressure of level is inbetween bg pressures of guessed levels kk and kk+1
               elseif (pr(ip,jp,k) .lt. prbg(i,j,kk) .and.
     +                 pr(ip,jp,k) .gt. prbg(i,j,kk+1)) then

                     if (bgdata(i,j,kk)   .lt. msngflag.and.
     .                   bgdata(i,j,kk+1) .lt. msngflag)then

                        fact=(pr(ip,jp,k)-prbg(i,j,kk))/
     .                       (prbg(i,j,kk+1)-prbg(i,j,kk))

                        bgdatavi(i,j,k)=bgdata(i,j,kk)
     .                       +(bgdata(i,j,kk+1)-bgdata(i,j,kk))*fact
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     nguess_int = nguess_int + 1
                     goto 20
               endif

               noguess = noguess + 1

c analysis pressure level is below lowest bg pressure level 
               if (pr(ip,jp,k) .gt. prbg(i,j,1)) then

                  if (bgdata(i,j,1) .lt. msngflag) then
                     bgdatavi(i,j,k)=bgdata(i,j,1)
                  else
                     bgdatavi(i,j,k)=msngflag
                  endif
                  goto 20


c analysis pressure level is above highest bg pressure level 
               elseif (pr(ip,jp,k) .lt. prbg(i,j,nzbg)) then

                  if (bgdata(i,j,nzbg) .lt. msngflag) then
                     bgdatavi(i,j,k)=bgdata(i,j,nzbg)
                  else
                     bgdatavi(i,j,k)=msngflag
                  endif
                  goto 20

               endif

               do kk=1,nzbg
c analysis pressure of level equals bg pressure of level
                  if (pr(ip,jp,k) .eq. prbg(i,j,kk)) then
                     if (bgdata(i,j,kk) .lt. msngflag) then
                        bgdatavi(i,j,k)=bgdata(i,j,kk)
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     kkguess=kk
                     goto 20

c analysis pressure of level is inbetween bg pressures of levels kk and kk+1
                  elseif (pr(ip,jp,k) .lt. prbg(i,j,kk) .and.
     +                    pr(ip,jp,k) .gt. prbg(i,j,kk+1)) then

                     if (bgdata(i,j,kk)   .lt. msngflag.and.
     .                   bgdata(i,j,kk+1) .lt. msngflag)then

                        fact=(pr(ip,jp,k)-prbg(i,j,kk))/
     .                       (prbg(i,j,kk+1)-prbg(i,j,kk))

                        bgdatavi(i,j,k)=bgdata(i,j,kk)
     .                       +(bgdata(i,j,kk+1)-bgdata(i,j,kk))*fact
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     kkguess=kk
                     goto 20
                  endif
               enddo ! kk
 20            continue
            enddo ! i
         enddo ! j
      enddo ! k

      write(6,*)' nguess_eq/nguess_int/noguess = ',nguess_eq,nguess_int
     1                                            ,noguess
c
      return
c
      end
