
      subroutine vinterp(nz_laps,nx,ny,nx_lp,ny_lp,
     .	nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,
     .  prlaps, prbght,prbgsh,prbguv,prbgww,
     .  htbg,tpbg,shbg,uwbg,vwbg,wwbg,
     .  htvi,tpvi,shvi,uwvi,vwvi,wwvi)
c
      implicit none
      include 'bgdata.inc'
c
      integer nx,ny,nx_lp,ny_lp
      integer ip,jp
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
      real   prlaps(nx_lp,ny_lp,nz_laps),prilaps,fact1,fact2
      real   datmsg,datmsg1,datmsg2
      integer i,j,k,kk,lencm

      interface
         subroutine vinterp_sub(msngflag
     .,nx,ny,nx_lp,ny_lp,nz,nzbg,pr,prbg,bgdata,bgdatavi)
         integer    nx,nx_lp,ny_lp,ny,nz
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

      if(nzbg_ht.ne.nzbg_tp)then
         print*,'vinterp requires nzbg_ht=nzbg_tp'
         print*,'no interp performed. terminating'
         stop
      endif

      nzbght=nzbg_ht
      nzbgtp=nzbg_tp
      nzbgsh=nzbg_sh
      nzbguv=nzbg_uv
      nzbgww=nzbg_ww

      datmsg = 0.
      do k=1,nz_laps
         do j=1,ny
            do i=1,nx
               ip = min(i,nx_lp)
               jp = min(j,ny_lp)
               prilaps=1./prlaps(ip,jp,k)
               do kk=1,nzbght

c lowest bg pressure level is above analysis lowest pressure levels
                  if (prlaps(ip,jp,k) .gt. prbght(i,j,1)) then
                     datmsg = max(htbg(i,j,1),tpbg(i,j,1))

                     if (datmsg .lt. missingflag) then
                        fact2=14.642857*alog(prbght(i,j,1)*prilaps)

                        tpvi(i,j,k)=tpbg(i,j,1)
     +                       +(prlaps(ip,jp,k)-prbght(i,j,1))*0.056

                        htvi(i,j,k)=htbg(i,j,1)
     +                       +(tpvi(i,j,k)+tpbg(i,j,1))*fact2

                     else
                        htvi(i,j,k)=missingflag
                        tpvi(i,j,k)=missingflag
                     endif
                     goto 10

c highest bg pressure level is below analysis highest pressure level
                  elseif (prlaps(ip,jp,k) .lt. prbght(i,j,nzbght)) then       

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
                  elseif (prlaps(ip,jp,k) .eq. prbght(i,j,kk)) then
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
      call vinterp_sub(missingflag,nx,ny,nx_lp,ny_lp,nz_laps,nzbgsh
     .                     ,prlaps,prbgsh,shbg,shvi)
      call vinterp_sub(missingflag,nx,ny,nx_lp,ny_lp,nz_laps,nzbguv
     .                     ,prlaps,prbguv,uwbg,uwvi)
      call vinterp_sub(missingflag,nx,ny,nx_lp,ny_lp,nz_laps,nzbguv
     .                     ,prlaps,prbguv,vwbg,vwvi)
      call vinterp_sub(missingflag,nx,ny,nx_lp,ny_lp,nz_laps,nzbgww
     .                     ,prlaps,prbgww,wwbg,wwvi)

      return
      end

      subroutine vinterp_sub(msngflag,nx,ny,nx_lp,ny_lp,nz,nzbg
     .          ,pr,prbg,bgdata,bgdatavi)

      implicit none

      integer  nx,ny,nx_lp,ny_lp,nz
      integer  ip,jp
      integer  nzbg

      real, intent(in)  ::   pr(nx_lp,ny_lp,nz)
      real, intent(in)  ::   prbg(nx,ny,nzbg)
      real, intent(in)  ::   bgdata(nx,ny,nzbg)
      real, intent(out) ::   bgdatavi(nx,ny,nz)

      real     msngflag
      real     fact

      integer  i,j,k,kk

      do k=1,nz
         do j=1,ny
            do i=1,nx
               ip = min(i,nx_lp)
               jp = min(j,ny_lp)
               do kk=1,nzbg


c lowest bg pressure level is above analysis lowest pressure levels
                  if (pr(ip,jp,k) .gt. prbg(i,j,1)) then

                     if (bgdata(i,j,1) .lt. msngflag) then
                        bgdatavi(i,j,k)=bgdata(i,j,1)
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     goto 20


c highest bg pressure level is below analysis highest pressure level
                  elseif (pr(ip,jp,k) .lt. prbg(i,j,nzbg)) then

                     if (bgdata(i,j,nzbg) .lt. msngflag) then
                        bgdatavi(i,j,k)=bgdata(i,j,nzbg)
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     goto 20

c analysis pressure of level equals bg pressure of level
                  elseif (pr(ip,jp,k) .eq. prbg(i,j,kk)) then
                     if (bgdata(i,j,kk) .lt. msngflag) then
                        bgdatavi(i,j,k)=bgdata(i,j,kk)
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
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
