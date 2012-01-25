
      subroutine vinterp_ht(nz_laps,nx,ny,
     .	nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,
     .  htlaps, prbght,prbgsh,prbguv,prbgww,
     .  htbg,tpbg,shbg,uwbg,vwbg,wwbg,
     .  prvi,tpvi,shvi,uwvi,vwvi,wwvi)
c
      implicit none
      include 'bgdata.inc'
      include 'constants.inc'
c
      integer nx,ny
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
     .       prvi(nx,ny,nz_laps), !pressure (mb)
     .       shvi(nx,ny,nz_laps), !specific humidity (kg/kg)
     .       uwvi(nx,ny,nz_laps), !u-wind (m/s)
     .       vwvi(nx,ny,nz_laps), !v-wind (m/s)
     .       wwvi(nx,ny,nz_laps)  !w-wind (omega [pa/s])
c
      real   htlaps(nx,ny,nz_laps),fact1,fact2,dlogp_dz,grinv2   
      real   datmsg,datmsg1,datmsg2
      integer i,j,k,kk,lencm

      interface
         subroutine vinterp_ht_sub(msngflag
     .,nx,ny,nz,nzbg,ht,htbg,bgdata,bgdatavi)
         integer    nx,ny,nz
         integer    nzbg
         real  ::   ht(nx,ny,nz)
         real  ::   htbg(nx,ny,nzbg)
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

      write(6,*)' Start vinterp_ht'

      if(nzbg_ht.ne.nzbg_tp)then
         print*,'vinterp_ht requires nzbg_ht=nzbg_tp'
         print*,'no interp performed. terminating'
         stop
      endif

      if(nzbg_ht.ne.nzbg_sh)then
         print*,'vinterp_ht requires nzbg_ht=nzbg_sh'
         print*,'no interp performed. terminating'
         stop
      endif

      if(nzbg_ht.ne.nzbg_uv)then
         print*,'vinterp_ht requires nzbg_ht=nzbg_uv'
         print*,'no interp performed. terminating'
         stop
      endif

      if(nzbg_ht.ne.nzbg_ww)then
         print*,'vinterp_ht requires nzbg_ht=nzbg_ww'
         print*,'no interp performed. terminating'
         stop
      endif

      nzbght=nzbg_ht
      nzbgtp=nzbg_tp
      nzbgsh=nzbg_sh
      nzbguv=nzbg_uv
      nzbgww=nzbg_ww

      datmsg = 0.

c     In 'vinterp' we derived the heights/temps by interpolating from the 
c     pressures done in log P space. Extrapolation is also done.

c     In 'vinterp_ht' we first get temperatures calling the 'vinterp_ht_sub' 
c     We then want to derive the (log) pressures by interpolating from the 
c     heights? Assume dlogp/dz is a approximately constant, though we can
c     refine this knowing the lapse rate. This is given by the hypsometric
c     equation in the from dlogp/dz = g/RTbar. We assume the model is 
c     hydrostatic when we interpolate between levels. We can change this to
c     non-hydrostatic by explicitly calculating dlogp/dz from the model P,Ht
c     fields instead of using Tbar.

      call vinterp_ht_sub(missingflag,nx,ny,nz_laps,nzbgtp
     .                   ,htlaps,htbg,tpbg,tpvi)

      grinv2 = 2. * (grav / r_d) 

      do k=1,nz_laps
         do j=1,ny
            do i=1,nx
               do kk=1,nzbght-1

c analysis height level is below lowest bg height level 
                  if (htlaps(i,j,k) .lt. htbg(i,j,1)) then
                     datmsg = max(prbght(i,j,1),tpbg(i,j,1))

                     if (datmsg .lt. missingflag) then
                        dlogp_dz = -grinv2 / (tpvi(i,j,k)+tpbg(i,j,1))
                        prvi(i,j,k) = prbght(i,j,1) 
     1                    * exp(dlogp_dz * (htlaps(i,j,k)-htbg(i,j,1)))       
                     else
                        prvi(i,j,k)=missingflag
                     endif

                     goto 10

c analysis height level is above highest bg height level 
                  elseif (htlaps(i,j,k) .gt. htbg(i,j,nzbght)) then

                     datmsg = max(htbg(i,j,nzbght),tpbg(i,j,nzbght))

                     if (datmsg .lt. missingflag) then
                        dlogp_dz =-grinv2/(tpvi(i,j,k)+tpbg(i,j,nzbght))              
                        prvi(i,j,k)=prbght(i,j,nzbght)
     1                   *exp(dlogp_dz*(htlaps(i,j,k)-htbg(i,j,nzbght))) 
                     else
                        prvi(i,j,k)=missingflag
                     endif

                     goto 10

c analysis height of level equals bg height of level
                  elseif (htlaps(i,j,k) .eq. htbg(i,j,kk)) then
                     datmsg = max(htbg(i,j,kk),tpbg(i,j,kk))

                     if (datmsg .lt. missingflag) then
                        prvi(i,j,k)=prbght(i,j,kk)
                     else
                        prvi(i,j,k)=missingflag
                     endif

                     goto 10

c analysis height of level is inbetween bg height of levels kk and kk+1
                  elseif (htlaps(i,j,k) .gt. htbg(i,j,kk) .and. 
     +                    htlaps(i,j,k) .lt. htbg(i,j,kk+1)) then

                     datmsg = max(htbg(i,j,kk),tpbg(i,j,kk))

                     if (datmsg .lt. missingflag)then
                        dlogp_dz = -grinv2 / (tpvi(i,j,k)+tpbg(i,j,kk))       
                        prvi(i,j,k) = prbght(i,j,kk) 
     1                    * exp(dlogp_dz * (htlaps(i,j,k)-htbg(i,j,kk)))
                     else
                        prvi(i,j,k)=missingflag
                     endif

                     goto 10

                  endif
               enddo ! kk
 10            continue 
            enddo ! i
         enddo ! j
      enddo ! k
c
c second loops for remaining variables
c
      call vinterp_ht_sub(missingflag,nx,ny,nz_laps,nzbgsh
     .                     ,htlaps,htbg,shbg,shvi)
      call vinterp_ht_sub(missingflag,nx,ny,nz_laps,nzbguv
     .                     ,htlaps,htbg,uwbg,uwvi)
      call vinterp_ht_sub(missingflag,nx,ny,nz_laps,nzbguv
     .                     ,htlaps,htbg,vwbg,vwvi)
      call vinterp_ht_sub(missingflag,nx,ny,nz_laps,nzbgww
     .                     ,htlaps,htbg,wwbg,wwvi)

      return
      end

      subroutine vinterp_ht_sub(msngflag,nx,ny,nz,nzbg
     .          ,ht,htbg,bgdata,bgdatavi)

      implicit none

      integer  nx,ny,nz
      integer  nzbg

      real, intent(in)  ::   ht(nx,ny,nz)
      real, intent(in)  ::   htbg(nx,ny,nzbg)
      real, intent(in)  ::   bgdata(nx,ny,nzbg)
      real, intent(out) ::   bgdatavi(nx,ny,nz)

      real     msngflag
      real     fact

      integer  i,j,k,kk

      do k=1,nz
         do j=1,ny
            do i=1,nx
               do kk=1,nzbg-1


c analysis height level is below lowest bg height level 
                  if (ht(i,j,k) .lt. htbg(i,j,1)) then

                     if (bgdata(i,j,1) .lt. msngflag) then
                        bgdatavi(i,j,k)=bgdata(i,j,1)
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     goto 20


c analysis height level is above highest bg height level 
                  elseif (ht(i,j,k) .gt. htbg(i,j,nzbg)) then

                     if (bgdata(i,j,nzbg) .lt. msngflag) then
                        bgdatavi(i,j,k)=bgdata(i,j,nzbg)
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     goto 20

c analysis height of level equals bg height of level
                  elseif (ht(i,j,k) .eq. htbg(i,j,kk)) then
                     if (bgdata(i,j,kk) .lt. msngflag) then
                        bgdatavi(i,j,k)=bgdata(i,j,kk)
                     else
                        bgdatavi(i,j,k)=msngflag
                     endif
                     goto 20

c analysis height of level is inbetween bg heights of levels kk and kk+1
                  elseif (ht(i,j,k) .gt. htbg(i,j,kk) .and.
     +                    ht(i,j,k) .lt. htbg(i,j,kk+1)) then

                     if (bgdata(i,j,kk)   .lt. msngflag.and.
     .                   bgdata(i,j,kk+1) .lt. msngflag)then

                        fact=(ht(i,j,k)-htbg(i,j,kk))/
     .                       (htbg(i,j,kk+1)-htbg(i,j,kk))

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
