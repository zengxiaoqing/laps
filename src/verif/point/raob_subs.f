      subroutine get_laps(nx,ny,nk,dir_in,i4time_syn,a9_time,
     1     uLapsGS,vLapsGS,tLapsGS,rhLapsGS,htLapsGS,
     1     uLapsGP,vLapsGP,tLapsGP,rhLapsGP,
     1     htLapsGP,htLgaGS,htLgaGP,
     1     balance,laps_levels,status)

      implicit none
      include 'netcdf.inc'

      integer       nx, ny, nk
      character*(*)   dir_in
      integer       i4time_syn,status(2,6)
      real          uLapsGS(nx,ny,nk), vLapsGS(nx,ny,nk), 
     1                tLapsGS(nx,ny,nk), rhLapsGS(nx,ny,nk),
     1                htLapsGS(nx,ny,nk),uLapsGP(nx,ny,nk), 
     1                vLapsGP(nx,ny,nk),tLapsGP(nx,ny,nk), 
     1                rhLapsGP(nx,ny,nk), htLapsGP(nx,ny,nk),
     1                htLgaGS(nx,ny,nk), htLgaGP(nx,ny,nk)
      real          laps_levels(nk) !laps pressure levels in mb
      real          make_rh

      character*4     lvl_coord(nk)
      character*3     var(nk), ext
      character*10    units(nk)
      character*20    cmdltype
      character*125   comment(nk)
      character*180   fdir_in
      character*9     fnameSyn,fnamePrev,fnamePPrev,a9_time
      integer       i, j, k, lvl(nk), fdir_len, i4time_Pprev,
     1                i4time_prev, len_dir_in, istatus,
     1                balance
      integer         i4time_init,i4time_fcst

! begin
      i4time_prev = i4time_syn - 3600
      call make_fnam_lp(i4time_syn,fnameSyn,istatus)
      call make_fnam_lp(i4time_prev,fnamePrev,istatus)

! status(1,x) = Syn  (2,x) = Prev
! status(x,1)=u (x,2)=v (x,3)=t (x,4)=ht (x,5)=rh (x,6)=LGAht

      status = 1

      call s_len(dir_in,len_dir_in)

      j = 6

      var='ht'
      lvl=laps_levels
c
c JS: this needs to consider nest7grid.parms fdda_model_source
c     since that can be our background too. Done: 10-11-02 with
c     mod to wind3d.pl
c     ext = 'lga'

      fdir_in = dir_in(1:len_dir_in)//'lga/'
      call s_len(fdir_in,fdir_len)

c     write(6,*) dir_in
c     write(6,*) fdir_in

c JS: need to determine what background was actually used in the
c     analysis and configure i4time_syn appropriately (likely will
c     need an i4time_syn_valid too as second argument to this sub)
C     read LGA HT valid i4time_syn

      call bkgd_wgi(a9_time,i4time_init,i4time_fcst
     +,ext,cmdltype,balance,istatus)
      if(istatus.ne.1)then
         print*,'Failure in bkgd_wgi to get model bkgd time'
         status=0
         return
      endif

      call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,htLgaGS,istatus)

c     call get_modelfg_3d(i4time_syn,var(1),nx,ny,nk,htLgaGS
c    1,istatus)

      if (istatus .ne. 1) then

        write (6,*) ' Error in readlapsdata for LGA HT Syn',
     1             fdir_in(1:fdir_len)//fnameSyn//'0000.'//ext(1:3)
      endif

C       try previous hr run, 1 hr fcst
c       call read_laps(i4time_prev,i4time_syn,fdir_in,ext,
c    1                      nx,ny,nk,nk,var,lvl,lvl_coord,
c    1                      units,comment,htLgaGS,istatus)
c       call get_modelfg_3d(i4time_prev,var(1),nx,ny,nk,htLgaGS
c    1,istatus)
c       if (istatus .ne. 1) then
c         write (6,*) ' Error in readlapsdata for LGA HT Syn',
c    1            fdir_in(1:fdir_len)//fnamePrev//'0100.'//ext(1:3)
c         status(1,j) = 0
c       else
c         write(6,*) 'Success reading LGA HT from ',
c    1            fdir_in(1:fdir_len)//fnamePrev//'0100.'//ext(1:3)
c       endif
c     else
c       write(6,*) 'Success reading LGA HT from ',
c    1             fdir_in(1:fdir_len)//fnameSyn//'0000.'//ext(1:3)
c     endif
C     read LGA HT valid i4time_prev
c JS: probably don't need to get the previous hour lga.  If the
c     current one doesn't exist, then we cannot verify. 
c     call read_laps(i4time_prev,i4time_prev,fdir_in,ext,
c    1                    nx,ny,nk,nk,var,lvl,lvl_coord,
c    1                    units,comment,htLgaGP,istatus)
c     call get_modelfg_3d(i4time_prev,var(1),nx,ny,nk,htLgaGP
c    1,istatus)
c     if (istatus .ne. 1) then
c       write (6,*) ' Error in readlapsdata for LGA HT Prev',
c    1             fdir_in(1:fdir_len)//fnamePrev//'0000.'//ext(1:3)
c       try previous hr run, 1 hr fcst
c       i4time_Pprev = i4time_prev-3600
c       call make_fnam_lp(i4time_Pprev,fnamePPrev,istatus)
c       call read_laps(i4time_Pprev,i4time_prev,fdir_in,ext,
c    1                      nx,ny,nk,nk,var,lvl,lvl_coord,
c    1                      units,comment,htLgaGP,istatus)
c       call get_modelfg_3d(i4time_Pprev,var(1),nx,ny,nk,htLgaGP
c    1,istatus)
c       if (istatus .ne. 1) then
c         write (6,*) ' Error in readlapsdata for LGA HT PPrev',
c    1            fdir_in(1:fdir_len)//fnamePPrev//'0100.'//ext(1:3)
c         status(2,j) = 0
c       else
c         write(6,*) 'Success reading LGA HT from ',
c    1            fdir_in(1:fdir_len)//fnamePPrev//'0100.'//ext(1:3)
c       endif
c     else
c       write(6,*) 'Success reading LGA HT from ',
c    1             fdir_in(1:fdir_len)//fnamePrev//'0000.'//ext(1:3)
c     endif

c JS: new code for verifying background with sounding or profiler
c   we'll use balance code = 2 to indicate verification of background.

      if(balance.eq.2)then

         ext = 'lga'
         var = 't3'
         call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,tLapsGS,istatus)
         if(istatus.ne.1)then
            write (6,*) ' Error in readlapsdata for LGA HT PPrev',
     1fdir_in(1:fdir_len)//fnamePPrev//'0100.'//ext(1:3)
            status(1,3)=0
            return
         endif
         var='ht'
         call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,htLapsGS,istatus)
         if(istatus.ne.1)then
            write (6,*) ' Error in readlapsdata for LGA HT PPrev',
     1fdir_in(1:fdir_len)//fnamePPrev//'0100.'//ext(1:3)
            status(1,4)=0
            return
         endif
         var='u3'
         call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,uLapsGS,istatus)
         if(istatus.ne.1)then
            write (6,*) ' Error in readlapsdata for LGA HT PPrev',
     1fdir_in(1:fdir_len)//fnamePPrev//'0100.'//ext(1:3)
            status(1,1)=0
            return
         endif
         var='v3'
         call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,vLapsGS,istatus)
         if(istatus.ne.1)then
            write (6,*) ' Error in readlapsdata for LGA HT PPrev',
     1fdir_in(1:fdir_len)//fnamePPrev//'0100.'//ext(1:3)
            status(1,2)=0
            return
         endif
         var='sh'
         call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,rhLapsGS,istatus)
         if(istatus.ne.1)then
            write (6,*) ' Error in readlapsdata for LGA HT PPrev',
     1fdir_in(1:fdir_len)//fnamePPrev//'0100.'//ext(1:3)
            status(1,5)=0
            return
         endif

         do k=1,nk
          do j=1,ny
           do i=1,nx
              rhLapsGS(i,j,k)=make_rh(laps_levels(k)
     1,tLapsGS(i,j,k)-273.16,rhLapsGS(i,j,k)*1000.,-132.)*100.
           enddo
          enddo
         enddo

      else

      j = 1
      var = 'u3'
      ext = 'lw3'

      if (balance .eq. 1) then
        fdir_in = dir_in(1:len_dir_in)//'balance/lw3/'
      elseif(balance .eq. 0)then
        fdir_in = dir_in(1:len_dir_in)//'lw3/'
      endif

      call s_len(fdir_in,fdir_len)

      call read_laps_data(i4time_syn,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,uLapsGS,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' Error in readlapsdata for LW3 U Syn'
        status(1,j) = 0
      else
        write(6,*) 'Success reading U3 from ',
     1             fdir_in(1:fdir_len)//fnameSyn//'.'//ext(1:3)
      endif

      call read_laps_data(i4time_prev,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,uLapsGP,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' Error in readlapsdata for LW3 U Prev'
        status(2,j) = 0
      else
        write(6,*) 'Success reading U3 from ',
     1             fdir_in(1:fdir_len)//fnamePrev//'.'//ext(1:3)
      endif

      j = 2
      do i = 1, nk
        var(i) = 'v3'
      enddo

      call read_laps_data(i4time_syn,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,vLapsGS,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' Error in readlapsdata for LW3 V Syn'
        status(1,j) = 0
      else
        write(6,*) 'Success reading V3 from ',
     1             fdir_in(1:fdir_len)//fnameSyn//'.'//ext(1:3)
      endif

      call read_laps_data(i4time_prev,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,vLapsGP,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' Error in readlapsdata for LW3 V Prev'
        status(2,j) = 0
      else
        write(6,*) 'Success reading V3 from ',
     1             fdir_in(1:fdir_len)//fnamePrev//'.'//ext(1:3)
      endif

      ext = 'lt1'

      j = 3
      do i = 1, nk
        var(i) = 't3'
      enddo
      if (balance .eq. 1) then
        fdir_in = dir_in(1:len_dir_in)//'balance/lt1/'
      else
        fdir_in = dir_in(1:len_dir_in)//'lt1/'
      endif
      call s_len(fdir_in,fdir_len)

      call read_laps_data(i4time_syn,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,tLapsGS,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' Error in readlapsdata for LT1 T3 Syn '
        status(1,j) = 0
        status = 0
      else
        write(6,*) 'Success reading T3 from ',
     1             fdir_in(1:fdir_len)//fnameSyn//'.'//ext(1:3)
      endif

      call read_laps_data(i4time_prev,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,tLapsGP,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' Error in readlapsdata for LT1 T3 Prev '
        status(2,j) = 0
      else
        write(6,*) 'Success reading T3 from ',
     1             fdir_in(1:fdir_len)//fnamePrev//'.'//ext(1:3)
      endif

      j = 4
      do i = 1, nk
        var(i) = 'ht'
      enddo

      call read_laps_data(i4time_syn,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,htLapsGS,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' Error in readlapsdata for LT1 HT Syn'
        status(1,j) = 0
      else
        write(6,*) 'Success reading HT from ',
     1             fdir_in(1:fdir_len)//fnameSyn//'.'//ext(1:3)
      endif

      call read_laps_data(i4time_prev,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,htLapsGP,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' Error in readlapsdata for LT1 HT Prev'
        status(2,j) = 0
      else
        write(6,*) 'Success reading HT from ',
     1             fdir_in(1:fdir_len)//fnamePrev//'.'//ext(1:3)
      endif

      ext = 'lh3'

      j = 5
c     do i = 1, nk
c       var(i) = 'rh3'
c     enddo

      var = 'rhl'

      if (balance .eq. 1) then
        fdir_in = dir_in(1:len_dir_in)//'balance/lh3/'
      else
        fdir_in = dir_in(1:len_dir_in)//'lh3/'
      endif
      call s_len(fdir_in,fdir_len)

      call read_laps_data(i4time_syn,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,rhLapsGS,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' Error in readlapsdata for LH3 RH3 Syn'
        status(1,j) = 0
      else
        write(6,*) 'Success reading RH from ',
     1             fdir_in(1:fdir_len)//fnameSyn//'.'//ext(1:3)
      endif

      call read_laps_data(i4time_prev,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,rhLapsGP,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' Error in readlapsdata for LH3 RH3 Prev'
        status(2,j) = 0
      else
        write(6,*) 'Success reading RH from ',
     1             fdir_in(1:fdir_len)//fnamePrev//'.'//ext(1:3)
      endif

      endif

      return
      end
!1.....................................................................................
      subroutIne ascend_w(timeSyn, timeRel, numSigW, numW,fileAvail,
     1                    wmoStaNum, staLat, staLon, staElev, maxW,
     1                    max_ht_m_proc,typeW,
     1                    maxRaob, nx, ny, nk, raobId, statusL,
     1                    htSigW, wdSigW, wsSigW,  !remember (maxW,maxRaob)
     1                    uLapsGS,vLapsGS,uLapsGP,vLapsGP,  !(nx,ny,nk)
     1                    htLapsGS, htLapsGP, raob_missing_data,
     1                    verif_missing_data,
     1                    lat, lon, laps_levels,  !(nx,ny)
!................................variables below returned...................
     1                    riW,rjW,rkW,latW,lonW,htW, uIW,vIW,
     1                    uPW,vPW,timeLapsW,status)

      include 'trigd.inc'

      implicit none

      real          height_to_zcoord3
      integer       timeSyn, timeRel, numSigW, numW, fileAvail,
     1                wmoStaNum
      real          staLat, staLon, staElev,max_ht_m_proc
      character*1     typeW(maxW,maxRaob),typeWM(maxW,maxRaob)
      integer       maxW, maxRaob,nx,ny,nk, raobId,statusL(2,6)
      real          htSigW(maxW,maxRaob),
     1                wdSigW(maxW,maxRaob), wsSigW(maxW,maxRaob),
     1                uLapsGS(nx,ny,nk), vLapsGS(nx,ny,nk),
     1                uLapsGP(nx,ny,nk), vLapsGP(nx,ny,nk),
     1                htLapsGS(nx,ny,nk), htLapsGP(nx,ny,nk),
     1                raob_missing_data, verif_missing_data,
     1                lat(nx,ny), lon(nx,ny),
     1                laps_levels(nk), 
     1                riW(maxW,maxRaob),rjW(maxW,maxRaob),
     1                rkW(maxW,maxRaob),latW(maxW,maxRaob),
     1                lonW(maxW,maxRaob),htW(maxW,maxRaob),
     1		      uIW(maxW,maxRaob),vIW(maxW,maxRaob),
     1		      uPW(maxW,maxRaob),vPW(maxW,maxRaob)
      integer       timeLapsW(maxW,maxRaob), status

      integer       j,lvl,timePrev,index, time_tot, istatus,
     1                timePrevComp, timeSynComp,int_ri,int_rj
      real          prSigW, uSigW(maxW),vSigW(maxW),
     1                u_disp, v_disp, delta_h,delta_t, delta_u,
     1                delta_v, rise_rate, maxHtProc, rilaps, 
     1                rjlaps, rklaps, u_laps, v_laps, delta_s

! BEGIN
      status = 1   !assume a good return
      timePrev = timeSyn - 3600
      timeSynComp = timeSyn
      timePrevComp = timePrev
      if (fileAvail .eq. 1) timePrevComp = 0
      if (fileAvail .eq. 2) timeSynComp = 0

      if ((timePrevComp .eq. 0) .and. (timeSynComp .eq. 0)) then
        write(6,*) ' No wind data to process '
        istatus = 0
        return
      endif

      maxHtProc = max_ht_m_proc  !max htSigW value to process (in meters)

!     setup for ascention
      index = 0
      u_disp = 0.
      v_disp = 0.
      time_tot = timeRel
      rise_rate = 300. / 60. ! m/s

!     loop thru levels
      numW = -1
      lvl = 1
      do while (lvl .le. numSigW)
        if ((htSigW(lvl,raobId) .gt. 0) .and.
     1      (htSigW(lvl,raobId) .le. maxHtProc) .and.
     1      (htSigW(lvl,raobId) .ne. raob_missing_data)) then
          index = index + 1

          if (index .ge. 1) then
            if(index .gt. 1)then ! Lower level is above 1st (sfc) level
              delta_h = htSigW(lvl,raobId) - htSigW(lvl-1,raobId)
            else                 ! Lower level is at the sfc
              delta_h = htSigW(lvl,raobId) - staElev
            endif

            htW(index,raobId) = htSigW(lvl,raobId)
            typeWM(index,raobId) = typeW(index,raobId)

            if ((wdSigW(lvl,raobId) .eq. raob_missing_data) .or.
     1          (wsSigW(lvl,raobId) .eq. raob_missing_data)) then
              index = index - 1  ! don't use this lvl ob...go to next
              lvl = lvl + 1
              goto 777
            else
              call disp_to_uv(wdSigW(lvl,raobId),wsSigW(lvl,raobId),
     1                        uPW(index,raobId),vPW(index,raobId))
            endif

            delta_t = delta_h / rise_rate
            time_tot = time_tot + nint(delta_t + 0.5)
            if (abs(time_tot - timeSynComp) .gt. 
     1          abs(time_tot - timePrevComp)) then
              timeLapsW(index,raobId) = timePrev
            else
              timeLapsW(index,raobId) = timeSyn
            endif
            delta_s = wsSigW(lvl,raobId) * delta_t
            delta_u = -sind(wdSigW(lvl,raobId)) * delta_s
            delta_v = -cosd(wdSigW(lvl,raobId)) * delta_s
            u_disp = u_disp + delta_u
            v_disp = v_disp + delta_v
            latW(index,raobId) = staLat + v_disp / 111000.
            lonW(index,raobId) = staLon + u_disp / 
     1          (111000.* cosd((staLat+latW(index,raobId))/2))

!           Determine LAPS i,j,k, calc deltaU and deltaV
!           index has number of elements to correlate

            call latlon_to_rlapsgrid(latW(index,raobId),
     1                               lonW(index,raobId),
     1                               lat,lon,nx,ny,
     1                               riW(index,raobId),
     1                               rjW(index,raobId),istatus)

            if (istatus .ne. 1) then  !if this level is out, rest probably will be too
              numW = index -1
              lvl = numSigW + 1
              if (index .gt. 1) then  !keep what have
                goto 777
              else
                write(6,*) 'Raob ',wmoStaNum,' not in Laps domain  ',
     1                   latW(index,raobId),' ',lonW(index,raobId)
                goto 777
              endif
            endif
             
C           get integral ri, rj for height_to_zcoord3
            int_ri = nint(riW(index,raobId))
            int_rj = nint(rjW(index,raobId))

            if (timeLapsW(index,raobId) .eq. timeSyn) then   !use htLapsGS
              rkW(index,raobId) = height_to_zcoord3(htSigW(lvl,raobId),
     1                                       htLapsGS,laps_levels,
     1                                       nx,ny,nk,int_ri,
     1                                       int_rj,istatus)
            else  !use htLapsGP
              rkW(index,raobId) = height_to_zcoord3(htSigW(lvl,raobId),
     1                                       htLapsGP,laps_levels,
     1                                       nx,ny,nk,int_ri,
     1                                       int_rj,istatus)
            endif

            if ((istatus .eq. 1) .and. (rkW(index,raobId) .le. nk)) then

!             determine which LAPS grid to use
              if (timeLapsW(index,raobId) .eq. timeSyn) then
                call trilinear_laps(riW(index,raobId),rjW(index,raobId),
     1                              rkW(index,raobId),nx,ny,nk,
     1                              uLapsGS, uIW(index,raobId))

                call trilinear_laps(riW(index,raobId),rjW(index,raobId),
     1                              rkW(index,raobId),nx,ny,nk,
     1                              vLapsGS, vIW(index,raobId))
              else
                call trilinear_laps(riW(index,raobId),rjW(index,raobId),
     1                              rkW(index,raobId),nx,ny,nk,
     1                              uLapsGP, uIW(index,raobId))

                call trilinear_laps(riW(index,raobId),rjW(index,raobId),
     1                              rkW(index,raobId),nx,ny,nk,
     1                              vLapsGP, vIW(index,raobId))
              endif

            else  !height is above Laps domain
              numW = index - 1
              lvl = numSigW + 1
              if (index .gt. 1) then  !save what have
                goto 777
              else
                write(6,*) 'Raob ',wmoStaNum, ' above Laps domain  ',
     1                     rkW(index,raobId)
                goto 777
              endif
            endif

          endif
        endif  !if data for lvl not missing

        lvl = lvl + 1
777     continue

      enddo

      if (numW .eq. -1) numW = index

C     pass back modified typeW
      do j = 1, numW
        typeW(j,raobId) = typeWM(j,raobId)
      enddo
 
      return
      end
!2.....................................................................................
      subroutine ascend_t(timeSyn, timeRel, numSigT,
     1                    numT,fileAvailUV,fileAvail,wmoStaNum,
     1                    staLat, staLon, staElev, typeT,
     1                    maxT, maxW, max_ht_m_proc,
     1                    maxRaob, nx, ny, nk, raobId, statusL,
     1                    prSigT, tSigT, tdSigT,htSigT,
     1                    tLapsGS,rhLapsGS,htLapsGS,tLapsGP,  !(nx,ny,nk)
     1                    rhLapsGP, htLapsGP,
     1                    htLgaGS, htLgaGP, raob_missing_data,
     1                    verif_missing_data, lat,lon, !(nx,ny,nk)
     1                    laps_levels, numSigW,
     1                    htSigW, wdSigW, wsSigW, !remember (maxW,maxRaob)
!................................variables below returned...................
     1                    riT,rjT,rkT,latT,lonT,htT,prIT,tIT,tdIT,
     1                    prPT,tPT,tdPT,timeLapsT,status)

      include 'trigd.inc'

      implicit none

      real          ztopsa, zcoord_of_logpressure, psatoz,
     1                make_ssh, k_to_c, make_td, c_to_k

      integer       timeSyn, timeRel, numSigT, numT,
     1                fileAvailUV, fileAvail,wmoStaNum
      real          staLat, staLon, staElev,max_ht_m_proc
      character*1     typeT(maxT,maxRaob),typeTM(maxT,maxRaob)
      integer       maxT, maxW, maxRaob,nx,ny,nk
      integer       raobId, statusL(2,6)
      real          prSigT(maxT,maxRaob),
     1                tSigT(maxT,maxRaob), 
     1                tdSigT(maxT,maxRaob), 
     1                htSigT(maxT,maxRaob), 
     1                tLapsGS(nx,ny,nk), htLapsGS(nx,ny,nk),
     1                tLapsGP(nx,ny,nk), htLapsGP(nx,ny,nk),
     1                rhLapsGS(nx,ny,nk), rhLapsGP(nx,ny,nk),
     1                htLgaGS(nx,ny,nk), htLgaGP(nx,ny,nk),
     1                raob_missing_data,verif_missing_data,
     1                lat(nx,ny), lon(nx,ny), laps_levels(nk)
      integer	      numSigW
      real          htSigW(maxW,maxRaob), wdSigW(maxW,maxRaob),
     1                wsSigW(maxW,maxRaob),
     1                riT(maxT,maxRaob),rjT(maxT,maxRaob),
     1                rkT(maxT,maxRaob),latT(maxT,maxRaob),
     1                lonT(maxT,maxRaob),htT(maxT,maxRaob),
     1		      prIT(maxT,maxRaob),tIT(maxT,maxRaob),
     1		      tdIT(maxT,maxRaob),prPT(maxT,maxRaob),
     1                tPT(maxT,maxRaob),tdPT(maxT,maxRaob)
      integer       timeLapsT(maxT,maxRaob), status

      integer       timePrevComp, timeSynComp,istatus,
     1                j,lvl,timePrev,index, time_tot,
     1                bkgPrev, bkgSyn,int_ri,int_rj
      real          wd, ws, rh, q, t_ref, tC,
     1                u_disp, v_disp,delta_h,delta_t,delta_u,
     1                delta_v, rise_rate, maxHtProc,htM,htZ, 
     1                t_laps, ht, htPrev,delta_s,pres_pa

! BEGIN
      status = 1   !assume a good return

      timePrev = timeSyn - 3600
      timeSynComp = timeSyn
      timePrevComp = timePrev
      if ((fileAvail .eq. 1) .or. (fileAvailUV .eq. 1)) timePrevComp = 0
      if ((fileAvail .eq. 2) .or. (fileAvailUV .eq. 2)) timeSynComp = 0
      if (statusL(1,6) .eq. 1) then
        bkgSyn = 1 
      else
        bkgSyn = 0 
      endif
      if (statusL(2,6) .eq. 1) then
        bkgPrev = 1 
      else
        bkgPrev = 0 
      endif

      if ((bkgSyn .eq. 0) .and. (bkgPrev .eq. 0)) then
        write(6,*) 'No LGA background available...',
     1             'using psatoz for pressure/height calcs for temp'
      endif
      if ((timePrevComp .eq. 0) .and. (timeSynComp .eq. 0)) then
        write(6,*) ' No temp/wind data to process '
        status = 0
        return
      endif

      maxHtProc = max_ht_m_proc  !max htSigW value to process (in meters)

!     setup for ascention
      index = 0
      u_disp = 0.
      v_disp = 0.
      time_tot = timeRel
      rise_rate = 300. / 60. ! m/s
      htPrev = staElev

!     loop thru levels
      numT = -1
      lvl = 1

      do while (lvl .le. numSigT)

c       write(6,*) 'prSigT(lvl,raobId)=',prSigT(lvl,raobId), lvl, raobId

        if (htSigT(lvl,raobId) .ne. verif_missing_data) then 
          write(6,*) 'Height from htMan used'
          ht = htSigT(lvl,raobId)
        else
          if (prSigT(lvl,raobId) .ne. raob_missing_data) then
C           quick and dirty for now...need ri,rj,timeLaps for pressure_to_height
            ht = psatoz(prSigT(lvl,raobId))
          else
            ht = raob_missing_data
          endif
        endif

c       write(6,*) 'height from ztopsa=',ht, prSigT(lvl,raobId)

C       stop below htSigW(max) because of ascend
        if ((ht .gt. 0) .and. (ht .le. maxHtProc) .and.
     1      (ht .ne. raob_missing_data) .and.
     1      (ht .le. htSigW(numSigW,raobId))) then

          call interpWind2T(maxW, maxRaob, numSigW, wdSigW, 
     1                      wsSigW, htSigW, raobId, maxHtProc, 
     1                      ht, ws, wd, istatus)
          if (istatus .eq.1) then
            index = index + 1

            if (index .ge. 1) then

C             fill in prPT, tPT, tdPT from sigT vars
              prPT(index,raobId) = prSigT(lvl,raobId)
              tPT(index,raobId) = tSigT(lvl,raobId)
              tdPT(index,raobId) = tdSigT(lvl,raobId)
              typeTM(index,raobId) = typeT(lvl,raobId)

              delta_h = ht - htPrev

              delta_t = delta_h / rise_rate
              time_tot = time_tot + nint(delta_t + 0.5)
              if (abs(time_tot - timeSynComp) .gt. 
     1            abs(time_tot - timePrevComp)) then
                timeLapsT(index,raobId) = timePrev
              else
                timeLapsT(index,raobId) = timeSyn
              endif

              delta_s = ws * delta_t
              delta_u = -sind(wd) * delta_s
              delta_v = -cosd(wd) * delta_s
              u_disp = u_disp + delta_u
              v_disp = v_disp + delta_v
              latT(index,raobId) = staLat + v_disp / 111000.
              lonT(index,raobId) = staLon + u_disp / 
     1                 (111000.* cosd((staLat+latT(index,raobId))/2))
              htPrev = ht

!             Determine LAPS i,j,k

              call latlon_to_rlapsgrid(latT(index,raobId),
     1                                 lonT(index,raobId),
     1                                 lat,lon,nx,ny,
     1                                 riT(index,raobId),
     1                                 rjT(index,raobId),
     1                                 istatus)

              if (istatus .ne. 1) then  !out of Laps domain...rest probably too
                 numT = index - 1
                 lvl = numSigT + 1
                if (index .gt. 1) then  !save what have
                  goto 777
                else
                  write(6,*) 'Raob ',wmoStaNum,' not in Laps domain  ',
     1                      latT(index,raobId),' ',lonT(index,raobId)
                  goto 777
                endif
              endif

C             re-calc height more accurately now that have ri, rj
              if (htSigT(lvl,raobId) .eq. verif_missing_data) then  !not Man height...recalc
              else
                htM = ht
              endif

                pres_pa = prSigT(lvl,raobId) * 100
                int_ri = nint(riT(index,raobId))
                int_rj = nint(rjT(index,raobId))

                if ((bkgSyn.eq.1).and.(bkgPrev.eq.1)) then  !use bkg closest to timeLaps(index,raobId)
                  if (timeLapsT(index,raobId) .eq. timeSyn) then
                    call pressure_to_height(pres_pa,htLgaGS,nx,ny,
     1                                      nk,int_ri,int_rj,ht,istatus)
                  else
                    call pressure_to_height(pres_pa,htLgaGP,nx,ny,
     1                                      nk,int_ri,int_rj,ht,istatus)
                  endif
                else
                  if (bkgSyn .eq. 1) then  ! use htLgaGS
                    call pressure_to_height(pres_pa,htLgaGS,nx,ny,
     1                                      nk,int_ri,int_rj,ht,istatus)
                  else
                    if (bkgPrev .eq. 1) then  ! use htLgaGP
                      call pressure_to_height(pres_pa,htLgaGP,nx,ny,
     1                                      nk,int_ri,int_rj,ht,istatus)
                    else
                      ! leave calc from psatoz
                    endif
                  endif
                endif
c             endif

              if (htSigT(lvl,raobId) .eq. verif_missing_data) then  !not Man height...recalc
              htT(index,raobId) = ht
              else
              htT(index,raobId) = htM
              htZ = psatoz(prSigT(lvl,raobId))
              write(6,*) 'ht=',ht,'  htM=',htM,'  htZ=',htZ

              endif
              rkT(index,raobId) = 
     1          zcoord_of_logpressure(prSigT(lvl,raobId)*100.)

              if (rkT(index,raobId) .le. nk) then

!               determine which LAPS grid to use and interpolate T in the vertical 
                if (timeLapsT(index,raobId) .eq. timeSyn) then

                  call trilinear_laps(riT(index,raobId),
     1                                rjT(index,raobId),
     1                                rkT(index,raobId),nx,ny,nk,
     1                                tLapsGS, tIT(index,raobId))

                  call trilinear_laps(riT(index,raobId),
     1                                rjT(index,raobId),
     1                                rkT(index,raobId),nx,ny,nk,
     1                                htLapsGS, ht)

                  call trilinear_laps(riT(index,raobId),
     1                                rjT(index,raobId),
     1                                rkT(index,raobId),nx,ny,nk,
     1                                rhLapsGS, rh)

                  rh = rh/100.0
                  prIT(index,raobId) = ztopsa(ht)
                  t_ref = -132.0
                  tC = k_to_c(tIT(index,raobId))
                  q = make_ssh(prIT(index,raobId),tC,rh,t_ref) 
                  tdIT(index,raobId) = 
     1              make_td(prIT(index,raobId),tC,q,t_ref)
                  tdIT(index,raobId) = c_to_k(tdIT(index,raobId))

                else

                  call trilinear_laps(riT(index,raobId),
     1                                rjT(index,raobId),
     1                                rkT(index,raobId),nx,ny,nk,
     1                                tLapsGP, tIT(index,raobId))

                  call trilinear_laps(riT(index,raobId),
     1                                rjT(index,raobId),
     1                                rkT(index,raobId),nx,ny,nk,
     1                                htLapsGS, ht)

                  call trilinear_laps(riT(index,raobId),
     1                                rjT(index,raobId),
     1                                rkT(index,raobId),nx,ny,nk,
     1                                rhLapsGP, rh)

                  prIT(index,raobId) = ztopsa(ht)
                  t_ref = -132.0
                  tC = k_to_c(tIT(index,raobId))
                  q = make_ssh(prIT(index,raobId),tC,rh,t_ref) 
                  tdIT(index,raobId) = 
     1               make_td(prIT(index,raobId),tC,q,t_ref)

                endif

!               check tIT(index,raobId) for validity
                if (tIT(index,raobId) .ge. 350) then
                  index = index - 1
                endif
              else
                numT = index - 1
                lvl = numSigT + 1
                if (index .gt. 1) then  !save what have
                  goto 777
                else
                 write(6,*) 'Raob ',wmoStaNum,' is above Laps domain  ',
     1                       rkT(index,raobId)
                  goto 777
                endif
              endif

            endif
          endif
        endif
 
        lvl = lvl + 1
777     continue

      enddo

      if (numT .eq. -1) numT = index

C     pass back modified typeT
      do j = 1, numT
        typeT(j,raobId) = typeTM(j,raobId)
      enddo
 
      return
      end
!3.....................................................................................
      subroutine interpWind2T(maxW, maxRaob, numSigW, wdSigW, 
     1                        wsSigW, htSigW, raobId, maxHtProc, 
     1                        ht, ws, wd, istatus)

      implicit none
      integer    maxW, maxRaob, numSigW, raobId, maxHtProc, 
     1             istatus, u, l, i, true, false, found, last
      real       htSigW(maxW,maxRaob), wdSigW(maxW,maxRaob),
     1             wsSigW(maxW,maxRaob), ht, ws, wd, deltaWS

! begin
      istatus = 1
      true = 1
      false = 0

!     We know that ht < maxHtProc and ht <= htSigW(numSigW,raobId)
        l = 1
      if (htSigW(1,raobId) .eq. 0.0) then
        if (numSigW .ge. 3) then
          l = 2
          u = 3
        else
          istatus = 0
          return
        endif
      elseif (numSigW .ge. 2) then
        l = 1
        u = 2
      endif

      if (ht .lt. htSigW(l,raobId)) then
        istatus = 0
        return
      endif

!     We know that ht >= htSigW(l,raobId) and ht <= htSigW(numSigW,raobId)

      found = false
      last = false
      do while ((u .le. numSigW) .and. (found .eq. false) .and.
     1          (last .eq. false)) 
        if ((ht .ge. htSigW(l,raobId)) .and.
     1      (ht .le. htSigW(u,raobId))) then
          found = true
        elseif (u .eq. numSigW) then
          last = true
        elseif (ht .gt. htSigW(u,raobId)) then
          l = l + 1
          u = u + 1
        endif
      enddo

      if (found .eq. true) then
        ws = wsSigW(l,raobId) + ((ht-htSigW(l,raobId)) *
     1       ((wsSigW(u,raobId)-wsSigW(l,raobId)) /
     1        (htSigW(u,raobId)-htSigW(l,raobId))))
        wd = wdSigW(l,raobId) + ((ht-htSigW(l,raobId)) *
     1       ((wdSigW(u,raobId)-wdSigW(l,raobId)) /
     1        (htSigW(u,raobId)-htSigW(l,raobId))))
        istatus = 1
      else
        istatus = 0
      endif

      return
      end
!4....................................................................................
      subroutine mergeTW(MAX_HTS,MAX_RAOBS,maxW,maxT,numW,numT,
     1                   nRaobs,use_raob,n_raobs_use,
     1                   ri,rj,rk,lat,lon,hts,type,uP,uI,vP,vI,
     1                   tP,tI,tdP,tdI,prP,prI,timeLaps,nHts,
     1                   riT,rjT,rkT,latT,lonT,htT,typeT,
     1                   riW,rjW,rkW,latW,lonW,htW,typeW,
     1                   uIW,vIW,uPW,vPW,timeLapsW,
     1                   prIT,tIT,tdIT,prPT,tPT,tdPT,timeLapsT,
     1                   verif_missing_data, raob_missing_data,
     1                   istatus)

      implicit none

      integer		MAX_HTS,MAX_RAOBS,maxW,maxT,
     1                  nRaobs,use_raob(MAX_RAOBS),
     1                  n_raobs_use

C     data written out
      real            ri(MAX_HTS,MAX_RAOBS),
     1                  rj(MAX_HTS,MAX_RAOBS),
     1                  rk(MAX_HTS,MAX_RAOBS),
     1                  lat(MAX_HTS,MAX_RAOBS),
     1                  lon(MAX_HTS,MAX_RAOBS),
     1                  hts(MAX_HTS,MAX_RAOBS)
      character*1       type(MAX_HTS,MAX_RAOBS)
      real            uP(MAX_HTS,MAX_RAOBS),
     1                  uI(MAX_HTS,MAX_RAOBS),
     1                  vP(MAX_HTS,MAX_RAOBS),
     1                  vI(MAX_HTS,MAX_RAOBS),
     1                  tP(MAX_HTS,MAX_RAOBS),
     1                  tI(MAX_HTS,MAX_RAOBS),
     1                  tdP(MAX_HTS,MAX_RAOBS),
     1                  tdI(MAX_HTS,MAX_RAOBS),
     1                  prP(MAX_HTS,MAX_RAOBS),
     1                  prI(MAX_HTS,MAX_RAOBS)
      integer         timeLaps(MAX_HTS,MAX_RAOBS),
     1                  nHts(MAX_RAOBS)

C     temp variables
      integer         numT(MAX_RAOBS),
     1			timeLapsT(maxT,MAX_RAOBS)
      character*1       typeT(maxT,MAX_RAOBS)
      real            riT(maxT,MAX_RAOBS),
     1                  rjT(maxT,MAX_RAOBS),
     1                  rkT(maxT,MAX_RAOBS),
     1                  latT(maxT,MAX_RAOBS),
     1                  lonT(maxT,MAX_RAOBS),
     1                  htT(maxT,MAX_RAOBS),
     1                  prIT(maxT,MAX_RAOBS),
     1                  tIT(maxT,MAX_RAOBS),
     1                  tdIT(maxT,MAX_RAOBS),
     1                  prPT(maxT,MAX_RAOBS),
     1                  tPT(maxT,MAX_RAOBS),
     1                  tdPT(maxT,MAX_RAOBS)

C     wind variables
      integer         numW(MAX_RAOBS),
     1 			timeLapsW(maxW,MAX_RAOBS)
      character*1       typeW(maxW,MAX_RAOBS)
      real            riW(maxW,MAX_RAOBS),
     1                  rjW(maxW,MAX_RAOBS),
     1                  rkW(maxW,MAX_RAOBS),
     1                  latW(maxW,MAX_RAOBS),
     1                  lonW(maxW,MAX_RAOBS),
     1                  htW(maxW,MAX_RAOBS),
     1                  uIW(maxW,MAX_RAOBS),
     1                  vIW(maxW,MAX_RAOBS),
     1                  uPW(maxW,MAX_RAOBS),
     1                  vPW(maxW,MAX_RAOBS)

      real		verif_missing_data,
     1                  raob_missing_data
      integer		istatus

C     local variables

      integer		tPtr,wPtr,index,i,j,nDiff
      real		sumDiff,avgDiff,pctDiff,
     1                  sumPctDiff,avgPctDiff,
     1                  rkDiff

C
C     BEGIN
C
      istatus = 1   !assume good return
      nDiff = 0
      sumDiff = 0

C     loop through raobs
      do i = 1, nRaobs
        if (use_raob(i) .eq. 1) then

          write(6,*) 'Raob ',i,' before merge'
          write(6,*) numW(i)
          write(6,*) 'htW, typeW, uPW, vPW'
          do j = 1, numW(i)
            write(6,*) j,htW(j,i),typeW(j,i),uPW(j,i),
     1                 vPW(j,i)
          enddo

          write(6,*) numT(i)
          write(6,*) 'htT, typeT, prPT, tPT, tdPT'
          do j = 1, numT(i)
            write(6,*) j,htT(j,i),typeT(j,i),prPT(j,i),
     1                 tPT(j,i),tdPT(j,i)
          enddo

C         setup pointers
          tPtr = 1
          wPtr = 1
          index = 1

C         both htT and htW should have no missing data...if it does, don't use it
C         drive by sigT because it has more levels

          do while (tPtr .le. numT(i))
            if (wPtr .le. numW(i)) then

C             verify ri,rj,rk,lat,lon,ht are valid

              if (htW(wPtr,i) .lt. htT(tPtr,i)) then  ! W ob
                hts(index,i) = htW(wPtr,i)
                ri(index,i) = riW(wPtr,i)
                rj(index,i) = rjW(wPtr,i)
                rk(index,i) = rkW(wPtr,i)
                lat(index,i) = latW(wPtr,i)
                lon(index,i) = lonW(wPtr,i)
                uP(index,i) = uPW(wPtr,i)
                uI(index,i) = uIW(wPtr,i)
                vP(index,i) = vPW(wPtr,i)
                vI(index,i) = vIW(wPtr,i)
                tP(index,i) = verif_missing_data
                tI(index,i) = verif_missing_data
                tdP(index,i) = verif_missing_data
                tdI(index,i) = verif_missing_data
                prP(index,i) = verif_missing_data
                prI(index,i) = verif_missing_data
                type(index,i) = 'W'
                if (typeW(wPtr,i) .eq. 'M') then 
                  write(6,*)'Raob ',i,' has mismatch Man htW.lt.htT',
     1                      wPtr,' ',tPtr,' ',htW(wPtr,i),' ',
     1                      htT(tPtr,i)
                endif
                wPtr = wPtr + 1
                index = index + 1
              else
                if (htT(tPtr,i) .lt. htW(wPtr,i)) then
                  hts(index,i) = htT(tPtr,i)
                  ri(index,i) = riT(tPtr,i)
                  rj(index,i) = rjT(tPtr,i)
                  rk(index,i) = rkT(tPtr,i)
                  lat(index,i) = latT(tPtr,i)
                  lon(index,i) = lonT(tPtr,i)
                  uP(index,i) = verif_missing_data
                  uI(index,i) = verif_missing_data
                  vP(index,i) = verif_missing_data
                  vI(index,i) = verif_missing_data
                  tP(index,i) = tPT(tPtr,i)
                  tI(index,i) = tIT(tPtr,i)
                  tdP(index,i) = tdPT(tPtr,i)
                  tdI(index,i) = tdIT(tPtr,i)
                  prP(index,i) = prPT(tPtr,i)
                  prI(index,i) = prIT(tPtr,i)
                  type(index,i) = 'T'
                  if (typeT(tPtr,i) .eq. 'M') then 
                    write(6,*)'Raob ',i,' has mismatch Man htT.lt.htW',
     1                        wPtr,' ',tPtr,' ',htW(wPtr,i),' ',
     1                        htT(tPtr,i)
                  endif
                  tPtr = tPtr + 1
                  index = index + 1
                else
                  if (htW(wPtr,i) .eq. htT(tPtr,i)) then
                    if ((typeW(wPtr,i) .eq. 'M') .and.  !find corresponding entry in T data
     1                  (riW(wPtr,i) .eq. riT(tPtr,i)) .and.
     1                  (rjW(wPtr,i) .eq. rjT(tPtr,i)) .and.
     1                  (latW(wPtr,i) .eq. latT(tPtr,i)) .and.
     1                  (lonW(wPtr,i) .eq. lonT(tPtr,i))) then
c                     if (abs(rkDiff) .lt. 0.02) then
                        hts(index,i) = htW(wPtr,i)
                        ri(index,i) = riW(wPtr,i)
                        rj(index,i) = rjW(wPtr,i)
                        rk(index,i) = rkW(wPtr,i)
                        lat(index,i) = latW(wPtr,i)
                        lon(index,i) = lonW(wPtr,i)
                        uP(index,i) = uPW(wPtr,i)
                        uI(index,i) = uIW(wPtr,i)
                        vP(index,i) = vPW(wPtr,i)
                        vI(index,i) = vIW(wPtr,i)
                        tP(index,i) = tPT(tPtr,i)
                        tI(index,i) = tIT(tPtr,i)
                        tdP(index,i) = tdPT(tPtr,i)
                        tdI(index,i) = tdIT(tPtr,i)
                        prP(index,i) = prPT(tPtr,i)
                        prI(index,i) = prIT(tPtr,i)
                        type(index,i) = 'M'
c                     else   !mismatch on nav...write out wind only
c                       hts(index,i) = htW(wPtr,i)
c                       ri(index,i) = riW(wPtr,i)
c                       rj(index,i) = rjW(wPtr,i)
c                       rk(index,i) = rkW(wPtr,i)
c                       lat(index,i) = latW(wPtr,i)
c                       lon(index,i) = lonW(wPtr,i)
c                       uP(index,i) = uPW(wPtr,i)
c                       uI(index,i) = uIW(wPtr,i)
c                       vP(index,i) = vPW(wPtr,i)
c                       vI(index,i) = vIW(wPtr,i)
c                       tP(index,i) = verif_missing_data
c                       tI(index,i) = verif_missing_data
c                       tdP(index,i) = verif_missing_data
c                       tdI(index,i) = verif_missing_data
c                       prP(index,i) = verif_missing_data
c                       prI(index,i) = verif_missing_data
c                       type(index,i) = 'M'

c                     endif
                    else
                      hts(index,i) = htW(wPtr,i)
                      ri(index,i) = riW(wPtr,i)
                      rj(index,i) = rjW(wPtr,i)
                      rk(index,i) = rkW(wPtr,i)
                      lat(index,i) = latW(wPtr,i)
                      lon(index,i) = lonW(wPtr,i)
                      uP(index,i) = uPW(wPtr,i)
                      uI(index,i) = uIW(wPtr,i)
                      vP(index,i) = vPW(wPtr,i)
                      vI(index,i) = vIW(wPtr,i)
                      tP(index,i) = verif_missing_data
                      tI(index,i) = verif_missing_data
                      tdP(index,i) = verif_missing_data
                      tdI(index,i) = verif_missing_data
                      prP(index,i) = verif_missing_data
                      prI(index,i) = verif_missing_data
                      type(index,i) = 'M'

                      rkDiff = rkW(wPtr,i) - rkT(tPtr,i)
                      write(6,*) 'rkDiff = ',rkDiff
  
                      pctDiff = ((rkW(wPtr,i)-rkT(tPtr,i))/
     1                         rkW(wPtr,i))*100
                      sumPctDiff = sumPctDiff + pctDiff
                      sumDiff = sumDiff + (rkW(wPtr,i)-rkT(tPtr,i))
                      nDiff = nDiff + 1
                      write(6,*)'Raob ',i, ' has Man data'
     1                         , ' rkW rkT pres %diff(w-t)/w ht'
                      write(6,*)rkW(wPtr,i),rkT(tPtr,i),prPT(tPtr,i),
     1                          pctDiff,hts(index,i)
                      write(6,*) 'latW/T lonW/T riW/T rjW/T'
                      write(6,*) latW(wPtr,i), latT(tPtr,i),
     1                           lonW(wPtr,i), lonT(tPtr,i),
     1                           riW(wPtr,i), riT(tPtr,i),
     1                           rjW(wPtr,i), rjT(tPtr,i)
                      write(6,*)
                    endif

                    wPtr = wPtr + 1
                    tPtr = tPtr + 1
                    index = index + 1

                  endif   !htW(wPtr,i) .eq. htT(tPtr,i)
                endif   !htT(tPtr,i) .lt. htW(wPtr,i)
              endif   !htW(wPtr,i) .lt. htT(tPtr,i) 

            else  !wPtr .gt numW(i)

              do while (tPtr .le. numT(i))
                hts(index,i) = htT(tPtr,i)
                ri(index,i) = riT(tPtr,i)
                rj(index,i) = rjT(tPtr,i)
                rk(index,i) = rkT(tPtr,i)
                lat(index,i) = latT(tPtr,i)
                lon(index,i) = lonT(tPtr,i)
                uP(index,i) = verif_missing_data
                uI(index,i) = verif_missing_data
                vP(index,i) = verif_missing_data
                vI(index,i) = verif_missing_data
                tP(index,i) = tPT(tPtr,i)
                tI(index,i) = tIT(tPtr,i)
                tdP(index,i) = tdPT(tPtr,i)
                tdI(index,i) = tdIT(tPtr,i)
                prP(index,i) = prPT(tPtr,i)
                prI(index,i) = prIT(tPtr,i)
                type(index,i) = 'T'

                tPtr = tPtr + 1
                index = index + 1
              enddo
            endif  !wPtr .le. numW(i)

          enddo    !tPtr .le. numT(i)

C         save number of heights for raob(i)
          nHts(i) = index - 1

C         do final sweep on data to make sure all missing values equal verif_missing_data
          do j = 1, nHts(i)
            if (uP(j,i).eq.raob_missing_data)  
     1          uP(j,i) = verif_missing_data
            if (uI(j,i).eq.raob_missing_data)  
     1          uI(j,i) = verif_missing_data
            if (vP(j,i).eq.raob_missing_data)  
     1          vP(j,i) = verif_missing_data
            if (vI(j,i).eq.raob_missing_data)  
     1          vI(j,i) = verif_missing_data
            if (tP(j,i).eq.raob_missing_data)  
     1          tP(j,i) = verif_missing_data
            if (tI(j,i).eq.raob_missing_data)  
     1          tI(j,i) = verif_missing_data
            if (tdP(j,i).eq.raob_missing_data)  
     1          tdP(j,i) = verif_missing_data
            if (tdI(j,i).eq.raob_missing_data)  
     1          tdI(j,i) = verif_missing_data
            if (prP(j,i).eq.raob_missing_data)  
     1          prP(j,i) = verif_missing_data
            if (prI(j,i).eq.raob_missing_data)  
     1          prI(j,i) = verif_missing_data
          enddo 

          write(6,*) 'Raob ',i,' after merge'
          write(6,*) nHts(i)
          write(6,*) 'hts, type, prP, tP, uP, vP'
          do j = 1, nHts(i)
            write(6,*) j,hts(j,i),type(j,i), prP(j,i),
     1                 tP(j,i), uP(j,i), vP(j,i)
          enddo

        endif    !use_raob(i) .eq. 1
      enddo    !loop through raobs

      if (nDiff .gt. 0) then
        avgPctDiff = sumPctDiff/nDiff
        avgDiff = sumDiff/nDiff
        write(6,*) 'Avg diff on rk = ',avgDiff,'  nDiff = ',nDiff
        write(6,*) 'Avg pct diff = ',avgPctDiff
        write(6,*)
      endif

      return
      end
!5....................................................................................
