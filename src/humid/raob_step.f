cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine raob_step (i4time, data, laps_pressure, 
     1     raob_lookback,
     1     lat,lon,  laps_t, ii,jj,kk)




c       this module modifies the input data array of specific humidity
c       with available raob data
c       this module will use .snd data if it is available.
c       in addition there is a way to turn off this module in the presence
c       of .snd data.   simply modify the static file
c       ../static/raob_moisture_switch.txt
c       by putting a 0 in the first record.  the advantage of this approach
c       is that it is external to the program and will not require a
c       recompile.
c
c
c
c       theory:
c
c       this code utilizes the second pass of a barnes analysis.
c       in this form, a background is already available i.e., from the
c       model forecast that the moisture analysis obtains from .lga
c       or in 4dda mode from a model run
c
c       the second pass barnes scales the background to the raob data
c       the scaling approach prevents illogical adjustments in large
c       domains, common for AFWA application.
c
c       this is a "one-pass" analysis.
c
c       the radius of influence of a raob is 600km.  At that distance,
c       the weight associated with a raob is 0.5.  it falls off rapidly
c       as distance increases
c
c       a note on the look_back_time parameter:
c
c       look_back_time was devised to allow the generation of snd files
c       only at the time of observation (not in the event of no obs)
c       hence this is really not required in the current way the "system"
c       as a whole is being coded as a group.  in the "group" approach,
c       a .snd file will be created each laps time and will automatically
c       contain the raob data for a predescribed "look back time".  
c       typically this will be 3 hours.  if there is no data, the .snd file 
c       will exist and will be empty.
c
c       if for some reason this group approach for lookback time is ever
c       dropped or modified, this code contains a switch to accommodate that
c       scenario if needed.
c
c       on a related note, if the .snd files exist with a more recent date
c       than the current runtime, the presence of these will be detected
c       and the current runtime .snd file will be sought instead.  if this
c       current file is not avialable, the code will report it is
c       unable to find .snd data for the current time.


      implicit none

c     include 'lapsparms.for'
c     include 'parmtrs.inc'

c     parameters for sounding max and mins

      integer snd_lev, snd_tot
      parameter (snd_lev = 300)
      parameter (snd_tot = 1000)



c     input parameters

      integer i4time, ii,jj,kk, raob_lookback
      real data(ii,jj,kk), laps_pressure (kk)
      real lat(ii, jj), lon(ii, jj)
      real laps_t (ii,jj,kk)  ! laps temp field for QC checking

c  dynamic dependent parameters                     

      real  q_r(kk, snd_tot)
      real diff(kk, snd_tot)
      real scale(kk, snd_tot)
      real weight (ii,jj, snd_tot)
      integer mask(ii,jj)
      real climo_w(kk)
      integer x_sum_save(kk)

c  normal internal parameters

      integer i,j,k,ks,is
      integer look_back_time
      integer raob_i4time       ! actual time of raob data used
      integer raob_i4time_chk   ! used to validate actual raob
      integer numoffiles, istatus, maxfiles
      character*256 c_filenames(200),pathname_in
      character*9 filename
      character*200 fname
      integer idx, len
      real lat_r( snd_tot), lon_r( snd_tot)
      integer  i_r( snd_tot), j_r( snd_tot)
      integer lev_r( snd_tot), dummy, isound
      real rdummy
      character*5 cdummy
      integer idummy
      character*9 r_filename(ii*jj)
      character*8 snd_type(ii*jj)
      real p_r(snd_lev, snd_tot), td_r(snd_lev, snd_tot), 
     1     t_r(snd_lev, snd_tot)
      real temt, temtd, ssh2
      real r, rmin
      real d2r,pi
      real tot_weight,tot_scale_weight,scale_used,weight_avg,counter
      real wt_ck(ii,jj) ! weight check to test for ample RAOB data
      real rspacing_dum
      real rmd
      real x(snd_tot)
      integer x_sum
      real ave(kk),adev(kk),sdev(kk),var(kk),skew(kk),curt(kk)

c     climate model (for QC)

      real center_lat
      integer jday
      real 
     1     standard_press(40),
     1     tempertur_guess(40),
     1     mixratio_guess(40)


        
c *** begin routine

      write(6,*) '1.21; allowing satsnd data in routine, = weight'

      pi = acos(-1.0)
      d2r = pi/180.
      look_back_time = raob_lookback
      maxfiles = 200
      do j = 1,jj
         do i = 1,ii
            wt_ck(i,j) = 0.0
         enddo
      enddo

      call get_r_missing_data(rmd, istatus)   

c     *** read in the raob data

c     +++ get latest raob file name

      call get_directory('snd',pathname_in,len)
      pathname_in = pathname_in(1:len)//'*'
c     pathname_in = '../lapsprd/snd/*'

      call get_file_names (pathname_in,numoffiles,c_filenames,
     1     maxfiles,istatus)

      if (istatus.ne.1) then    ! no data to read
         return

      elseif(numoffiles.eq.0) then ! no data to read
         return

      endif

      idx = index(c_filenames(numoffiles), ' ')

      filename = c_filenames(numoffiles)(idx-13:idx-4)

      call i4time_fname_lp (filename, raob_i4time, istatus)

c     +++ validate the raob time

      if (raob_i4time + look_back_time .lt. i4time) then ! too old
         write (6,*) 'raob data found is too old to be used'
         return
      endif

      if (raob_i4time .gt. i4time) then ! too new
         write (6,*) 'warning raob data found is too new to be used'
         write (6,*) 'assigning current file time for .snd'
         raob_i4time = i4time
         call make_fnam_lp (raob_i4time, filename, istatus)
      endif


c     +++ read raob file
      call get_directory('snd',fname,len)

      open (12, file = fname(1:len)//filename//'.snd',
     1     form='formatted',status='old',err=18)

      isound = 0  ! initialize isound (enable time rejection) 5/24/99 db

      do is = 1, snd_tot

         isound = isound + 1
                

 15      read (12,511,end=16) idummy,idx, lat_r(isound),lon_r (isound),
     1        rdummy, cdummy,r_filename(isound),snd_type(isound)
 511     format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)
         print*, idx
         lev_r(isound) = 1

         do i = 1, idx

            read(12,*) rspacing_dum,
     1           p_r (lev_r(isound),isound),
     1           t_r (lev_r(isound),isound),
     1           td_r(lev_r(isound),isound)

            if (
     1           (p_r (lev_r(isound),isound).ne.rmd)
     1           .and.
     1           (t_r (lev_r(isound),isound).ne.rmd)
     1           .and.
     1           (td_r (lev_r(isound),isound).ne.rmd)
     1           ) then

               lev_r(isound) = lev_r(isound) + 1

            else

               continue

            endif

         enddo !  i

         lev_r(isound) = lev_r(isound) - 1
c     reject on time condition (one hour lookback)
         if(r_filename(isound) .ne. filename) then ! examine closer
            call i4time_fname_lp (r_filename(isound), raob_i4time_chk, 
     1           istatus)
            if (abs(raob_i4time - raob_i4time_chk) .gt. 3600) then ! over 1hr
               write(6,*) 'rejecting on time bounds', r_filename(isound)
               isound = isound -1 !reject -- out of time bounds
            else
               write(6,*) 'accepting.. ', r_filename(isound)
                                !COMMENTING OUT RAOB SPECIFIC
                                !1.21 MOD  8/6/99 BEGIN TEST D.B.
c     test for non-RAOB type sounding
c               if (snd_type(isound) .ne. "RAOB    ") then ! reject (wrong type)
c                  write(6,*) 'rejecting ', snd_type(isound)
c                  isound = isound -1
c               endif
            endif
         else                   ! accept implicitly
            write(6,*) 'accepting.. time exact ', snd_type(isound)
c     test for non-RAOB type sounding
c            if (snd_type(isound) .ne. "RAOB    ") then ! reject (wrong type)
c               write(6,*) 'rejecting ', snd_type(isound)
c               isound = isound -1
c            endif

         endif


      enddo  ! dummy is loop

 16   close (12)

      isound = isound -1

      if (isound .le. 0 ) then ! no data found
         write(6,*) 'No usable RAOB data in database'
         return
      else
         write(6,*) isound, ' Number of RAOBs considered in analysis'
      endif

c     +++ compute the nearest i,j for all sounding locations
c     also compute weighting function for all locations for given raob
c     this approach used to accommodate any map projection or sized 
c     domain


      do k = 1,isound

         rmin = 1.e30

         do j = 1,jj
            do i = 1,ii

               r = sin(lat_r(k)*d2r)*sin(lat(i,j)*d2r)
     1              +cos(lat_r(k)*d2r) *cos(lat(i,j)*d2r)
     1              * cos( abs( lon_r(k) - lon(i,j)) *d2r )

               r = acos (r)

               call check_nan (r,istatus)
               if(istatus.eq.0) then ! nan detected bail
                  write(6,*) 'nan detected in reading raobs (r)'
                  write(6,*) 'i = ', i, 'k = ', k
                  return
               endif

               weight(i,j,k) = exp ( -1.*(r * 6419)**2/519370.)
               wt_ck(i,j) = wt_ck(i,j) + weight(i,j,k)

               call check_nan(weight(i,j,k),istatus)
               if(istatus.eq.0) then ! nan detected bail
                  write(6,*) 'nan detected in reading raobs (weight)'
                  write(6,*) i,j,k, 'i,j,k'
                  return
               endif 
                 
c     r is in radians, convert to degrees
c     r = r/d2r
c     r is in degrees convert to nautical miles
c     r = r * 60.  ! one minute = 1 nautical mile
c     r is in n. miles, convert to km
c     r = r * 1.852  ! km
c     note the factor 6419 converts radians to km.                  
c     used value of 519370 came from r=600km ---> 0.5 weight
c     using 600 km as "cutoff" distance, it was assumed that a reasonable
c     distance for haveing enough raob influence was about 1800 miles
c     this relates to R=1800*1.6 km or R = 2880 km.  at this distance
c     the associated weight is 1.159e-7.  The algorithm will accumulate
c     the weight for each gridpoint and see if there is enough influence
c     (sum >= 1.159e-7). this is the criteria used for determining enough
c     raob data to continue.
c     this step was inserted at this point in the code to avoid further
c     computations if there were not enough raobs.

c     rmin is the minimum distance to the raob (i,j of raob essentially)
               rmin = min (rmin,r)
               if (rmin.eq.r) then
                  i_r(k) = i
                  j_r(k) = j
               endif

            enddo
         enddo

      enddo
c     now check for enough raob data

      rdummy = 1000000.

      do j = 1,jj
         do i = 1,ii
            rdummy = min(wt_ck(i,j),rdummy)
         enddo
      enddo

      write(6,*) rdummy/float(isound), 'value of min weight in domian'

      if (rdummy/float(isound) .lt. 1.159e-7) then
         write(6,*) 'Not enough raob weighting in the domain'
         write(6,*) 'raob step aborting with no mods to field'
         return
      endif


c     *** interpolate q in the vertical at each raob location to the laps
c     pressure levels

      do is = 1,isound

         do k = 1,kk            ! each laps pressure level (ks)

            call locate(p_r(1,is),lev_r(is),laps_pressure(k),ks)

            if(ks.eq.0  .or. ks.eq.lev_r(is)) then

               q_r(k,is) = rmd

            else                ! interpolate to the desired level

               call interp (log(laps_pressure(k)),
     1              log(p_r(ks,is)),log(p_r(ks+1,is)),
     1              t_r(ks,is),t_r(ks+1,is),temt)

               call interp (log(laps_pressure(k)),
     1              log(p_r(ks,is)),log(p_r(ks+1,is)),
     1              td_r(ks,is),td_r(ks+1,is),temtd)

c     adjustment to t_ref for sounding data assuming now Td WRT liq at
c     temps below freezing, DB  1.20 change

               q_r(k,is) = ssh2 (laps_pressure(k),temt,temtd,-129.)/1000.

c     accept only valid raob data (temp must be > -40 c
c     and +/- 3c temp check on laps background to raob temp)  DB 1.18 change
               if(  temt.lt.-40. ! raob is not effective .lt. -40c
     1              .or. 
     1              abs(temt+273.-laps_t(i_r(is),j_r(is),k)).ge. 3.0
     1              ) then  !qc failure on temp
                  write(6,*) 'qc failure on temp, raob_step 1.18 mod'
                  write(6,*) laps_pressure(k), temt,
     1                 abs(temt+273.-laps_t(i_r(is),j_r(is),k)), is
                  q_r(k,is) = rmd
               endif
 

               call check_nan(q_r(k,is),istatus)
               if(istatus.eq.0) then ! nan detected bail
                  write(6,*) 'nan detected in processing raobs (q_r)'
                  write(6,*) 'k, is ', k,is
                  return
               endif 
       
            endif

         enddo                  ! k

      enddo                     ! is

c *** difference the raob data at each gridpoint location and height

      do is = 1,isound

         do k = 1,kk

            if (
     1           (data(i_r(is),j_r(is),k) .eq. rmd)
     1           .or.
     1           (q_r(k,is) .eq. rmd)
     1           .or.
     1           (data(i_r(is),j_r(is),k) .le. 0.0)
     1           ) then

               diff(k,is) = rmd

            else

               diff(k,is) = q_r(k,is) - data(i_r(is),j_r(is),k)
               scale(k,is) = q_r(k,is)/data(i_r(is),j_r(is),k)
            endif

         enddo

      enddo

c     compute level mean and sigma  modified from original... now examines
c     bias data (diff variable)

      write(6,*) 'compting layer statistics'
      
      do k = 1,kk               ! do for each layer
         x_sum = 0
         do is = 1,isound       ! do for all differences each layer
            if(diff(k,is) .ne. rmd ) then
               x_sum=x_sum +1
               x(x_sum) = diff(k,is)
            endif

c     all data assembled for given layer, compute stats for that layer.

         enddo                  ! is
      

 220     call moment_b (x,x_sum,ave(k),adev(k),sdev(k),
     1        var(k),skew(k),curt(k), istatus)

         if(istatus.ne.1) go to 222

c     this loop put in to recursivly recompute sigma after changing
c     the population.  here use 3*sigma criteria  (normal full population)
         
         do is = 1, x_sum
            if(x(is).gt.ave(k)+3.*sdev(k) .or. 
     1           x(is).lt.ave(k)-3.*sdev(k)) then ! reject
               do i = is,x_sum-1 ! discard bad x(is), shift
                  x(is) = x(is + 1)
               enddo
               x_sum = x_sum -1  !decrement x_sum
               go to 220        !recompute mean
            endif
         enddo



 222     if (istatus.ne.1) then
            write(6,*) 'not enough data for computing moments level ', k
            write(6,*) 'this will result in NO raob use for this level'
            ave(k) = rmd
            sdev(k) = rmd
         else
            x_sum_save (k) = x_sum
            write (6,*) x_sum,ave(k),sdev(k),' level ',k
         endif

      enddo                     ! k

c     new   QC check, data adjustment against climo values
c     vsn 1.19

c     call climo routine for middle latitude

      center_lat = lat(ii/2,jj/2)
      read (filename(3:5),44) jday
 44   format (i3)

      call climate_sm(center_lat, jday, standard_press, tempertur_guess,
     1    mixratio_guess, istatus)

c     interpolate from climo coordinates to laps coords for w

      do k = 1, kk  ! fill all laps levels from climo data
         call locate(standard_press,40,laps_pressure(k),ks)
         if(ks.eq.0  .or. ks.eq.40) then
            climo_w(k) = rmd
         else                   !interpolate
            call interp(log(laps_pressure(k)),
     1              log(standard_press(ks)),log(standard_press(ks+1)),
     1              mixratio_guess(ks),mixratio_guess(ks+1),climo_w(k))
         endif
      enddo ! k

c     mixing ratio and SH are close so we won't differentiate for this QC
c     step
     
c     compute all ratios of clim_w(k) to diffs

      do k = 1,kk
         if( ave(k) .ne. rmd) then
            write(6,*) (abs( ave(k))+3.*sdev(k)) / climo_w (k)*1.e-3
     1           /x_sum_save(k) ,
     1           laps_pressure(k)
            if ((abs( ave(k))+3.*sdev(k)) / climo_w (k)*1.e-3
     1           /x_sum_save(k) .gt. 0.5e-6) then ! reject criteria
               write(6,*) 'Climo Q QC reject criteria met 1.19 ', k
               ave(k) = rmd
               sdev(k) = rmd
            endif
         endif
      enddo




c     compare computed diff values with QC thresholding... assign rmd if bad
c     here use 2*sigma cutoff for addnl constraint.
      
      do is = 1,isound
         do k = 1,kk
            
            if(sdev(k).ne.rmd) then
 
               if   (
     1              diff(k,is) .le. ave(k)+2.*sdev(k) 
     1              .and. 
     1              diff (k,is) .ge. ave(k)-2.*sdev(k)
     1              ) then      ! in accept range
                  write(6,*) 'accepting ',k,is,diff(k,is)
               else
                  write(6,*) 'tossing (exceeds limits)',k,is, diff(k,is)
                  diff(k,is) = rmd
               endif
            else
               write(6,*) 'tossing level, too few stats level ',k
               diff(k,is) = rmd
            endif
            
         enddo
      enddo


c     *** analyze the field (barnes second scaling pass to background)

      do k = 1,kk
         do j = 1,jj
            do i = 1,ii

               tot_weight = 0.0
               tot_scale_weight = 0.0
               scale_used = 0.0
               weight_avg = 0.0
               counter = 0.0

               do is = 1,isound

                  if (diff(k,is).eq.rmd) then ! dont modify
                     continue
                  elseif(data(i,j,k).eq.rmd) then ! dont modify
                     continue
                  else
                     tot_scale_weight = weight(i,j,is)*scale(k,is) 
     1                    + tot_scale_weight
                     tot_weight = tot_weight + weight(i,j,is)
                     counter = counter + 1.0
                  endif

               enddo

               if   (
     1              (data(i,j,k) .ne. rmd)
     1              .and.
     1              (tot_weight.ne.0.0) 
     1              .and.
     1              (counter.ne.0.0)
     1              )   then

c     compute the appropriate scale factor to apply (weighted average)

                  scale_used = tot_scale_weight/tot_weight

c     compute the weight average for this location, this introduces 
c     possiblity for circles, but on the other hand reduces them due to
c     exaggerated saturations going into regions that are near saturation
c     already

                  weight_avg = tot_weight/counter

c     apply the computed weights and averages to the data point in question
c     the first approach (commented out) was dropped after the second 
c     method was found to work well after assurance of a critical mass of
c     data points.

c             data(i,j,k) = data(i,j,k)*(scale_used-1.)*weight_avg 
c     1                 + data (i,j,k)

                  data(i,j,k) = data(i,j,k)*(scale_used-1.)
     1                 + data (i,j,k)

                  call check_nan(data(i,j,k),istatus)
                  if(istatus.eq.0) then ! nan detect
                     write(6,*) 'Nan Detected during raob_step'
                     write(6,*) 'data(i,j,k), i,j,k'
                     write(6,*) data(i,j,k), i,j,k
                     write(6,*) 'diff(k,is),is = 1,isound'
                     write(6,*) (diff(k,is),is = 1,isound)
                  endif
                

                  if (data(i,j,k) .lt. 0.0) then
                     write(6,*) 'data zero scale_used, weight_avg, i,j'
                     write(6,*) scale_used, weight_avg, i,j 
                     data(i,j,k) = 0.0
                  endif
               endif
c     if conditions are not met, data is unchanged

            enddo
         enddo
      enddo

c     *** return modified data array

      write(6,*) 'end RAOB step'

      return

 18   write(6,*) 'error opening requested .snd file, raob step skipped'
      return

 24   write(6,*)  'raob moisture switch missing!, ignoring raobs'

      return

      end

