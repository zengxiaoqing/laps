      subroutine adjust_geog(nnxp,nnyp,ncat,ctype
     &,istat_dat,lat,topt_out           !istattmp,istatslp,lat,topt_out,path_to_soiltemp
     &,landmask,geog_data,istatus)      !soiltemp_1deg,greenfrac,islope,istatus)
c
c This routine uses landmask to make the course resolution
c soil temp and green fraction data conform to water-land mask.
c It also fills small islands or isolated land
c bodies with appropriate geog values as necessary
c
c J. Smart NOAA/FSL
c     "    06-22-01 original
c     "    12-01-03 put subroutine into gridgen_utils, make arrangements for soil
c                   type data, and allow separate calls for each geog type to reduce
c                   memory requirement.
c
      implicit none

      integer nnxp,nnyp
      integer i,j,l,ii,jj
      integer is,js
      integer istatus
      integer istat_dat    !Input. =0 indicates input geog data was NOT processed properly; =1 otherwise.
      integer ijsthresh    !search distance to look for representative vlues to fill inconsistency.
      integer isum
      integer l1,l2
      integer lent
      integer ncat
      integer ifixw,ifixl
      integer ifixws,ifixwt,ifixls

      integer isc
      integer ic(ncat) 

      character*132 path_to_soiltemp
      character*(*) ctype

      logical endsearch

      real,intent(inout)  ::    lat(nnxp,nnyp)
      real,intent(inout)  ::    landmask(nnxp,nnyp)
      real,intent(inout)  ::    topt_out(nnxp,nnyp)
      real,intent(inout)  ::    geog_data(nnxp,nnyp,ncat)

c     real,intent(inout)  ::    soiltemp_1deg(nnxp,nnyp)
c     real,intent(inout)  ::    greenfrac(nnxp,nnyp,ncat)
c     real,intent(inout)  ::    islope(nnxp,nnyp)

      real,   allocatable ::    geog_tmp      (:,:,:) !temporary holder of input geog data.

c     real,   allocatable ::    grnfrctmp     (:,:,:)
c     real,   allocatable ::    soiltmp       (:,:)
c     real,   allocatable ::    islopetmp     (:,:)

      real,   allocatable ::    rmeanlattemp  (:)

      real    avgtmp
      real    avggrn(ncat)
      real    avgcat
      real    tatlmx,tatlmn
      real    rlatmx,rlatmn
      real    r_missing_data
      real    rmngrn
      real    islp
      real    sumt
      real    sumg
      real    sum(ncat)
      real    tslp

      istatus=0
      if(istat_dat.eq.0)then  !.and.istattmp.eq.0.and.istatslp.eq.0)then
	 print*,'Unable to process geog data in adjust_geog ...'
     &,' processsing of data failed prior to subroutine call.'
	 return
      endif


c use moist adiabatic laps rate (6.5 deg/km) to get new temp
 
      allocate (geog_tmp(nnxp,nnyp,ncat))

      call get_r_missing_data(r_missing_data,istatus)

      geog_tmp=geog_data

      where(geog_tmp.eq.0.0)geog_tmp=r_missing_data

      if(ctype.eq.'greenfrac')then
         rmngrn=minval(geog_tmp(1:nnxp,1:nnyp,1:1))
         print*,'minimum green fraction = ',rmngrn
      endif

c determine average soiltemp and greenfrac in domain
      sumt=0.0
      sum=0.0
      isc=0
      ic=0
      islp=0
      isum=0
      if(ctype.eq.'soiltemp')then
         do j = 1,nnyp
         do i = 1,nnxp
            if(geog_tmp(i,j,1).ne. r_missing_data)then
               sumt=geog_tmp(i,j,1)+sumt
               isc=isc+1
            endif
         enddo
         enddo
      elseif(ctype.eq.'greenfrac')then
c greenfrac is assumed to be continuous globally; only use land points
         do j = 1,nnyp
         do i = 1,nnxp
            do l=1,ncat
               if(landmask(i,j) .ne. 0 .and.
     &geog_tmp(i,j,l).lt.r_missing_data)then
                  sum(l)=geog_tmp(i,j,l)+sum(l)
                  ic(l)=ic(l)+1
               endif
            enddo
         enddo
         enddo
      elseif(ctype.eq.'islope'.or.ctype.eq.'soiltype')then
         do j = 1,nnyp
         do i = 1,nnxp
            if(geog_tmp(i,j,1).ne. r_missing_data)then
               isum=geog_tmp(i,j,1)+isum
               islp=islp+1
            endif
         enddo
         enddo
      endif
c
c get average information for filling as necessary
c ----------------------------------------------------
c
c 1. deep soil temp.
c
      if(isc.gt.0)then
         avgtmp=sumt/float(isc)
         print*,'Domain average annual mean temp = ',avgtmp
      elseif(ctype.eq.'soiltemp')then
c
c this section uses mean latitudinally averaged temps derived from 
c the raw 1 deg temp data. File in raw geog annual mean deep soil
c temp directory (see wrfsi.nl variable soiltemp_1deg).
c
         call get_path_to_soiltemp_1deg(path_to_soiltemp,istatus)
         sumt=0.0
         print*,'*** Using mean lat temps -> LATMEANTEMP ***'
         allocate (rmeanlattemp(180))
         call  get_directory_length(path_to_soiltemp,lent)
         call  get_meanlattemp(path_to_soiltemp(1:lent-1)
     &,rmeanlattemp,istatus)
         if(istatus.ne.1)then
            print*,'Error returned: get_meanlattemp'
            return
         endif
         l1=90-nint(minval(lat(:,1)))+1
         l2=90-nint(maxval(lat(:,nnyp)))+1
         if(l1.gt.180)l1=180
         if(l2.gt.180)l2=180

         rlatmn=minval(lat(:,1))
         rlatmx=maxval(lat(:,nnyp))
         tatlmx=rmeanlattemp(l2)
         tatlmn=rmeanlattemp(l1)
         if(tatlmx.ne.0.and.tatlmn.ne.0.and.
     .      tatlmx.ne.tatlmn)then
            tslp=(rlatmx-rlatmn)/(tatlmx-tatlmn)
            do j=1,nnyp
            do i=1,nnxp
               geog_tmp(i,j,1)=tatlmn+(rlatmx-lat(i,j))/tslp
               sumt=geog_tmp(i,j,1)+sumt
            enddo
            enddo
            avgtmp=sumt/(nnyp*nnxp)
         elseif(tatlmx.gt.0)then
            avgtmp=tatlmx
         elseif(tatlmn.gt.0)then
            avgtmp=tatlmn
         else
            print*,'Unusual condition with meanlattemp'
            istatus = 0
            return
         endif
         deallocate (rmeanlattemp)
         print*,'Domain average annual mean temp = ',avgtmp
      endif
c
c 2. green fraction.
c
      if(ctype.eq.'greenfrac')then
         avggrn=r_missing_data
         do l=1,ncat
            if(ic(l).gt.0)then
               avggrn(l)=sum(l)/float(ic(l))
            else
               avggrn(l)=0.0
            endif
            print*,'Domain average greenfrac = ',l,avggrn(l)
         enddo
c
c 3. Terrain slope index or soil type category..
c
      elseif(ctype.eq.'islope'.or.ctype.eq.'soiltype')then
        avgcat=r_missing_data
        if(islp.gt.0)then
           avgcat=float(isum)/islp
           print*,'Average category: ',avgcat
        else
           print*,'Could not compute average category ?'
           print*,'Average category: ',avgcat
        endif
      endif

c extend search to a fraction of the domain size. Could improve this for
c ratio geog-data-res/domain-res (possibly) to avoid unreasonable
c search distance.

      ijsthresh = int(nnxp/2)
      ifixw=0
      ifixl=0

      if(ctype.eq.'soiltemp')then

         do j = 1,nnyp
         do i = 1,nnxp

            if(landmask(i,j).eq.1)then                   !a land point
               if(geog_tmp(i,j,1).eq.r_missing_data)then !inconsistent because it is water
 
                  is=1
                  js=1
                  endsearch = .false.

                  sumt=0.0
                  isc=0

                  do while (.not.endsearch)

                   do jj=j-js,j+js
                   do ii=i-is,i+is

                      if((ii.ge.1) .and. (ii.le.nnxp) .and.
     &                   (jj.ge.1) .and. (jj.le.nnyp)) then

                       if(geog_tmp(ii,jj,1).ne.r_missing_data)then
                          sumt=sumt+geog_tmp(ii,jj,1)
                          isc=isc+1
                       endif

                      endif

                   enddo
                   enddo

                   if(isc.gt.0)then
                      geog_data(i,j,1)=-0.0065*topt_out(i,j)+sumt/isc
                      ifixw=ifixw+1    !count the # of inconsistent water points fixed with rep land value
                      endsearch=.true.
                   else
                      is=is+1
                      js=js+1
                      if(is.gt.ijsthresh)endsearch=.true.
                   endif

                  enddo

               else  !no inconsistency

                  geog_data(i,j,1)=-0.0065*topt_out(i,j)
     &+geog_data(i,j,1)

               endif

            else     !landmask says this is a water point

               if(geog_tmp(i,j,1).ne.r_missing_data)then
                  geog_data(i,j,1)=r_missing_data
                  ifixl=ifixl+1            !count the # of fixed land points
               endif

            endif

         enddo
         enddo

      elseif(ctype.eq.'greenfrac')then

         do j = 1,nnyp
         do i = 1,nnxp

            if(landmask(i,j).eq.1)then                  !a land point

               if(geog_tmp(i,j,1).eq.r_missing_data
     &.and.rmngrn.lt.r_missing_data)then  !inconsistency and there is valid data to search for.

                  endsearch = .false.

                  sum=0.0
                  ic=0
                  is=1
                  js=1
        
                  do while (.not.endsearch)

                   do ii=i-is,i+is
                   do jj=j-js,j+js

                    if( (ii.ge.1) .and. (ii.le.nnxp)
     &             .and.(jj.ge.1) .and. (jj.le.nnyp)) then

                       if(landmask(ii,jj).eq.1.and.
     &         geog_tmp(ii,jj,1).lt.r_missing_data)then

                          do l=1,12
                             sum(l)=sum(l)+geog_tmp(ii,jj,l)
                             ic(l)=ic(l)+1
                          enddo
                       endif

                    endif

                   enddo
                   enddo

                   if(ic(1).gt.0)then
                      do l=1,12
                         geog_data(i,j,l)=sum(l)/float(ic(l))
                      enddo
                      ifixw=ifixw+1
                      endsearch=.true.
                   else
                      is=is+1
                      js=js+1
                      if(is.gt.ijsthresh.or.js.gt.ijsthresh)
     &                   endsearch=.true.
                   endif

                  enddo

               else

                  geog_data(i,j,:)=geog_tmp(i,j,:)

               endif

            else  !this is a water point

c              do l=1,12
c                 if(geog_tmp(i,j,l).ne. 0.0 .and.
c    &               geog_tmp(i,j,l).lt.r_missing_data)then
c                    geog_data(i,j,l)=0.0
c                    ifixl=ifixl+1  !count # of land points that should be water
c                 endif
c              enddo

               if(geog_tmp(i,j,1).ne. 0.0 .and.
     &               geog_tmp(i,j,1).lt.r_missing_data)then
                     geog_data(i,j,1:12)=0.0
                     ifixl=ifixl+1  !count # of land points that should be water
               endif

            endif

         enddo
         enddo
c
c -------------------------------------------------------
c this section for categories (terrain slope or soiltype)
c
      elseif(ctype.eq.'islope'.or.ctype.eq.'soiltype')then

         do j = 1,nnyp
         do i = 1,nnxp

            if(landmask(i,j).eq.1)then                  !a land point

               if(geog_tmp(i,j,1).eq.r_missing_data)then !an inconsistency

                  endsearch = .false.

                  sumt=0.0
                  isc=0
                  islp=-9
                  is=1
                  js=1

                  do while (.not.endsearch)
                     do ii=i-is,i+is
                     do jj=j-js,j+js

                        if( (ii.ge.1) .and. (ii.le.nnxp)
     &                 .and.(jj.ge.1) .and. (jj.le.nnyp)) then

                           if(landmask(ii,jj).eq.1.and.
     &                        geog_tmp(ii,jj,1).lt.r_missing_data)then
                              islp=geog_tmp(ii,jj,1)
                              sumt=sumt+islp
                              isc=isc+1
                           endif

                        endif

                     enddo
                     enddo
 
                     if(ctype.eq.'islope')then
                        if(isc.gt.0)then
                           if(islp.gt.7.0)islp=13.0 !force mean to glacial ice as necessary
                           if(sumt.lt.isc)then
                              islp=1.0
                           else
                              islp=float(nint(sumt/float(isc)))
                           endif
                           geog_data(i,j,1)=islp
                           endsearch=.true.
                           ifixw=ifixw+1  !water point fixed to be land point
                        else
                           is=is+1
                           js=js+1
                           if(is.gt.ijsthresh.or.js.gt.ijsthresh)
     &                        endsearch=.true.
                        endif
                     else
                        if(isc.gt.0)then
                           if(sumt.lt.isc)then
                              islp=1.0
                           else
                              islp=float(nint(sumt/float(isc)))
                           endif
                           geog_data(i,j,1)=islp
                           endsearch=.true.
                           ifixw=ifixw+1  !water point fixed to be land point
                        else
                           is=is+1
                           js=js+1
                           if(is.gt.ijsthresh.or.js.gt.ijsthresh)
     &                        endsearch=.true.
                        endif
                     endif

                  enddo

               else

                  geog_data(i,j,1)=geog_tmp(i,j,1)

               endif

            else     !this is a water point

               if(geog_tmp(i,j,1).ne.r_missing_data)then
                  geog_data(i,j,1)=r_missing_data
                  if(ctype.eq.'islope')
     &            ifixl=ifixl+1   !land point fixed to be water point
               endif

            endif

         enddo
         enddo

      else

         print*,'Unknown type input to adjust_geog!'
         return

      endif ! ctype switch

      deallocate (geog_tmp)

c if the above search failed to find nearby soil temp or greenness frac
c then use average value

      if(ctype.eq.'soiltemp')then
         do j = 1,nnyp
         do i = 1,nnxp
            if(landmask(i,j).eq.1)then  !a land point
               if(geog_data(i,j,1).eq.r_missing_data)then
                  geog_data(i,j,1) = avgtmp-0.0065*topt_out(i,j)
                  ifixw=ifixw+1
               endif
            endif
         enddo
         enddo
         print*
         print*,'Soil Temp stats'

      elseif(ctype.eq.'greenfrac')then

         do j = 1,nnyp
         do i = 1,nnxp
            if(landmask(i,j).eq.1)then  !a land point
               do l = 1,12
                  if(geog_data(i,j,l).eq.0.0)then
                     geog_data(i,j,l)=float(nint(avggrn(l)))
                     if(l==1)ifixw=ifixw+1
                  endif
               enddo
            endif
         enddo
         enddo
         print*
         print*,'Green fraction stats'

      elseif(ctype.eq.'islope'.or.ctype.eq.'soiltype')then

         do j = 1,nnyp
         do i = 1,nnxp

            if(landmask(i,j).eq.1)then  !a land point
               if(geog_data(i,j,1).eq.r_missing_data)then
                  geog_data(i,j,1)=float(nint(avgcat))
                  if(ctype.eq.'islope')then
                     ifixw=ifixw+1
                  else
                     ifixl=ifixl+1
                  endif
               endif
            endif

        enddo
        enddo
        print*
        if(ctype.eq.'islope')then
           print*,'Terrain Slope Index stats'
        else
           print*,'Soil Type Category stats'
        endif

        where(geog_data.eq.r_missing_data)geog_data=0.0  !put water category back to original value

      else

        print*,'unknown geog data from ctype ',ctype

      endif

      print*,'--------------------------'
      if(ctype.ne.'islope'.and.ctype.ne.'soiltype')then
         if(ctype.eq.'greenfrac')then
            print*,'Fixed ',ifixw,'  points with rep land value'
            print*,'Fixed ',ifixl,'  points with missing value = water'
            print*
            do j = 1,nnyp
            do i = 1,nnxp
            do l = 1,ncat
               if(landmask(i,j) .ne. 0 .and.
     &geog_data(i,j,l).gt.0.0)then
                  sum(l)=geog_data(i,j,l)+sum(l)
                  ic(l)=ic(l)+1
               endif
            enddo
            enddo
            enddo
            do l=1,ncat
               if(ic(l).gt.0)then
                  avggrn(l)=sum(l)/float(ic(l))
               else
                  avggrn(l)=0.0
               endif
               print*,'Domain average greenfrac = ',l,avggrn(l)
            enddo
         else
            print*,'Fixed ',ifixw,'  points with rep land value'
            print*,'Fixed ',ifixl,'  points with missing value = water'
         endif
      elseif(ctype.eq.'islope')then
         print*,'Fixed ',ifixw,'  points with rep land value'
         print*,'Fixed ',ifixl,'  points with missing value = water'
      else
         print*,'Fixed ',ifixw,'  points with rep land value'
         print*,'Fixed ',ifixl,'  points with missing value = water'
      endif

      istatus=1

      return
      end
