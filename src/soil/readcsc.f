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
         subroutine readcsc(i4time_cur,imax,jmax,cscnew)
C
         integer  imax,jmax

         real csc(imax,jmax),lcv(imax,jmax),csctot(imax,jmax)
         real csc2nd(imax,jmax)
         real snow_total(imax,jmax)
         real snow_accum(imax,jmax)
         real missval,cscnew(imax,jmax)
         integer i4time,i4time_cur
         integer i4time_csc(imax,jmax),i4time_snow(imax,jmax)
         CHARACTER*24 TIME
c ************************************************************
c
         do j=1,jmax
            do i=1,imax
               snow_total(i,j)=0.0
            enddo
         enddo
        call get_laps_cycle_time(laps_cycle_time,istatus)
        if(istatus.ne.1)then
           write(6,*)'Error - get_laps_cycle_time'
           return
        endif
        call get_r_missing_data(missval,istatus)
        if(istatus.ne.1)then
           write(6,*)'Error - get_r_missing_data'
           return
        endif

        i4time=i4time_cur-48*3600 !force 48 hours 

c ************************************************************
c initialize arrays                                    
c
        csctot=missval
        csc2nd=missval
        i4time_csc = -999
        i4time_snow = -999
c
c now loop through past 48 hours, looking only at csc files
c
        rmxlcv=0.0
        rmnlcv=100.
        icsc=0
        ncycle_times=48*int(3600./float(laps_cycle_time))
        do itime=1,ncycle_times
c
           call cv_i4tim_asc_lp(i4time,time,istatus)
c           if(time(13:14).eq.'13')goto 52
           write(6,101)time
 101       format(1x,'First csc loop data search time: ',a17)
c
	   CALL GETLAPSLCV(I4TIME,LCV,CSC,IMAX,JMAX,ISTATUS)
	   IF (ISTATUS .NE. 1)go to 522

           icsc=icsc+1
           nspts=0
           nrpts=0
           do i=1,imax
           do j=1,jmax
             if(csc(i,j).le.1.1)then
                csctot(i,j)=csc(i,j)
                nspts=nspts+1
             endif
             if(csc(i,j).le.0.1)then
                nrpts=nrpts+1
                csctot(i,j)=0.0
             endif
             rmxlcv=max(lcv(i,j),rmxlcv)
             rmnlcv=min(lcv(i,j),rmnlcv)
           enddo
           enddo
           nmpts=(imax*jmax)-nspts-nrpts
           nspts=nspts-nrpts
           write(6,501)itime,nmpts,nspts,nrpts
 501       format(' Cycle ',i4,' nmiss/nsnow/nreject ',3i9)
 522       continue
c
           i4time=i4time+laps_cycle_time
        enddo
c
        if(icsc.eq.0)then
        write(6,*) '**** No lcv data available over 48 hours ****'
          do i=1,imax
          do j=1,jmax
             csctot(i,j)=missval
          enddo
          enddo
        endif
c
        call analyze(csctot,cscnew,imax,jmax,missval)
c
c now loop through past 48 hours, looking for latest snow_total and
c  latest csc obs. snow_total only applies to a point - it is not
c  spread horizontally. If there are then any later csc obs at that 
c  point, then the csctot point is set to the csc ob.
c
        i4time=i4time_cur-48*3600
        ncycle_times=48*int(3600./float(laps_cycle_time))
        write(6,*)'---------------------------------------'
        write(6,*)

        do itime=1,ncycle_times
c
           call cv_i4tim_asc_lp(i4time,time,istatus)
c          if(time(13:14).eq.'13')goto 52
           write(6,102)time
 102       format(1x,'Second loop data search time: ',a17)
c
c First get csc field again. Will not analyze - just use this to
c   write over any previous time's snow_total obs.
	   CALL GETLAPSLCV(I4TIME,LCV,CSC,IMAX,JMAX,ISTATUS)
	   IF (ISTATUS .NE. 1)goto 52

           do i=1,imax
           do j=1,jmax
             if(csc(i,j).le.1.1)then
                 cscnew(i,j)=csc(i,j)  
                 csc2nd(i,j)=csc(i,j)  
                 i4time_csc(i,j) = i4time
             endif
	     if(csc(i,j).le.0.1)then
                 cscnew(i,j)=0.0
                 csc2nd(i,j)=0.0
                 i4time_csc(i,j) = i4time
             endif
           enddo
           enddo
 52        continue
           i4time=i4time+laps_cycle_time
        enddo
c
c now get snow_total obs
c J.Smart (11-3-97). Compute 48 hour snow accumulation instead of
c                    snow total to avoid false snow cover from "old"
c                    snow total accumulation.
c Recent snow is given more weight
c
        i4time=i4time_cur-48*3600
        ncycle_times=48*int(3600./float(laps_cycle_time))
        write(6,*)
        write(6,*)'Accumulate hrly Snow for 48 hour total '
        write(6,*)'---------------------------------------'

        do itime=1,ncycle_times

           call cv_i4tim_asc_lp(i4time,time,istatus)
           write(6,103)time
 103       format(1x,'Snow Accum data search time: ',a17)

           CALL GETLAPSL1S(I4TIME,snow_accum,imax,jmax,1,ISTATUS)
           IF (ISTATUS .EQ. 1)then   
             do j=1,jmax
             do i=1,imax
               snow_total(i,j)=snow_total(i,j)
     &       +(snow_accum(i,j)*(float(itime)/float(ncycle_times)))
               if(snow_accum(i,j) .gt. 0.)then
                 i4time_snow(i,j) = i4time
               endif
             enddo
             enddo
             call array_range(snow_accum,imax,jmax,rmin1,rmax1,missval)          
             call array_range(snow_total,imax,jmax,rmin2,rmax2,missval)          
             write(6,111)rmin1,rmax1,rmin2,rmax2 
111          format(' L1S/S01, snow total ranges ',4f9.5)
           else
             write(6,*)'WARNING - L1S/S01 unavailable'
           endif     

           i4time=i4time+laps_cycle_time
        enddo
c
c Adjust snow cover field based upon snow total. ------
c   if snowfall total >= 2cm, then set snow cover to 1.
c   if snowfall total  = 1cm, then set snow cover to 0.5
c   if snowfall total between 1 and 2cm, then ramp snow cover 0.5 to 1.
c
c We added an extra condition that snowfall is more recent than satellite
c data.                                                                   

        nrecent_snow = 0
        nrecent_csc2nd = 0
        n_neither = 0
        do i=1,imax
        do j=1,jmax
c         if(csc2nd(i,j) .eq. missval)then
c         if(i4time_csc(i,j) .lt. i4time_cur-86400)then
          if(i4time_snow(i,j) .gt. i4time_csc(i,j))then
            nrecent_snow = nrecent_snow + 1

            if(snow_total(i,j) .ge. 0.02)cscnew(i,j)=1.0
            if((snow_total(i,j).ge.0.01).and.(snow_total(i,j).lt.0.02))
     1         cscnew(i,j)=snow_total(i,j)/0.02

c           if(snow_total(i,j).ge.0.02.and.cscnew(i,j).gt.0.25)
c    1         cscnew(i,j)=1.0
c           if((snow_total(i,j).ge.0.01).and.(snow_total(i,j).lt.0.02))
c    1         cscnew(i,j)=snow_total(i,j)/0.02

          elseif(csc2nd(i,j) .ne. missval)then ! test csc2nd is recent   
              nrecent_csc2nd = nrecent_csc2nd + 1
          else
              n_neither = n_neither + 1
          endif

        enddo
        enddo     

        write(6,*)' recent snow(S01)/csc(sat)/none = ',nrecent_snow
     1                                                ,nrecent_csc2nd
     1                                                ,n_neither
c
c       do j=jmax,1,-1
c          write(6,76)j,(cscnew(i,j),i=1,18)
c 76       format(1x,i2,1x,18f4.1)
c       enddo
c
        return
        end
C-------------------------------------------------------------------------------
C
      SUBROUTINE GETLAPSLCV(I4TIME,LCV,CSC,IMAX,JMAX,ISTATUS)
C
      Integer*4 imax,jmax
      Integer*4 kmax
      parameter(kmax = 2)
c
      INTEGER*4 I4TIME,LVL(KMAX),I,J,ERROR(2),ISTATUS
C
      REAL*4 lcv(imax,jmax),csc(imax,jmax)
      Real*4 readv(imax,jmax,kmax)
C
      CHARACTER*150 LDIR
      CHARACTER*31 EXT
      CHARACTER*3 VAR(KMAX)
      CHARACTER*4 LVL_COORD(KMAX)
      CHARACTER*10 UNITS(KMAX)
      CHARACTER*125 COMMENT(KMAX)
C
C-------------------------------------------------------------------------------
C
	ERROR(1)=1
	ERROR(2)=0
C
        call get_directory('lcv',ldir,len)

	EXT='lcv'
	VAR(1)='LCV'
	VAR(2)='CSC'
	LVL(1)=0
	LVL(2)=0
C
C ****  Read LAPS lcv and csc.
C
	CALL READ_LAPS_DATA(I4TIME,LDIR,EXT,IMAX,JMAX,KMAX,KMAX,VAR,LVL,
     1      LVL_COORD,UNITS,COMMENT,readv,ISTATUS)
C
	IF (ISTATUS .NE. 1) THEN
		ISTATUS=ERROR(2)
		RETURN
	ENDIF
C
        do j=1,jmax
        do i=1,imax
          csc(i,j)=readv(i,j,2) !Cloud analysis implied snow cover
          lcv(i,j)=readv(i,j,1) !LAPS cloud cover
        enddo
        enddo
c
c        print *,' '
c        print *,'lcv field'
c        do j=1,jmax
c          write(6,20)j,(lcv(i,j),i=21,40)
c        enddo
c        print *,' '
c        print *,'csc field'
c        do j=1,jmax
c          write(6,20)j,(csc(i,j),i=21,40)
c        enddo
c        print *,'csc(21,55)=',csc(21,55)
c 20     format(1x,i2,1x,20f4.1)
c
	ISTATUS=ERROR(1)
	RETURN
C
	END
C
C-------------------------------------------------------------------------------
C
      SUBROUTINE GETLAPSL1S(I4TIME,snow_accum,imax,jmax,kmax,
     &ISTATUS)
C
      integer*4 imax,jmax,kmax
c
      INTEGER*4 I4TIME,LVL(KMAX),I,J,ERROR(2),ISTATUS
C
      REAL*4 snow_accum(imax,jmax),readv(imax,jmax,kmax)		
C
      CHARACTER*150 LDIR
      CHARACTER*31 EXT
      CHARACTER*3 VAR(KMAX)
      CHARACTER*4 LVL_COORD(KMAX)
      CHARACTER*10 UNITS(KMAX)
      CHARACTER*125 COMMENT(KMAX)
C
C-------------------------------------------------------------------------------
C
	ERROR(1)=1
	ERROR(2)=0
C
        call get_directory('l1s',ldir,len)

	EXT='L1S'
	VAR(1)='S01'
	LVL(1)=0
C
C ****  Read LAPS l1s
C
	CALL READ_LAPS_DATA(I4TIME,LDIR,EXT,IMAX,JMAX,KMAX,KMAX,VAR,LVL,
     1      LVL_COORD,UNITS,COMMENT,readv,ISTATUS)
C
	IF (ISTATUS .NE. 1) THEN
		ISTATUS=ERROR(2)
		RETURN
	ENDIF
C
        do j=1,jmax
        do i=1,imax
          snow_accum(i,j)=readv(i,j,1)
        enddo
        enddo
c
c
	ISTATUS=ERROR(1)
	RETURN
C
	END
C
C
C-------------------------------------------------------------------------------
         subroutine analyze(csctot,cscnew,imax,jmax,missval)
C
         integer imax,jmax
         integer nboxes, nruns
         parameter(nboxes=100,nruns=5)
c
         real csctot(imax,jmax),cscnew(imax,jmax)
         dimension iimin(nruns,nboxes),
     1       iimax(nruns,nboxes),jjmin(nruns,nboxes),
     2       jjmax(nruns,nboxes),nbtot(nruns)
         real*4 missval
c
c set up boundaries for boxes to do averaging over
c
         do i=1,imax
         do j=1,jmax
            cscnew(i,j)=csctot(i,j)
         end do
         end do
         nbtot(1)=1
         nbtot(2)=4
         nbtot(3)=9
         nbtot(4)=36
         nbtot(5)=100
         call setboxes(imax,jmax,iimin,iimax,jjmin,jjmax,
     1    nbtot,nruns,nboxes)   
c  
c **************************************************************
c
c Don't know ahead of time how much missing data, so:
c   1) Take average for whole grid, apply to missing points (create
c        cscnew grid).
c   2) Divide grid into 4 boxes. Take average from original grid 
c        in each box. If there are some non-missing points in the
c        box, apply average to missing points (modify cscnew grid). 
c   3) Divide grid into 9 boxes. Take average from original grid 
c        in each box. If there are some non-missing points in the
c        box, apply average to missing points (modify cscnew grid).
c     -> Actually do this again for 36 and 100 box subdivisions.
c        Thus the cscnew grid will contain the average values from
c        the smallest boxes (from the 100 box subdivision). The
c        missing grid points receive this value while the non-missing
c        grid points get overwritten with the original csc values.
c   4) Apply original grid non-missing values to cscnew grid.
c
       do nrun=1,nruns
c 
       do nn=1,nbtot(nrun)
         total=0.
         sum=0.
         do i=iimin(nrun,nn),iimax(nrun,nn)
         do j=jjmin(nrun,nn),jjmax(nrun,nn)
           if(csctot(i,j).le.1.1)then
             sum=sum+csctot(i,j)
             total=total+1.
           endif
         enddo
         enddo
         if(total.eq.0.)then
c          print *,'all values missing for this box'
c          print *,'nrun,nbox=',nrun,nn
c          if(nrun.eq.1)then
c             print *,'ALL VALUES MISSING - CANNOT CONTINUE'
c             stop
c          endif
          go to 50
       else
          avrg=sum/total
c          print *,'nrun,nbox=',nrun,nn
c          print *,'box average=',avrg,' npts=',total
          do i=iimin(nrun,nn),iimax(nrun,nn)
          do j=jjmin(nrun,nn),jjmax(nrun,nn)
             cscnew(i,j)=avrg
          enddo
          enddo
       endif
 50    continue
c
       enddo
c
c       print *,'nrun,nbox=',nrun,nn
c       do j=jmax,1,-1
c       write(6,75)j,(cscnew(i,j),i=1,imax)
c       enddo
c 75    format(1x,i2,61f2.1)
c
       enddo
c
c now overwrite non-missing points
c
       do i=1,imax
       do j=1,jmax
         if(csctot(i,j).le.1.1)cscnew(i,j)=csctot(i,j)
       enddo
       enddo
c
        return
        end
c
c ***********************************************************
c
         subroutine setboxes(imax,jmax,iimin,iimax,jjmin,jjmax,
     1    nbtot,nruns,nboxes)   
c
c set up boundaries for boxes to do averaging over
c   
         dimension iimin(nruns,nboxes),
     1       iimax(nruns,nboxes),jjmin(nruns,nboxes),
     2       jjmax(nruns,nboxes),nbtot(nruns)
c
         do nn=1,nruns
           idir=sqrt(real(nbtot(nn)))
           inc=imax/idir+0.99
c           print *,'nn,idir,inc=',nn,idir,inc
           nbox=1
           do i=1,idir
             do j=1,idir
               iimin(nn,nbox)=1+(i-1)*inc
               iimax(nn,nbox)=min(imax,1+i*inc)
               jjmin(nn,nbox)=1+(j-1)*inc
               jjmax(nn,nbox)=min(jmax,1+j*inc)
               nbox=nbox+1
             enddo
           enddo
c           do nbox=1,nbtot(nn)
c             write(6,10)nn,nbox,inc,iimin(nn,nbox),iimax(nn,nbox),
c    1        jjmin(nn,nbox),jjmax(nn,nbox)
c10          format(1x,'nn,nbox,inc=',3i3,'imin,max,jmin,max=',4i4)
c           enddo
         enddo
c
         return
         end
