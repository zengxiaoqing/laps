c
c
       subroutine read_metar(data_path,nsta,maxsta,
     &                      stanam,lat,long,elev,
     &                      timeobs,rptype,atype,
     &                      cvr,ht,vis,weather,
     &                      slp,t,tt,td,ttd,dd,ff,ffg,alt,
     &                      dpchar,dp,precip1,status)
c
c.......................................................................
c
c     Routine to read the the netcdf metar data.  Based on a routine
c     by T. Smith, FSL.
c
c     Original:  2 Oct 1996  Peter A. Stamus, FSL
c     Changes:   9 Jun 1997  Changes for new Metar CDL
c
c.......................................................................
c
       include 'netcdf.inc'
c      
       integer*4 lenstr
       parameter (lenstr=MAXNCNAM)
c
c.....  Variables for the data.
c
       real*8 timeobs(maxsta)
       real*4 lat(maxsta), long(maxsta), elev(maxsta)
c
       real*4 t(maxsta), td(maxsta), tt(maxsta), ttd(maxsta)
       real*4 dd(maxsta), ff(maxsta), ffg(maxsta)
       real*4 slp(maxsta), alt(maxsta), vis(maxsta)
       real*4 ht(6,maxsta), precip1(maxsta), dp(maxsta)
c
       integer*2 dpchar(maxsta)
       character cvr(6,maxsta)*8, weather(maxsta)*25
       character stanam(maxsta)*5, rptype(maxsta)*6, atype(maxsta)*6
c
c.....  NETCDF variables
c
       integer*4 cdfid, rcode, recdim
       integer*4 ndims, natts, nvars
       integer*4 dimsiz(MAXNCDIM)
       integer*4 vdtype, adtype(MAXNCATT)
       integer*4 attlen(MAXNCATT), namlen
       integer*4 nvdims, vdims(MAXVDIMS), nvatts
       character*(MAXNCNAM) dimnam(MAXNCDIM),varnam,
     &                      attnam(MAXNCATT),
     &                      string(MAXNCATT)
       integer*4 varid
       real*8 value(MAXNCATT)
c
       integer*4 START(MAXVDIMS),COUNT(MAXVDIMS)
       integer*4 status, cnt
       character data_path*80 
       integer min_stations minutes_to_wait_for_metars
       parameter(minutes_to_wait_for_metars=10)
       parameter(min_stations=1000)
       data cnt/0/
c
c.....  Start the program
c      call ncpopt(0)
c       
c.....  Open the desired netCDF data file
c
       do while(nsta.lt.min_stations.and.
     +           cnt.lt.minutes_to_wait_for_metars)

          cdfid = ncopn(data_path,ncnowrit,rcode)
          if(rcode.ne.0 .or. cdfid.lt.0) then
             status = rcode
             if(status .eq. 0) status = 0
             return
          else
             status = 1
          endif
c 
c.....  Find out the number of dimensions, variables, attributes
c.....     and record dimensions
c
          call ncinq(cdfid,ndims,nvars,natts,recdim,rcode)

c 
c.....  Inquire about dimension names and sizes
c
          do k=1,ndims
             call ncdinq(cdfid,k,dimnam(k),dimsiz(k),rcode)
             if(dimnam(k).eq.'recNum') nsta=dimsiz(k)
          enddo                 !k

          print *,'Total records found: ',nsta
          if(nsta.lt.min_stations.and.
     +        cnt.lt.minutes_to_wait_for_metars-1) then
             print*,'Waiting for more stations'
             call ncclos(cdfid,rcode)
             call waiting_c(60)
             cnt = cnt+1
          endif
       enddo
       if(nsta.le.0) then
         print*,'Error, No stations in file ',data_path
         stop 'read_metar'
       else if(nsta.lt.min_stations) then
          print*,'Proceeding with less than ',min_stations,' stations'
       endif
       



c.....  Inquire about variables- datatype, name, number of dimensions,
c.....     dimensions, number of attributes
c
       do 20 k=1,nvars
         call ncvinq(cdfid,k,varnam,vdtype,nvdims,vdims,
     &               nvatts,rcode)
c      
c.....  Get the attributes' name, length, and value for each variable
c
         do 25 j=1,nvatts       
           call ncanam (cdfid,k,j,attnam(j),rcode)        
           call ncainq (cdfid,k,attnam(j),adtype(j),attlen(j),rcode)
c
c.....  Make sure you call the right routine to get the attribute string/value
c
           If(adtype(j).ne.NCCHAR) then
             call ncagt (cdfid,k,attnam(j),value(j),rcode)
           Else
             call ncagtc (cdfid,k,attnam(j),string(j),lenstr,rcode)
           Endif
c           
c.....  Find out which index matches the desired variable
c
           If(attnam(j)(1:9).eq.'long_name') then  !loop over all variables
c
c.....  Check for the desired variable names, get the proper id number for
c.....  that variable and the dimensions, then get the data.
c
c
             If(string(j)(1:attlen(j)).eq.
     &                     'alphanumeric station identification') then
               varid=ncvid(cdfid,varnam,rcode)
               namlen=1
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
                 namlen=namlen*dimsiz(vdims(i))
               enddo !I
               call ncvgtc(cdfid,varid,start,count,stanam,namlen,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                          'latitude') then
               varid=ncvid(cdfid,varnam,rcode)  !get the id # for this var.
               Do I=1,nvdims  !get dimentions to use for hyperslab of data
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,lat,rcode) !get the data
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                          'longitude') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,long,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                          'elevation') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,elev,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                          'time of observation') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,timeobs,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                                        'report type') then
               varid=ncvid(cdfid,varnam,rcode)
               namlen=1
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
                 namlen=namlen*dimsiz(vdims(i))
               enddo !I
               call ncvgtc(cdfid,varid,start,count,rptype,namlen,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                              'automated station type') then
               varid=ncvid(cdfid,varnam,rcode)
               namlen=1
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
                 namlen=namlen*dimsiz(vdims(i))
               enddo !I
               call ncvgtc(cdfid,varid,start,count,atype,namlen,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq. 'sky cover') then
               varid=ncvid(cdfid,varnam,rcode)
               namlen=1
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
                 namlen=namlen*dimsiz(vdims(i))
               enddo !I
               call ncvgtc(cdfid,varid,start,count,cvr,namlen,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                               'sky cover layer base') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,ht,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                                         'visibility') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,vis,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                                    'present weather') then
               varid=ncvid(cdfid,varnam,rcode)
               namlen=1
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
                 namlen=namlen*dimsiz(vdims(i))
               enddo !I
               call ncvgtc(cdfid,
     &                     varid,start,count,weather,namlen,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                                 'sea level pressure') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,slp,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                                        'temperature') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,t,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)).eq.
     &       'temperature from tenths of a degree Celsius') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,tt,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq. 'dewpoint') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,td,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)).eq.
     &       'dewpoint from tenths of a degree Celsius') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,ttd,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                                   'wind direction') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,dd,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                                       'wind speed') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,ff,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                                       'wind gust') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,ffg,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                                'altimeter setting') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,alt,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                                '1 hour precipitation') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,precip1,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                     'character of pressure change') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,dpchar,rcode)
c
c
             ElseIf(string(j)(1:attlen(j)) .eq.
     &                           '3 hour pressure change') then
               varid=ncvid(cdfid,varnam,rcode)
               Do I=1,nvdims
                 START(I)=1
                 COUNT(I)=dimsiz(vdims(i))
               enddo !I
               call ncvgt(cdfid,varid,start,count,dp,rcode)
c
c
             EndIf
           Endif
c           
25       continue
20     continue
c
c.....  Close the netCDF file
c
       call ncclos (cdfid,rcode)
c
c.....  Set the return status to be the return code
c
       if(rcode .eq. 0) then
          status = 1
       else
          status = -1
       endif
c        
       return
       end
