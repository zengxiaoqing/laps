       subroutine GSI_radarbufr
!***********************************************************************
! Program name : GSI_radarbufr(GSI_radarbufr.f)
!
! Description: To read in LAPS's level II polar radar data in NetCDF format,
!              and write out NCEP's GSI level III radar BUFR data.
!
! Called Function:
!    get_systime(${LAPS_SRC_ROOT}/src/lib/get_systime.f)
!    get_radarbufr_parms
!    get_rfile_time
!    write_radvel_bufr
!    write_radref_bufr
!
! Date :
!   Original     -- May. 04, 2006 (Shiow-Ming, Deng)
!***********************************************************************

       parameter (MAXRADARS=99,MAXTIMES=999)
       character*9 a9_time
       character*150 path_lvl2(MAXRADARS),path_table,path_output
       character*150 file_table,file_output,path
       character*9 cdf_times(MAXTIMES),cdfdate
       integer num_elevs(MAXTIMES),ielev(30,MAXTIMES)
       integer jelev(MAXTIMES)

!-----------------------------------------------------------------------
!c  To get LAPS's system time.

       call get_systime(i4time,a9_time,istatus)
       if( istatus.ne.1 )go to 901

!-----------------------------------------------------------------------
!c  To get paths of radar bufr table, LAPS netcdf level 2 radar data.

       call get_radarbufr_parms(nradars,path_table,path_output   &
                               ,path_lvl2,itimeb,itimea,istatus)
       if( istatus.ne.1 )go to 902
       print*,'Get paths of radar bufr table, and '
       print*,'    the paths of LAPS netcdf level 2 radar data.'
       print*,'Total number of radar: ',nradars

!-----------------------------------------------------------------------
!c  To open radar BUFR table file: "radar.bufrtable"

       lent=index(path_table,' ')
       file_table=path_table(1:lent-1)//'/radar.bufrtable'
       open(12,file=file_table,err=903)

!-----------------------------------------------------------------------
!c  To get output radar level III BUFR data file name: "radarbufr"

       len=index(path_output,' ')
       file_output(1:len-1)=path_output(1:len-1)
       file_output(len:len+9)='/radarbufr'
       file_output(len+10:lev+10)=char(0)

       icheck=0
       do 100 i=1,nradars

!-----------------------------------------------------------------------
!c  To get LAPS radar polar NetCDF data time.

          print*,' '
          path=path_lvl2(i)
          print*,'The input radar path: ',path
          call get_rfile_time(a9_time,itimeb,itimea,path,num_times  &
                             ,cdf_times,num_elevs,ielev,istatus)
          if( istatus.ne.1 )go to 100
          icheck=icheck+1

!-----------------------------------------------------------------------
!c  To open radar BUFR data.

          if( icheck.eq.1 )then
!             open(11,file=file_output,form='unformatted')
              open(11,file='radarbufr',form='unformatted')
              call openbf(11,'OUT',12)
              lunin=11
          endif

!-----------------------------------------------------------------------
!c  To write out BUFR data for radial velocity field.

          len=index(path,' ')
          len=len-1
          do k=1,num_times
             cdfdate(1:9)=cdf_times(k)(1:9)
             numelev=num_elevs(k)
             do j=1,numelev
                jelev(j)=ielev(j,k) 
             enddo
             call write_radvel_bufr(lunin,path,len,a9_time,cdfdate  &
                                   ,numelev,jelev,istatus)
          enddo
 100   continue
       if( icheck.eq.0 )go to 904

!-----------------------------------------------------------------------
!c  To force BUFR data to be closed.

!      call writsa(-11,ibfmsg,libf)
       call closbf(11)

       return
 901   continue
       print*,'Cannot get system time.'
       return
 902   continue
       print*,'Cannot read in parameter of namelist file: '
       print*,'Please check file: static/radar_bufr.nl'
       return
 903   continue
       print*,'Cannot open radar BUFR table.'
       print*,'radar BUFR table file: ',file_table(1:lent+15)
       return
 904   continue
       print*,'No radar data can be written to BUFR.'
       return

       end

       subroutine write_radvel_bufr(lunin,directory,len,a9_time,cdfdate  &
                                   ,num_elev,ielev,istatus)
!***********************************************************************
! Function Name : write_radvel_bufr
!
! Usage : call write_radvel_bufr(lunin,directory,len,a9_time,cdfdate
!        1                      ,num_elev,ielev,istatus)
!
! Description      : To write out level III radar's radial velocity data
!                    in NCEP BUFR format.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     lunin      integer     absolute value is Fortran logical unit
!                                 number for BUFR file.
!    I     directory  C*150       input directory of NetCDF level II radar data.
!    I     len        integer     the character's length of directory.
!    I     a9_time    C*9         system time. (= yyJJJhhmm )
!    I     cdfdate    C*9         the date: "yyJJJhhmm".
!                                 yy: Year, JJJ: Julian day, hh: Hour, mm: Minute.
!    I     num_elev   integer     number of elevation angle.
!    I ielev(num_elev) int array  elevation angle ID.
!    O     istatus    integer     work or no work message,
!                                 =0, the subroutine doesnot work.
!                                 =1, the subroutine does work.
!
! Sub./Fun. Called : 
!    verjulian, 
!    read_radar_polar
!    int_ppi,
!    openmb,  ufbint,  writsa  (NCEP BUFR API)
!
! Date :
!   Original     -- Apr. 12, 2006 (Shiow-Ming, Deng)
!***********************************************************************
       parameter( ngates=1600,nbeams=420,ix=401 )
       parameter( mxbf=16000 )

       character directory*150,a9_time*9,cdfdate*9
       integer ielev(num_elev)

       character radarfile*150,ca2*2
       character radarName*5
       real azim(nbeams),elev(nbeams)
       integer*8 itim(nbeams)

       integer iz(ngates*nbeams),iv(ngates*nbeams),iw(ngates*nbeams)
       integer it(ngates*nbeams) 
       integer*8 istarttime,isec

       real xd(ix)
       integer ndatv(ix,ix),ndatw(ix,ix),ndatt(ix,ix)

       real*8 hdr(10),rwnd(7,ngates)
       character*4 sstn
       equivalence (sstn,hdr(1))
       character*8 subset
       character cout*49,time*12
       integer ibfmsg(mxbf/4)

       istatus=0
       cout(1:49)='RPID CLAT CLON SELV ANEL YEAR MNTH DAYS HOUR MINU'
       i_missing=32767
       pi=acos(-1.)
       earth=6371.25
       rmax=200.
       nx=ix
       do i=1,nx
          xd(i)=-rmax+(i-1.)
       enddo

!-----------------------------------------------------------------------
!c  To get NetCDF file names.

       if( directory(len:len).eq.'/' )then
           len=len-1
       endif 

!-----------------------------------------------------------------------
!c  To get idate from year,month,day,and hour. 

       read(a9_time,'(i2,i3,2i2)')iyear,julian,ihour,minute
       iyear=iyear+2000
       call verjulian(iyear,julian,month,iday,ier)
       if( ier.ne.0 )return
       idate=iyear*1000000+month*10000+iday*100+ihour
       read(cdfdate,'(i2,i3,2i2)')iyear,julian,ihour,minute
       iyear=iyear+2000
       call verjulian(iyear,julian,month,iday,ier)
       if( ier.ne.0 )return
       write(time,'(i4,4i2)')iyear,month,iday,ihour,minute

!-----------------------------------------------------------------------
!c  To open BUFR message.

       subset(1:8)='NC006001'
       call openmb(lunin,subset,idate)

       icheck=0
       do 100 k=1,num_elev

!-----------------------------------------------------------------------
!c  To read in radar's data from NetCDF polor-coordinated file.

          write(ca2,'(i2)')ielev(k)
          if( k.lt.10 )ca2(1:1)='0'
          radarfile(1:len)=directory(1:len)
          radarfile(len+1:len+1)='/'
          radarfile(len+2:len+10)=cdfdate(1:9)
          radarfile(len+11:len+15)='_elev'
          radarfile(len+16:len+17)=ca2(1:2)
          radarfile(len+18:len+18)=char(0)
          call read_radar_polar(radarfile,i_missing           &
                   ,radarName,slat,slon,salt,istarttime       &
                   ,inumber,elevation,nobeam,nogateZ,nogateV  &
                   ,rangeZ,rangeV,spaceZ,spaceV,snyquist      &
                   ,azim,elev,itim,iz,iv,iw,istat)
          if( istat.eq.0 )go to 100

!-----------------------------------------------------------------------
!c  To interpolate data for PPI scan.

          range=1000.*rangeV
          space=1000.*spaceV
          call int_ppi(nogateV,nobeam,iv,i_missing,rmax,range,space &
                      ,elevation,azim,nx,xd,ndatv)
          ickw=0
          do i=1,nobeam*nogateV
             if( iw(i).ne.i_missing )ickw=1
          enddo
          if( ickw.eq.1 )then
              call int_ppi(nogateV,nobeam,iw,i_missing,rmax,range,space  &
                      ,elevation,azim,nx,xd,ndatw)
          else
              do j=1,nx
              do i=1,nx
                 ndatw(i,j)=i_missing
                 if( ndatv(i,j).ne.i_missing )ndatw(i,j)=250
              enddo
              enddo
          endif
          ij=0
          do j=1,nobeam
          do i=1,nogateV
             ij=ij+1
             it(ij)=(itim(j)-istarttime)*10/6.
          enddo
          enddo
          call int_ppi(nogateV,nobeam,it,i_missing,rmax,range,space  &
                      ,elevation,azim,nx,xd,ndatt)

          idy=istarttime/86400
          isec=istarttime-86400*idy
          ii=isec/60
          isecst=isec-60*ii
          ih=ii/60
          imst=ii-60*ih 

          icheck=1
          sstn(1:4)=radarName(1:4) 
          hdr(2)=slat
          hdr(3)=slon
          hdr(4)=salt
          hdr(5)=elevation
          hdr(6)=iyear
          hdr(7)=month
          hdr(8)=iday
          hdr(9)=ih
          hdr(10)=imst

          ct=tan(elevation*pi/180.)
          nn=0
          nnn=0
          do 50 j=1,nx
          do 50 i=1,nx
             if( ndatv(i,j).eq.i_missing )go to 50
             x=xd(i)
             y=xd(j)
             r=sqrt(x*x+y*y)
             if( r.lt.0.1 )go to 50
             plat=y/earth*180./pi
             platm=slat+0.5*plat
             earth1=earth*cos(platm*pi/180.)
             ht=1000*r*ct+r**2/17.
             th=acos( y/r )*180./pi
             if( x.lt.0. )th=360.-th
             nn=nn+1
             nnn=nnn+1
             rwnd(1,nn)=0.01*ndatt(i,j)
             rwnd(2,nn)=slat+plat
             rwnd(3,nn)=slon+x/earth1*180./pi
             rwnd(4,nn)=salt+ht
             rwnd(5,nn)=0.01*ndatv(i,j)
             rwnd(6,nn)=th
             rwnd(7,nn)=0.01*ndatw(i,j)
             if( nn.eq.810 )then
                 call ufbint(lunin,hdr,10,1,levs,cout)
                 call ufbint(lunin,rwnd,7,nn,levs                  &
                        ,'STDM SUPLAT SUPLON HEIT RWND RWAZ RSTD')
                 call writsa(lunin,ibfmsg,ibf)
                 nn=0
             endif
 50       continue
          if( nn.gt.4 )then
              call ufbint(lunin,hdr,10,1,levs,cout)
              call ufbint(lunin,rwnd,7,nn,levs                    &
                        ,'STDM SUPLAT SUPLON HEIT RWND RWAZ RSTD')
              call writsa(lunin,ibfmsg,ibf)
          endif
          if( nnn.gt.4 )then
              print*,'output time: ',time(1:12) &
                    ,'   sweeep number of input data: ',ca2(1:2) &
                    ,'   total number of output data: ',nnn
          endif
 100   continue
       if( icheck.eq.0 )return

       istatus=1
       return
       end

       subroutine int_ppi(ng,nb,idat,miss,rmax,range,space,ang,azim  &
                         ,nx,xd,ndat)
!***********************************************************************
! Function Name: int_ppi
!
! Usage :
!    call int_ppi(ng,nb,idat,miss,rmax,range,space,ang,azim
!   1            ,nx,xd,ndat)
!
! Description      : Using the nearest point method intepolate  
!                    the PPI dBZ data  to grid point of  RADARs
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     ng         integer     gate number in one ray.
!    I     nb         integer     total ray number of input sweep scan layer.
!    I    idat(ng,nb) int array   input data.
!    I     miss       integer     missing or bad value.
!    I     rmax       real        the max. radius in kilometers.
!    I     range      real        range to first gate. (meters)
!    I     space      real        space between gates. (meters)
!    I     ang        real        sweep angle of input sweep scan layer,
!                                 (degrees)
!    I     azim(nb)   real array  azimuthal angle of a ray in one sweep.
!                                 (degrees)
!    I     nx         integer     the dimension of domain.
!    I     xd(nx)     real array  the x- or y-grid in unit of km.
!    O    ndat(nx,nx) int array   the output data.
!***********************************************************************

       dimension idat(ng,nb)
       dimension azim(nb),xd(nx)
       dimension ndat(nx,nx)
       dimension iwk(720),wk(720)

       rspace=0.001*space
       rmin=0.001*range
       rmin1=max(rmin,0.1)
       pi=acos(-1.)
       rcos=cos(ang*pi/180.)
       angg=max(ang,0.1)
       rrmax=0.001*(range+(ng-1)*space)

       do i=1,nb
          wk(i)=azim(i)
       enddo
       wk(nb+1)=azim(1)
       do 10 i=1,720
          iwk(i)=-1
          th=0.5*i
          do k=1,nb
             diff=abs( wk(k+1)-wk(k) )
             diff1=wk(k+1)-th
             diff2=wk(k)-th
             if( diff1*diff2.le.0. )then
                 if( diff.lt.10. )then
                     if( abs(diff1).lt.abs(diff2) )then
                         iwk(i)=k+1
                         if( iwk(i).gt.nb )iwk(i)=1
                     else
                         iwk(i)=k
                     endif
                     go to 10 
                 endif
             else
                 if( diff.gt.355. )then
                     dd1=abs(diff1)
                     dd2=abs(diff2)
                     if( dd1.gt.dd2 )then
                         dd1=abs(360.-dd1)
                     else
                         dd2=abs(360.-dd2)
                     endif
                     if( dd1.lt.dd2 )then
                         iwk(i)=k
                     else
                         iwk(i)=k+1
                         if( iwk(i).gt.nb )iwk(i)=1
                     endif
                     go to 10
                 endif
             endif
          enddo
 10    continue

       do 100 j=1,nx
       do 100 i=1,nx
          ndat(i,j)=miss
          rr=sqrt( xd(i)**2+xd(j)**2 )
          rr1=rr/rcos
          if( (rr1.lt.rmin1).or.(rr1.ge.rrmax) )go to 100
          irad=(rr1-rmin)/rspace+1
          if( irad.lt.1 )go to 100
          th=acos( xd(j)/rr )*180./pi
          if( xd(i).lt.0. )th=360.-th
          ith=2.*th
          if( ith.lt.1 )then
              ith=720
              if( th.gt.0.25 )ith=1
          endif 
          iazi=iwk(ith)
          if( iazi.gt.0 )ndat(i,j)=idat(irad,iazi)
 100   continue

       return
       end

       subroutine verjulian(iyear,julian,month,iday,ier)
!***********************************************************************
! Function Name : verjulian
!
! Usage : call verjulian(iyear,julian,month,iday,ier)
!
! Description      : To get month, day from year, and julian day.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     iyear      integer     year, ex: 1998
!    I     julian     integer     the Julian day. 
!                                 =1, for Jan. 1st.
!    O     month      integer     month.
!    O     iday       integer     day
!    O     ier        integer     =0, success.
!                                 =1, failure.
!
! Sub./Fun. Called : none
!
! Date :
!   Original     -- Apr. 05, 2006 (Shiow-Ming, Deng)
!***********************************************************************

       dimension idate1(12),idate2(12)
       data idate1/31,59,90,120,151,181,212,243,273,304,334,365/
       data idate2/31,60,91,121,152,182,213,244,274,305,335,366/

       ier=0
       month=-1
       iday=-1
       if( (iyear.le.0).or.(julian.lt.1).or.(julian.gt.366) )then
            ier=1
            return
       endif

       ind=0
       i1=iyear/4
       i4=iyear-4*i1
       if( i4.eq.0 )ind=1
       i1=iyear/400
       i2=iyear-400*i1
       if( i2.eq.0 )ind=1
       i1=iyear/100
       i3=iyear-100*i1
       if( (i3.eq.0).and.(i2.ne.0) )ind=0
       if( (ind.eq.0).and.(julian.gt.365) )then
           ier=1
           return
       endif

       if( ind.eq.0 )then
           do i=1,12
              if( julian.le.idate1(i) )then
                  month=i
                  if( month.eq.1 )then
                      iday=julian
                  else
                      iday=julian-idate1(i-1)
                  endif
                  return
              endif
           enddo
       else
           do i=1,12
              if( julian.le.idate2(i) )then
                  month=i
                  if( month.eq.1 )then
                      iday=julian
                  else
                      iday=julian-idate2(i-1)
                  endif
                  return
              endif
           enddo
       endif
       return
       end

       subroutine get_rfile_time(a9_time,itimeb,itimea,path,num_times  &
                                ,cdf_times,num_elevs,ielev,istatus)
!***********************************************************************
! Subroutine/Function : get_rfile_time
!
! Usage :
!      call get_rfile_time(a9_time,itimeb,itimea,path,num_times
!     1                   ,cdf_times,num_elevs,ielev,istatus)
!
! Description      : To get NetCDF level II radar file names.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     a9_time    C*9         system time. (= yyJJJhhmm )
!    I     itimeb     integer     time window before.
!    I     itimea     integer     time window after. 
!                                 Set these variables for time filtering the data.
!                                 Units are seconds.  For example, if you want data with an
!                                 observation time from 60 min before to 30 min after the analysis
!                                 time to be included in the BUFR file, use 3600 and 1800 for
!                                 itimeb and itimea, respectively.
!    I     path       C*150       full path to directory containing a set of input
!                                 radar tilts/volumes.
!    O     num_times  integer     number of time during system_time - itimeb and
!                                                       system_time + itimea.
!    O  cdf_times(999) C*9 array  time to write out BUFR table. ( = yyJJJhhmm )
!    O  num_elevs(999) int array  number of elevation in a volue w.r.t. cdf_time.
!    O  ielevs(30,999) int array  elevation ID in a volue w.r.t. cdf_time.
!    O     istatus    integer     the work or not message.
!                                 =1, read in data.
!                                 =0, didnot read any data.
!
! Modules Called :
!   getfilenames(getiofile.c)
!   difference_time
!
! Date :
!   Original     -- Apr. 12, 2006 (Shiow-Ming, Deng)
!***********************************************************************
       parameter (MAXTIMES=999)
       character a9_time*9,path*150
       character*9 cdf_times(MAXTIMES)
       integer num_elevs(MAXTIMES),ielev(30,MAXTIMES)

       character directory*150
       character*900000 cdat
       character*16 file(9900),ftemp
       integer iwork(30)
       integer isyear,isjulian,ishour,isminute
       integer iyear,julian,ihour,minute,jsec,ier

       istatus=0
       num_times=0
       do j=1,MAXTIMES
          cdf_times(j)=' '
          num_elevs(j)=0
          do i=1,30
             ielev(i,j)=0
          enddo
       enddo

       read(a9_time,'(i2,i3,2i2)',err=901)isyear,isjulian,ishour  &
                                         ,isminute
       len=index(path,' ')
       if( len.lt.1 )go to 902
       directory(1:len-1)=path(1:len-1)
       directory(len:len)=char(0)

       isize=9000
       do i=1,isize
          cdat(i:i)=char(0)
       enddo 
       call getfilenames(directory,isize,cdat,ier)
       if( ier.ne.0 )go to 903

       numfile=0
       do i=isize,1,-1
          if( cdat(i:i).ne.char(0) )then
              ii=i
              go to 10 
          endif
       enddo
       return
 10    continue
       i1=1
       j=1
 20    continue
       j=j+1
       if( j.gt.ii )go to 30
       if( cdat(j:j).eq.' ' )then
           i2=j-1
           jj=i2-i1+1
           if( (jj.eq.16).and.(cdat(i1+9:i1+13).eq.'_elev') )then
               numfile=numfile+1
               file(numfile)(1:16)=cdat(i1:i2)
           endif
           i1=j+1
       else
           if( j.eq.ii )then
               i2=ii
               jj=i2-i1+1
               if( (jj.eq.16).and.(cdat(i1+9:i1+13).eq.'_elev') )then
                    numfile=numfile+1
                    file(numfile)(1:16)=cdat(i1:i2)
               endif
               go to 30
           endif 
       endif
       go to 20
 30    continue
       if( numfile.eq.0 )go to 904

       do 50 k=1,numfile
          ftemp(1:16)=file(k)(1:16)
          read(ftemp,'(i2,i3,2i2)',err=50)iyear,julian,ihour,minute
          call difference_time(isyear,isjulian,ishour,isminute        &
                              ,iyear,julian,ihour,minute,jsec,ier)
          if( ier.ne.0 )go to 50
          if( jsec.gt.itimea )go to 50
          if( abs(jsec).gt.itimeb )go to 50
          if( num_times.eq.0 )then
              num_times=num_times+1
              cdf_times(num_times)(1:9)=ftemp(1:9)
              go to 50
          endif
          do i=1,num_times
             if( cdf_times(i)(1:9).eq.ftemp(1:9) )go to 50
          enddo
          num_times=num_times+1
          cdf_times(num_times)(1:9)=ftemp(1:9)
 50    continue
       if( num_times.eq.0 )go to 905

       do 100 k=1,num_times
          do i=1,30
             iwork(i)=0
          enddo
          id=0
          do 60 i=1,numfile
             ftemp(1:16)=file(i)(1:16)
             if( cdf_times(k)(1:9).eq.ftemp(1:9) )then
                 read(ftemp,'(14x,i2)',err=60)iva
                 num_elevs(k)=num_elevs(k)+1
                 id=id+1
                 iwork(id)=iva
             endif
             if( num_elevs(k).ge.30 )go to 70
 60       continue
 70       continue
          if( id.eq.0 )go to 100
          if( id.eq.1 )then
              ielev(1,k)=iwork(1)
              go to 100              
          endif
          do i=1,id-1
             do j=i+1,id
                if( iwork(i).gt.iwork(j) )then
                    iw1=iwork(i)
                    iwork(i)=iwork(j)
                    iwork(j)=iw1
                endif
             enddo
          enddo
          do i=1,id
             ielev(i,k)=iwork(i)
          enddo
 100   continue

       istatus=1
       return
 901   continue
       print*,'Cannot read system time a9_time(yyjjjhhmm): '
       print*,'a9_time: ',a9_time(1:9)
       return
 902   continue
       print*,'The path of radar level II file does not exit.'
       print*,'path: ',path(1:150)
       return
 903   continue
       print*,'Cannot read in the file names in directory: '  &
             ,path(1:len)
       return
 904   continue
       print*,'No any radar level II file can be found.'
       print*,'path: ',path(1:len)
       return
 905   continue
       print*,'No any time data can be found.'
       print*,'path: ',path(1:len)
       print*,'system time: ',a9_time(1:9)
       print*,'itimeb: ',itimeb,' itimea: ',itimea
       return 

       end

       subroutine difference_time(isyear,isjulian,ishour,isminute    &
                                 ,iyear,julian,ihour,minute,jsec,ier)
!***********************************************************************
! Subroutine/Function : difference_time
!
! Usage :
!      call difference_time(isyear,isjulian,ishour,isminute
!     1                    ,iyear,julian,ihour,minute,jsec,ier)
!
! Description      : To compute diference second between
!                    target time and system time.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     isyear     integer     system year.
!    I     isjulian   integer     system julian day.
!    I     ishour     integer     system hour.
!    I     isminute   integer     system minute.
!    I     iyear      integer     target year.
!    I     julian     integer     target julian day.
!    I     ihour      integer     target hour.
!    I     minute     integer     target minute.
!    O     jsec       integer     second between target time and system time.
!                                 =999999 for too large second. (ier=1)
!    O     ier        integer     error message.
!                                 =1, failure.
!                                 =0, success.
!
! Modules Called : none
!
! Date :
!   Original     -- Apr. 12, 2006 (Shiow-Ming, Deng)
!***********************************************************************

       integer isyear,isjulian,ishour,isminute
       integer iyear,julian,ihour,minute,jsec,ier

       ier=1
       jsec=999999
       jyr=iyear-isyear
       if( abs(jyr).gt.1 )return
       if( jyr.eq.0 )then
           jday=julian-isjulian
       else
           if( jyr.eq.1 )then
               if( julian.ge.isjulian )return
               ind=0
               i1=isyear/4
               i4=isyear-4*i1
               if( i4.eq.0 )ind=1
               i1=isyear/400
               i2=isyear-400*i1
               if( i2.eq.0 )ind=1
               i1=isyear/100
               i3=isyear-100*i1
               if( (i3.eq.0).and.(i2.ne.0) )ind=0
               if( ind.eq.0 )then
                   jday=julian+365-isjulian
               else
                   jday=julian+366-isjulian
               endif
           else
               if( isjulian.ge.julian )return
               ind=0
               i1=iyear/4
               i4=iyear-4*i1
               if( i4.eq.0 )ind=1
               i1=iyear/400
               i2=iyear-400*i1
               if( i2.eq.0 )ind=1
               i1=iyear/100
               i3=iyear-100*i1
               if( (i3.eq.0).and.(i2.ne.0) )ind=0
               if( ind.eq.0 )then
                   jday=julian-365-isjulian
               else
                   jday=julian-366-isjulian
               endif
           endif
       endif
       jsec=86400*jday+3600*(ihour-ishour)+60*(minute-isminute)
       ier=0
       return
       end

       subroutine get_radarbufr_parms(nradars,path_table,path_output   &
                                     ,path_lvl2,itimeb,itimea,istatus)
!***********************************************************************
! Subroutine/Function : get_radarbufr_parms
!
! Usage :
!      call get_radarbufr_parms(nradars,path_table,path_output
!     1                        ,path_lvl2,itimeb,itimea,istatus)
!
! Description      : To get parameters of namelist file: 
!                    'static/radar_bufr.nl' for radar GSI BUFR data.
!
! Arguments :
!  I/O/W   name,      type,       description
!    O     nradars    integer     number of radars (and/or radar types) to loop through
!                                 and process.
!    O     path_table C*150       directory for input BUFR table.
!    O    path_output C*150       directory for output BUFR data.
!    O  path_lvl2(99) C*150 array full path to each directory containing a set of input
!                                 radar tilts/volumes.
!    O     itimeb     integer     time window before.
!    O     itimea     integer     time window after. 
!                                 Set these variables for time filtering the data.
!                                 Units are seconds.  For example, if you want data with an
!                                 observation time from 60 min before to 30 min after the analysis
!                                 time to be included in the BUFR file, use 3600 and 1800 for
!                                 itimeb and itimea, respectively.
!    O     istatus    integer     the work or not message.
!                                 =1, read in data.
!                                 =0, didnot read any data.
!
! Modules Called :
!   get_directory($LAPS_SRC_ROOT/src/lib/get_directory.f)
!
! Date :
!   Original     -- May. 05, 2006 (Shiow-Ming, Deng)
!***********************************************************************
       integer MAXRADARS
       parameter (MAXRADARS=99)
       character*150 path_table,path_output,path_lvl2(MAXRADARS)

       character*150 path_to_radar_a(MAXRADARS),path_to_vrc_nl
       character*4 laps_radar_ext_a(MAXRADARS) 
       logical l_line_ref_qc,l_hybrid_first_gate,l_unfold
       character*3 ext

       namelist /remap_nl/ n_radars_remap,max_times,path_to_radar_a      &
                          ,laps_radar_ext_a,path_to_vrc_nl               &
                          ,ref_min,min_ref_samples,min_vel_samples,dgr   &
                          ,abs_vel_min,l_line_ref_qc,l_hybrid_first_gate &
                          ,l_unfold 
       character*150 dir,filename

       nradars=0 
       istatus=0

       call get_directory('nest7grid',dir,len)
       if( dir(len:len).eq.'/' )then
           path_table(1:len-1)=dir(1:len-1)
           path_table(len:len)=' '
       else
           path_table(1:len)=dir(1:len)
       endif
       itimeb=3600
       itimea=3600

       filename = dir(1:len)//'/remap.nl'
       open(1,file=filename,status='old',err=900)
       read(1,remap_nl,err=901)
       close(1)

       nradars=0
       do i=1,n_radars_remap
          ext(1:3)=laps_radar_ext_a(i)(1:3)
          read(ext,'(1x,i2)',err=10)num
          nradars=nradars+1
          if( nradars.gt.MAXRADARS )then
      print*,'stop get radar bufr data -- increase parameter MAXRADARS.'
              return
          endif
          path_lvl2(nradars)=path_to_radar_a(i)
 10       continue
       enddo
       if( nradars.le.0 )return

       call get_directory('log',dir,len)
       if( dir(len:len).eq.'/' )then
           path_output(1:len-1)=dir(1:len-1)
           path_output(len:len)=' '
       else
           path_output(1:len)=dir(1:len)
       endif

       istatus=1 
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading remap_in ',filename
       write(*,remap_nl)
       istatus = 0
       return

       end

       subroutine read_radar_polar(radarfile,i_missing          &
                      ,radarName,slat,slon,salt,istarttime      &
                      ,inumber,elevation,nobeam,nogateZ,nogateV &
                      ,rangeZ,rangeV,spaceZ,spaceV,snyquist     &
                      ,azim,elev,itim,iz,iv,iw,istatus)
!***********************************************************************
! Subroutine/Function : read_radar_polar
!
! Usage :
!      call read_radar_polar(radarfile,i_missing
!    1                ,radarName,slat,slon,salt,istarttime
!    2                ,inumber,elevation,nobeam,nogateZ,nogateV
!    3                ,rangeZ,rangeV,spaceZ,spaceV,snyquist
!    4                ,azim,elev,itim,iz,iv,iw,istatus)
!
! Description      : To read radar polor-coordinate NetCDF-formatted
!                    data for LAPS's radar importting system.
!
! Parameters:
!    ngates: the maximum gate number. (=944) 
!    nbeams: the maximum radial ray number. (=420)
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     radarfile  C*150       input radar file name.
!                                 ='yyJJJhhmm_elev01','yyJJJhhmm_elev02',...
!                                 = or 'yyJJJhhmm_elevnn'.
!                                   yy: Year, JJJ: Julian day, hh: Hour, mm: Minute,
!                                   nn: the elevation number.
!    I     i_missing  integer     integer missing value.
!    O     radarName  C*5         radar name. (='RCWF ','RCKT ',... )
!    O     slat       real        latitude of radar site. (degrees)
!    O     slon       real        longitude of radar site. (degrees)
!    O     salt       real        height of radar site. (meters)
!    O     istarttime integer*8      starting time, seconds since 1970-1-1 00:00:00.00
!    O     inumber    integer     number of constant PPI-slope.
!    O     elevation  real        elevation angle of PPI slope. (degrees)
!    O     nobeam     integer     total radial ray number of sweep scan.
!    O     nogateZ    integer     gate number in one ray for reflectivity (Z) data.
!    O     nogateV    integer     gate number in one ray for wind (V) data.
!                                 (for radial velocity or spectral width)
!    O     rangeZ     real        range to first gate for Z field. (km)
!    O     rangeV     real        range to first gate for V field. (km)
!    O     spaceZ     real        space between gates for Z field. (km)
!    O     spaceV     real        space between gates for V field. (km)
!    O     snyquist   real        Nyquist velocity. (m/s)
!    O   azim(nbeams) real array  azimuthal angle of a ray in one sweep. (degrees)
!    O   elev(nbeams) real array  elevation angle of a ray in one sweep. (degrees)
!    O   itim(nbeams) int*8 array time, seconds since 1970-1-1 00:00:00.00 
!    O iz(ngates*nbeams) int array the reflectivity data. (100*dBZ)
!    O iv(ngates*nbeams) int array the radial velocity data. (100*m/s)
!    O iw(ngates*nbeams) int array the spectral width data. (100*m/s)
!    O     istatus    integer     work or no work message,
!                                 =0, the subroutine doesnot work.
!                                 =1, the subroutine does work.
!
! Modules Called :
!   netcdf fortran API.
!
! Date :
!   Original     -- Apr. 05, 2006 (Shiow-Ming, Deng)
!***********************************************************************
       parameter( ngates=944,nbeams=420 )
       character*150 radarfile
       integer i_missing

       character radarName*5
       real slat,slon,salt
       real*8 starttime
       integer*8 istarttime
       integer inumber,nobeam,nogateZ,nogateV
       real elevation,rangeZ,rangeV,spaceZ,spaceV,snyquist 
       real azim(nbeams),elev(nbeams)
       integer*8 itim(nbeams)
       integer istatus
       integer iz(ngates*nbeams),iv(ngates*nbeams),iw(ngates*nbeams)

       integer*2 elevationNumber,numRadials,numGatesZ,numGatesV
       real*8 radialTime(nbeams)
       integer*1 zbyte(ngates*nbeams),vbyte(ngates*nbeams)
       integer*1 wbyte(ngates*nbeams) 

       include 'netcdf.inc'

       istatus=0

!-----------------------------------------------------------------------
!c  To open NetCDF file.

       ier=nf_open(radarfile,nf_nowrite,ncid)
       if( ier.ne.0 )return 

!-----------------------------------------------------------------------
!c  To get dimensions: Z_bin, V_bin, and unlimit.

       ier=nf_inq_dimid(ncid,'Z_bin',id_z_bin)
       if( ier.ne.0 )return 
       ier=nf_inq_dimlen(ncid,id_z_bin,len_z_bin)
       if( ier.ne.0 )return 
       ier=nf_inq_dimid(ncid,'V_bin',id_v_bin)
       if( ier.ne.0 )return 
       ier=nf_inq_dimlen(ncid,id_v_bin,len_v_bin)
       if( ier.ne.0 )return 
       ier=nf_inq_unlimdim(ncid,id_unlimit)
       if( ier.ne.0 )return 
       ier=nf_inq_dimlen(ncid,id_unlimit,len_unlimit)
       if( ier.ne.0 )return 

!-----------------------------------------------------------------------
!c  To get radar name, latitude, longitude, and altitude.

       ier=nf_inq_varid(ncid,'radarName',id_radarName)
       if( ier.ne.0 )return
       ier=nf_get_var_text(ncid,id_radarName,radarName)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'siteLat',id_siteLat)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_siteLat,slat)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'siteLon',id_siteLon)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_siteLon,slon)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'siteAlt',id_siteAlt)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_siteAlt,salt)
       if( ier.ne.0 )return

!-----------------------------------------------------------------------
!c  To get elevationNumber, elevationAngle, numRadials, esStartTime,
!c         firstGateRangeZ, firstGateRangeV, gateSizeZ, gateSizeV,
!c         numGatesZ, numGatesV, and nyquist

       ier=nf_inq_varid(ncid,'elevationNumber',id_elevationNumber)
       if( ier.ne.0 )return
       ier=nf_get_var_int2(ncid,id_elevationNumber,elevationNumber)
       if( ier.ne.0 )return
       inumber=elevationNumber
       ier=nf_inq_varid(ncid,'elevationAngle',id_elevationAngle)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_elevationAngle,elevation)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'numRadials',id_numRadials)
       if( ier.ne.0 )return
       ier=nf_get_var_int2(ncid,id_numRadials,numRadials)
       if( ier.ne.0 )return
       nobeam=numRadials
       ier=nf_inq_varid(ncid,'esStartTime',id_esStartTime)
       if( ier.ne.0 )return
       ier=nf_get_var_double(ncid,id_esStartTime,starttime)
       if( ier.ne.0 )return
       istarttime=starttime
       ier=nf_inq_varid(ncid,'firstGateRangeZ',id_firstGateRangeZ)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_firstGateRangeZ,rangeZ)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'firstGateRangeV',id_firstGateRangeV)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_firstGateRangeV,rangeV)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'gateSizeZ',id_gateSizeZ)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_gateSizeZ,spaceZ)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'gateSizeV',id_gateSizeV)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_gateSizeV,spaceV)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'numGatesZ',id_numGatesZ)
       if( ier.ne.0 )return
       ier=nf_get_var_int2(ncid,id_numGatesZ,numGatesZ)
       if( ier.ne.0 )return
       nogateZ=numGatesZ
       ier=nf_inq_varid(ncid,'numGatesV',id_numGatesV)
       if( ier.ne.0 )return
       ier=nf_get_var_int2(ncid,id_numGatesV,numGatesV)
       if( ier.ne.0 )return
       nogateV=numGatesV
       ier=nf_inq_varid(ncid,'nyquist',id_nyquist)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_nyquist,snyquist)
       if( ier.ne.0 )return

!-----------------------------------------------------------------------
!c  To get check (len_z_bin,numGatesZ), (len_v_bin,numGatesV), and
!c               (len_unlimit,numRadials)

       if( len_z_bin.ne.nogateZ )then
!          print*,'Z_bin: ',len_z_bin,' numGatesZ: ',numGatesZ
           return
       endif
       if( len_v_bin.ne.nogateV )then
!          print*,'V_bin: ',len_v_bin,' numGatesV: ',numGatesV
           return
       endif
       if( len_unlimit.ne.nobeam )then
!          print*,'unlimit: ',len_unlimit,' numRadials: ',numRadials
           return
       endif

!-----------------------------------------------------------------------
!c  To get radialAzim, radialElev, and radialTime.

       ier=nf_inq_varid(ncid,'radialAzim',id_radialAzim)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_radialAzim,azim)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'radialElev',id_radialElev)
       if( ier.ne.0 )return
       ier=nf_get_var_real(ncid,id_radialElev,elev)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'radialTime',id_radialTime)
       if( ier.ne.0 )return
       ier=nf_get_var_double(ncid,id_radialTime,radialTime)
       if( ier.ne.0 )return
       do i=1,nobeam
          itim(i)=radialTime(i)
       enddo

!-----------------------------------------------------------------------
!c  To get Z, V, and W.

       ier=nf_inq_varid(ncid,'Z',id_zbyte)
       if( ier.ne.0 )return
       ier=nf_get_var_int1(ncid,id_zbyte,zbyte)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'V',id_vbyte)
       if( ier.ne.0 )return
       ier=nf_get_var_int1(ncid,id_vbyte,vbyte)
       if( ier.ne.0 )return
       ier=nf_inq_varid(ncid,'W',id_wbyte)
       if( ier.ne.0 )return
       ier=nf_get_var_int1(ncid,id_wbyte,wbyte)
       if( ier.ne.0 )return
       do i=1,nobeam*nogateZ
          if( (zbyte(i).eq.-1).or.(zbyte(i).eq.0) )then
              iz(i)=i_missing
          else
              ii=zbyte(i)
              if( ii.lt.0 )ii=ii+256
              iz(i)=50*(ii-2)-3200
          endif
       enddo
       do i=1,nobeam*nogateV
          if( (vbyte(i).eq.-1).or.(vbyte(i).eq.0) )then
              iv(i)=i_missing
          else
              ii=vbyte(i)
              if( ii.lt.0 )ii=ii+256
              iv(i)=50*(ii-129)
          endif
          if( (wbyte(i).eq.-1).or.(wbyte(i).eq.0) )then
              iw(i)=i_missing
          else
              ii=wbyte(i)
              if( ii.lt.0 )ii=ii+256
              iw(i)=50*(ii-129)
          endif
       enddo

!-----------------------------------------------------------------------
!c  To close NetCDF file.

       ier=nf_close(ncid)
       if( ier.ne.0 )return
       istatus=1

       return
       end
