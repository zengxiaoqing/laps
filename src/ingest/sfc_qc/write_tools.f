c
c
      subroutine writev(a,imax,jmax,m,cs,atime,onoff,offset)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      parameter (badflag = -99.9)
      real a(m,m)
      integer ii(m,m)   ,onoff,m 
      real smax,smin,s,sum,rms
      character*12 cs,atime*24
      integer i,j
c
      smin=1.e25
      smax=-1.e25
      sum=0.
      rms=0.
      cnt=0.
      do j=1,jmax
      do i=1,imax
         ii(i,j)=0 
         if(a(i,j).ne.badflag) then
            sum=sum+a(i,j)-offset
            if((a(i,j)-offset).gt.smax) smax=a(i,j)-offset
            if((a(i,j)-offset).lt.smin) smin=a(i,j)-offset
            cnt=cnt+1.
         endif
      enddo !i
      enddo !j
      sum=sum/cnt              
      s=amax1(abs(smax),abs(smin)) 
      do i=1,25
         if((s/10.**(i-12)).lt.1000.) then
            iii=i-12
            go to 1
         endif
      enddo !i
    1 write(6,999) cs
 999  format(//1x,a12)
      write(6,1000) atime,iii, smax,smin
 1000 format(13x,a24,' x 10**',i3,' Max: ',e10.4,'  Min: ',e10.4)
      do j=1,jmax
      do i=1,imax
         if(a(i,j).ne.badflag) then
            rms=rms+((a(i,j)-offset)-sum)**2 
            ii(i,j)=(a(i,j)-offset)*10.**(-iii)+.5
         endif
      enddo !i
      enddo !j
      rms=sqrt(rms/cnt)                
      write(6,1004) rms,sum,offset
 1004 format(1x,'RMS  = ',e12.4,' MEAN = ',e12.4,' OFFSET = ',e12.4)
 1001 format(1x,'transposed vector')
      if (jmax.eq.1) then
         write(6,1001)
         if(onoff.eq.1) write(6,1002) (ii(i,1),i=1,imax)
      else
         continue
         if(onoff.eq.1) then
            do i=1,imax
               write(6,1005) i
 1005          format(/i3)
               write(6,1002) (j,ii(i,j),j=1,jmax)
            enddo !i
         endif
      endif
 1002 format(10i7)
c
      return
      end 
c
c
      subroutine writei(a,imax,m,cs,atime)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      integer a(m)
      character*72 cs,atime*24
      integer i
c
      write(6,1000) cs,atime
 1000 format(//1x,a72)
      write(6,1002) (a(i)  ,i=1,imax)
 1002 format(10i7)
c
      return
      end 
c
c
      subroutine writec(a,imax,m,cs,atime)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      character*5   a(m)
      character*14 cs,atime*24
c
      write(6,1000) cs,atime
 1000 format(//1x,a14)
      write(6,1002) (a(i)(1:5)  ,i=1,imax)
 1002 format(10(2x,a5))
c
      return
      end 
c
c
       Subroutine writewvv(wvr,maxsta,m,ncycles,wvv,it) 

c
c*********************************************************************
c     Subroutine puts the latest sets of errors into wvv     for
c      error averaging over ncycles
c
c     Original: John McGinley, NOAA/FSL  Summer 2004
c
c*********************************************************************
c
       real wvr(m)
       real wvv(ncycles,m)      
c
          do l=ncycles-1,1,-1
             do i=1,maxsta
                wvv(l+1,i)=wvv(l,i)       
             enddo !i
          enddo !l
          do i=1,maxsta
             wvv(1,i)=wvr(i)       
          enddo !i
c
       return
       end
c
       Subroutine writemon(ta,tda,ua,va,pmsla,alta,
     &   nvar,maxsta,m,ncycles,monster,it) 
c
c*********************************************************************
c     Subroutine puts the latest sets of obs into monster for
c     fourier processing.
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          New version.
c
c*********************************************************************
c
       real ta(m),tda(m),ua(m),va(m),pmsla(m),alta(m)
       real monster(m,ncycles,nvar)
c
       do k=1,nvar
          do l=ncycles-1,1,-1
             do i=1,maxsta
                monster(i,l+1,k)=monster(i,l,k)
             enddo !i
          enddo !l
       enddo !k
          do i=1,maxsta
             monster(i,1,1)=ta(i)
             monster(i,1,2)=tda(i)
             monster(i,1,3)=ua(i)
             monster(i,1,4)=va(i)  
             monster(i,1,5)=pmsla(i)
             monster(i,1,6)=alta(i)
          enddo !i
c
       return
       end
c
       Subroutine writeqc(qta,qtda,qua,qva,qpmsla,qalta,
     &   nvar,maxsta,m,ncycles,qmonster,it) 
c
c*********************************************************************
c     Subroutine puts the latest sets of obs into monster for
c     fourier processing.
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          New version.
c
c*********************************************************************
c
       integer qta(m),qtda(m),qua(m),qva(m),qpmsla(m),qalta(m)
       integer qmonster(m,ncycles,nvar)
c
       do k=1,nvar
          do l=ncycles-1,1,-1
             do i=1,maxsta
                qmonster(i,l+1,k)=qmonster(i,l,k)
             enddo !i
          enddo !l
       enddo !k
          do i=1,maxsta
             qmonster(i,1,1)=qta(i)
             qmonster(i,1,2)=qtda(i)
             qmonster(i,1,3)=qua(i)
             qmonster(i,1,4)=qva(i)  
             qmonster(i,1,5)=qpmsla(i)
             qmonster(i,1,6)=qalta(i)
          enddo !i
c
       return
       end
c
c
      subroutine write_qc_cdf(i4time_file, dir, ext, m, 
     1   num_sta, 
     1   stations, provider, reptype, lat, lon, elev,
     1   qcstat,  t,    tb,    ta,    tc,    te,    tf,  
     1   qcstatd, td,   tdb,   tda,   tdc,   tde,   tdf,  
     1   qcstauv, u,    ub,    ua,    uc,    ue,    uf,  
     1            v,    vb,    va,    vc,    ve,    vf,  
     1   qcstapm, pmsl, pmslb, pmsla, pmslc, pmsle, pmslf, 
     1   qcstal,  alt,  altb,  alta,  altc, alte,   altf, 
     1   status)

c*********************************************************************
c     Original: Linda Wharton, NOAA/FSL  05 Oct 1998
c       Subroutine write_qc_cdf to write output to a netCDF
c         file rather than the "temploc" file.
c
c*********************************************************************

	integer i4time_file  ! valid time of file
	character dir*(*)      ! directory for output
	character ext*(*)      ! ext for output
	integer   m            ! dimension of data arrays
        integer   num_sta      ! number of valid stations in data arrays
        character*(*) stations(m) !should be length 20 
	character*(*) provider(m) !should be length 11
        character*(*) reptype(m)  ! should be length 6
	real    lat(m), lon(m), elev(m)
        integer   qcstat(m),qcstatd(m),qcstauv(m),qcstapm(m),qcstal(m)
        real    t(m), td(m), u(m), v(m),  pmsl(m),  alt(m)
        real    ta(m),tda(m),ua(m),va(m), pmsla(m), alta(m)
        real    tb(m),tdb(m),ub(m),vb(m), pmslb(m), altb(m)
        real    tc(m),tdc(m),uc(m),vc(m), pmslc(m), altc(m)
        real    te(m),tde(m),ue(m),ve(m), pmsle(m), alte(m)
        real    tf(m),tdf(m),uf(m),vf(m), pmslf(m), altf(m)
	integer   status

c Local variables

        character filename*200, cdl_path*190, gtime*9, fcst_hh_mm*4
	integer   fn_len, ext_len, sta_len, pro_len, type_len
	integer   istatus, error(2) 
	integer i_reftime, i_valtime
c
c begin
c
        error(1)=1 ! success
        error(2)=0 ! failure


c make output file name from i4time_file 

        call make_fnam_lp(i4time_file,gtime,istatus)
        if (istatus .ne. 1) then
          write (6,*)
     1'Error converting i4time to file name... writing QC file aborted.'
          status=error(2)
          return
        endif

        call s_len(ext, ext_len)

        fcst_hh_mm = '0000'

        call cvt_fname_v3(dir, gtime, fcst_hh_mm, ext, ext_len,
     1                    filename, fn_len, istatus)

c get cdl_path
	call get_directory('cdl',cdl_path, cdl_len)

c setup lengths of character variables for passing into C code
        sta_len = len(stations(1))
        pro_len = len(provider(1))
        type_len = len(reptype(1))

c setup valtime and reftime for writing into output file
	i_reftime = i4time_file - 315619200
	i_valtime = i_reftime

c pass data into C code for generating netCDF output file
	call write_qc( filename, cdl_path, num_sta, m, i_reftime, 
     1        i_valtime, stations, provider, reptype, sta_len, 
     1        pro_len, type_len, fn_len, cdl_len, lat, lon, elev, 
     1        qcstat,  t,    tb,    ta,    tc,    te,    tf,  
     1        qcstatd, td,   tdb,   tda,   tdc,   tde,   tdf,  
     1        qcstauv, u,    ub,    ua,    uc,    ue,    uf,  
     1                 v,    vb,    va,    vc,    ve,    vf,  
     1        qcstapm, pmsl, pmslb, pmsla, pmslc, pmsle, pmslf, 
     1        qcstal,  alt,  altb,  alta,  altc, alte,   altf, 
     1        istatus)

	if (istatus .eq. error(1)) then   !successful completion of 
	  status = error(1)
	  goto 99
	else
          if (istatus .eq. 2)  goto 96  !error writing file
          if (istatus .eq. -1)  goto 95  !error reading dimension data
          if (istatus .eq. -2)  goto 94  !error creating netCDF file
	endif

94	write (6,*) 'Error creating netCDF file...write aborted.'
	status=error(2)
	goto 99

95	write (6,*) 
     1'Error reading dimensions from netCDF file...write aborted.'
	status=error(2)
	goto 99

96	write (6,*) 'Error writing text to  netCDF file...write aborted.'
	status=error(2)
	goto 99

99      return
	end
