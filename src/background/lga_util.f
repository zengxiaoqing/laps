      subroutine get_laps_lat_lon(dir,ext,nx,ny,lat,lon,istatus)
c
c *** Modified from code generated by gennet (B. Schwartz, FSL)  
c
      
c

      include 'netcdf.inc'

c
      integer nx,ny
      real*4 lat(nx,ny),
     .       lon(nx,ny)
c
      integer len,elen,i,istatus
      integer start(10),count(10),
     .          ncid,ntp,nvdim,nvs,lenstr,ndsize,rcode
      integer   vdims(10) 
c
      character*(*) dir
      character*500 ldir
      character*(*) ext
      character*500 lext
      character*31  dummy
c
      logical exists
c_______________________________________________________________________________
c
c *** Create netcdf file name and check for existence.
c
      ldir=dir
c      len=index(ldir,' ')-1

      call s_len(ldir,len)
      if (ldir(len:len) .ne. '/') then
         ldir(len+1:len+1)='/'
         len=len+1
      endif
      lext=ext
c      elen=index(lext,' ')-1

      call s_len(lext,elen)
      inquire(file=ldir(1:len)//'static.'//lext(1:elen),exist=exists)
      if (.not. exists) then
         print *,' *** LAPS static file does not exist:  ',
     .           ldir(1:len)//'static.'//lext(1:elen)
         istatus=0
         return
      endif
c
c *** Open netcdf file.
c
      rcode=NF_OPEN(ldir(1:len)//'static.'//lext(1:elen)
     +     ,NF_NOWRITE,ncid)
c
c *** Statements to fill lat, lon.                             
c
      call NCVINQ(ncid,1,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do i=1,nvdim
         call NCDINQ(ncid,vdims(i),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(i)=1
         count(i)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,1,start,count,lat)
c
      call NCVINQ(ncid,2,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do i=1,nvdim
         call NCDINQ(ncid,vdims(i),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(i)=1
         count(i)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,2,start,count,lon)
c
c *** Close netcdf file.
c
      rcode= NF_CLOSE(ncid)
c
      istatus=1
      return
      end
c
c===============================================================================
c
      subroutine gdtost(ab,ix,iy,stax,stay,staval)
c
c *** Subroutine to return stations back-interpolated values(staval)
c        from uniform grid points using overlapping-quadratics.
c        gridded values of input array a dimensioned ab(ix,iy), where
c        ix = grid points in x, iy = grid points in y.  Station
c        location given in terms of grid relative station x (stax)
c        and station column.
c *** Values greater than 1.0e30 indicate missing data.
c
      dimension ab(ix,iy),r(4),scr(4)
c_______________________________________________________________________________
c
      iy1=int(stay)-1
      iy2=iy1+3
      ix1=int(stax)-1
      ix2=ix1+3
      staval=1e30
      fiym2=float(iy1)-1
      fixm2=float(ix1)-1
      ii=0
      do i=ix1,ix2
         ii=ii+1
         if (i .ge. 1 .and. i .le. ix) then 
            jj=0
            do j=iy1,iy2
               jj=jj+1
               if (j .ge. 1 .and. j .le. iy) then
                  r(jj)=ab(i,j)
               else
                  r(jj)=1e30
               endif
            enddo
            yy=stay-fiym2
            if (yy .eq. 2.0) then
               scr(ii)=r(2)
            else
               call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii))
            endif
         else 
            scr(ii)=1e30
         endif
      enddo
      xx=stax-fixm2
      if (xx .eq. 2.0) then
         staval=scr(2)
      else
         call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval)
      endif
c
      return
      end
c
c===============================================================================
c
      subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
c
      yyy=1e30
      if (x2 .gt. 1.e19 .or. x3 .gt. 1.e19 .or.
     .    y2 .gt. 1.e19 .or. y3 .gt. 1.e19) return
c
      wt1=(xxx-x3)/(x2-x3)
      wt2=1.0-wt1
c
      if (y4 .lt. 1.e19 .and. x4 .lt. 1.e19) then
c        yz22=(xxx-x3)*(xxx-x4)/((x2-x3)*(x2-x4))
         yz22=wt1*(xxx-x4)/(x2-x4)
c        yz23=(xxx-x2)*(xxx-x4)/((x3-x2)*(x3-x4))
         yz23=wt2*(xxx-x4)/(x3-x4)
         yz24=(xxx-x2)*(xxx-x3)/((x4-x2)*(x4-x3))
      else
         yz22=wt1
         yz23=wt2
         yz24=0.0
      endif
c
      if (y1 .lt. 1.e19 .and. x1 .lt. 1.e19) then
         yz11=(xxx-x2)*(xxx-x3)/((x1-x2)*(x1-x3))
c        yz12=(xxx-x1)*(xxx-x3)/((x2-x1)*(x2-x3))
         yz12=wt1*(xxx-x1)/(x2-x1)
c        yz13=(xxx-x1)*(xxx-x2)/((x3-x1)*(x3-x2))
         yz13=wt2*(xxx-x1)/(x3-x1)
      else
         yz11=0.0
         yz12=wt1
         yz13=wt2
      endif
c
      if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
         yyy=wt1*y2+wt2*y3
      else
         yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)+wt2*(yz22*y2+yz23*y3+yz24*y4)
      endif
c
      return
      end
c
c===============================================================================
c
      subroutine uvgrid_to_uvtrue_a(u,v,lon,std_lon,nx,ny,nz,angle)
c
c *** Convert grid north winds to true north winds.
c
      
c
      integer nx,ny,nz,i,j
c
      real*4    u(nx,ny,nz),
     .          v(nx,ny,nz),
     .          lon(nx,ny),
     .          std_lon,
     .          angle(nx,ny)
c_______________________________________________________________________________
c
      do j=1,ny
      do i=1,nx
         angle(i,j)=lon(i,j)-std_lon
      enddo
      enddo
c
      call rotate_vec_a(u,v,angle,nx,ny,nz)
c
      return
      end
c
c===============================================================================
c
      subroutine uvtrue_to_uvgrid_a(u,v,lon,std_lon,nx,ny,nz,angle)
c
c *** Convert true north winds to grid north winds.
c
c
      
c
      integer nx,ny,nz,i,j
c
      real*4    u(nx,ny,nz),
     .          v(nx,ny,nz),
     .          lon(nx,ny),
     .          std_lon,
     .          angle(nx,ny)
c_______________________________________________________________________________
c
      do j=1,ny
      do i=1,nx
         angle(i,j)=std_lon-lon(i,j)
      enddo
      enddo
c
      call rotate_vec_a(u,v,angle,nx,ny,nz)
c
      return
      end
c
c===============================================================================
c
      subroutine rotate_vec_a(u,v,angle,nx,ny,nz)
c
      
c
      integer nx,ny,nz,i,j,k
c
      real*4    u(nx,ny,nz),
     .          v(nx,ny,nz),
     .          angle(nx,ny),
     .          utmp
c_______________________________________________________________________________
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
         if (u(i,j,k) .gt. -400. .and. u(i,j,k) .lt. 400. .and.
     .       v(i,j,k) .gt. -400. .and. v(i,j,k) .lt. 400.) then
            utmp    =+u(i,j,k)*cosd(angle(i,j))
     .               +v(i,j,k)*sind(angle(i,j))
            v(i,j,k)=-u(i,j,k)*sind(angle(i,j))
     .               +v(i,j,k)*cosd(angle(i,j))
            u(i,j,k)=utmp
         endif
      enddo
      enddo
      enddo
c
      return
      end
c
c===============================================================================
c
      subroutine thvpc2tq(thv,pc,p,t,q)
c
c *** Subprogram:  thvpc2tq - Calculates temperature (K) and specific 
c                             humidity (kg/kg) given virtual potential 
c                             temperature (K), pressure (mb), and
c                             condensation pressure (mb).
c
c *** Program history log:
c        93-12-20  S. Benjamin - Original version 
c        96-09-17  J. Snook    - esw calculated in a table
c
c *** Usage:  call thvpc2tq(thv,pc,p,t,q)
c
c *** Input argument list:
c        thv    - real  virtual potential temperature (K)
c        pc     - real  condensation pressure (mb)
c        p      - real  pressure (mb)
c
c *** Output argument list:
c        t      - real  temperature (K)
c        q      - real  specific humidity (kg/kg)
c
c *** Subprograms called:
c        tv2tq  - calculate temp and spec. hum. from virtual
c                    temp and relative humidity
c        esw    - calculate saturation vapor pressure (from a table)
c_______________________________________________________________________________
c
      
c
      real tv,rh,p,t,q,thv,kappa,templcl,x,x1,pc
      integer it
      data kappa/0.285714/
c
      real*4 esat,esw
      common /estab/esat(15000:45000),esw(15000:45000)
c_______________________________________________________________________________
c
      tv=thv*(p*0.001)**kappa
      templcl=thv*(pc*0.001)**kappa
      it=tv*100
      it=min(45000,max(15000,it))
      x =esw(it)
      it=templcl*100
      it=min(45000,max(15000,it))
      x1=esw(it)
      rh=x1/x * (p-x) / (pc-x1)
      call tv2tq(tv,rh,p,t,q)
      return
      end
c
c===============================================================================
c
      subroutine tv2tq(tv,rh,p,t,q)
c
c *** Subprogram:  tv2tq - Calculates temperature (K) and specific 
c                          humidity (kg/kg) given virtual temperature (K),
c                          pressure (mb), and relative humidity.
c
c *** Program history log:
c        93-01-12  S. Benjamin - Original version
c
c *** Usage:  call tv2tq(tv,rh,p,t,q)
c
c *** Input argument list:
c        tv     - real  virtual temperature (K)
c        rh     - real  relative humidity (range 0.0-1.0)
c        p      - real  pressure (mb)
c
c *** Output argument list:
c        t      - real  temperature (K) 
c        q      - real  specific humidity (kg/kg)
c
c *** Reamrks:
c        It uses an iterative newton-raphson technique.  Four iterations are
c        generally adequate to provide convergence to 5 decimal places.
c        the wobus function for saturation vapor pressure over liquid water
c        is used.
c_______________________________________________________________________________
c
      
c
      real tv,rh,p,t,q,t1,estv1,etv,t2,estv2,dt,dum
c
      integer j,it
c
      real*4 esat,esw
      common /estab/esat(15000:45000),esw(15000:45000)
c_______________________________________________________________________________
c
      t1 = tv
c
c *** estv = saturation vapor pressure (mb) for tv.
c
      it=t1*100
      it=min(45000,max(15000,it))
      estv1=esw(it)
c
      do j=1,3
c
c ****** etv = vapor pressure (mb) for tv and rh*
c
         etv=estv1*rh
c
c ****** q = mixing ratio for tv (kg/kg).
c
         q=0.62197*etv/(p-etv)
         t2=tv/(1.+0.608*q)
         if (abs(t2-t1) .lt. 0.001) goto 77
c
c ****** estv2 = saturation vapor pressure (mb) for estimated t (=t2).
c
         it=t2*100
         it=min(45000,max(15000,it))
         estv2=esw(it)
c
c ****** etv = vapor pressure (mb) for estimated t and rh.
c
         etv=estv2*rh
c
c ****** q = mixing ratio for estimated t and rh (kg/kg).
c
         q=0.62197*etv/(p-etv)
         t=t2*(1.+0.608*q)
c
c ****** Recalc. tv.
c
         dt=tv-t
         dum=(estv2-estv1)/(t2-t1)
         etv=estv2+dum*dt*rh
c
c ****** Reset t1 and estv1 before next iteration.
c
         t1=t2
         estv1=estv2
c
      enddo
77    continue
      t=t2
c
      return
      end
c
c===============================================================================
c
      subroutine esat_init
c
      common /estab/esat(15000:45000),esw(15000:45000)
c
c *** Create tables of the saturation vapour pressure with up to
c        two decimal figures of accuraccy:
c
      do it=15000,45000
         t=it*0.01
         p1 = 11.344-0.0303998*t
         p2 = 3.49149-1302.8844/t
         c1 = 23.832241-5.02808*alog10(t)
         esat(it) = 10.**(c1-1.3816E-7*10.**p1+
     .               8.1328E-3*10.**p2-2949.076/t)
c
         t=t-273.15
         pol=   0.99999683     + t*(-0.90826951E-02 +
     1       t*(0.78736169E-04 + t*(-0.61117958E-06 +
     2       t*(0.43884187E-08 + t*(-0.29883885E-10 +
     3       t*(0.21874425E-12 + t*(-0.17892321E-14 +
     4       t*(0.11112018E-16 + t*(-0.30994571E-19)))))))))
         esw(it)=6.1078/pol**8
c
      enddo
c
      return
      end
c
c===============================================================================
c
      subroutine filter_2dx(field,ix,iy,iz,smth)
c
c *** Subprogram:  smooth - Smooth a meteorological field.
c     Author:  Stan Benjamin 
c     Date  :  90-06-15
c
c *** Abstract:  Shapiro smoother. 
c 
c *** Program history log: 
c        85-12-09  S. Benjamin - Original version
c        96-06-16  J. Snook    - Modified to do 3d RAMS fields
c                              - hold array is dynamically allocated
c 
c *** Usage:  call smooth(field,ix,iy,iz,smth) 
c
c *** Input argument list: 
c        field    - real array  field(ix,iy,iz)
c                               Meteorological field
c        ix       - integer     x coordinates of field
c        iy       - integer     y coordinates of field
c        iz       - integer     z coordinates of field
c        smth     - real      
c
c *** Output argument list:   
c        field    - real array  field(ix,iy,iz)
c                               Smoothed meteorological field
c 
c *** Remarks:  Reference:  Shapiro, 1970: "Smoothing, filtering, and
c        boundary effects", Rev. Geophys. Sp. Phys., 359-387.
c
c     This filter is of the type 
c        z(i) = (1-s)z(i) + s(z(i+1)+z(i-1))/2
c     for a filter which is supposed to damp 2dx waves completely
c     but leave 4dx and longer with little damping,
c     it should be run with 2 passes using smth (or s) of 0.5
c     and -0.5.
c_______________________________________________________________________________
c   
      
c
      integer ix,iy,iz,i,j,k,i1,i2,it
c
      real field(ix,iy,iz),
     .     hold(ix,2),
     .     smth,smth1,smth2,smth3,smth4,smth5,
     .     sum1,sum2
c_______________________________________________________________________________
c
      smth1=0.25*smth*smth
      smth2=0.50*smth*(1.-smth)
      smth3=(1.-smth)*(1.-smth)
      smth4=(1.-smth)
      smth5=0.5*smth
c
      do k=1,iz
c
         do j=1,2
         do i=1,ix
            hold(i,j)=0.
         enddo
         enddo
c
         i1=2
         i2=1
         do j=2,iy-1
            it=i1
            i1=i2
            i2=it
            do i=2,ix-1
               sum1=field(i-1,j+1,k)+field(i-1,j-1,k)
     .             +field(i+1,j+1,k)+field(i+1,j-1,k)
               sum2=field(i  ,j+1,k)+field(i+1,j  ,k)
     .             +field(i  ,j-1,k)+field(i-1,j  ,k)
               hold(i,i1)=smth1*sum1+smth2*sum2+smth3*field(i,j,k)
            enddo
            if (j .eq. 2) goto 200
            do i=2,ix-1
               field(i,j-1,k)=hold(i,i2)
            enddo
200         continue
         enddo
c
         do i=2,ix-1
            field(i,iy-1,k)=hold(i,i1)
         enddo
c
         do i=2,ix-1
            field(i,1,k)=smth4*field(i,1,k) 
     .                  +smth5*(field(i-1,1,k)+field(i+1,1,k))
            field(i,iy,k)=smth4*field(i,iy,k) 
     .                   +smth5*(field(i-1,iy,k)+field(i+1,iy,k))
         enddo
c
         do j=2,iy-1
            field(1,j,k)=smth4*field(1,j,k) 
     .                  +smth5*(field(1,j-1,k)+field(1,j+1,k))
            field(ix,j,k)=smth4*field(ix,j,k) 
     .                   +smth5*(field(ix,j-1,k)+field(ix,j+1,k))
         enddo
c
      enddo
c
      return
      end
