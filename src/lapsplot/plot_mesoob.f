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
c
c
        subroutine plot_mesoob(dir,spd,gust,t,td,p,ri,rj,lat,lon,
     &                   imax,jmax,relsize_in,zoom,icol_in,du2,iflag)

        include 'lapsparms.cmn'

        real*4 lat(imax,jmax),lon(imax,jmax)
        character*3 t1,td1,p1

        call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)
!       write(6,1234) mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype
 1234   format(1x,4i5,4e12.4,i4)

        zoom_eff = max((zoom / 3.0),1.0)
        relsize = relsize_in / zoom_eff

        du_b=(imax)/300. * relsize

        jsize = nint(0.4 * relsize) - 1

        write(6,*)' relsize,du_b,jsize,zoom = ',relsize,du_b,jsize,zoom       

        call get_border(imax,jmax,x_1,x_2,y_1,y_2)
        call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax))

        rot = (standard_longitude - lon(nint(ri),nint(rj))) / 57.295

!       Convert ri and rj to x1 and y1 (U and V)
!       call supcon(alat,alon,x1,y1)
        x1 = umin + (umax - umin) * (ri-1.) / float(imax-1)
        y1 = vmin + (vmax - vmin) * (rj-1.) / float(jmax-1)

        xsta=ri
        ysta=rj

        if(dir .ge. 0.  .and. spd .ge. 0. .and.
     1     dir .le. 360 .and. spd .le. 200.       )then
            call barbs(spd,dir,ri,rj,du_b,rot,-1e10,+1e10,-1e10,+1e10)
            if(spd .ge. 1.0)then
                call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
                call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)
            endif
        else
            call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
            call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)
        endif

        u = ri
        v = rj

        if(iflag .eq. 3)then ! Plot on top of station location for 'tmg'
            du = 0.
        else
            du = du_b
        endif

        dv   = 1.2 * du
        du_t = 3.0 * du
        du_p = 3.0 * du

        charsize = .0040 / zoom_eff

!       Plot Temperature       
        if(t.gt.-75. .and. t.lt.140.) then 
           write(t1,100,err=20) nint(t)
!          call pwrity(u-du_t,v+dv,t1,3,jsize,0,0)
           CALL PCLOQU(u-du_t,v+dv,t1,charsize,ANGD,CNTR)
        endif
 100    format(i3)
c

!       Plot Dew Point
 20     if(td.gt.-75. .and. td.lt.100.) then
           write(td1,100,err=30) nint(td)
!          call pwrity(u-du_t,v-dv,td1,3,jsize,0,0)
           CALL PCLOQU(u-du_t,v-dv,td1,charsize,ANGD,CNTR)
        endif
c
 30     if(p .gt. 0. .and. p .lt. 10000.) then
           if(p .gt. 1000.) p = p - 1000.
           ip = ifix( p )
           write(p1,101,err=40) ip
 101       format(i3.3)
!          call pwrity(u+du_p,v+dv,p1,3,jsize,0,0)
           CALL PCLOQU(u+du_p,v+dv,p1,charsize,ANGD,CNTR)
        endif

        if(iflag .eq. 1)then ! Plot Gusts (FSL WWW)
           if(gust .gt. 40)then
               ig = int(gust)
               write(p1,102,err=40) ig
               call setusv_dum(2HIN,4)
               dg = 3.0 * du
!              call pwrity(u,v+dg,p1,3,jsize,0,0)          ! On Top
!              call pwrity(u+du_p,v-dv,p1,3,jsize,0,0)     ! Lower Right
 102           format('G',i2)
               call setusv_dum(2HIN,icol_in)
           endif
        endif           

c
 40     continue
        return
        end
