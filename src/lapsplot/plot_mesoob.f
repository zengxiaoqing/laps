cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
c
c
        subroutine plot_mesoob(dir,spd,gust,t,td,p,ri,rj
     1                        ,lat,lon,imax,jmax,relsize_in,zoom
     1                        ,icol_in,du2,iflag,iflag_cv)

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

        if(iflag_cv .eq. 0)then ! Normal obs plot
            if(dir .ge. 0.  .and. spd .ge. 0. .and.
     1         dir .le. 360 .and. spd .le. 200.       )then
                call barbs(spd,dir,ri,rj,du_b,rot
     1                    ,-1e10,+1e10,-1e10,+1e10)
                if(spd .ge. 1.0)then
                    call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
                    call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)
                endif
            else
                call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
                call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)
            endif

!           Plot Temperature       
            if(t.gt.-75. .and. t.lt.140.) then 
               write(t1,100,err=20) nint(t)
!              call pwrity(u-du_t,v+dv,t1,3,jsize,0,0)
               CALL PCLOQU(u-du_t,v+dv,t1,charsize,ANGD,CNTR)
            endif
 100        format(i3)
 20         continue

!           Plot Dew Point
            if(td.gt.-75. .and. td.lt.100.) then
               write(td1,100,err=30) nint(td)
               CALL PCLOQU(u-du_t,v-dv,td1,charsize,ANGD,CNTR)
            endif
 30         continue
 
!           Plot Pressure
            if(p .gt. 0. .and. p .lt. 10000.) then
               if(p .gt. 1000.) p = p - 1000.
               ip = ifix( p )
               write(p1,101,err=40) ip
 101           format(i3.3)
!              call pwrity(u+du_p,v+dv,p1,3,jsize,0,0)
               CALL PCLOQU(u+du_p,v+dv,p1,charsize,ANGD,CNTR)
            endif

!           Plot Gusts (FSL WWW)
            if(iflag .eq. 1)then 
               if(gust .gt. 40)then
                   ig = int(gust)
                   write(p1,102,err=40) ig
                   call setusv_dum(2HIN,4)
                   dg = 3.0 * du
!                  call pwrity(u,v+dg,p1,3,jsize,0,0)          ! On Top
!                  call pwrity(u+du_p,v-dv,p1,3,jsize,0,0)     ! Lower Right
 102               format('G',i2)
                   call setusv_dum(2HIN,icol_in)
               endif
            endif           

        else ! C&V plot
            call plot_circle(u,v,du*0.8)

!           Plot Visibility
            if(td.gt.-75. .and. td.lt.100.) then
               write(td1,100,err=31) nint(td)
               call left_justify(td1)
               CALL PCLOQU(u+du_t,v-dv,td1,charsize,ANGD,CNTR)
            endif
 31         continue

        endif

c
 40     continue
        return
        end
