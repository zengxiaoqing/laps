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
        subroutine plot_types_2d(cldpcp_type_2d,interval,size,c2_field,l
     1_meta
     1                                  ,imax,jmax,lat,lon,ifield_2d)

        byte cldpcp_type_2d(imax,jmax)
        integer*4 ifield_2d(imax,jmax)

        logical l_meta

        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)

        byte barg
        integer*4 iarg,byte_to_i4

        character*2 c2_field,c2_type

!       character*2 c2_cloud_types(0:10)
!       1       /'  ','St','Sc','Cu','Ns','Ac','As','Cs','Ci','Cc','Cb'/

        character*2 c2_cloud_types(0:10)
     1  /'  ','ST','SC','CU','NS','AC','AS','CS','CI','CC','CB'/

        character*2 c2_precip_types(0:10)
     1  /'  ','Rn','Sn','Zr','Sl','Ha','L ','Zl','  ','  ','  '/

        character*1 c1_precip_types(0:10)
     1  /' ','R','*','Z','I','H','L','F',' ',' ',' '/


!       Pull out relavant bits
        do i = 1,imax
        do j = 1,jmax
           barg = cldpcp_type_2d(i,j)
           iarg = byte_to_i4(barg)
           if(c2_field .eq. 'tc')then
               ifield_2d(i,j) = iarg - iarg/16*16
           elseif(c2_field .eq. 'tp')then
               ifield_2d(i,j) = iarg/16
           endif
        enddo ! j
        enddo ! i

        do j=jmax,1,-3
            if(c2_field .eq. 'tc')then
                write(6,500)(c2_cloud_types(ifield_2d(i,j)),i=1,imax,2)
            elseif(c2_field .eq. 'tp')then
                write(6,500)(c2_precip_types(ifield_2d(i,j)),i=1,imax,2)
            endif
500         format(1x,29(1x,a2))
        enddo ! j

        if(.not. l_meta)return

        if(c2_field .eq. 'tc')then
            nc = 2

            if(imax .gt. 110 .or. jmax .gt. 110)then
                isize = 0
                iskip = 3
            elseif(imax .gt. 75 .or. jmax .gt. 75)then
                isize = 0
                iskip = 2
            else
                isize = 1
                iskip = 2
            endif

        elseif(c2_field .eq. 'tp')then
            nc = 1

            if(imax .gt. 110 .or. jmax .gt. 110)then
                iskip = 2
                isize = 0
            elseif(imax .gt. 90 .or. jmax .gt. 90)then
                iskip = 2
                isize = 1
            elseif(imax .gt. 75 .or. jmax .gt. 75)then
                iskip = 1
                isize = 0
            else
                iskip = 1
                isize = 1
            endif

        endif

        write(6,*)' isize,iskip',isize,iskip

        do j = 1,jmax,iskip
        do i = 1,imax,iskip

            alat = lat(i,j)
            alon = lon(i,j)

            if(ifield_2d(i,j) .ne. 0)then
                if(c2_field .eq. 'tc')then
                    c2_type = c2_cloud_types(ifield_2d(i,j))
                    nc = 2
                elseif(c2_field .eq. 'tp')then
                    c2_type = c1_precip_types(ifield_2d(i,j))
                    nc = 1
                endif

!               if(i .eq. j)then
!                   c2_type = '*'
!               else
!                   c2_type = ' '
!               endif
     
                ri = i
                rj = j
                call plot_typeob(c2_type,nc,ri,rj,isize,imax,jmax)
            endif

        enddo ! i
        enddo ! j

        return
        end

        subroutine plot_typeob(c_type,nc,ri,rj,isize,imax,jmax)

        include 'lapsparms.cmn'

        character*(*) c_type

        call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)

!       write(6,1234) mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype
 1234   format(1x,4i5,4e12.4,i4)

        call get_border(imax,jmax,x_1,x_2,y_1,y_2)
        call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax))

        call pwrity (ri, rj, c_type, nc, isize, 0, 0)

    1   continue
        return
        end
