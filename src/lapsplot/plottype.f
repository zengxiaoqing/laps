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
        subroutine plot_types_2d(icldpcp_type_2d,interval,size
     1                          ,c2_field,l_meta,plot_parms
     1                          ,imax,jmax,lat,lon,ifield_2d)

        integer icldpcp_type_2d(imax,jmax)
        integer ifield_2d(imax,jmax)

        logical l_meta,l_new

        real lat(imax,jmax)
        real lon(imax,jmax)

        character barg
        integer iarg,byte_to_i4

        include 'lapsplot.inc'

        character*2 c2_field,c2_type

        character*2 c2_cloud_types(0:10)
!       data c2_cloud_types
!    1  /'  ','St','Sc','Cu','Ns','Ac','As','Cs','Ci','Cc','Cb'/
        data c2_cloud_types
     1  /'  ','ST','SC','CU','NS','AC','AS','CS','CI','CC','CB'/

        character*2 c2_precip_types(0:10)
        data c2_precip_types
     1  /'  ','Rn','Sn','Zr','Sl','Ha','L ','Zl','  ','  ','  '/

        character*1 c1_precip_types(0:10)
        data c1_precip_types
     1  /' ','R','*','Z','I','H','L','F','M','M',' '/

        ANGD = 0.

!       Pull out relavant bits
        do i = 1,imax
        do j = 1,jmax
!          barg = cldpcp_type_2d(i,j)
!          iarg = byte_to_i4(barg)
           iarg = icldpcp_type_2d(i,j)
           if(c2_field .eq. 'tc' .or. c2_field .eq. 'cy')then
               ifield_2d(i,j) = iarg ! - iarg/16*16
           elseif(c2_field .eq. 'tp' .or. c2_field .eq. 'py')then
               ifield_2d(i,j) = iarg ! /16
           endif
        enddo ! j
        enddo ! i

        do j=jmax,1,-3
            if(c2_field .eq. 'tc' .or. c2_field .eq. 'cy')then
                write(6,500)(c2_cloud_types(ifield_2d(i,j)),i=1,imax,2)
            elseif(c2_field .eq. 'tp' .or. c2_field .eq. 'py')then
                write(6,500)(c2_precip_types(ifield_2d(i,j)),i=1,imax,2)
            endif
500         format(1x,29(1x,a2))
        enddo ! j

        if(.not. l_meta)return

        if(c2_field .eq. 'tc' .or. c2_field .eq. 'cy')then
            l_new = .true. ! Testing as T/F for cloud type

            nc = 2

            iskip = max( nint( (jmax / 30) / plot_parms%zoom), 1)

        elseif(c2_field .eq. 'tp' .or. c2_field .eq. 'py')then
            l_new = .true.

            nc = 1

            iskip = max( nint( (jmax / 60) / plot_parms%zoom), 1)

        endif

        if(l_new)then
            size_z = size / plot_parms%zoom
            write(6,*)' size_z, iskip = ',size_z,iskip
            call get_border(imax,jmax,x_1,x_2,y_1,y_2)
            call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax),1)
        else
            write(6,*)' isize,iskip',isize,iskip
        endif

        do j = 1,jmax,iskip
        do i = 1,imax,iskip

            alat = lat(i,j)
            alon = lon(i,j)

            if(ifield_2d(i,j) .ne. 0)then
                if(c2_field .eq. 'tc' .or. c2_field .eq. 'cy')then
                    c2_type = c2_cloud_types(ifield_2d(i,j))
                    nc = 2
                elseif(c2_field .eq. 'tp' .or. c2_field .eq. 'py')then
                    c2_type = c1_precip_types(ifield_2d(i,j))
                    nc = 1
                endif

!               if(i .eq. j)then
!                   c2_type = '*'
!               else
!                   c2_type = ' '
!               endif
     
                if(.not. l_new)then ! old way
!                   write(6,*)'Plot type: ',c2_type,i,j
                    ri = i
                    rj = j
                    call plot_typeob(c2_type,nc,ri,rj,isize,imax,jmax)
                else
!                   write(6,*)'Plot type: ',c2_type,i,j
                    ri = i
                    rj = j
!                   CALL PCLOQU(ri, rj, c2_type, size_z, ANGD, CNTR)
                    CALL PCMEQU(ri, rj, c2_type, size_z, ANGD, CNTR)       
                endif
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
        call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax),ltype)       

        call pwrity (ri, rj, c_type, nc, isize, 0, 0)

    1   continue
        return
        end
