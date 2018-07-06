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
        subroutine array_diagnosis(a,imax,jmax,name)
c
c.....  Routine to calculate stats on array 'a', that has a 10-character
c.....  name 'name'.
c
        dimension a(imax,jmax)
        character name*10
c
        write(6,1001) name,imax,jmax
1001    format(1x,'ARRAY DIAGNOSIS ROUTINE FOR ARRAY ',A10/
     1 ' WITH DIMENSION ',2I6)
c
        suma=0
        std=0
        aae=0
        amin=1.e30
        amax=-1.e30
        cnt=0
        cnt0=0
        do j=1,jmax
        do i=1,imax
          if(a(i,j).eq.0.0) then
            cnt0=cnt0+1.
          else
            suma=suma+a(i,j)
            cnt=cnt+1.
            if(a(i,j).gt.amax) then
              imxsv=i
              jmxsv=j
              amax=a(i,j)
            endif
            if(a(i,j).lt.amin) then
              imnsv=i
              jmnsv=j
              amin=a(i,j)
            endif
          endif
        enddo
        enddo
c
        if(cnt.eq.0.) then
          write(6,1002)
1002      format(1x,' The array is all zeros')
          return
        endif
c
        amean=suma/cnt
        do j=1,jmax
        do i=1,imax
          if(a(i,j).ne.0.) then
            df=a(i,j)-amean
            std=(df)**2
            aae=abs(df)+aae
          endif
        enddo
        enddo
c
        std=sqrt(std/cnt)
        aae=aae/cnt
c
        write(6,1000) cnt,cnt0,amean,std,aae,amax,imxsv,jmxsv,
     1  amin,imnsv,jmnsv
1000    format(1x,'There are ',f8.0,'  grid points with data '/
     1 1x,' There are ',f8.0,' points with zero '/
     2 1x,' The data mean is ',e12.3,' ; with std dev of ',e12.4/
     3 1x,' The abs deviation from the mean is ',e12.4/
     4 1x,' The maximum is ',e12.4, ' at I,J ',2i6/
     5 1x,' The minimum is ',e12.4, ' at I,J ',2i6)
c
        return
        end
c
c------------------------------------------------------------------
c
        subroutine get_mxmn_2d(nx,ny,array_2d,rmx2d,rmn2d
     &,imx,jmx,imn,jmn)

        implicit none
        integer  nx,ny
        real     array_2d(nx,ny)
        integer  i,j,istatus
        integer  imx,jmx,imn,jmn
        real     rmx2d,rmn2d
        real     r_missing_data

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus.ne.1)then
           print*,'get_mxmn_2d: Error getting r_missing_data'
           stop
        endif

        rmx2d=(-1.)*r_missing_data
        rmn2d=      r_missing_data
        do j=1,ny
        do i=1,nx
           if(array_2d(i,j).gt.rmx2d.and.
     .        array_2d(i,j).ne.r_missing_data)then
              rmx2d=array_2d(i,j)
              imx=i
              jmx=j
           endif
           if(array_2d(i,j).lt.rmn2d)then
              rmn2d=array_2d(i,j)
              imn=i
              jmn=j
           endif
        enddo
        enddo
        return
        end


       subroutine array_range(a,ni,nj,rmin,rmax,r_missing_data)

       real a(ni,nj)

       rmin =  abs(r_missing_data)
       rmax = -abs(r_missing_data)

       nmsg = 0

       do i = 1,ni
       do j = 1,nj
           if(a(i,j) .ne. r_missing_data)then
               rmin = min(rmin,a(i,j))
               rmax = max(rmax,a(i,j))
           else
               if(nmsg .eq. 0)then
                   write(6,*)' Missing data detected in array_range'
               endif
               nmsg = nmsg + 1
           endif
       enddo ! j
       enddo ! i

       frac_msg = float(nmsg) / (float(ni*nj))
       write(6,*)' Fraction of array with data is ',1.-frac_msg

       return
       end

