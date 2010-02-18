
        subroutine radarhist(NX_L,NY_L,NZ_L,dbz_3d
     1                      ,ilow,ihigh,jlow,jhigh
     1                      ,lmask_3d,lun_out)
 
        real dbz_3d(NX_L,NY_L,NZ_L)
        logical lmask_3d(NX_L,NY_L,NZ_L)

        parameter (nbins=10)

        integer ihist(nbins)
        integer ihist_2d(NZ_L,nbins)

!       Initialize histograms
        ihist = 0
        ihist_2d = 0

!       Define area of interest
        klow = 1
        khigh = NZ_L

        hist_low = 0.
        hist_intvl = 10.

!       Note bin #1 is from 0-10 dbz and so on at 10dbz increments
        do k = klow,khigh
        do i = ilow,ihigh
        do j = jlow,jhigh
            if(lmask_3d(i,j,k))then
                ibin = int( (dbz_3d(i,j,k)-hist_low) / hist_intvl ) + 1
                if(ibin .ge. 1 .and. ibin .le. 10)then
                   ihist(ibin) = ihist(ibin) + 1
                   ihist_2d(k,ibin) = ihist_2d(k,ibin) + 1
                endif
            endif
        enddo ! j
        enddo ! i
        enddo ! k

!       Write full volume histogram
        do ibin = 1,nbins
            v1 = hist_low + hist_intvl * float(ibin-1)
            v2 = hist_low + hist_intvl * float(ibin  )
            write(lun_out,1)ibin,v1,v2,ihist(ibin)
 1          format(' bin/range/count = ',i10,f6.1,f6.1,i10)
        enddo ! ibin

!       Write histogram by level
        write(lun_out,*)
        write(lun_out,*)' Histogram by level'
        write(lun_out,*)'  Level    bin1      bin2      bin3      bin4'
     1                 ,'      bin5      bin6      bin7      bin8  '     
     1                 ,'    bin9     bin10'   
        do k = klow,khigh
            write(lun_out,2)k,(ihist_2d(k,ibin),ibin=1,nbins)
 2          format(i4,2x,10i10)
        enddo ! k

        return
        end
