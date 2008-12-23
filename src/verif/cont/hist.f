
        subroutine radarhist(NX_L,NY_L,NZ_L,dbz_3d,lmask_3d,lun_out)
 
        real dbz_3d(NX_L,NY_L,NZ_L)
        logical lmask_3d(NX_L,NY_L,NZ_L)

        parameter (nbins=10)

        integer ihist(nbins)

        ihist = 0

!       Define area of interest
        ilow = 1
        ihigh = NX_L

        jlow = 1
        jhigh = NY_L

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
                endif
            endif
        enddo ! j
        enddo ! i
        enddo ! k

        do ibin = 1,nbins
            v1 = hist_low + hist_intvl * float(ibin-1)
            v2 = hist_low + hist_intvl * float(ibin  )
            write(lun_out,1)ibin,v1,v2,ihist(ibin)
 1          format(' bin/range/count = ',i10,f6.1,f6.1,i10)
        enddo ! ibin

        return
        end
