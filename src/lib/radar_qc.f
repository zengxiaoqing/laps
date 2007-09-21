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

        subroutine radar_qc(ni,nj,nk,ref_3d,istatus_qc)

        common /laps_diag/ no_laps_diag

        real ref_3d(ni,nj,nk)

        logical l_qc

        if(no_laps_diag .eq. 0)
     1  write(6,*)' Performing Radar QC'

        n_0_qc = 0
        n_10_qc = 0
        n_20_qc = 0
        n_30_qc = 0
        n_40_qc = 0
        n_50_qc = 0
        n_60_qc = 0
        n_70_qc = 0

        do j = 1,nj
        do i = 1,ni
          if(ref_3d(i,j,nk) .ge. 0.)then
            n_0_qc = n_0_qc + 1

            if(ref_3d(i,j,nk) .ge. 10.)then
              n_10_qc = n_10_qc + 1

              if(ref_3d(i,j,nk) .ge. 20.)then
                n_20_qc = n_20_qc + 1

                if(ref_3d(i,j,nk) .ge. 30.)then
                  n_30_qc = n_30_qc + 1

                  if(ref_3d(i,j,nk) .ge. 40.)then
                    n_40_qc = n_40_qc + 1

                    if(ref_3d(i,j,nk) .ge. 50.)then
                      n_50_qc = n_50_qc + 1

                      if(ref_3d(i,j,nk) .ge. 60.)then
                        n_60_qc = n_60_qc + 1

                        if(ref_3d(i,j,nk) .ge. 70.)then
                          n_70_qc = n_70_qc + 1

                        endif ! 70
                      endif ! 60
                    endif ! 50
                  endif ! 40
                endif ! 30
              endif ! 20
            endif ! 10
          endif ! 0
        enddo ! i
        enddo ! j

        npts = ni*nj

        frac0 = float(n_0_qc) / float(npts)
        frac10 = float(n_10_qc) / float(npts)
        frac20 = float(n_20_qc) / float(npts)
        frac30 = float(n_30_qc) / float(npts)
        frac40 = float(n_40_qc) / float(npts)
        frac50 = float(n_50_qc) / float(npts)
        frac60 = float(n_60_qc) / float(npts)
        frac70 = float(n_70_qc) / float(npts)

        if(frac40 .gt. .15)then
            l_qc = .false.
            istatus_qc = 0
        else
            l_qc = .true.
        endif

        if(no_laps_diag .eq. 0)
     1  write(6,1)frac0,frac10,frac20,frac30,frac40,frac50,frac60,frac70
     1,l_qc
1       format(1x,' %0-70',8f7.3,' QC=',l2)

        return
        end
