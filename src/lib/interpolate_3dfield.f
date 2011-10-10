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

        subroutine interpolate_3dfield(field_3d_in,nx,ny,
     1          nlvl_in,rlvls_in,nlvl_out,rlvls_out,field_3d_out)

        real field_3d_in(nx,ny,nlvl_in)
        real field_3d_out(nx,ny,nlvl_out)
        real rlvls_in(nlvl_in)
        real rlvls_out(nlvl_out)

        do kout = 1,nlvl_out

!          Default if out of bounds
           kin_low = 1
           kin_high = 1
           frac_low = 1.0
           frac_high = 0.0

           do kk = 1,nlvl_in-1
             if(rlvls_in(kk  ) .ge. rlvls_out(kout) .and.
     1          rlvls_in(kk+1) .le. rlvls_out(kout)            )then
                  kin_low = kk
                  kin_high = kin_low + 1
                  frac_high = (rlvls_out(kout)  - rlvls_in(kk)) /
     1                        (rlvls_in(kk+1)   - rlvls_in(kk))
                  frac_low  = 1.0 - frac_high
             endif
           enddo ! kk

c         write(6,1)kout,kin_low,kin_high,frac_low,frac_high
c       1            ,rlvls_in(kin_low),rlvls_out(kout),rlvls_in(kin_high)
c1        format(1x,i5,3x,2i5,3x,2f8.3,2x,3f8.2)

          do j = 1,ny
          do i = 1,nx
           field_3d_out(i,j,kout) = field_3d_in(i,j,kin_low ) * frac_low
     1                           + field_3d_in(i,j,kin_high) * frac_high
          enddo ! i
          enddo ! j

        enddo ! kout

        return
        end
