!dis   
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis    
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis   
!dis 

 subroutine snowfall (tsfc, prcpinc, preciptype, imax, jmax, &
                      snowinc, snowtot) 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--
!--   Name: Snow Accumulation Algorithm
!--
!--   Purpose
!--   =======
!--   This algorithm calculates incremental snow accumulation and
!--   total snow accumulation.
!--
!--   Method
!--   ======
!--   - If precip type is snow
!--      - calculate surface temperature on dot
!--      - calculate incremental precip on dot
!--      - if surface temp is >= 10 F
!--        - use a 10:1 snow/liquid ratio
!--      - if surface temp is < 10 F
!--        - use a 15:1 snow/liquid ratio
!--   - If precip type is not snow
!--      - snow accumulation is zero
!--
!--   Variables
!--   =========
!--   Name             Type      I/O     Definition
!--   ----             ----      ---     ----------
!--   imax             integer  input   grid dimension, i direction
!--   jmax             integer  input   grid dimension, j direction    
!--   prcpinc          real     input   3hr precip accum
!--                                     array(2-D)
!--   preciptype       integer  input   0 - is no precip
!--                                     1 - is rain
!--                                     2 - is freezing rain
!--                                     3 - is ice/mixed
!--                                     4 - is snow
!--   snowinc          real     output  3hr snow accum
!--                                         array(2-D)
!--   snowtot          real     output  total snow accum
!--                                        array(2-D)
!--   tsfc             real     input   surface temp array (2-D)
!--   tsfcf            real     local   surface temp in F
!--
!--   Updates
!--   =======
!--   20 Feb 98  Initial version..................Capt. John Lewis/DNXT
!--   10 Nov 98  Changed preciptype flag for snow 
!--              from 4 to 5; this corresponds 
!--              with changes made to the precip 
!--              type algorithm (wintprec.f)...Capt David Beberwyk/DNXT   
!--    5 Jan 99  Removed interpolations to dot grid as mmpost now
!--              operates on the cross grid........................DNXM
!--    4 Jan 01  Adapted by FSL for use with LAPS.. B. Shaw, NOAA/FSL
!--
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      use constants
      implicit none

      real,       external        :: fahren
      integer                     :: i
      integer,    intent(in)      :: imax
      integer                     :: j
      integer,    intent(in)      :: jmax
      real,       intent(in)      :: prcpinc    ( imax , jmax )
      integer,    intent(in)      :: preciptype ( imax , jmax )
      real,       intent(out)     :: snowinc    ( imax , jmax )
      real,       intent(inout)   :: snowtot    ( imax , jmax )
      real,       intent(in)      :: tsfc       ( imax , jmax )
      real                        :: tsfcf
      
!--------------------------------------------------------------------------
!     -- Begin the main double-do-loop
!-------------------------------------------------------------------------
        
      do j = 1, jmax
        do i = 1, imax    
            
!-----------------------------------------------------------------------
!--       Check if precipitation type is snow.
!-----------------------------------------------------------------------

          if (preciptype(i,j) == 5) then

            tsfcf = fahren(tsfc(i,j) - t0)

!-----------------------------------------------------------------------
!--         Liquid equivalent of snow depends of surface temperature.
!-----------------------------------------------------------------------

            if (tsfcf >= 10.0) then

              snowinc(i,j) = 10.0 * prcpinc(i,j)

            else

              snowinc(i,j) = 15.0 * prcpinc(i,j)

            end if

!-----------------------------------------------------------------------
!--         If precip type is not snow then snow accum is zero.
!-----------------------------------------------------------------------

          else

            snowinc(i,j) = 0.0

          end if
            
        end do
      end do

!-------------------------------------------------------------------
!--   Update snow total
!-------------------------------------------------------------------


      snowtot = snowtot + snowinc

      end subroutine snowfall
