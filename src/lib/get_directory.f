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

        subroutine get_directory(ext_in,directory,len_dir)

!       1993         S. Albers

!       This routine is a multipurpose routine that accepts an extension.
!       If the extension is a domain name, the location of the static file
!       directory is returned.

!       If the extension is a product file, the domain is tested, the extension
!       is modified for the domain if necessary, and the directory containing
!       the products is returned. This routine is designed for a UNIX
!       environment.

        character*(*) ext_in     ! input
        character*80 ext         ! local
        character*(*) directory  ! output
        integer*4 len_dir        ! output

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common
        character *200 laps_data_root

        call GETENV('LAPS_DATA_ROOT',laps_data_root) 

!       Test if the extension is the domain name. If so, then return the
!       directory containing the static files.

        call s_len(ext_in,len_ext_in)



        call s_len(grid_fnam_common,len_grid_fnam_common)
        call s_len(laps_data_root,len_lapsroot)

c        print *, laps_data_root(1:len_lapsroot),len_lapsroot

        if(laps_data_root(len_lapsroot:len_lapsroot).ne.'/') then
           len_lapsroot=len_lapsroot+1
           laps_data_root(len_lapsroot:len_lapsroot)='/'
        endif
        if(len_ext_in .eq. len_grid_fnam_common)then
          if(len_ext_in .ne. 0)then
            if(ext_in(1:len_ext_in) .eq.
     1         grid_fnam_common(1:len_ext_in))then
                directory = laps_data_root(1:len_lapsroot)//'static/'
                goto 999
            endif ! ext_in .eq. grid_fnam_common
          endif ! len_ext .ne. 0
        endif ! lens are =

        call downcase(ext_in,ext)

!       Test if the extension is for the thermo directory.
!       If so, then return the directory containing the static thermo files.
!       added direct test for static dir
        if(ext(1:3) .eq. 'dat' .or. ext(1:6) .eq. 'static')then
           directory = laps_data_root(1:len_lapsroot)//'static/'
           goto 999
        endif

        if(ext(1:3) .eq. 'etc' .or. ext(1:4).eq.'time')then
c           directory = laps_data_root(1:len_lapsroot)//'etc/'
           directory = laps_data_root(1:len_lapsroot)//'time/'
           goto 999
        endif

        if(ext(1:3) .eq. 'cdl')then
           directory = laps_data_root(1:len_lapsroot)//'data/cdl/'
           goto 999
        endif

        if(ext(1:3) .eq. 'log')then
           directory = laps_data_root(1:len_lapsroot)//'log/'
           goto 999
        endif

        if(ext(1:8) .eq. 'radioprd') then
           directory = laps_data_root(1:len_lapsroot)//'radioprd/'
           goto 999
        endif

        if(ext(1:7) .eq. 'nowrad6') then
           directory = laps_data_root(1:len_lapsroot)//'nowrad6/'
           goto 999
        endif


!       In this section, we assume the extension points to a product file
        call s_len(ext,len_ext)


        directory = laps_data_root(1:len_lapsroot)
     +       //'lapsprd/'//ext(1:len_ext)//'/'

!       Get Length of directory

 999    call s_len(directory,len_dir)

        return
        end
