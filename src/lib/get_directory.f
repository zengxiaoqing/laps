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

        subroutine get_directory(ext_in,directory_out,len_dir)

!       1993         S. Albers

!       This routine is a multipurpose routine that accepts an extension.
!       If the extension is a domain name, the location of the static file
!       directory is returned.

!       Added wrfsi as a possible domain name (namelist) in addition to nes7grid.
!       Thus, made the data root a bit more generic. Retain functionality for
!       calls independent of get_laps_config; ie.,  still gets the env or command
!       line input
!       2000         J. Smart

!       If the extension is a product file, the domain is tested, the extension
!       is modified for the domain if necessary, and the directory containing
!       the products is returned. This routine is designed for a UNIX
!       environment.

        character*(*) ext_in         ! input
        character*80 ext             ! local
        character*(*) directory_out  ! output
        integer len_dir            ! output

        include 'grid_fname.cmn'

        character*201 directory

        call s_len(grid_fnam_common,len_dir)

cc        len_dir = index(grid_fnam_common,'/',.true.)
        if(len_dir.eq.0) then
           call get_config(istatus)
           if(istatus.ne.1)then
              print*,'ERROR: get_config failed'
              return
           endif
        endif

!       Test if the extension is the domain name. If so, then return the
!       directory containing the static files.

        call s_len(grid_fnam_common,len_grid_fnam)
        call s_len(ext_in,len_ext_in)
        call s_len(generic_data_root,len_root)

        if(len_root .gt. 200) then
          write(6,*)'get_directory ERROR: The env var for the DATA ROOT'
     1             ,' is too long; shorten to <200'
          stop
        endif 
c        print *, generic_data_root(1:len_root),len_root

        if(generic_data_root(len_root:len_root).ne.'/') then
           len_root=len_root+1
           generic_data_root(len_root:len_root)='/'
        endif
c
c both laps and wrfsi have static files in 'static' at the moment.
        if(len_ext_in .eq. len_grid_fnam)then
          if(len_ext_in .ne. 0)then
            if(ext_in(1:len_ext_in).eq.grid_fnam_common(1:len_grid_fnam)
     1)then
              if(ext_in(1:len_ext_in) .eq. 'nest7grid'
     1 .or.      ext_in(1:len_ext_in) .eq. 'wrfsi')then
                directory = generic_data_root(1:len_root)//'static/'
                goto 999
              endif
            endif ! ext_in .eq. grid_fnam_common
          endif ! len_ext .ne. 0
        endif ! lens are = ?

        call downcase(ext_in,ext)

!       Test if the extension is for the thermo directory.
!       If so, then return the directory containing the static thermo files.
!       added direct test for static dir
        if(ext(1:3) .eq. 'dat' .or. ext(1:6) .eq. 'static')then
           directory = generic_data_root(1:len_root)//'static/'
           goto 999
        endif

        if(ext(1:3) .eq. 'etc' .or. ext(1:4).eq.'time')then
c           directory = generic_data_root(1:len_root)//'etc/'
           directory = generic_data_root(1:len_root)//'time/'
           goto 999
        endif

        if(ext(1:3) .eq. 'cdl')then
           directory = generic_data_root(1:len_root)//'cdl/'
           goto 999
        endif

        if(ext(1:3) .eq. 'log')then
           directory = generic_data_root(1:len_root)//'log/'
           goto 999
        endif

        if(ext(1:4) .eq. 'root' .OR. ext(1:8) .eq. 'dataroot')then
           directory = generic_data_root(1:len_root)//'/'
           goto 999
        endif

!       In this section, we assume the extension points to a product file
        call s_len(ext,len_ext)


        directory = generic_data_root(1:len_root)
     +       //'lapsprd/'//ext(1:len_ext)//'/'

!       Get Length of directory

 999    call s_len(directory,len_dir)

!       Test length of input directory
        if(len_dir .gt. len(directory_out))then
            write(6,*)' get_directory ERROR: input directory is too'     
            write(6,*)' short, lengthen from ',len(directory_out)
     1               ,' to ',len_dir     
            stop

        elseif(len_dir .eq. 201)then
            write(6,*)' get_directory ERROR: LAPS_DATA_ROOT is too long'       
     1               ,' total directory length >= ',len_dir
            stop ! consider increasing max string lengths in this routine

        else ! len_dir fits all constraints
            directory_out = directory(1:len_dir)

        endif ! len_dir
        return
        end
