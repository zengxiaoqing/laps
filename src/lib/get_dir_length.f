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

C
        subroutine get_directory_length(c_fname,lenf)
C
C********************************************************************
C
C       This routine takes as input a full path filename and returns an
C       index that points at the end of the directory portion of the pathname.
C       Simple-minded algorithm just searches backwards for the first
C       occurance of a `/' (for UNIX) or a `]' (for VMS).
C
C       Input/Output:
C
C       Name            Type    I/O     Description
C       ----            ---     --      -----------
C       c_fnames        char    I       file name.
C       lenf            I       O       index to end of directory
C
C********************************************************************
C
        character c_fname*(*)
        integer*4 lenf
C
        integer*4 i, strlen
C
C****************************
C
        strlen = len(c_fname)
C
        i = strlen
        do while (i .gt. 0)
        if( (c_fname(i:i) .ne. ']')
     1.and. (c_fname(i:i) .ne. '/') )then
           i = i-1
        else
           goto 100
        endif
        enddo
C
100     lenf = i
C
        return
        end

C
        subroutine get_time_length(c_fname,lenf)
C
C********************************************************************
C
C       This routine takes as input a full path filename and returns an
C       index that points at the end of the filetime portion of the pathname.
C       Simple-minded algorithm just searches backwards for the first
C       occurance of a `.'.
C
C       Input/Output:
C
C       Name            Type    I/O     Description
C       ----            ---     --      -----------
C       c_fnames        char    I       file name.
C       lenf            I       O       index to end of directory
C
C********************************************************************
C
        character c_fname*(*)
        integer*4 lenf
C
        integer*4 i, strlen
C
C****************************
C
        strlen = len(c_fname)
C
        i = strlen
        do while (i .gt. 0)
        if (c_fname(i:i) .ne. '.')then
           i = i-1
        else
           goto 100
        endif
        enddo
C
100     lenf = i - 1
C
        return
        end

C
C#####################################################
C
        subroutine s_len(string,s_length)

C*********************************************************************
C
C       This routines receives a fortran string, and
C       returns the number of characters in the string
C       before the first "space" is encountered.  It
C       considers ascii characters 33 to 126 to be valid
C       characters, and ascii 0 to 32, and 127 to be "space"
C       characters.
C
C       Name            Type      I/O     Description
C       ----            ---       --      -----------
C       string          char       I       string
C       s_length        integer*4  O       valid number characters
C                                            in string

        implicit none

        character*(*)   string
        integer*4       s_length, i, len_str, aval
        logical         space

        space = .false.
        i = 1
        len_str = len(string)
        s_length = len_str      !default, used if string
                                !  is full, with no "spaces"
        do while ((i .le. len_str) .and. (.not. space))
          aval = ichar(string(i:i))
          if ((aval .lt. 33) .or. (aval .gt. 126)) then
            s_length = i - 1
            space = .true.
          endif
          i = i + 1
        enddo

        return
        end

C####################################################################
C
        subroutine get_fname_length(c_fname,lenf)
C
C********************************************************************
C
C       This routine takes as input a full path filename and returns
C       the length of the filename portion of the path.
C       Simple-minded algorithm just searches backwards for the first
C       occurance of a `/'. The number of characters searched up to
C       that point indicates the filename length.
C
C       Input/Output:
C
C       Name            Type    I/O     Description
C       ----            ---     --      -----------
C       c_fnames        char    I       file name.
C       lenf            I       O       length of filename part of c_fnames
C
C********************************************************************
C
        character c_fname*(*)
        integer*4 lenf
C
        integer*4 i, strlen, i_char_len
C
C****************************
C
        i_char_len = len(c_fname)
        call s_len(c_fname,strlen)
C
        i = i_char_len
        do while (i .gt. 0)
        if(c_fname(i:i) .ne. '/')then
           i = i-1
        else
           goto 100
        endif
        enddo
C
100     lenf = strlen - i
C
        return
        end
C
C#####################################################################
C
        subroutine get_filetime_length(c_fname,lent)
C
C*********************************************************************
C
C       This routine takes as input a full path filename and returns
C       the length of the file time length portion of the path.
C       Simple-minded algorithm uses other simple minded algorithms
C       found in this file to determine the length of the time portion
C       of the character string.
C

        character*(*) c_fname
        integer*4 lend,lent,lenf

        call get_directory_length(c_fname,lend)

        lenf = len(c_fname) - lend

        if(lenf .ge. 13 .and. c_fname(lend+9:lend+9).eq.'_')then
           lent=13

        elseif(lenf .ge. 9)then !assume 9 chars for time portion of filename
           lent=9

        else
           lent = lenf

        endif

c       write(6,*)'Time portion of string = ',c_fname(lend+1:lend+lent)       

        return
        end

C
C#####################################################
C
        subroutine filter_string(string)

C       It considers ascii characters 33 to 126 to be valid
C       characters, and ascii 0 to 32, and 127 to be "space"
C       characters.

C*********************************************************************
C
C       This routines filters a fortran string of unprintable characters
C
C       Name            Type      I/O     Description
C       ----            ---       --      -----------
C       string          char       I       string

        implicit none

        character*(*)   string
        integer*4       i, len_str, aval
        logical         space

        space = .false.
        len_str = len(string)

        do i = 1,len_str
          aval = ichar(string(i:i))
          if ((aval .lt. 33) .or. (aval .gt. 126)) then
            space = .true.
            string(i:i) = ' '
          endif
        enddo

        return
        end

