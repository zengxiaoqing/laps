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
cdoc    This routine takes as input a full path filename and returns an
cdoc    index that points at the end of the directory portion of the pathname.
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
        integer lenf
C
        integer i, strlen
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
cdoc    This routine takes as input a full path filename and returns an
cdoc    index that points at the end of the filetime portion of the pathname.
C       Simple-minded algorithm just searches backwards for the first
C       occurance of a `.'.
C
C       Input/Output:
C
C       Name            Type    I/O     Description
C       ----            ---     --      -----------
C       c_fnames        char    I       file name.
C       lenf            I       O       index to end of filetime
C
C********************************************************************
C
        character c_fname*(*)
        integer lenf
C
        integer i, strlen
C
C****************************
C
        call s_len(c_fname,strlen)
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
cdoc    This routine receives a fortran string, and
cdoc    returns the number of characters in the string
cdoc    before the first "space" is encountered.  
C       It considers ascii characters 33 to 126 to be valid
C       characters, and ascii 0 to 32, and 127 to be "space"
C       characters.
C
C       Name            Type      I/O     Description
C       ----            ---       --      -----------
C       string          char       I       string
C       s_length        integer  O       valid number characters
C                                            in string

        implicit none

        character*(*)   string
        integer       s_length, i, len_str, aval
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
cdoc    This routine takes as input a full path filename and returns
cdoc    the length of the filename portion of the path.
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
        integer lenf
C
        integer i, strlen, i_char_len
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
cdoc    This routine takes as input a full path filename and returns
cdoc    the length of the file time length portion of the path.
C       Simple-minded algorithm uses other simple minded algorithms
C       found in this file to determine the length of the time portion
C       of the character string.
C

        character*(*) c_fname
        integer lend,lent,lenf

        call get_directory_length(c_fname,lend)
        call s_len(c_fname,len_fname)

        lenf = len_fname - lend

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
C#####################################################################
C
        subroutine get_filetime_type(c_fname,c20_type,leni,lent)
C
C*********************************************************************
C
cdoc    This routine takes as input a full path filename and returns
cdoc    the length and type of the file time length portion of the path.
C       Simple-minded algorithm uses other simple minded algorithms
C       found in this file to determine the length and type of the time portion
C       of the character string.
C
C       1998          Steve Albers

        character*(*) c_fname
        character*20 c20_type

        integer lend ! Directory length including the last 'slash'. 
        integer lenf ! Length of the filename excluding the path.
        integer lent ! Length of the filetime portion.
        integer leni ! Length of the initial non-filetime portion.
                       ! This is often but not always the directory length.

        call filter_non_numeric_fnames(c_fname,1,num_out,1
     1                                    ,istatus)

c initialize
        c20_type = 'unknown'

        if(num_out.eq.1)then

           call s_len(c_fname,len_fname)
           call get_directory_length(c_fname,lend)
           lenf = len_fname - lend

           if(lenf .eq. 20)then
            if(c_fname(lend+14:lend+14) .eq. '_')then
               c20_type = 'yyyyjjjhhmmss'                         ! RSA radar type
               leni = lend
               lent=13
               return
            endif
           endif

           if(lenf .eq. 13)then
            if(c_fname(lend+1:lend+5) .eq. 'raob.')then
               c20_type = 'yymmddhh'                              ! AFWA raob
               leni = lend+5
               lent=8
               return
c
c - repositioned to "lenf .ge. 13" switch below. J.Smart 3-20-00 
c           elseif(c_fname(lend+1:lend+2) .eq. 'nf'  )            !.or.
c     +             c_fname(lend+1:lend+2) .eq. 're' )then         ! Taiwan FA model
c              c20_type = 'ymmddhh'
c              leni = lend+2
c              lent = 7
c              return

            endif
           endif

           if(lenf .ge. 13)then
            if(c_fname(lend+9:lend+9) .eq. '_')then
               c20_type = 'yyyymmdd_hhmm'                         ! WFO type
               leni = lend
               lent=13
               return
            elseif(c_fname(lend+1:lend+2) .eq. 'nf'.or.
     +             c_fname(lend+1:lend+2) .eq. 're'.or.
     +             c_fname(lend+1:lend+2) .eq. 'gb'.or.
     +             c_fname(lend+1:lend+2) .eq. 'sb'.or.
     +             c_fname(lend+1:lend+2) .eq. 'gs')then
               c20_type = 'yyyymmddhh'                            !Taiwan/CWB nfs, gfs, or sb (tropical cyclone) model
               leni = lend+2
               lent = 10
               return
            endif
           endif

           if(lenf .eq. 16)then
            if(c_fname(lend+1:lend+4) .eq. 'temp')then
               c20_type = 'yymmddhh'                              ! CWB raob
               leni = lend+4
               lent=8
               return
            endif
           endif

           if(lenf .ge. 9)then ! assume 9 chars for time portion of filename
              c20_type = 'yyjjjhhmm'                             ! NIMBUS/LAPS
              leni = lend
              lent=9
              return
           endif

           c20_type = 'unknown'
           leni = lend
           lent = lenf

        endif

!       write(6,*)'Time portion of string = ',c20_type
!    1                                       ,c_fname(leni+1:leni+lent)
        
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
cdoc    This routine filters a fortran string of unprintable characters
C
C       Name            Type      I/O     Description
C       ----            ---       --      -----------
C       string          char       I       string

        implicit none

        character*(*)   string
        integer       i, len_str, aval
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


        function l_string_contains(string,substring,istatus)

cdoc    Returns boolean on whether 'string' contains 'substring'
cdoc    The 'string' and 'substring' should be free of blanks.

        logical l_string_contains

        character*(*) string
        character*(*) substring
 
        l_string_contains = .false.

        call s_len(string,len1)

        call s_len(substring,len2)

        if(len1 .ge. len2)then
            i_search_end = len1 - len2 + 1

            do i = 1,i_search_end
                if(string(i:i+len2-1) .eq. substring(1:len2))then
                    l_string_contains = .true.
                endif
            enddo ! i

            istatus = 1
            return

        else
            istatus = 0
            return

        endif

        end


        function l_parse(string1,string2)

cdoc    Returns boolean on whether 'string1' contains 'string2'
cdoc    Similar to 'l_string_contains' except that blanks are allowed

        logical l_parse

        character*(*) string1,string2

!       integer slen1,slen2

        len1 = len(string1)
        len2 = len(string2)

!       call s_len(string1,slen1)
!       call s_len(string2,slen2)

        l_parse = .false.

        if(len2 .gt. len1)return

        i_offset_max = len1-len2

        do i = 0,i_offset_max
            if(string1(i+1:i+len2) .eq. string2(1:len2))then
                l_parse = .true.
            endif ! match is found
        enddo ! i             

        return
        end        


        subroutine s_len2(string,len_string)

cdoc    This routine finds the length of the string counting intermediate
cdoc    blanks

        character*(*) string

        len1 = len(string)

        len_string = 0

        do i = 1,len1
            call s_len(string(i:i),len2)
            if(len2 .eq. 1)len_string = i
        enddo ! i

        return
        end

