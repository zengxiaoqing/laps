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

        subroutine      downcase(input,output)

!       This routine only handles strings up to 500 characters long.

        character*(*)   input,
     1          output
cc        character*500   string

        integer*4       nchar,
!       1               lnblnk,
     1          l,
     1          len,
     1          chr,
     1          i

cc        string=input

!       nchar=lnblnk(string)

cc        if(string(500:500) .ne. ' ')then
cc            write(6,*)'String truncated to 500 characters.'
cc        endif

cc        do i = 500,1,-1
cc            if(string(i:i) .ne. ' ')then
cc                nchar = i
cc                go to 10
cc            endif
cc        enddo

10      continue

        l=len(output)
        do i=1,l
                output(i:i)=' '
        enddo

        call s_len(input,nchar)

        do i=1,nchar
cc                chr=ichar(string(i:i))
                chr=ichar(input(i:i))
                if (chr .ge. 65 .and. chr .le. 90) chr=chr+32
                output(i:i)=char(chr)
        enddo

        return

        end

