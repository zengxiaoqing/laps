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
#ifndef solaris
        function i4_to_byte(i4_in)

        character i4_to_byte

        integer*4 i4_in,i4

        byte barg(4)
        equivalence (barg,i4)

        i4 = i4_in
        i4_to_byte = barg(4)

        return

        end
#endif

        function byte_to_i4(b_in)

        integer*4 i4,byte_to_i4
        character b_in

        character barg(4)
        equivalence (barg,i4)

        i4 = 0

        barg(4) = b_in
#ifndef hpux
        if(i4 .ge. 128)i4 = i4 - 256
#endif

        byte_to_i4 = i4

!       write(6,*)i4,byte_to_i4

        return

        end
#ifdef solaris
      real function trigd(val)
      real val
      real pi_180,invpi_180
      parameter (invpi_180=57.29577951308)
      parameter (pi_180=0.017453292519943)

      entry sind(val)
        trigd = sin(val*pi_180)
        return
      entry cosd(val)
        trigd = cos(val*pi_180)
        return
      entry tand(val)
        trigd = tan(val*pi_180)
        return
      entry atand(val)
        trigd = atan(val)*invpi_180
        return
      entry acosd(val)
        trigd = acos(val)*invpi_180
        return
      entry asind(val)
        trigd = asin(val)*invpi_180
        return
      end

      real function atan2d(val1,val2)
      real val1,val2
      real invpi_180
      parameter (invpi_180=57.29577951308)
      atan2d = atan2(val1,val2)*invpi_180
      return
      end

#endif
        subroutine open_lapsprd_file_append(lun,i4time,ext,istatus)

!       1997   Steve Albers

        character*(*)    ext
        character*150    directory
        character*13 filename13

!       Test for proper length of extension
        call s_len(ext,len_ext)
        if(len_ext .ne. 3)then
           write(6,*)' Error in open_lapsprd_file: ext has length'
     1               ,len_ext
           istatus = 0
           return
        endif

        call get_directory(ext,directory,len_dir)
#ifdef hpux
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))       
     1          ,status='unknown',ACCESS='APPEND',err=998)

#else
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))       
     1          ,status='unknown',POSITION='APPEND',err=998)
#endif
        go to 999

 998    write(6,*)
     1  ' Error in open_lapsprd_file_append: cannot open the file'       
        istatus = 0
        return

 999    istatus = 1
        return
        end ! open_lapsprd_file_append

