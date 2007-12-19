
      subroutine write_orb_att(path,sat_name_in,num_atts,orb_att)
C
C**********************************************************************
C
        implicit  none
C
        character*(*)  path
        character*(*)  sat_name_in
        character*15   sat_name_down

        integer          num_atts
        real*8         orb_att(num_atts)
C
        integer        sat_name_len,
     1                 end_path,
     1                 noerror, error,
     1                 istatus,
     1                 cdl_len,
     1                 filename_len
C
        character*150  filename
        character*150  cdlfile
C
C ****  Begin
C
        noerror=1   !SUCCESS
        error=0     !ERROR

C ****  downcase sat_name_in for making file name
        call downcase(sat_name_in, sat_name_down)
        call s_len(sat_name_down,sat_name_len)

C ****  Make filename: path/satname_orbatt.dat

        call s_len(path,end_path)
        if (path(end_path:end_path) .eq. '/') then
        else
          path = path//'/'
          end_path = end_path + 1
        endif
        filename = 
     1path(1:end_path)//sat_name_down(1:sat_name_len)//'_orbatt.dat'
        call s_len(filename,filename_len) 
C
C **** make full cdl file name        
C
        call get_directory('cdl',cdlfile, cdl_len)
        if (cdlfile(cdl_len:cdl_len) .eq. '/') then
        else
          cdlfile = cdlfile(1:cdl_len)//'/'
          cdl_len = cdl_len + 1
        endif
        cdlfile = cdlfile(1:cdl_len)//'orb.cdl'
        cdl_len = cdl_len + 7
C
C **** call c program write_att_c to write file
C
        print*,'calling write_att_c'
        call write_att_c(filename, filename_len, cdlfile, cdl_len, 
     1                   sat_name_in, sat_name_len, num_atts, orb_att, 
     1                   istatus)

        if (istatus .eq. noerror)  goto 999  !SUCCESS
        if (istatus .eq. -1) goto 950 
        if (istatus .eq. -2) goto 960
        if (istatus .eq. -3) goto 970
        if (istatus .eq. 2) goto 980
C
C ****  Return normally.
C
        istatus=noerror
999     return
C
C ****  Error trapping.
C
950     write (6,*) 'Error opening netCDF file...write aborted.'
        istatus=error
        goto 999
C
960     write (6,*) 'Error creating netCDF file...write aborted.'
        istatus=error
        goto 999
C
970     write (6,*) 'Error writing netCDF file...write aborted.'
        istatus=error
        goto 999
C
980     write (6,*) 'Error in array dimensioning...write aborted.'
        istatus=error
        goto 999
C
992     continue
        istatus=error
        goto 999
C
        END

