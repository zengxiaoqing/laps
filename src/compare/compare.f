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

      program diff_test

      integer imax,jmax,kmax,nlvl

      include 'lapsparms.cmn'
c
       parameter (imax=125)
       parameter (jmax=105)
c      parameter (imax=NX_L_MAX)
c      parameter (jmax=NY_L_MAX)
      parameter (kmax=200)
      parameter (nlvl=200)
      
      INTEGER	I4TIME,		!I4time of data
     1		KDIM,		!K dimension of DATA array
     1		LVL(nlvl),        !Level of each field (4 digit max)
     1		LVL_AVAIL(nlvl),        !Level of each field (4 digit max)
     1		FLAG,
     1		ERROR(2),
     1		INDX,
     1		I,J,K,start,
     1		ISTATUS
C
      REAL	DATA1(imax,jmax,kmax)	!Raw data to be written
      REAL	DATA2(imax,jmax,kmax)	!Raw data to be written
C
      CHARACTER*200	DIR_in  !Directory to be written to
      CHARACTER*200	DIR_out !Directory to be written to
      CHARACTER*31	EXT		!File name ext (up to 31 chars)
      CHARACTER*3	VAR(nlvl) 	        !3 letter ID of each field
      CHARACTER*3	LAPS_VAR_AVAIL(nlvl) 	!3 letter ID of each field
      character*19      VAR_AVAIL(nlvl)
      CHARACTER*4	LVL_COORD(nlvl)	!Vertical coordinate for each field
      CHARACTER*10	UNITS(nlvl)	!units of each field
      CHARACTER*125	COMMENT1(nlvl)	!Comments for each field
      CHARACTER*125	COMMENT2(nlvl)	!Comments for each field
      CHARACTER*9	GTIME
      CHARACTER*91	FILE_NAME
      CHARACTER*17	ASCTIME
      
      CHARACTER*4	VERSION
      CHARACTER*131	MODEL 		!Meteorological model in file
      CHARACTER*131	ORIGIN		!Location where file was created
      CHARACTER*11    LAPS_DOM_FILE   !Name of domain file e.g. NEST7GRID
      LOGICAL         l_packed_data

      logical l_pass

      character*20 filename
      character*9 a9_time
      character*3 var_last


      var_last = '   '
      idiff_msg_flag = 0
      diff_max_all = 0.
      diff_max_all_rel = 0.
      n_files = 0
      ndiff_all = 0
      l_pass = .true.

      write(6,*)' Enter 1 if comparing different machines, 0 if same'
      read(5,*)machine_factor

      write(6,*)' filename?'
5     read(5,1)filename
1     format(a)

      if(filename(1:3) .eq. 'end')then
!         write(6,*)' end'
          goto 999
      endif

      read(5,1,err=999)dir_in
      read(5,1,err=999)dir_out

      if(filename(1:6) .eq. 'static')then
          khmax = 4
          var(1) = 'LAT'
          var(2) = 'LON'
          var(3) = 'AVG'
          var(4) = 'LDF'
          ext = 'nest7grid'
          call rd_laps_static(dir_in,ext,imax,jmax,khmax,var,units
     1                     ,comment,data1,grid_spacing_m,istatus)

          if(istatus .ne. 1)then
              write(6,*)' Bad status reading 1st static file'
              stop
          endif

          call rd_laps_static(dir_out,ext,imax,jmax,khmax,var,units
     1                     ,comment,data2,grid_spacing_m,istatus)

          if(istatus .ne. 1)then
              write(6,*)' Bad status reading 2nd static file'
              stop
          endif

          thresh_write_pair = .01
          thresh_count_diff = .01

          ihmax = imax
          jhmax = jmax
 
      else
          n1 = index(filename,'.')
          n2 = index(filename,' ')

          a9_time = filename(1:9)
          ext = filename(n1+1:n2-1)

          call downcase(ext,ext)
          call cv_asc_i4time(a9_time,i4time)


          if(dir_in(1:3) .eq. 'end')then    
!             write(6,*)' end'
              goto 999
          endif

          do i = 1,nlvl
              var(i) = '   '
          enddo

          call READ_LAPS_HEADER(I4TIME,DIR_IN,EXT,IHMAX,JHMAX,KHMAX,
     1                         LAPS_DOM_FILE,ASCTIME,VERSION,
     1                         MODEL,ORIGIN,VAR,LVL,NUM_VARIABLES,
     1                         VAR_AVAIL,LAPS_VAR_AVAIL,NUM_LEVELS,
     1                         LVL_AVAIL,LVL_COORD,UNITS,
     1                         COMMENT1,L_PACKED_DATA,ISTATUS)


 
!         For this extension, set default values of:
!         thresh_write_pair, thresh_count_diff, or num_diff_field_thresh

          thresh_write_pair = 1e-05
          thresh_count_diff = 0.
          num_diff_field_thresh = 10

!         For this particular extension, set values of: 
!         thresh_write_pair, thresh_count_diff, or num_diff_field_thresh

          if(ext(1:3) .eq. 'lc3')then
              do i = 1,42
                  var(i) = 'LC3'
                  lvl(i) = i
              enddo
              thresh_write_pair = .01
              thresh_count_diff = .01

          elseif(ext(1:3) .eq. 'lcp')then
              thresh_write_pair = .01
              thresh_count_diff = .01

          elseif(ext(1:3) .eq. 'lcv')then
              thresh_write_pair = .01
              thresh_count_diff = .01

          elseif(ext(1:3) .eq. 'lwc')then
              thresh_write_pair = .01
              thresh_count_diff = .01

          elseif(ext(1:3) .eq. 'lw3')then
              thresh_write_pair = .01
              thresh_count_diff = .1

          elseif(ext(1:3) .eq. 'lwm')then
              thresh_write_pair = .01
              thresh_count_diff = .1

          elseif(ext(1:3) .eq. 'lps')then
              thresh_write_pair = .01
              thresh_count_diff = .1

          elseif(ext(1:3) .eq. 'lvd')then
              thresh_write_pair = 1.
              thresh_count_diff = 1.

          elseif(ext(1:3) .eq. 'vrc')then
              thresh_write_pair = 0.1
              thresh_count_diff = 0.1

          elseif(ext(1:3) .eq. 'l1s')then
              do i = 1,4
                  lvl(i) = 0
              enddo
              var(1) = 'R01'
              var(2) = 'RTO'
              var(3) = 'S01'
              var(4) = 'STO'
              thresh_write_pair = .0001
              thresh_count_diff = .0001

          elseif(ext(1:3) .eq. 'lps')then
              do i = 1,21
                  lvl(i) = 1150 - 50 * i            
              enddo

          elseif(ext(1:3) .eq. 'lsx')then
              thresh_write_pair = .01
              thresh_count_diff = .01
              num_diff_field_thresh = imax * jmax

          endif

!         Adjust values for this extension (dependent on machine) of: 
!         thresh_count_diff and num_diff_field_thresh

          if(machine_factor .eq. 0)then
              thresh_count_diff = 0
              num_diff_field_thresh = 0
          endif

          if(var(1) .eq. 'RH')var(1) = 'LHE'

          write(6,*)' Reading: ',dir_in,a9_time,'.',ext(1:3)

      
          call read_laps_data(i4time,dir_in,ext,ihmax,jhmax,khmax,
     1             khmax,var,lvl,lvl_coord,units,comment1,data1,
     1             istatus)

          if(istatus .ne. 1)then
            if(ext(1:3) .ne. 'lsx' .and. ext(1:3) .ne. 'liw')then
              if(l_pass)then
                  write(6,*)' READ ERROR: OVERALL CRITERIA FAILURE'
                  l_pass = .false.     
              endif
            endif
          endif

          write(6,*)' Reading: ',dir_out,a9_time,'.',ext(1:3)

      
          call read_laps_data(i4time,dir_out,ext,ihmax,jhmax,khmax,
     1             khmax,var,lvl,lvl_coord,units,comment2,data2,
     1             istatus)

          if(istatus .ne. 1)then
            if(ext(1:3) .ne. 'lsx' .and. ext(1:3) .ne. 'liw')then
              if(l_pass)then
                  write(6,*)' READ ERROR: OVERALL CRITERIA FAILURE'
                  l_pass = .false.     
              endif
            endif
          endif

      endif ! static file

!       write(6,*)' Hit RETURN to CONTINUE'
!	read(5,*)

        thresh_write_pair = thresh_write_pair * machine_factor

        diff_max_file_rel = 0.
        diff_max_file = 0.
        diff_max_var = 0.
        ndiff_file = 0
        nvar = 1
        n_levels = 0

        do k = 1,khmax

!       Test whether we should switch variables within this file
        if(k .gt. 1)then
            if(var(k) .ne. var(k-1))then
                if(n_levels .gt. 1)then
                    write(6,*)' Max diff for variable ',var(k-1)(1:3)
     1                      ,' =',diff_max_var
                    write(6,*)
                    diff_max_var = 0.
                    nvar = nvar + 1
                endif
                n_levels = 0
            endif
        endif

!       Test if new variable
        if(var(k)(1:3) .ne. var_last)then

            write(6,*)
            write(6,*)' New var = ',var(k)(1:3)

!           For this variable, set values of:
!           thresh_write_pair, thresh_count_diff, or num_diff_field_thresh

            if    (ext(1:3) .eq. 'lvd')thresh_write_pair = 1.0

            if    (ext(1:3) .eq. 'lt1' .and. var(k)(1:2) .eq. 'HT')then       
                thresh_count_diff = 2.0 * machine_factor
            elseif(ext(1:3) .eq. 'lt1' .and. var(k)(1:2) .eq. 'T3')then
                thresh_count_diff = .01 * machine_factor
            elseif(ext(1:3) .eq. 'lmt' .and. var(k)(1:3) .eq. 'LMT')then
                thresh_count_diff = 2.0 * machine_factor
            elseif(ext(1:3) .eq. 'lmt' .and. var(k)(1:3) .eq. 'LLR')then
                thresh_count_diff = 2.0 * machine_factor
            elseif(ext(1:3) .eq. 'lvd' .and. var(k)(1:3) .eq. 'SVN')then
                thresh_count_diff = 10.0 * machine_factor
            elseif(ext(1:3) .eq. 'lvd' .and. var(k)(1:3) .eq. 'ALB')then
                thresh_count_diff = 0.1 * machine_factor
                thresh_write_pair = 0.1                                   !JSmart addition 11-1-96
!           elseif(ext(1:3) .eq. 'lw3' .and. var(k)(1:2) .ne. 'OM')then
!               thresh_count_diff = .1 * machine_factor
            endif

            write(6,*)
     1      ' Threshold to write (first ten) grid point pairs = '
     1                                               ,thresh_write_pair
            write(6,*)
     1      ' Threshold to count lvl grid point differences   = '
     1                                               ,thresh_count_diff
            write(6,*)
     1      ' Max allowed count of lvl grid point differences = '
     1                                           ,num_diff_field_thresh
            write(6,*)


        endif

        n_levels = n_levels + 1

        diff_max_field = 0.
        diff_max_field_rel = 0.
        abs_value_max = 0.
        imaxd = 0
        jmaxd = 0
        iwrite = 0
        ndiff = 0
        inan = 0
        ndiff_msg = 0

	do i = 1,ihmax
           do j = 1,jhmax

c          if(
c!    1       data1(i,j,k) .le. r_min_normal()   .or.
c     1       data1(i,j,k) .ge. r_max_normal()   .or.         
c!    1       data2(i,j,k) .le. r_min_normal()   .or.         
c     1       data2(i,j,k) .ge. r_max_normal()            
c     1                                                    )then
              if(isnan(data1(i,j,k)).ne.0 .or. isnan(data2(i,j,k)).ne.0)
     +             then
                 iwrite = iwrite + 1
                 if(iwrite .le. 10)then
                    write(6,21)i,j,k,' Nan'
                 endif
                 inan = inan + 1
              else
                 diff     = abs(data1(i,j,k)-data2(i,j,k))

!           Test if one of the points is missing and the other isn't
            if(   (data1(i,j,k) .eq. r_missing_data .or.
     1             data2(i,j,k) .eq. r_missing_data      )        
     1                          .AND.
     1                     diff .gt. 0.                     )then

                ndiff_msg = ndiff_msg + 1
                idiff_msg_flag = 1

            else ! Both data points are non-missing

                diff_max_file = max(diff_max_file,diff)
                diff_max_var = max(diff_max_var,diff)
                if(diff .gt. diff_max_field)then
                    diff_max_field = diff
                    imaxd = i
                    jmaxd = j
                endif

            endif

            if(data1(i,j,k) .ne. r_missing_data)then
                abs_value_max = max(abs_value_max,abs(data1(i,j,k)))
            endif

            if(data2(i,j,k) .ne. r_missing_data)then
                abs_value_max = max(abs_value_max,abs(data2(i,j,k)))
            endif


            if(diff .gt. thresh_count_diff)then
                ndiff = ndiff + 1
                ndiff_file = ndiff_file + 1
                ndiff_all = ndiff_all + 1
            endif

            if(diff .gt. thresh_write_pair)then
                iwrite = iwrite + 1 
                if(iwrite .le. 10)then
                    write(6,21,err=22)i,j,k,data1(i,j,k),data2(i,j,k)
     1                                                  ,diff
21                  format(1x,3i5,2f12.6,f12.6)
22              endif
            endif
          endif ! Nan test
        enddo ! j
        enddo ! i

        if(comment1(k) .ne. comment2(k))then
            write(6,*)comment1(k)(1:80)
            write(6,*)comment2(k)(1:80)
        endif
        if(inan .gt. 0)write(6,*)' # of Nans = ',inan

        if(abs_value_max .gt. 0.)then
!          if(ndiff_msg .eq. 0)then
               diff_max_field_rel = diff_max_field / abs_value_max
!          endif
        else
           diff_max_field_rel = 0.
        endif

        diff_max_file_rel  = max(diff_max_file_rel,diff_max_field_rel)

 	write(6,*)' df_mx - fld #',k,' ',var(k)
     1  ,lvl(k),' abs/rel/#',diff_max_field,diff_max_field_rel
     1  ,ndiff,imaxd,jmaxd

        if(k .eq. khmax .and. nvar .gt. 1 
     1                  .and. n_levels .gt. 1)then
                write(6,*)
                write(6,*)' Max diff for variable ',var(k)(1:3),' ='
     1                      ,diff_max_var
                nvar = nvar + 1
        endif

        if(ndiff_msg .gt. 0)then
            write(6,*)
            write(6,*)' WARNING: # OF POINTS DIFFERING '
     1                ,'WRT MISSING DATA = ',
     1                  ndiff_msg
        endif

        if(ndiff + ndiff_msg .gt. num_diff_field_thresh)then
          if(l_pass)then
            write(6,*)' OVERALL CRITERIA FAILURE'
     1                           ,ndiff,ndiff_msg,num_diff_field_thresh
            l_pass = .false.     
          endif
        endif

        write(6,*)

        var_last = var(k)(1:3)

        enddo ! k


      if (istatus .ne. 1) write (6,*)'Error in readlapsdata'

	   write(6,*)' OVERALL FILE diff_max (',ext(1:3)
     1      ,') [abs/rel/#] = ',diff_max_file,diff_max_file_rel
     1               ,ndiff_file
           write(6,*)
           n_files = n_files + 1
           diff_max_all     = max(diff_max_all,diff_max_file)
           diff_max_all_rel = max(diff_max_all_rel,diff_max_file_rel)
c           goto 5
      
999     continue

        if(n_files .gt. 1)then
          if(l_pass)then
             write(6,*)' MAX difference (all files)  [abs/rel/#] = '
     1              ,diff_max_all,diff_max_all_rel
     1              ,ndiff_all,' PASSED'
          else
             write(6,*)' MAX difference (all files)  [abs/rel/#] = '
     1              ,diff_max_all,diff_max_all_rel
     1              ,ndiff_all,' FAILED'
          endif
           write(6,*)
        else
          if(l_pass)then
            write(6,*)' PASSED'
            write(6,*)
          else
            write(6,*)' FAILED'
            write(6,*)
          endif
        endif
 
 
        if(idiff_msg_flag .eq. 1)then
            write(6,*)
            write(6,*)' WARNING: DIFFERENCES WRT MISSING DATA DETECTED'
            write(6,*)
        endif

        stop

      end
      
