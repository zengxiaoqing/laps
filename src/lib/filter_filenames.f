
        subroutine    Filter_non_numeric_fnames(c_fnames,
     1                           i_nbr_files_in,
     1                           i_nbr_files_out,
     1                           max_files,
     1                           istatus)
c
c JSmart  9-9-96
cdoc    Routine takes the result of 'get_file_names' and filters
cdoc    it for certain filename types. It tests the first two characters
c       to determine if they are allowed in i4time_fname_lp. A
c       filename of allowable types qualifies. If not then
c       i_nbr_files_in is reduced and this filename type is discarded from
c       the list. The number of
c       files returned (i_nbr_files_out) are files allowed in i4time_fname_lp.
c
c       This routine is used in the source get_file_time.f. LAPS can
c       only deal with specific filename types. See i4time_fname_lp for
c       more details of the types allowed.
c   
c      File type starting with the two characters 'nf' is now allowed since
c      this represents the FA (Taiwan) filename type for the the FA model.
c
c      File type starting with 'gb' and 'gs' now allowed since this is Taiwan
c      CWB global and tropical cyclone model backgrounds for lga (JS: 07/04)
c
       implicit none

       include 'lapsparms.for'

       integer   max_numeric_char
       parameter(max_numeric_char = 10)
       integer   max_2letter_strings
       parameter(max_2letter_strings = 5)
       integer   max_files_filtered
       parameter(max_files_filtered=MAX_RADAR_FILES)

       integer   max_files
       integer   i_nbr_files_in
       integer   i_nbr_files_out
       integer   i_fnames_filtered
       integer   i_files_qualifying
       integer   lend
       integer   istatus
       integer   i,j,k,n

       logical found_qualifying

       include 'filter_fnames.inc'   !contains lwant_filter_output logical

       character*(*) c_fnames(max_files)
       character*255 c_files_qualifying(max_files_filtered)
       character*1   c_qualifying_numeric_char(max_numeric_char)
       character*2   c_qualifying_2letter_string(max_2letter_strings)
       character*255 c_fnames_filtered(max_files_filtered)

       data c_qualifying_numeric_char/'0','1','2','3','4','5'
     &,'6','7','8','9'/
       data c_qualifying_2letter_string/'nf','re','te','gb','sb'/ 

       i_fnames_filtered = 0
       i_files_qualifying= 0
       istatus = 1

       do i=1,i_nbr_files_in

          call get_directory_length(c_fnames(i),lend)

          found_qualifying=.false.

          do k=1,max_numeric_char
             if(c_fnames(i)(lend+1:lend+1).eq.c_qualifying_
     1numeric_char(k))found_qualifying=.true.
          enddo

          do k=1,max_2letter_strings
             if(c_fnames(i)(lend+1:lend+2).eq.c_qualifying_
     12letter_string(k)(1:2))then
                do j=1,max_numeric_char
                   if(c_fnames(i)(lend+3:lend+3).eq.
     +c_qualifying_numeric_char(j))found_qualifying=.true.
                enddo
             endif
          enddo

          if(found_qualifying)then
             i_files_qualifying = i_files_qualifying+1
             c_files_qualifying(i_files_qualifying)=c_fnames(i)
          else
             i_fnames_filtered=i_fnames_filtered+1
             c_fnames_filtered(i_fnames_filtered)=c_fnames(i)
          endif

       enddo

       do i=1,i_nbr_files_in
          c_fnames(i)=' '
       enddo

       if(i_fnames_filtered.gt.0.and.lwant_filter_output)then
        do i=1,i_fnames_filtered
           n=index(c_fnames_filtered(i),' ')
           write(6,*)' filtered filename = ',c_fnames_filtered(i)(1:n)
        enddo
       endif
       if(i_fnames_filtered.gt.0)then
          do i=1,i_files_qualifying
             c_fnames(i)=c_files_qualifying(i)
          enddo
       endif

       do i=1,i_files_qualifying
          c_fnames(i)=c_files_qualifying(i)
       enddo
       i_nbr_files_out = i_files_qualifying

       return
       end
