
        subroutine    Filter_non_numeric_fnames(c_fnames,
     1                           i_nbr_files_in,
     1                           i_nbr_files_out,
     1                           max_files,
     1                           istatus)
c
c JSmart  9-9-96
c       Routine takes the result of get_file_names and filters
c       it for non-numeric filenames. Test the first two characters
c       to determine if they are characters. A filename of this type
c       is disqualified.
c       - i_nbr_files_in is reduced. The number of
c       files returned (i_nbr_files_out) are files with only numerals
c       in the filename.
c
c       This routine is used in the source get_file_time.f that can
c       only deal with numeric filenames (particularly of the form
c       yyjjjhhmm and yyyymmdd_hhmm).
c

       implicit none

       integer*4 max_char
       parameter(max_char = 9)
       integer*4 max_files_filtered
       parameter(max_files_filtered=3000)

       integer*4 max_files
       integer*4 i_nbr_files_in
       integer*4 i_nbr_files_out
       integer*4 i_fnames_filtered
       integer*4 i_files_qualifying
       integer*4 lend
       integer*4 istatus
       integer*4 i,k,n

       logical found_qualifying_char

       character*(*) c_fnames(max_files)
       character*255 c_files_qualifying(max_files_filtered)
       character*1   c_qualifying_char(max_char)
       character*255 c_fnames_filtered(max_files_filtered)

       data c_qualifying_char/'1','2','3','4','5','6','7','8','9'/

       i_fnames_filtered = 0
       i_files_qualifying= 0

       do i=1,i_nbr_files_in

          call get_directory_length(c_fnames(i),lend)

          found_qualifying_char=.false.
          do k=1,max_char

             if(c_fnames(i)(lend+1:lend+1).eq.c_qualifying_char(k))
     1found_qualifying_char=.true.

          enddo

          if(found_qualifying_char)then
             i_files_qualifying = i_files_qualifying+1
             c_files_qualifying(i_files_qualifying)=c_fnames(i)
          else
             i_fnames_filtered=i_fnames_filtered+1
             c_fnames_filtered(i_fnames_filtered)=c_fnames(i)
          endif

       enddo

       if(i_fnames_filtered.gt.0)then
          write(6,*)'Filenames have been filtered'
          do i=1,i_fnames_filtered
             n=index(c_fnames_filtered(i),' ')
             write(6,*)i,' filtered fname = ',c_fnames_filtered(i)(1:n)
          enddo
          do i=1,i_files_qualifying
             c_fnames(i)=c_files_qualifying(i)
          enddo
c      else
c         write(6,*)'No filenames have been filtered'
       endif

       do i=1,i_files_qualifying
          c_fnames(i)=c_files_qualifying(i)
       enddo
       i_nbr_files_out = i_files_qualifying

       istatus=0
       return
       end
