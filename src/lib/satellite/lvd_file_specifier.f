      subroutine lvd_file_specifier(c_type,
     &lvd_index,istatus)

      implicit none

      integer lvd_index
      integer istatus
      integer n
      character*3 c_type

      istatus = 0
      lvd_index=0
      n=index(c_type,' ')-1
      if(n.le.0)n=3

      if(c_type(1:n).eq.'vis')lvd_index=1

      if( (c_type(1:n).eq.'39').or.(c_type(1:n).eq.'i39').or.
     &    (c_type(1:n).eq.'4u') )
     &    lvd_index =2

      if( (c_type(1:n).eq.'wv').or.(c_type(1:n).eq.'iwv').or.
     &    (c_type(1:n).eq.'wvp') )
     &    lvd_index=3

      if( (c_type(1:n).eq.'ir').or.
     &    (c_type(1:n).eq.'ir ') )
     &    lvd_index=4

      if( (c_type(1:n).eq.'i11').or.
     &    (c_type(1:n).eq.'11u').or.
     &    (c_type(1:n).eq.'10p') )
     &    lvd_index=4

      if( (c_type(1:n).eq.'12').or.(c_type(1:n).eq.'i12').or.
     &    (c_type(1:n).eq.'12u') )
     &    lvd_index=5

      if(lvd_index.eq.0)then
         write(6,*)'Error getting lvd_index ',c_type
         istatus = -1
         stop ! for debugging
      endif

      return
      end
