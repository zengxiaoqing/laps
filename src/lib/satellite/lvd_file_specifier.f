      subroutine lvd_file_specifier(c_type,
     &lvd_index,istatus)

      implicit none

      integer lvd_index
      integer istatus
      integer n
      character*3 c_type
      character*3 ctype

      istatus = 0
      n=index(c_type,' ')
      if(n.le.0)then
         ctype=c_type
         n=3
      else
         ctype=c_type(1:2)//' '
      endif

      if(ctype(1:n).eq.'vis')lvd_index=1

      if( (ctype(1:n).eq.'39 ').or.(ctype(1:n).eq.'i39').or.
     &(ctype(1:n).eq.'4u ') )lvd_index =2

      if( (ctype(1:n).eq.'wv ').or.(ctype(1:n).eq.'iwv'))
     &lvd_index=3

      if( (ctype(1:n).eq.'ir').or.
     &(ctype(1:n).eq.'ir '))lvd_index=4

c      if((ct.eq.'ir').or.(ct.eq.'ir '))lvd_index=4

      if((ctype(1:n).eq.'i11').or.
     &(ctype(1:n).eq.'11u'))lvd_index=4

      if( (ctype(1:n).eq.'12 ').or.(ctype(1:n).eq.'i12')
     &.or.(ctype(1:n).eq.'12u')   )lvd_index =5

      if(lvd_index.eq.0)then
         write(6,*)'Error getting lvd_index'
         istatus = -1
      endif

      return
      end
