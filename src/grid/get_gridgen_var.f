
      subroutine get_gridgen_var(nf,ngrids,var)

      implicit none
      integer nf,ngrids
      character*(*)  var(nf)

      if(ngrids.eq.6)then

         var(1)    = 'LAT'
         var(2)    = 'LON'
         var(3)    = 'AVG'
         var(4)    = 'LDF'
         var(5)    = 'USE'
         var(6)    = 'ZIN'

      elseif(ngrids.eq.10)then

         var(1)    = 'LAT'
         var(2)    = 'LON'
         var(3)    = 'LAB'
         var(4)    = 'LOB'
         var(5)    = 'LAC'
         var(6)    = 'LOC'
         var(7)    = 'AVG'
         var(8)    = 'LDF'
         var(9)    = 'USE' 
         var(10)   = 'ZIN'

      endif
      return
      end
