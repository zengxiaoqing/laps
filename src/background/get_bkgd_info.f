      subroutine get_background_info(len,bgpaths,bgmodels
     +     ,oldest_forecast,max_forecast_delta,use_analysis)
      implicit none
      integer maxbgmodels,len
      parameter (maxbgmodels=10)
      character*150 nest7grid
      character*150 bgpaths(maxbgmodels)
      integer bgmodels(maxbgmodels), len_dir
      integer oldest_forecast, max_forecast_delta
      logical use_analysis
      namelist /background_nl/bgpaths,bgmodels
     +         ,oldest_forecast,max_forecast_delta,use_analysis

      max_forecast_delta=6
      oldest_forecast=18
      use_analysis=.false.
      call get_directory('nest7grid',nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif
      nest7grid = nest7grid(1:len_dir)//'background.nl'

      open(1,file=nest7grid(1:len_dir+13),status='old',err=900)
      read(1,background_nl,err=901)
      close(1)
      return
 900  print*,'error opening file ',nest7grid
      stop
 901  print*,'error reading background_nl in ',nest7grid
      write(*,background_nl)
      stop
      end
