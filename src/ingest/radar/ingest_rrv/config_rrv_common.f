      subroutine config_ingest_rrv_parms(nxv01,nyv01,nzv01,istatus)

      character nest7grid*150
      include 'ingest_rrv_dims.inc'
      include 'ingest_rrv_common.inc'

      namelist /ingest_rrv_nl/ path_to_raw_rrv,nxv01,nyv01,nzv01,
     +dxv01,dyv01

      istatus = 0

      call get_directory('nest7grid',nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/ingest_rrv.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,ingest_rrv_nl,err=901)
      close(1)
      istatus = 1
      return

 900  print*,'error opening file ',nest7grid
      return
 901  print*,'error reading ingest_rrv_nl in ',nest7grid
      write(*,ingest_rrv_nl)
      return
      end

