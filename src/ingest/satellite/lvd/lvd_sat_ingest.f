      Program lvd_sat_ingest
c
c	3-6-97	J. Smart	New main driver for laps satellite image ingest process.
c				Purpose for having this new top-level driver is to accomdate
c				dynamic memory requirements.
c				1. Include remapping LUT as subroutine.
c				2. Include acquisition of domain parameters from
c				   static/ nest7grid.parms.
c				3. lvd_driver now a subroutine called by this main routine.
c
c       9-12-97 J. Smart        Dynamic array development. Renamed this to sub(routine)
c       2-19-98 J. Smart        Incorporate satellite_master.nl. This eliminates the need for
c                               all the separate files containing nav info for each sat and
c                               each format type.
c                               Made this the main program once again.
      Implicit None

      Integer nx_l
      Integer ny_l
      Integer nlinesir,nelemir
      Integer nlineswv,nelemwv
      Integer nlinesvis,nelemvis
      Integer i,j,k
      Integer ispec
      Integer nchannels
      Integer istatus

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      Character*3 channeltypes(maxchannel)
c
c ========================== START ==============================
c 
      call get_grid_dim_xy(nx_l,ny_l,istatus)
      if(istatus.ne.1)then
         write(6,*)'Error getting nx_l/ny_l'
         goto 1000
      else
         write(6,*)'LAPS nx_l and ny_l obtained'
      endif
c
      call config_satellite_lvd(istatus)
      if(istatus.ne.1)then
         write(*,*)'Error - Cannot continue'
         stop
      endif
c
c---------------------------------------------------------------
c Compute array dimensions for ir, vis, and wv.
c
      do k=1,maxsat
       if(isats(k).eq.1)then

       do j=1,maxtype
        if(itypes(j,k).eq.1)then

        nchannels=0
        do 4 i=1,maxchannel
         if(ichannels(i,j,k).eq.1)then
          nchannels=nchannels+1
          channeltypes(nchannels)=c_channel_types(i,j,k)
          call lvd_file_specifier(c_channel_types(i,j,k),ispec,istatus)
          if(istatus.ne.0)then
             write(6,*)'Error status returned from lvd_file_specifier'
             goto 1000
          endif

          goto(1,2,3,2,2)ispec

1         nelemvis=i_end_vis(j,k)-i_start_vis(j,k)+1
          nlinesvis=j_end_vis(j,k)-j_start_vis(j,k)+1
          goto 4
2         nelemir=i_end_ir(j,k)-i_start_ir(j,k)+1
          nlinesir=j_end_ir(j,k)-j_start_ir(j,k)+1
          goto 4
3         nelemwv=i_end_wv(j,k)-i_start_wv(j,k)+1
          nlineswv=j_end_wv(j,k)-j_start_wv(j,k)+1
         endif
4       enddo
 
        write(6,*)'lvd process information'
        write(6,*)'==============================='
        write(6,*)'Satellite ID: ',c_sat_id(k)
        write(6,*)'Satellite TYPE: ',c_sat_types(j,k)
        write(6,39)(channeltypes(i),i=1,nchannels)
39      format(1x,'Satellite CHANNELS:',5(1x,a3))

        write(6,*)'line/elem dimensions: '
        write(6,*)'VIS: ',nlinesvis,nelemvis
        write(6,*)'IR:  ',nlinesir,nelemir
        write(6,*)'WV:  ',nlineswv,nelemwv
c
c ================================================================
c
        call lvd_driver_sub(nx_l,ny_l,k,j,n_images,
     &                      nlinesir,nelemir,
     &                      nlineswv,nelemwv,
     &                      nlinesvis,nelemvis,
     &                      channeltypes,maxchannel,nchannels,
     &                      istatus)

        if(istatus.ne.1)then
           write(6,*)'Error status returned from lvd_driver_sub'
        else
           write(6,*)'Finished - lvd_driver_sub'
        endif

c =================================================================
        endif 
       enddo
       endif
      enddo

1000  stop
      end
