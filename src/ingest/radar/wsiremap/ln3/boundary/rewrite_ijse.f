      subroutine rewrite_ijse(istrtin,jstrtin,iendin,jendin,
     +istatus)
c
c routine rewrites the ln3.nl namelist file with new values
c for istart,jstart,iend,jend.
c
      implicit none

      integer istrtin,jstrtin,iendin,jendin
      integer istart,jstart,iend,jend
      integer istatus,iostat,lend
      integer lun
      integer msng_radar

      character filename*200
      character dir*150

      istatus=1   !return failed

      call get_ln3_parameters(msng_radar,
     +     istart,jstart,iend,jend,iostat)

      if(iostat.ne.0)then
         print*,'error reading ln3 parameters'
         return
      endif

      call get_directory('static',dir,lend)
      if(dir(lend:lend).ne.'/') then
        lend=lend+1
        dir(lend:lend)='/'
      endif


      filename=dir(1:lend)//'ln3.nl'
      lun=11
      open(lun,file=filename,status='old',form='formatted',
     +     err=99)
      rewind(lun)

      write(lun,*,err=950)'&ln3_nl'
      write(lun,10,err=950)msng_radar
      write(lun,11,err=950)istrtin
      write(lun,12,err=950)jstrtin
      write(lun,13,err=950)iendin
      write(lun,14,err=950)jendin
      write(lun,*,err=950)'/'

      close(lun)


10    format(1x,'MSNG_RADAR = ',i3,',')
11    format(1x,'ISTART = ',i4,',')
12    format(1x,'JSTART = ',i4,',')
13    format(1x,'IEND = ',i4,',')
14    format(1x,'JEND = ',i4,',')

      istatus=0
      return

99    print*,'Error opening ln3.nl'
      return
950   print*,'Error writing to ln3.nl'
      return
      end
