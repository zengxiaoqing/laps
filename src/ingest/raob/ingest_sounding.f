
!      Driver for sounding ingest

!      Steve Albers      May-1999       Original Version

       character*150 path_to_raw_raob,path_to_raw_satsnd

       call get_laps_config('nest7grid',istatus)

       call get_snd_parms(path_to_raw_raob,path_to_raw_satsnd
     1                   ,istatus)
       if(istatus .ne. 1)goto999
 
       write(6,*)
       write(6,*)' Call ingest_raob'
       call ingest_raob(path_to_raw_raob)

       write(6,*)
       write(6,*)' Call ingest_satsnd'
       call ingest_satsnd(path_to_raw_satsnd)

 999   end


       subroutine get_snd_parms(path_to_raw_raob,path_to_raw_satsnd
     1                         ,istatus)

       character*150 path_to_raw_raob,path_to_raw_satsnd
       namelist /snd_nl/ path_to_raw_raob,path_to_raw_satsnd
 
       character*150 static_dir,filename
 
       call get_directory('nest7grid',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/snd.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,snd_nl,err=901)
       close(1)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading snd_nl in ',filename
       istatus = 0
       return

       end
