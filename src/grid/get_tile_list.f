      subroutine get_tile_list(min_lat,max_lat,min_lon,max_lon
     1,maxtiles,isbego,iwbego,iblksizo,ctiletype,itilesize
     1,num_tiles_needed,ctile_name_list,istatus)

      implicit  none

      real      max_lat,min_lat
      real      max_lon,min_lon

      integer   num_tiles_needed
      integer   isbego,iwbego
      integer   iblksizo
      integer   itilesize
      integer   itile_ns
      integer   imxlat
      integer   imnlat
      integer   imxlon
      integer   imnlon
      integer   iwoc,iwoc1,iwoc2
      integer   isoc,isoc1,isoc2
      integer   isocpt,isocpo
      integer   iwocpt,iwocpo,iwocph
      integer   istatus
      integer   maxtiles


      character*1 ctiletype
      character*3 nstitle
      character*4 ewtitle
      character*8 ctilenamelist(500)
      character*(*) ctile_name_list(maxtiles)

      double precision r8term

      istatus = 1

      r8term=(min_lat-float(isbego))/float(iblksizo)
      ISOC1=(INT(r8term+200.)-200)*IBLKSIZO+ISBEGO
      r8term=(min_lon-float(iwbego))/float(iblksizo)
      IWOC1=(INT(r8term+400.)-400)*IBLKSIZO+IWBEGO
      r8term=(max_lat-float(isbego))/float(iblksizo)
      ISOC2=(INT(r8term+200.)-200)*IBLKSIZO+ISBEGO
      r8term=(max_lon-float(iwbego))/float(iblksizo)
      IWOC2=(INT(r8term+400.)-400)*IBLKSIZO+IWBEGO

      num_tiles_needed=0

      do IWOC = IWOC1,IWOC2,itilesize

         IWOCPH=ABS(IWOC)/100
         IWOCPT=(ABS(IWOC)-IWOCPH*100)/10
         IWOCPO=ABS(IWOC)-IWOCPH*100-IWOCPT*10
         IF(IWOC.GE.0
     1               .and. IWOC .ne. 180 
     1                                      )THEN
            WRITE(EWTITLE,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'E'
         ELSE
            WRITE(EWTITLE,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'W'
         ENDIF

         if(ewtitle(1:1).eq.' ')ewtitle(1:1)='0'
         if(ewtitle(2:2).eq.' ')ewtitle(2:2)='0'

c        ewtitle=ewtitle2//ewtitle1

         do ISOC = ISOC1,ISOC2,itilesize

            ISOCPT=ABS(ISOC)/10
            ISOCPO=ABS(ISOC)-ISOCPT*10

            IF(ISOC.GE.0)THEN
              WRITE(NSTITLE,'(2I1,A1)')ISOCPT,ISOCPO,'N'
            ELSE
              WRITE(NSTITLE,'(2I1,A1)')ISOCPT,ISOCPO,'S'
            ENDIF

            num_tiles_needed=num_tiles_needed+1
            ctilenamelist(num_tiles_needed)=nstitle//ewtitle

         enddo
      enddo

      if(num_tiles_needed.le.maxtiles)then
         do itile_ns=1,num_tiles_needed
            ctile_name_list(itile_ns)=ctilenamelist(itile_ns)
         enddo
      else
         print*,'more tiles than array allocation'
         istatus = 0
      endif

      return
      end
