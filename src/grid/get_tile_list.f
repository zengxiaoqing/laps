      subroutine get_tile_list(min_lat,max_lat,min_lon,max_lon
     1,maxtiles,isbego,iwbego,itilesize,ctiletype
     1,num_tiles_needed,ctile_name_list,iwoc1,iwoc2
     1,isoc1,isoc2,istatus)

      implicit  none

      real      max_lat,min_lat
      real      max_lon,min_lon

      integer   num_tiles_needed
      integer   isbego,iwbego
      integer   itilesize
      integer   itile_ns
      integer   iwoc0
      integer   iwoc,iwoc1,iwoc2
      integer   isoc,isoc1,isoc2
      integer   isocpt,isocpo
      integer   iwocpt,iwocpo,iwocph
      integer   iwocdif
      integer   istatus
      integer   maxtiles
      integer   maxtiles_loc

      parameter (maxtiles_loc = 1000)

      character*1 ctiletype
      character*3 nstitle
      character*4 ewtitle
      character*8 ctilenamelist(maxtiles_loc)
      character*(*) ctile_name_list(maxtiles)

      double precision r8term

      istatus = 1

      r8term=(min_lat-float(isbego))/float(itilesize)
      ISOC1=(INT(r8term+200.)-200)*itilesize+ISBEGO
      r8term=(min_lon-float(iwbego))/float(itilesize)
      IWOC1=(INT(r8term+400.)-400)*itilesize+IWBEGO
      r8term=(max_lat-float(isbego))/float(itilesize)
      ISOC2=(INT(r8term+200.)-200)*itilesize+ISBEGO
      r8term=(max_lon-float(iwbego))/float(itilesize)
      IWOC2=(INT(r8term+400.)-400)*itilesize+IWBEGO

      num_tiles_needed=0

      if(IWOC1.lt.-180)then
         if(itilesize.eq.180)then
            iwocdif=iwoc2-iwoc1
            IWOC1=IWOC1+IWOCDIF
            IWOC2=IWOC2+abs(IWOCDIF)
c        else
c           IWOC1=360+IWOC1
         endif
c     elseif(IWOC1.gt.180)then
c        IWOC1=360-IWOC1
      endif
      print*,'Noddy IWOC1, IWOC2 ',IWOC1,IWOC2
      do IWOC = IWOC1,IWOC2,itilesize

         IWOC0 = IWOC
         IF(IWOC.LT.-180)IWOC0=360+IWOC0
c        IF(IWOC.GT.+180)IWOC0=360-IWOC0

         IWOCPH=ABS(IWOC0)/100
         IWOCPT=(ABS(IWOC0)-IWOCPH*100)/10
         IWOCPO=ABS(IWOC0)-IWOCPH*100-IWOCPT*10
!
! TH: 8 Aug 2002 We now allow 180E longitudes (and greater). The only 
! time we want to assign W is when the longitude is less than 0.
!
         IF(IWOC0.GE.0.and.IWOC0.LT.180) THEN
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

            if(num_tiles_needed .gt. maxtiles_loc)then
                print*,'more tiles needed than array allocation A'
     1                ,num_tiles_needed,maxtiles_loc
                istatus = 0
                return
            endif

            ctilenamelist(num_tiles_needed)=nstitle//ewtitle

         enddo
      enddo

      if(num_tiles_needed.eq.3.and.maxtiles.eq.2)then
         num_tiles_needed = 2
      endif

      if(num_tiles_needed.le.maxtiles)then
         do itile_ns=1,num_tiles_needed
            ctile_name_list(itile_ns)=ctilenamelist(itile_ns)
         enddo
      else
         print*,'more tiles than array allocation B'
     1         ,num_tiles_needed,maxtiles
         istatus = 0
      endif

      return
      end
