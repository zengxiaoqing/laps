       subroutine GSI_noaa1d2bufr
!***********************************************************************
! Program name : GSI_noaa1d2bufr(GSI_noaa1d2bufr.f90)
! 
! Description: To read in NOAA-NN satellite 1d data (AMSU-A,AMSU-B,AIRS3,AIRS4,or MHS) 
!              and transform to NCEP bufr format.
! 
! Environmental variable:
!      POLARSAT : directory of polar satellite input data.
!
! Files:
!      satellite.bufrtable : satellite bufr table (input file)
!      amsual1d_noaaNN_yyyymmdd_hhmm*.l1d: NOAA-NN AMSU-A 1D format file. (input file)
!      amsuabufrNN_yyyymmdd_hhmm or amsuabufr: NOAA-NN AMSU-A NCEP bufr data file. (output file)
!
!  or  amsubl1d_noaaNN_yyyymmdd_hhmm*.l1d: NOAA-NN AMSU-B 1D format file. (input file)
!      amsubbufrNN_yyyymmdd_hhmm or amsubbufr: NOAA-NN AMSU-B NCEP bufr data file. (output file)
!
!  or  mhsl1d_noaa18_yyyymmdd_hhmm*.l1d: NOAA-18 MHS 1D format file. (input file)
!      amsubbufr18_yyyymmdd_hhmm or amsubbufr: NOAA-18 MHS NCEP bufr data file. (output file)
!
!  or  hirsl1d_noaaNN_yyyymmdd_hhmm*.l1d: NOAA-NN HIRS 1D format file. (input file)
!      hirs3bufrNN_yyyymmdd_hhmm or hirs3bufr: NOAA-NN HIRS NCEP bufr data file. (output file)
!      (yyyymmdd_hhmm: year, month, day, hour and minute; NN is satellite id.)
!
! 1-D format Structure:
!          general information:        4*24 bytes
!          1-D Data Structure: depending on different instruments
! AMSUA-1D Data Structure: 
!          header record: 4*606 bytes
!          spare record: 2088 bytes
!          data record scan line 1: 4*1152 bytes
!          data record scan line 2: 4*1152 bytes
!          .....
!          data record scan line n: 4*1152 bytes
! AMSUB-1D or MHS-D Data Structure:
!          header record: 4*2278 bytes
!          spare record: 3080 bytes
!          data record scan line 1: 4*3072 bytes
!          data record scan line 2: 4*3072 bytes
!          .....
!          data record scan line n: 4*3072 bytes
! HIRS-1D Data Structure:
!          header record: 4*3823 bytes
!          spare record: 484 bytes
!          data record scan line 1: 4*3968 bytes
!          data record scan line 2: 4*3968 bytes
!          .....
!          data record scan line n: 4*3968 bytes
!
! NCEP NOAA satellite BUFR table-A:
!     NC021021: MSG TYPE 021-021 processed hirs-2 1B data (NOAA 14)
!     NC021022: MSG TYPE 021-022 processed msu-2  1B data (NOAA 14)
!     NC021023: MSG TYPE 021-023 processed amsu-a 1B data (NOAA 15-18)
!     NC021024: MSG TYPE 021-024 processed amsu-b 1B data (NOAA 15-17)
!     NC021025: MSG TYPE 021-025 processed hirs-3 1B data (NOAA 15-17)
!     NC021027: MSG TYPE 021-027 processed mhs    1B data (NOAA 18)
!     NC021028: MSG TYPE 021-028 processed hirs-4 1B data (NOAA 18)
!
! NCEP BUFR FORMAT:
!   amsu-a, amsu-b, mhs:
!-----------------------------------------------------------------------
!  NC021sss  | YEAR  MNTH  DAYS  HOUR  MINU  SECO  CLAT  CLON  SAID
!  NC021sss  | SIID  FOVN  LSQL  SAZA  SOZA  HOLS  HMSL  SOLAZI  BEARAZ
!  NC021sss  | "BRIT"xx
!  BRIT      | CHNM  TMBR
!-----------------------------------------------------------------------
!  where xx=15 for amsu-a,  =5 for amsu-b/mhs
!
!   hirs-3, hirs4:
!-----------------------------------------------------------------------
!  NC021sss  | YEAR  MNTH  DAYS  HOUR  MINU  SECO  CLAT  CLON  SAID
!  NC021sss  | SIID  FOVN  LSQL  SAZA  SOZA  HOLS  HMSL  SOLAZI  BEARAZ
!  NC021sss  | "BRIT"20
!  BRIT      | CHNM  TMBR
!-----------------------------------------------------------------------
! 
!   hirs-2, msu:
!-----------------------------------------------------------------------
!  NC021sss  | YEAR  MNTH  DAYS  HOUR  MINU  SECO  CLAT  CLON  SAID
!  NC021sss  | SIID  FOVN  LSQL  SAZA  SOZA  HOLS  HMSL  
!  NC021sss  | "BRIT"xx
!  BRIT      | CHNM  TMBR
!-----------------------------------------------------------------------
!  where xx=20 for hirs-2,  =4 for msu
!  for SIID:
!      =605, for hirs-2;    =606, for hirs-3;    =607, for hirs-4;
!      =623, for msu;       =570, for amsu-a:    =574, for amsu-b;
!      =203, for mhs
! 
! Called Function:  
!    rioreadbyte(getiofile.c),
!    getnoaafile,  dc_gen_inf,  dc_hirs,  dc_amsua,  dc_mhs
!    openbf, openmb, ufbseq, writsb, closbf (NCEP BUFR routine)
!    
! Date : 
!   Original     -- May  02, 2007 (Shiow-Ming Deng)  
!***********************************************************************

       implicit none

       character*80 path_satdat
       integer ipath,numfile,numcha(900)
       character*180 file(900)

       character*180 dir,sat_table
       integer len,len_sat_table
       logical sattab

       character filein*80,fileout*80
       integer isize,machine
       parameter( isize=90000000,machine=1 )
       integer ifile,ier,isize1
       byte buf(isize),genbyte(96)
       byte scanbyte(15872)

       integer i,j,k,kk,l,ll
       integer sateid,instrument,nline
       integer year,month,day,hour,minute,icheck
       integer siid,nbyte,xtrack,nchannel

       integer iscan,isqc
       real second,hmsl

       integer ifqc(90),isty(90)
       real rlat(90),rlon(90),ht(90),stzn(90),sozn(90)
       real staz(90),soaz(90)
       real brtmp_hirs(20,56),brtmp_amsua(15,30),brtmp_mhs(5,90)

       integer lnbufr,idate,iret,n
       character subset*8
       integer ndat
       parameter( ndat=100 )
       real*8 bufrf(ndat)

!-----------------------------------------------------------------------
!c  To get directory of satellite input data.
!c  Environmental variable: POLARSAT

       call getenv('POLARSAT',path_satdat)
       ipath=index(path_satdat,' ')
       if( ipath.le.1 )then
           print*,'No setting environmental variable: POLARSAT'
           return
       endif
       print*,'directory of input data: ',path_satdat(1:ipath-1)
       path_satdat(ipath:ipath)=char(0)

!-----------------------------------------------------------------------
!c  To get NOAA satellite 1D format input files.

       call getnoaafile(path_satdat,numfile,numcha,file,ier)
       if( (ier.ne.0).or.(numfile.lt.1) )then
           print*,'No any NOAA satellite hdf format file.'
           print*,'directory of input data: ',path_satdat(1:ipath-1)
           return
       endif

!-----------------------------------------------------------------------
!c  To get satellite bufr table: satellite.bufrtable

       call getenv('LAPS_DATA_ROOT',dir)
       len=index(dir,' ')
       if( len.le.1 )then
           print*,'Cannot get LAPS_DATA_ROOT directory.'
           return
       endif
       sat_table(1:len-1)=dir(1:len-1)
       sat_table(len:len+23)='/log/satellite.bufrtable'
       len_sat_table=len+23
       inquire(file=sat_table(1:len_sat_table),exist=sattab)
       if( .not.sattab )then
           print*,'satellite bufr table is not exist.'
           print*,'satellite bufr table: ',sat_table(1:len_sat_table)
           return
       endif 
       sat_table(len_sat_table+1:len_sat_table+1)=char(0)
       print*,'satellite bufr table: ',sat_table(1:len_sat_table)

      do 100 l=1,numfile
    
         filein(1:ipath-1)=path_satdat(1:ipath-1)
         filein(ipath:ipath)='/'
         ll=numcha(l)
         filein(ipath+1:ipath+ll)=file(l)(1:ll)
         filein(ipath+ll+1:ipath+ll+1)=char(0)
         print*,' '
         print*,'process file: ',filein(1:ipath+ll)
         print*,' '

!-----------------------------------------------------------------------
!c  To read in NOAA-NN 1D data bytes.

       isize1=isize
       call rioreadfile(filein,buf,isize1,ier)
       if( ier.ne.0 )then
           print*,'Cannot read in byte data.'
           print*,'file: ',filein(1:ipath+ll)
           go to 100
       endif
       if( isize1.lt.4608 )then
           print*,'This is not noaa-NN level 1d file.'
           print*,'isize1: ',isize1
           print*,'file: ',filein(1:ipath+ll)
           go to 100
       endif

       do k=1,ndat
          bufrf(k)=10.0e10
       enddo

!-----------------------------------------------------------------------
!c  To decode header record.

       do i=1,96
          genbyte(i)=buf(i)
       enddo 
       call dc_gen_inf(genbyte,machine,sateid,instrument,nline &
                      ,year,month,day,hour,minute,ier)
       if( ier.ne.0 )then
           print*,'Cannot decode 1D general information.'
           print*,'Input file may be not a noaa-NN 1D data.'
           print*,'input file: ',filein(1:ipath+ll)
           go to 100
       endif 
       if( instrument.eq.10 )icheck=4608*(nline+1) !amsu-a
       if( instrument.eq.11 )icheck=12288*(nline+1) !amsu-b
       if( instrument.eq.12 )icheck=12288*(nline+1) !mhs
       if( instrument.eq.5 )icheck=15872*(nline+1) !hirs
       if( icheck.ne.isize1 )then
           print*,'The size of input 1D data has question.'
           print*,'nline: ',nline,' icheck: ',icheck
           print*,'size of read in data: ',isize1,' (bytes)'
           print*,'input file: ',filein(1:ipath+ll)
           go to 100
       endif

!-----------------------------------------------------------------------
!c  To get output bufr format file name.

       if( instrument.eq.5 )then  !hirs
           write(fileout,'(8x,i2.2,1x,i4.4,2i2.2,1x,2i2.2)')sateid,year &
                                               ,month,day,hour,minute
           fileout(1:8)='hir3bufr'
           fileout(11:11)='_'
           subset(1:8)='NC021025'
           siid=606
           if( sateid.ge.18 )then
               subset(1:8)='NC021028'
               siid=607
           endif
           if( sateid.le.15 )then
               subset(1:8)='NC021021'
               siid=605
               fileout(4:4)='2'
            endif
           fileout(20:20)='_'
           fileout(25:25)=char(0)
           nbyte=15872
           xtrack=56
           nchannel=20
       endif
       if( instrument.eq.10 )then  !amsu-a
           write(fileout,'(9x,i2.2,1x,i4.4,2i2.2,1x,2i2.2)')sateid,year &
                                               ,month,day,hour,minute
           fileout(1:9)='amsuabufr'
           fileout(12:12)='_'
           fileout(21:21)='_'
           fileout(26:26)=char(0)
           siid=570
           subset(1:8)='NC021023'
           nbyte=4608
           xtrack=30
           nchannel=15
       endif
       if( instrument.eq.11 )then  !amsu-b
           write(fileout,'(9x,i2.2,1x,i4.4,2i2.2,1x,2i2.2)')sateid,year &
                                               ,month,day,hour,minute
           fileout(1:9)='amsubbufr'
           fileout(12:12)='_'
           fileout(21:21)='_'
           fileout(26:26)=char(0)
           siid=574
           subset(1:8)='NC021024'
           nbyte=12288
           xtrack=90
           nchannel=5
       endif
       if( instrument.eq.12 )then  !mhs
           write(fileout,'(9x,i2.2,1x,i4.4,2i2.2,1x,2i2.2)')sateid,year &
                                               ,month,day,hour,minute
           fileout(1:9)='amsubbufr'
           fileout(12:12)='_'
           fileout(21:21)='_'
           fileout(26:26)=char(0)
           siid=203
           subset(1:8)='NC021027'
           nbyte=12288
           xtrack=90
           nchannel=5
       endif


!-----------------------------------------------------------------------
!c  To write out bufr.

       print*,'output BUFR file name: ',fileout
       print*,' '

       open(11,file=fileout,form='unformatted')
       open(12,file=sat_table)
       lnbufr=11
       call openbf(lnbufr,'OUT',12)

       n=0
       do k=1,nline

!-----------------------------------------------------------------------
!c  To decode scan line data.

          do i=1,nbyte
             j=nbyte+nbyte*(k-1)+i
             scanbyte(i)=buf(j)
          enddo

!-----------------------------------------------------------------------
!c  To decode HIRS line data.

          if( instrument.eq.5 )then
              call dc_hirs(nbyte,scanbyte,machine,xtrack,iscan,isqc &
                     ,year,month,day,hour,minute,second,hmsl,ifqc,isty &
                     ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp_hirs,ier)
          endif 

!-----------------------------------------------------------------------
!c  To decode AMSU-A line data.

          if( instrument.eq.10 )then
              call dc_amsua(nbyte,scanbyte,machine,xtrack,iscan,isqc &
                     ,year,month,day,hour,minute,second,hmsl,ifqc,isty &
                     ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp_amsua,ier)
          endif 

!-----------------------------------------------------------------------
!c  To decode AMSU-B or MHS line data.

          if( (instrument.eq.11).or.(instrument.eq.12) )then
              call dc_mhs(nbyte,scanbyte,machine,xtrack,iscan,isqc &
                     ,year,month,day,hour,minute,second,hmsl,ifqc,isty &
                     ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp_mhs,ier)
          endif 

!-----------------------------------------------------------------------
!c  To write data into BUFR

          bufrf(1)=year
          bufrf(2)=month
          bufrf(3)=day
          bufrf(4)=hour
          bufrf(5)=minute
          bufrf(6)=second
          bufrf(9)=sateid+191
          bufrf(10)=siid
          bufrf(16)=hmsl
          idate=year*1000000+month*10000+day*100+hour
          do i=1,xtrack
             if( ifqc(i).eq.0 )then
                 bufrf(7)=rlat(i)
                 bufrf(8)=rlon(i)
                 bufrf(11)=i
                 bufrf(12)=isty(i)
                 bufrf(13)=stzn(i)
                 bufrf(14)=sozn(i)
                 bufrf(15)=ht(i)
                 bufrf(17)=soaz(i)
                 bufrf(18)=staz(i)
                 kk=18
                 if( instrument.eq.5 )then
                     do j=1,nchannel
                        kk=kk+1
                        bufrf(kk)=j
                        kk=kk+1
                        bufrf(kk)=brtmp_hirs(j,i)
                     enddo 
                 endif
                 if( instrument.eq.10 )then
                     do j=1,nchannel
                        kk=kk+1
                        bufrf(kk)=j
                        kk=kk+1
                        bufrf(kk)=brtmp_amsua(j,i)
                     enddo 
                 endif
                 if( (instrument.eq.11).or.(instrument.eq.12) )then
                     do j=1,nchannel
                        kk=kk+1
                        bufrf(kk)=j
                        kk=kk+1
                        bufrf(kk)=brtmp_mhs(j,i)
                     enddo 
                 endif
                 call openmb(lnbufr,subset,idate)
                 call ufbseq(lnbufr,bufrf,ndat,1,iret,subset)
                 call writsb(lnbufr)
                 n=n+1
             endif
          enddo
       enddo

       call closbf(lnbufr)
       close(lnbufr)
       close(12)  
       print*,'number of output BUFR: ',n
       print*,' '

!-----------------------------------------------------------------------
!c  To write out bufr again.

       if( instrument.eq.5 )then  !hirs
           fileout(1:8)='hir3bufr'
           fileout(9:9)=char(0)
       endif
       if( instrument.eq.10 )then  !amsu-a
           fileout(1:9)='amsuabufr'
           fileout(10:10)=char(0)
       endif
       if( (instrument.eq.11).or.(instrument.eq.12) )then  !amsu-b
           fileout(1:9)='amsubbufr'
           fileout(10:10)=char(0)
       endif
       print*,'output BUFR file name: ',fileout
       print*,' '

       open(11,file=fileout,form='unformatted')
       open(12,file=sat_table)
       lnbufr=11
       call openbf(lnbufr,'OUT',12)

       n=0
       do k=1,nline

!-----------------------------------------------------------------------
!c  To decode scan line data.

          do i=1,nbyte
             j=nbyte+nbyte*(k-1)+i
             scanbyte(i)=buf(j)
          enddo

!-----------------------------------------------------------------------
!c  To decode HIRS line data.

          if( instrument.eq.5 )then
              call dc_hirs(nbyte,scanbyte,machine,xtrack,iscan,isqc &
                     ,year,month,day,hour,minute,second,hmsl,ifqc,isty &
                     ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp_hirs,ier)
          endif 

!-----------------------------------------------------------------------
!c  To decode AMSU-A line data.

          if( instrument.eq.10 )then
              call dc_amsua(nbyte,scanbyte,machine,xtrack,iscan,isqc &
                     ,year,month,day,hour,minute,second,hmsl,ifqc,isty &
                     ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp_amsua,ier)
          endif 

!-----------------------------------------------------------------------
!c  To decode AMSU-B or MHS line data.

          if( (instrument.eq.11).or.(instrument.eq.12) )then
              call dc_mhs(nbyte,scanbyte,machine,xtrack,iscan,isqc &
                     ,year,month,day,hour,minute,second,hmsl,ifqc,isty &
                     ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp_mhs,ier)
          endif 

!-----------------------------------------------------------------------
!c  To write data into BUFR

          bufrf(1)=year
          bufrf(2)=month
          bufrf(3)=day
          bufrf(4)=hour
          bufrf(5)=minute
          bufrf(6)=second
          bufrf(9)=sateid+191
          bufrf(10)=siid
          bufrf(16)=hmsl
          idate=year*1000000+month*10000+day*100+hour
          do i=1,xtrack
             if( ifqc(i).eq.0 )then
                 bufrf(7)=rlat(i)
                 bufrf(8)=rlon(i)
                 bufrf(11)=i
                 bufrf(12)=isty(i)
                 bufrf(13)=stzn(i)
                 bufrf(14)=sozn(i)
                 bufrf(15)=ht(i)
                 bufrf(17)=soaz(i)
                 bufrf(18)=staz(i)
                 kk=18
                 if( instrument.eq.5 )then
                     do j=1,nchannel
                        kk=kk+1
                        bufrf(kk)=j
                        kk=kk+1
                        bufrf(kk)=brtmp_hirs(j,i)
                     enddo 
                 endif
                 if( instrument.eq.10 )then
                     do j=1,nchannel
                        kk=kk+1
                        bufrf(kk)=j
                        kk=kk+1
                        bufrf(kk)=brtmp_amsua(j,i)
                     enddo 
                 endif
                 if( (instrument.eq.11).or.(instrument.eq.12) )then
                     do j=1,nchannel
                        kk=kk+1
                        bufrf(kk)=j
                        kk=kk+1
                        bufrf(kk)=brtmp_mhs(j,i)
                     enddo 
                 endif
                 call openmb(lnbufr,subset,idate)
                 call ufbseq(lnbufr,bufrf,ndat,1,iret,subset)
                 call writsb(lnbufr)
                 n=n+1
             endif
          enddo
       enddo

       call closbf(lnbufr)
       close(lnbufr)
       close(12)  
       print*,'number of output BUFR: ',n
       print*,' '

 100   continue

       return
       end

       subroutine getnoaafile(path_satdat,numfile,numcha,file,ier)
!***********************************************************************
! Subroutine/Function : getnoaafile
!
! Usage :
!    call getnoaafile(path_satdat,numfile,numcha,file,ier)
!
! Description      : To get NOAA satellinte 1D format input file names.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I    path_satdat C*80        directory of NOAA satellite input data.
!    O     numfile    integer     total number of input file names.
!    O     file(900)  C*180       NOAA satellite 1-D file names.
!    O    numcha(900) int array
!    O     ier        integer     error message.
!                                 =0, success;  =1, failure.
!
! Modules Called : 
!   getfilenames (getiofile.c)
!***********************************************************************

       implicit none
       character*80 path_satdat
       character*180 file(900)
       integer numfile,numcha(900),ier

       character*900000 cdat
       integer i,isize,ii,j,i1,i2,jj

       ier=0
       numfile=0
       isize=900000
       do i=1,isize
          cdat(i:i)=char(0)
       enddo
       call getfilenames(path_satdat,isize,cdat,ier)
       if( ier.ne.0 )return

       ii=1
       do i=isize,1,-1
          if( cdat(i:i).ne.char(0) )then
              ii=i
              go to 10
          endif
       enddo
       return
 10    continue
       if( ii.lt.5 )then
           ier=1
           return
       endif

       j=1
 20    continue
       j=j+1
       if( j.gt.ii )go to 30
       if( cdat(j:j).eq.' ' )then
           i2=j-1
           jj=i2-i1+1
           if( cdat(i2-3:i2).eq.'.l1d' )then
               numfile=numfile+1
               if( numfile.gt.900 )then
                   print*,'Error in subroutine: getnoaafile'
                   print*,'Please modify parameter file(900) and numcha(900).'
                   print*,'numfile: ',numfile
                   ier=1
                   return
               endif
               file(numfile)(1:jj)=cdat(i1:i2)
               numcha(numfile)=jj
           endif
           i1=j+1
       else
           if( j.eq.ii )then
               i2=ii
               jj=i2-i1+1
               if( cdat(i2-3:i2).eq.'.l1d' )then
                   numfile=numfile+1
                   if( numfile.gt.900 )then
                       print*,'Error in subroutine: getnoaafile'
                       print*,'Please modify parameter file(900) and numcha(900).'
                       print*,'numfile: ',numfile
                       ier=1
                   return
                   endif
                   file(numfile)(1:jj)=cdat(i1:i2)
                   numcha(numfile)=jj
               endif
               go to 30
           endif
       endif
       go to 20
 30    continue

       return
       end

       subroutine dc_mhs(nbt,scan,machine,xtrack,iscan,isqc &
                    ,year,month,day,hour,minute,second,hmsl,ifqc,isty &
                    ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp,ier)
!***********************************************************************
! Subroutine/Function : dc_mhs
!
! Usage
!    call dc_mhs(nbt,head,machine,xtrack,iscan,isqc
!   1       ,year,month,day,hour,minute,second,hmsl,ifqc,isty
!   2       ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp,ier)
!
! Description      : To decode scan line data of AMSUB-1d or MHS data.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     nbt        integer     total number of input byte. (=12288)
!    I     scan(nbt)  byte        scan line data array.
!    I     machine    integer     machine index,
!                                 =0, for HP, SGI machine.
!                                 .not. 0, for DEC, PC machine.
!    I     xtrack     integer     number of FOV. (=90)
!    O     iscan      integer     scan line number.
!    O     isqc       integer     scan line quality flags.
!                                 =0, good
!                                 =1, bad
!    O     year       integer     scan line year.
!    O     month      integer     scan line month.
!    O     day        integer     scan line day.
!    O     hour       integer     scan line hour.
!    O     minute     integer     scan line minute.
!    O     second     real        scan line second.
!    O     hmsl       real        satellite altitude above refernce ellipsoid. (meters)
!    O   ifqc(xtrack) int array   FOV quality flag.
!                                 =0, good
!                                 =1, bad
!    O   isty(xtrack) int array   surface type of FOV.
!                                 =0, land
!                                 =1, sea
!    O   rlat(xtrack) real array  latitude in degrees of FOV.
!    O   rlon(xtrack) real array  longitude in degrees of FOV.
!    O     ht(xtrack) real array  surface height of FOV in metres.
!    O   stzn(xtrack) real array  satellite zenith angle in degrees of FOV.
!    O   sozn(xtrack) real array  solar zenith angle in degrees of FOV.
!    O   staz(xtrack) real array  satellite azimuth angle in degrees of FOV.
!    O   soaz(xtrack) real array  solar azimuth angle in degrees of FOV.
!    O    brtmp(5,90) real array  scene brightness temperature in K ch.16 - ch.20
!    O     ier        integer     error message.
!                                 =0, success.
!                                 =1, failure.
!
! Modules Called : by2int4,  jumonday1
! 
! Note: the line scan data structure
!
!       scan line information:  4*9 bytes
!       navigation:             4*723 bytes
!       earth observations:     4*1800 bytes
!       pre-processing output:  4*540 bytes
!
!  data    field  
!  type    description
!  ....    .................................................
!                     SCAN LINE INFORMATION (4*9 bytes)
!  I  4    scan line number
!  I  4    scan line year
!  I  4    scan line day of year
!  I  4    scan line UTC time of day in milliseconds
!  I  4    quality indicator bit field -- in all of the follwong,
!          if the bit is on (i.e., if it is set to 1) then the statement is true.
!          Otherwise it is false.
!          bit 31:  do not use scan for product generation
!          bit 30:  time sequence error detected with this scan (see below)
!          bit 29:  data gap precedes this scan
!          bit 28:  no calibration (see below)
!          bit 27:  no earth location (see below)
!          bit 26:  first good time following a clock update
!          bit 25:  instrument status changed with this scan
!          bit 24-0: spare <zero fill>
!  I  4    scan line quality flags -- if bit is on (=1) then true
!          Time Problem Code: (All bits off implies the scan time is as expected.)
!          bit 31-24: spare <zero fill>
!          bit 23:  time field is bad but can probaly be inferred from the previous
!                   good time.
!          bit 22:  time field is bad and cannot be inferred from the previous good time.
!          bit 21:  this record starts a sequence that is inconsistent with previous times
!                   (i.e., there is a time discontinuity).  This may or may not be
!                   associated with a spacecraft clock update.
!          bit 20:  start of a sequence that apparently repeates scan times that have
!                   been prevously accepted.
!          bit 19-16: spare <zero fill>
!          Calibration Problem Code: (Note these bits compliment the channel indicators;
!            all bits set to 0 indicates normal calibration.)
!          bit 15:  Scan line was not calibrated because of bad time.
!          bit 14:  Scan line was calibrated using fewer than the preferred number of scan
!                   lines because of proximity to start or end of data set or to a data gap.
!          bit 13:  Scan line was not calibrated because of bad or insufficient PRT data.
!          bit 12:  Scan line was calibrated but with marginal PRT data.
!          bit 11:  Some uncalibrated channels on this scan. (See channel indicators.)
!          bit 10:  Uncalibrated due to instrument mode.
!          bit 09:  Questionable calibration because of antenna position error of space view.
!          bit 08:  Questionable calibration because of antenna position error of blackbody.
!          Earth Location Problem Code: (all bits set to 0 implies the earth location
!            was normal)
!          bit 07:  Not earth located because of bad time.
!          bit 06:  Earth location questionable because of questionable time code.
!                   (See time problem flags above.)
!          bit 05:  Earth location questionable -- only marginal agreement with
!                   reasonableness check.
!          bit 04:  Earth location questionable -- fails reasonableness check.
!          bit 03:  Earth location questionable because of antenna position check.
!          bit 02-0: spare <zero fill>
!  I  4    mixer chan 18-20 instrument temperature (K*100)
!  I  4*2  spare
!  ....    .................................................
!                       NAVIGATION (4*723 bytes)
!  I  4    10000*(latitude in degrees of position 1)
!  I  4    10000*(longitude in degrees of position 1)
!  I  4    10000*(latitude in degrees of position 2)
!  I  4    10000*(longitude in degrees of position 2)
!          ....
!  I  4    10000*(latitude in degrees of position 90)
!  I  4    10000*(longitude in degrees of position 90)
!  I  4    surface height of position 1 in metres
!  I  4    surface type of position 1 (0=sea, 1=mixed, 2=land)
!  I  4    surface height of position 2 in metres
!  I  4    surface type of position 2 (0=sea, 1=mixed, 2=land)
!          ....
!  I  4    surface height of position 90 in metres
!  I  4    surface type of position 90 (0=sea, 1=mixed, 2=land)
!  I  4    100*(local zenith angle in degrees of position 1)
!  I  4    100*(local azimuth angle in degrees of position 1)
!  I  4    100*(solar zenith angle in degrees of position 1)
!  I  4    100*(solar azimuth angle in degrees of position 1)
!  I  4    100*(local zenith angle in degrees of position 2)
!  I  4    100*(local azimuth angle in degrees of position 2)
!  I  4    100*(solar zenith angle in degrees of position 2)
!  I  4    100*(solar azimuth angle in degrees of position 2)
!          ....
!  I  4    100*(local zenith angle in degrees of position 90)
!  I  4    100*(local azimuth angle in degrees of position 90)
!  I  4    100*(solar zenith angle in degrees of position 90)
!  I  4    100*(solar azimuth angle in degrees of position 90)
!  I  4    10*(satellite altitude above refernce ellipsoid. km)
!  I  4*2  spare
!  ....    .................................................
!                    EARTH OBSERVATIONS (4*1800 bytes)
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-A ch.1
!          (missing data indicator is -999999)
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-A ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-A ch.15
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-B ch.16
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-B ch.17
!          ....
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-B ch.20
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-A ch.1
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-A ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-A ch.15
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-B ch.16
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-B ch.17
!          ....
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-B ch.20
!          ........
!  I  4    100* scene brightness temp. in K. FOV 90, AMSU-A ch.1
!  I  4    100* scene brightness temp. in K. FOV 90, AMSU-A ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 90, AMSU-A ch.15
!  I  4    100* scene brightness temp. in K. FOV 90, AMSU-B ch.16
!  I  4    100* scene brightness temp. in K. FOV 90, AMSU-B ch.17
!          ....
!  I  4    100* scene brightness temp. in K. FOV 90, AMSU-B ch.20
!  ....    .................................................
!                   PRE-PROCESSING OUTPUT (4*540 bytes)
!  I  4    Pre-processing quality flags for FOV 1: (all bits off implies acceptable data)
!          bit 31:  spare <zero fill>
!          bit 30:  set if AMSU-B used secondary calibration from nearest neighbour AMSU-A:
!          bit 29:  set if AMSU-A used secondary calibration
!          bit 28:  set if AMSU-A data missing
!          bit 27:  maximum probability scheme cloud flag
!          bit 26:  scattering test (only set over the sea)
!          bit 25:  logistic precipitation probability test
!          bit 24:  Grody light rainfall test
!          bit 23:  mismatch between AMSU-A/B 89GHz values
!          bit 22:  mismatch between surface type from topography dataset and from 
!                   pre-processing
!          bit 21:  spare
!          bit 20:  set if 89GHz channel greatly different from that in surrounding fovs
!          bit 19:  scattering test (only set over the sea) - using AMSU-B 89GHz channel
!          bit 18:  mismatch between AMSU-A/B 89GHz values
!          bit 17-1: spare <zero fill>
!          bit 00:  set if AMSU-B data missing
!  I  4    nearest-neightbour AMSU-A estimated surface type for FOV 1:
!          1 = Bare young ice (i.e. new ice, no snow)
!          2 = Dry land (i.e. dry. with or without vegetation)
!          3 = Dry snow (i.e. snow with water less than 2%, over land)
!          4 = Multi-year ice (i.e. old ice with dry snow cover)
!          5 = Sea (i.e. open water, no islands, ice-free, WS < 14m/s)
!          6 = Wet forest (i.e. established forest with wet canopy)
!          7 = Wet land (i.e. non-forestedland with a wet suface)
!          8 = Wet snow (i.e. water content > 2%, over land or ice)
!          9 = Desert
!  I  4    cost function from PPASURF surface identification for FOV 1
!  I  4    scattering index (recalculated with AMSU-B 89GHz) for FOV 1
!  I  4*2  spare
!  I  4    Pre-processing quality flags for FOV 2
!  I  4    nearest-neightbour AMSU-A estimated surface type for FOV 2
!  I  4    cost function from PPASURF surface identification for FOV 2
!  I  4    scattering index (recalculated with AMSU-B 89GHz) for FOV 2
!  I  4*2  spare
!          ....
!  I  4    Pre-processing quality flags for FOV 90
!  I  4    nearest-neightbour AMSU-A estimated surface type for FOV 90
!  I  4    cost function from PPASURF surface identification for FOV 90
!  I  4    scattering index (recalculated with AMSU-B 89GHz) for FOV 90
!  I  4*2  spare
!
!***********************************************************************
       implicit none
       integer nbt,xtrack
       byte scan(nbt),by4(4)
       integer machine
       integer iscan,isqc,year,month,day,hour,minute,ier
       real second,hmsl
       integer ifqc(xtrack),isty(xtrack)
       real rlat(xtrack),rlon(xtrack),ht(xtrack)
       real stzn(xtrack),sozn(xtrack),staz(xtrack),soaz(xtrack)
       real brtmp(5,90)
       integer i,j,k,kk,ij,jday,mn,ihour
       integer in4

       ier=0
       if( nbt.ne.12288 )then
           print*,'Error in subroutine: dc_mhs.'
           print*,'input total byte number is error.'
           print*,'nbt(=12288): ',nbt
           ier=1
           return
       endif
       if( xtrack.ne.90 )then
           print*,'Error in subroutine: dc_mhs.'
           print*,'input total xtrack number is error.'
           print*,'xtrack(=90): ',xtrack
           ier=1
           return
       endif

!-----------------------------------------------------------------------
!c  To decode scan line information.

       do i=1,4
          by4(i)=scan(i)
       enddo
       call by2int4(machine,by4,in4)
       iscan=in4

       do i=1,4
          j=i+4
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       year=in4

       do i=1,4
          j=i+8
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       jday=in4
       call jumonday1(year,jday,month,day)

       do i=1,4
          j=i+12
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       mn=in4/60000
       second=0.001*(in4-60000*mn)
       hour=mn/60
       minute=mn-60*hour

       do i=1,4
          j=i+16
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       i=in4/2
       isqc=in4-2*i

!-----------------------------------------------------------------------
!c  To decode navigation.

       do k=1,90
          do i=1,4
             j=36+8*(k-1)+i
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          rlat(k)=0.0001*in4
          do i=1,4
             j=36+8*(k-1)+i+4
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          rlon(k)=0.0001*in4
       enddo
       do k=1,90
          do i=1,4
             j=756+8*(k-1)+i
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          ht(k)=in4
          do i=1,4
             j=756+8*(k-1)+i+4
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          isty(k)=-9
          if( in4.eq.2 )isty(k)=0
          if( in4.eq.0 )isty(k)=1
          if( in4.eq.1 )isty(k)=2
          if( (isty(k).ne.0).and.(isty(k).ne.1) )then
              isty(k)=1
              if( ht(k).gt.0 )isty(k)=0
          endif
       enddo
       do k=1,90
          do i=1,4
             j=1476+16*(k-1)+i
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          stzn(k)=0.01*in4
          do i=1,4
             j=1476+16*(k-1)+i+4
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          staz(k)=0.01*in4
          do i=1,4
             j=1476+16*(k-1)+i+8
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          sozn(k)=0.01*in4
          do i=1,4
             j=1476+16*(k-1)+i+12
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          soaz(k)=0.01*in4
       enddo
       do i=1,4
          j=2916+i
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       hmsl=100.*in4

!-----------------------------------------------------------------------
!c  To decode earth observations.

       do k=1,90
          ij=2928+80*(k-1)
          ifqc(k)=0
          do kk=1,5
             do i=1,4
                j=ij+60+4*(kk-1)+i
                by4(i)=scan(j)
             enddo
             call by2int4(machine,by4,in4)
             brtmp(kk,k)=0.01*in4
             if( in4.lt.0 )ifqc(k)=1
          enddo
       enddo
       if( isqc.eq.0 )then
           ij=0
           do k=1,90
              ij=ij+ifqc(k)
           enddo 
           if( ij.eq.90 )isqc=1
       endif

       return
       end

       subroutine dc_amsua(nbt,scan,machine,xtrack,iscan,isqc &
                    ,year,month,day,hour,minute,second,hmsl,ifqc,isty &
                    ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp,ier)
!***********************************************************************
! Subroutine/Function : dc_amsua
!
! Usage
!    call dc_amsua(nbt,scan,machine,xtrack,iscan,isqc
!   1       ,year,month,day,hour,minute,second,hmsl,ifqc,isty
!   2       ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp,ier)
!
! Description      : To decode scan line data of AMSUA-1d data.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     nbt        integer     total number of input byte. (=4608) 
!    I     scan(nbt)  byte        scan line data array.
!    I     machine    integer     machine index,
!                                 =0, for HP, SGI machine.
!                                 .not. 0, for DEC, PC machine.
!    I     xtrack     integer     number of FOV. (=30)
!    O     iscan      integer     scan line number.
!    O     isqc       integer     scan line quality flags.
!                                 =0, good
!                                 =1, bad
!    O     year       integer     scan line year.
!    O     month      integer     scan line month.
!    O     day        integer     scan line day.
!    O     hour       integer     scan line hour.
!    O     minute     integer     scan line minute.
!    O     second     real        scan line second.
!    O     hmsl       real        satellite altitude above refernce ellipsoid. (meters)
!    O   ifqc(xtrack) int array   FOV quality flag.
!                                 =0, good
!                                 =1, bad
!    O   isty(xtrack) int array   surface type of FOV.
!                                 =0, land
!                                 =1, sea
!                                 =2, mixed
!    O   rlat(xtrack) real array  latitude in degrees of FOV.
!    O   rlon(xtrack) real array  longitude in degrees of FOV.
!    O     ht(xtrack) real array  surface height of FOV in metres.
!    O   stzn(xtrack) real array  satellite zenith angle in degrees of FOV.
!    O   sozn(xtrack) real array  solar zenith angle in degrees of FOV.
!    O   staz(xtrack) real array  satellite azimuth angle in degrees of FOV.
!    O   soaz(xtrack) real array  solar azimuth angle in degrees of FOV.
!    O   brtmp(15,30) real array  scene brightness temperature in K ch.1 - ch.15
!    O     ier        integer     error message.
!                                 =0, success.
!                                 =1, failure.
!
! Modules Called : by2int4,  jumonday1
!
! Note: the scan line data structure
!
!       scan line information:  4*11 bytes
!       navigation:             4*243 bytes
!       earth observations:     4*600 bytes
!       pre-processing output:  4*210 bytes
!       spare bytes:            4*88 bytes 
!
!  data    field  
!  type    description
!  ....    .................................................
!                      SCAN LINE INFORMATION (4*11 bytes)
!  I  4    scan line number
!  I  4    scan line year
!  I  4    scan line day of year
!  I  4    scan line UTC time of day in milliseconds
!  I  4    quality indicator bit field -- in all of the follwong,
!          if the bit is on (i.e., if it is set to 1) then the statement is true.
!          Otherwise it is false.
!          bit 31:  do not use scan for product generation
!          bit 30:  time sequence error detected with this scan (see below)
!          bit 29:  data gap precedes this scan
!          bit 28:  no calibration (see below)
!          bit 27:  no earth location (see below)
!          bit 26:  first good time following a clock update
!          bit 25:  instrument status changed with this scan
!          bit 24-0: spare <zero fill>
!  I  4    scan line quality flags -- if bit is on (=1) then true
!          Time Problem Code: (All bits off implies the scan time is as expected.)
!          bit 31-24: spare <zero fill>
!          bit 23:  time field is bad but can probaly be inferred from the previous
!                   good time.
!          bit 22:  time field is bad and cannot be inferred from the previous good time.
!          bit 21:  this record starts a sequence that is inconsistent with previous times
!                   (i.e., there is a time discontinuity).  This may or may not be 
!                   associated with a spacecraft clock update.
!          bit 20:  start of a sequence that apparently repeates scan times that have
!                   been prevously accepted.
!          bit 19-16: spare <zero fill>
!          Calibration Problem Code: (Note these bits compliment the channel indicators;
!            all bits set to 0 indicates normal calibration.)
!          bit 15:  Scan line was not calibrated because of bad time.
!          bit 14:  Scan line was calibrated using fewer than the preferred number of scan
!                   lines because of proximity to start or end of data set or to a data gap.
!          bit 13:  Scan line was not calibrated because of bad or insufficient PRT data.
!          bit 12:  Scan line was calibrated but with marginal PRT data.
!          bit 11:  Some uncalibrated channels on this scan. (See channel indicators.)
!          bit 10:  Uncalibrated due to instrument mode.
!          bit 09:  Questionable calibration because of antenna position error of space view.
!          bit 08:  Questionable calibration because of antenna position error of blackbody.
!          Earth Location Problem Code: (all bits set to 0 implies the earth location 
!            was normal)
!          bit 07:  Not earth located because of bad time.
!          bit 06:  Earth location questionable because of questionable time code.
!                   (See time problem flags above.)
!          bit 05:  Earth location questionable -- only marginal agreement with 
!                   reasonableness check.
!          bit 04:  Earth location questionable -- fails reasonableness check.
!          bit 03:  Earth location questionable because of antenna position check.
!          bit 02-0: spare <zero fill>
!  I  4    AMSU-A1 instrument RF shelf temperature (K*100)
!  I  4    AMSU-A2 instrument RF shelf temperature (K*100)
!  I  4    AMSU-B instrument mixer scan 18-20 temperature (K*100)
!  I  4*2  spare
!  ....    .................................................
!                       NAVIGATION (4*243 bytes)
!  I  4    10000*(latitude in degrees of position 1)
!  I  4    10000*(longitude in degrees of position 1)
!  I  4    10000*(latitude in degrees of position 2)
!  I  4    10000*(longitude in degrees of position 2)
!          ....
!  I  4    10000*(latitude in degrees of position 30)
!  I  4    10000*(longitude in degrees of position 30)
!  I  4    surface height of position 1 in metres
!  I  4    surface type of position 1 (0=sea, 1=mixed, 2=land)
!  I  4    surface height of position 2 in metres
!  I  4    surface type of position 2 (0=sea, 1=mixed, 2=land)
!          ....
!  I  4    surface height of position 30 in metres
!  I  4    surface type of position 30 (0=sea, 1=mixed, 2=land)
!  I  4    100*(local zenith angle in degrees of position 1)
!  I  4    100*(local azimuth angle in degrees of position 1)
!  I  4    100*(solar zenith angle in degrees of position 1)
!  I  4    100*(solar azimuth angle in degrees of position 1)
!  I  4    100*(local zenith angle in degrees of position 2)
!  I  4    100*(local azimuth angle in degrees of position 2)
!  I  4    100*(solar zenith angle in degrees of position 2)
!  I  4    100*(solar azimuth angle in degrees of position 2)
!          ....
!  I  4    100*(local zenith angle in degrees of position 30)
!  I  4    100*(local azimuth angle in degrees of position 30)
!  I  4    100*(solar zenith angle in degrees of position 30)
!  I  4    100*(solar azimuth angle in degrees of position 30)
!  I  4    10*(satellite altitude above refernce ellipsoid. km)
!  I  4*2  spare
!  ....    .................................................
!                   EARTH OBSERVATIONS (4*600 bytes)
!  I  4    100* scene brightness temp. in K. FOV 1, ch.1
!          (missing data indicator is -999999)
!  I  4    100* scene brightness temp. in K. FOV 1, ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 1, ch.15
!  I  4    100* scene brightness temp. in K. FOV 1. AMSU-B ch.16
!  I  4    100* scene brightness temp. in K. FOV 1. AMSU-B ch.17
!          ....
!  I  4    100* scene brightness temp. in K. FOV 1. AMSU-B ch.20
!  I  4    100* scene brightness temp. in K. FOV 2, ch.1
!  I  4    100* scene brightness temp. in K. FOV 2, ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 2. AMSU-B ch.20
!          ........
!  I  4    100* scene brightness temp. in K. FOV 30, ch.1
!  I  4    100* scene brightness temp. in K. FOV 30, ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 30. AMSU-B ch.20
!  ....    .................................................
!                   PRE-PROCESSING OUTPUT (4*210 bytes)
!  I  4    Pre-processing quality flags for FOV 1: (all bits off implies acceptable data)
!          bit 31:  AMSU BTs considered contaminated. due e.g. to precip or surface type
!          bit 30:  set if AMSU-A used secondary calibration
!          bit 29:  set if AMSU-B used secondary calibration
!          bit 28:  set if AMSU-B missing due to insufficient good data
!          bit 27:  maximum probability scheme cloud flag
!          bit 26:  scattering test (only set over the sea)
!          bit 25:  logistic precipitation probability test
!          bit 24:  Grady light rainfall test
!          bit 23:  mismatch between AMSU-A/B 89GHz values
!          bit 22:  mismatch between surface type from topography data and from pre-processing
!          bit 21-1: space <zero fill>
!          bit 00:  set if AMSU-A data missing
!  I  4    estimated AMSU-A surface type, for FOV 1:
!          1 = Bare young ice (i.e. new ice, no snow)
!          2 = Dry land (i.e. dry. with or without vegetation)
!          3 = Dry snow (i.e. snow with water less than 2%, over land)
!          4 = Multi-year ice (i.e. old ice with dry snow cover)
!          5 = Sea (i.e. open water, no islands, ice-free, WS < 14m/s)
!          6 = Wet forest (i.e. established forest with wet canopy)
!          7 = Wet land (i.e. non-forestedland with a wet suface)
!          8 = Wet snow (i.e. water content > 2%, over land or ice)
!          9 = Desert
!  I  4    cost function from surface identification (PPASURF), for FOV 1
!  I  4    scattering index (PPASCAT), for FOV1
!  I  4    logistic precipitation probability (PPCROSBY), for FOV 1
!  I  4*2  spare
!  I  4    Pre-processing quality flags for FOV 2
!  I  4    estimated AMSU-A surface type, for FOV 2
!  I  4    cost function from surface identification for FOV 2
!  I  4    scattering index for FOV 2
!  I  4    logistic precipitation probability for FOV 2
!  I  4*2  spare
!          ....
!  I  4    Pre-processing quality flags for FOV 30
!  I  4    estimated AMSU-A surface type, for FOV 30
!  I  4    cost function from surface identification for FOV 30
!  I  4    scattering index for FOV 30
!  I  4    logistic precipitation probability for FOV 30
!  I  4*2  spare
!  ....    .................................................
!                           SPARE (4*88 bytes)
!
!***********************************************************************
       implicit none
       integer nbt,xtrack
       byte scan(nbt),by4(4)
       integer machine
       integer iscan,isqc,year,month,day,hour,minute,ier
       real second,hmsl
       integer ifqc(xtrack),isty(xtrack)
       real rlat(xtrack),rlon(xtrack),ht(xtrack)
       real stzn(xtrack),sozn(xtrack),staz(xtrack),soaz(xtrack)
       real brtmp(15,30)
       integer i,j,k,kk,ij,jday,mn,ihour
       integer in4

       ier=0
       if( nbt.ne.4608 )then
           print*,'Error in subroutine: dc_amsua.'
           print*,'input total byte number is error.'
           print*,'nbt(=4608): ',nbt
           ier=1
           return
       endif
       if( xtrack.ne.30 )then
           print*,'Error in subroutine: dc_amsua.'
           print*,'input total xtrack number is error.'
           print*,'xtrack(=30): ',xtrack
           ier=1
           return
       endif

!-----------------------------------------------------------------------
!c  To decode scan line information.

       do i=1,4
          by4(i)=scan(i)
       enddo
       call by2int4(machine,by4,in4)
       iscan=in4

       do i=1,4
          j=i+4
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       year=in4

       do i=1,4
          j=i+8
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       jday=in4
       call jumonday1(year,jday,month,day) 

       do i=1,4
          j=i+12
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       mn=in4/60000
       second=0.001*(in4-60000*mn)
       hour=mn/60
       minute=mn-60*hour

       do i=1,4
          j=i+16
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       i=in4/2
       isqc=in4-2*i

!-----------------------------------------------------------------------
!c  To decode navigation.

       do k=1,30
          do i=1,4
             j=44+8*(k-1)+i
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          rlat(k)=0.0001*in4
          do i=1,4
             j=44+8*(k-1)+i+4
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          rlon(k)=0.0001*in4
       enddo
       do k=1,30
          do i=1,4
             j=284+8*(k-1)+i
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          ht(k)=in4
          do i=1,4
             j=284+8*(k-1)+i+4
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          isty(k)=-9
          if( in4.eq.2 )isty(k)=0
          if( in4.eq.0 )isty(k)=1
          if( in4.eq.1 )isty(k)=2
          if( (isty(k).ne.0).and.(isty(k).ne.1) )then
              isty(k)=1
              if( ht(k).gt.0 )isty(k)=0
          endif
       enddo
       do k=1,30
          do i=1,4
             j=524+16*(k-1)+i
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          stzn(k)=0.01*in4
          do i=1,4
             j=524+16*(k-1)+i+4
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          staz(k)=0.01*in4
          do i=1,4
             j=524+16*(k-1)+i+8
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          sozn(k)=0.01*in4
          do i=1,4
             j=524+16*(k-1)+i+12
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          soaz(k)=0.01*in4
       enddo
       do i=1,4
          j=1004+i
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       hmsl=100.*in4

!-----------------------------------------------------------------------
!c  To decode earth observations.

       do k=1,30
          ij=1016+80*(k-1)
          ifqc(k)=0
          do kk=1,15
             do i=1,4
                j=ij+4*(kk-1)+i
                by4(i)=scan(j)
             enddo
             call by2int4(machine,by4,in4)
             brtmp(kk,k)=0.01*in4
             if( in4.lt.0 )ifqc(k)=1
          enddo
       enddo
       if( isqc.eq.0 )then
           ij=0
           do k=1,30
              ij=ij+ifqc(k)
           enddo
           if( ij.eq.30 )isqc=1
       endif
       
       return
       end

       subroutine dc_hirs(nbt,scan,machine,xtrack,iscan,isqc &
                    ,year,month,day,hour,minute,second,hmsl,ifqc,isty &
                    ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp,ier)
!***********************************************************************
! Subroutine/Function : dc_hirs
!
! Usage
!    call dc_hirs(nbt,scan,machine,xtrack,iscan,isqc
!   1       ,year,month,day,hour,minute,second,hmsl,ifqc,isty
!   2       ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp,ier)
!
! Description      : To decode scan line data of HIRS-1d data.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     nbt        integer     total number of input byte. (=15872)
!    I     scan(nbt)  byte        scan line data array.
!    I     machine    integer     machine index,
!                                 =0, for HP, SGI machine.
!                                 .not. 0, for DEC, PC machine.
!    I     xtrack     integer     number of FOV. (=56)
!    O     iscan      integer     scan line number.
!    O     isqc       integer     scan line quality flags.
!                                 =0, good
!                                 =1, bad
!    O     year       integer     scan line year.
!    O     month      integer     scan line month.
!    O     day        integer     scan line day.
!    O     hour       integer     scan line hour.
!    O     minute     integer     scan line minute.
!    O     second     real        scan line second.
!    O     hmsl       real        satellite altitude above refernce ellipsoid. (meters)
!    O   ifqc(xtrack) int array   FOV quality flag.
!                                 =0, good
!                                 =1, bad
!    O   isty(xtrack) int array   surface type of FOV.
!                                 =0, land
!                                 =1, sea
!                                 =2, mixed
!    O   rlat(xtrack) real array  latitude in degrees of FOV.
!    O   rlon(xtrack) real array  longitude in degrees of FOV.
!    O     ht(xtrack) real array  surface height of FOV in metres.
!    O   stzn(xtrack) real array  satellite zenith angle in degrees of FOV.
!    O   sozn(xtrack) real array  solar zenith angle in degrees of FOV.
!    O   staz(xtrack) real array  satellite azimuth angle in degrees of FOV.
!    O   soaz(xtrack) real array  solar azimuth angle in degrees of FOV.
!    O   brtmp(20,56) real array  scene brightness temperature in K for ch.1 - ch.19
!                                 scene radiance in Wm-2sr-1(cm-1)-1 for ch.20
!    O     ier        integer     error message.
!                                 =0, success.
!                                 =1, failure.
!
! Modules Called : by2int4,  jumonday1
! 
! Note: the line scan data structure
!
!       scan line information:  4*12 bytes
!       navigation:             4*451 bytes
!       earth observations:     4*2240 bytes
!       AVHRR:                  4*728 bytes
!       pre-processing output:  4*448 bytes
!       spare:                  4*89 bytes
!
!  data    field  
!  type    description
!  ....    .................................................
!                    SCAN LINE INFORMATION (4*12 bytes)
!  I  4    scan line number
!  I  4    scan line year
!  I  4    scan line day of year
!  I  4    scan line UTC time of day in milliseconds
!  I  4    quality indicator bit field -- in all of the follwong,
!          if the bit is on (i.e., if it is set to 1) then the statement is true.
!          Otherwise it is false.
!          bit 31:  do not use scan for product generation
!          bit 30:  time sequence error detected with this scan (see below)
!          bit 29:  data gap precedes this scan
!          bit 28:  no calibration (see below)
!          bit 27:  no earth location (see below)
!          bit 26:  first good time following a clock update
!          bit 25:  instrument status changed with this scan
!          bit 24-0: spare <zero fill>
!  I  4    scan line quality flags -- if bit is on (=1) then true
!          Time Problem Code: (All bits off implies the scan time is as expected.)
!          bit 31-24: spare <zero fill>
!          bit 23:  time field is bad but can probaly be inferred from the previous
!                   good time.
!          bit 22:  time field is bad and cannot be inferred from the previous good time.
!          bit 21:  this record starts a sequence that is inconsistent with previous times
!                   (i.e., there is a time discontinuity).  This may or may not be
!                   associated with a spacecraft clock update.
!          bit 20:  start of a sequence that apparently repeates scan times that have
!                   been prevously accepted.
!          bit 19-16: spare <zero fill>
!          Calibration Problem Code: (Note these bits compliment the channel indicators;
!            all bits set to 0 indicates normal calibration.)
!          bit 15:  Scan line was not calibrated because of bad time.
!          bit 14:  Scan line was calibrated using fewer than the preferred number of scan
!                   lines because of proximity to start or end of data set or to a data gap.
!          bit 13:  Scan line was not calibrated because of bad or insufficient PRT data.
!          bit 12:  Scan line was calibrated but with marginal PRT data.
!          bit 11:  Some uncalibrated channels on this scan. (See channel indicators.)
!          bit 10:  Uncalibrated due to instrument mode.
!          bit 09 and 08:  spare <zero fill>
!          Earth Location Problem Code: (all bits set to 0 implies the earth location
!            was normal)
!          bit 07:  Not earth located because of bad time.
!          bit 06:  Earth location questionable because of questionable time code.
!                   (See time problem flags above.)
!          bit 05:  Earth location questionable -- only marginal agreement with
!                   reasonableness check.
!          bit 04:  Earth location questionable -- fails reasonableness check.
!          bit 03-0:  spare <zero fill>
!  I  4    HIRS instrument baseplate temperature (K*100)
!  I  4    AMSU-A1 RF shelf instrument temperature (K*100)
!  I  4    AMSU-A2 RF shelf instrument temperature (K*100)
!  I  4    AMSU-B mixer chan 18-20 instrument temperature (K*100)
!  I  4*2  spare
!  ....    .................................................
!                     NAVIGATION (4*451 bytes)
!  I  4    10000*(latitude in degrees of position 1)
!  I  4    10000*(longitude in degrees of position 1)
!  I  4    10000*(latitude in degrees of position 2)
!  I  4    10000*(longitude in degrees of position 2)
!          ....
!  I  4    10000*(latitude in degrees of position 56)
!  I  4    10000*(longitude in degrees of position 56)
!  I  4    surface height of position 1 in metres
!  I  4    surface type of position 1 (0=sea, 1=mixed, 2=land)
!  I  4    surface height of position 2 in metres
!  I  4    surface type of position 2 (0=sea, 1=mixed, 2=land)
!          ....
!  I  4    surface height of position 56 in metres
!  I  4    surface type of position 56 (0=sea, 1=mixed, 2=land)
!  I  4    100*(local zenith angle in degrees of position 1)
!  I  4    100*(local azimuth angle in degrees of position 1)
!  I  4    100*(solar zenith angle in degrees of position 1)
!  I  4    100*(solar azimuth angle in degrees of position 1)
!  I  4    100*(local zenith angle in degrees of position 2)
!  I  4    100*(local azimuth angle in degrees of position 2)
!  I  4    100*(solar zenith angle in degrees of position 2)
!  I  4    100*(solar azimuth angle in degrees of position 2)
!          ....
!  I  4    100*(local zenith angle in degrees of position 56)
!  I  4    100*(local azimuth angle in degrees of position 56)
!  I  4    100*(solar zenith angle in degrees of position 56)
!  I  4    100*(solar azimuth angle in degrees of position 56)
!  I  4    10*(satellite altitude above refernce ellipsoid. km)
!  I  4*2  spare
!  ....    .................................................
!                    EARTH OBSERVATIONS (4*2240 bytes)
!  I  4    100* scene brightness temp. in K. FOV 1, HIRS ch.1
!          (missing data indicator is -999999)
!  I  4    100* scene brightness temp. in K. FOV 1, HIRS ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 1, HIRS ch.19
!  I  4    1000* scene radiance in W*m**2sr**-1*cm**-1 for FOV1, ch.20
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-A ch.1
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-A ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-A ch.15
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-B ch.16
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-B ch.17
!          ....
!  I  4    100* scene brightness temp. in K. FOV 1, AMSU-B ch.20
!  I  4    100* scene brightness temp. in K. FOV 2, HIRS ch.1
!  I  4    100* scene brightness temp. in K. FOV 2, HIRS ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 2, HIRS ch.19
!  I  4    1000* scene radiance in W*m**2sr**-1*cm**-1 for FOV2, ch.20
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-A ch.1
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-A ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-A ch.15
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-B ch.16
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-B ch.17
!          ....
!  I  4    100* scene brightness temp. in K. FOV 2, AMSU-B ch.20
!          ........
!  I  4    100* scene brightness temp. in K. FOV 56, HIRS ch.1
!  I  4    100* scene brightness temp. in K. FOV 56, HIRS ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 56, HIRS ch.19
!  I  4    1000* scene radiance in W*m**2sr**-1*cm**-1 for FOV56, ch.20
!  I  4    100* scene brightness temp. in K. FOV 56, AMSU-A ch.1
!  I  4    100* scene brightness temp. in K. FOV 56, AMSU-A ch.2
!          ....
!  I  4    100* scene brightness temp. in K. FOV 56, AMSU-A ch.15
!  I  4    100* scene brightness temp. in K. FOV 56, AMSU-B ch.16
!  I  4    100* scene brightness temp. in K. FOV 56, AMSU-B ch.17
!          ....
!  I  4    100* scene brightness temp. in K. FOV 56, AMSU-B ch.20
!  ....    .................................................
!                      AVHRR (4*728 bytes)
!  I  4*12*56  13 words for each HIRS FOV, 56 FOVS 
!  ....    .................................................
!                  PRE-PROCESSING OUTPUT (4*448 bytes)
!  I  4    FOV quality flags (HIRS) for FOV 1: (all bits off implies acceptable data)
!          bit 31:  spare <zero fill>
!          bit 30:  set if secondary calibration used
!          bit 29-22:  spare <zero fill>
!          bit 21:  HIRS cloud test (TBD)
!          bit 20-1:  bit n set to 1 if brightness temperature in HIRS channel n is
!                     missing or unreasonable
!          bit 0:  bad or missing data (in any or all channels)
!  I  4    FOV quality flags (HIRS) for FOV 2
!          ....
!  I  4    FOV quality flags (HIRS) for FOV 56
!  I  4    Pre-processing quality flags for FOV 1: (all bits off implies acceptable data)
!          bit 31:  set if AMSU-A surface types not all the same
!          bit 30:  set if AMSU-A used secondary calibration
!          bit 29:  set if AMSU-B used secondary calibration
!          bit 28:  set if AMSU-B missing
!          bit 27:  flag for cloud cost set for any AMSU-A
!          bit 26:  scattering tflag set for any AMSU-A (only set over the set)
!          bit 25:  logistic precipitation probability test
!                   calculated from AMSU-A data mapped to HIRS grid
!          bit 24:  Grody light rainfall test calculated on HIRS grid
!          bit 23:  mismatch between AMSU-A/B 89GHz values for any AMSU-A
!          bit 22:  mismatch between surface type from topography dataset and from
!                   pre-processing (any AMSU-A)
!          bit 21-4:  spare <zero fill>
!          bit 03:  set when AVHRR chan 3 is albedo, not bright. temp.
!          bit 02:  cloud cost flag (recalculated on HIRS grid)
!          bit 01:  scattering flag (recalculated on HIRS grid)
!          bit 00:  set if AMSU-A and AMSU-B data missing
!  I  4    estimated nearest AMSU surface type for FOV 1:
!          1 = Bare young ice (i.e. new ice, no snow)
!          2 = Dry land (i.e. dry, with or without vegetation)
!          3 = Dry snow (i.e. snow with water less than 2%, over land)
!          4 = Multi-year ice (i.e. old ice with dry snow cover)
!          5 = Sea (i.e. open water, no islands, ice-free, WS < 14 m/s)
!          6 = Wet forest (i.e. established forest with wt canopy)
!          7 = Wet land (i.e. non-forested land with a wet surface)
!          8 = Wet snow (i.e. water content > 2%, over land or ice)
!          9 = Desert
!  I  4    cost function from PPASURF surface identification for FOV 1 
!  I  4    scattering index for FOV 1
!  I  4    logistic precipitation probability for FOV 1
!  I  4*2  spare 
!  I  4    Pre-processing quality flags for FOV 2
!  I  4    estimated nearest AMSU surface type for FOV 2
!  I  4    cost function from PPASURF surface identification for FOV 2
!  I  4    scattering index for FOV 2
!  I  4    logistic precipitation probability for FOV 2
!  I  4*2  spare 
!          ....
!  I  4    Pre-processing quality flags for FOV 56
!  I  4    estimated nearest AMSU surface type for FOV 56
!  I  4    cost function from PPASURF surface identification for FOV 56
!  I  4    scattering index for FOV 56
!  I  4    logistic precipitation probability for FOV 56
!  I  4*2  spare 
!  ....    .................................................
!                      SPARE (4*89 bytes)
!
!***********************************************************************

       implicit none
       integer nbt,xtrack
       byte scan(nbt),by4(4)
       integer machine
       integer iscan,isqc,year,month,day,hour,minute,ier
       real second,hmsl
       integer ifqc(xtrack),isty(xtrack)
       real rlat(xtrack),rlon(xtrack),ht(xtrack)
       real stzn(xtrack),sozn(xtrack),staz(xtrack),soaz(xtrack)
       real brtmp(20,56),radiate
       integer i,j,k,kk,ij,jday,mn,ihour
       integer in4

       ier=0
       if( nbt.ne.15872 )then
           print*,'Error in subroutine: dc_hirs.'
           print*,'input total byte number is error.'
           print*,'nbt(=15872): ',nbt
           ier=1
           return
       endif
       if( xtrack.ne.56 )then
           print*,'Error in subroutine: dc_hirs.'
           print*,'input total xtrack number is error.'
           print*,'xtrack(=56): ',xtrack
           ier=1
           return
       endif

!-----------------------------------------------------------------------
!c  To decode scan line information.

       do i=1,4
          by4(i)=scan(i)
       enddo
       call by2int4(machine,by4,in4)
       iscan=in4

       do i=1,4
          j=i+4
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       year=in4

       do i=1,4
          j=i+8
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       jday=in4
       call jumonday1(year,jday,month,day)

       do i=1,4
          j=i+12
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       mn=in4/60000
       second=0.001*(in4-60000*mn)
       hour=mn/60
       minute=mn-60*hour

       do i=1,4
          j=i+16
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       i=in4/2
       isqc=in4-2*i

!-----------------------------------------------------------------------
!c  To decode navigation.

       do k=1,56
          do i=1,4
             j=48+8*(k-1)+i
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          rlat(k)=0.0001*in4
          do i=1,4
             j=48+8*(k-1)+i+4
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          rlon(k)=0.0001*in4
       enddo
       do k=1,56
          do i=1,4
             j=496+8*(k-1)+i
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          ht(k)=in4
          do i=1,4
             j=496+8*(k-1)+i+4
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          isty(k)=-9
          if( in4.eq.2 )isty(k)=0
          if( in4.eq.0 )isty(k)=1
          if( in4.eq.1 )isty(k)=2
          if( (isty(k).ne.0).and.(isty(k).ne.1) )then
              isty(k)=1
              if( ht(k).gt.0 )isty(k)=0
          endif
       enddo
       do k=1,56
          do i=1,4
             j=944+16*(k-1)+i
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          stzn(k)=0.01*in4
          do i=1,4
             j=944+16*(k-1)+i+4
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          staz(k)=0.01*in4
          do i=1,4
             j=944+16*(k-1)+i+8
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          sozn(k)=0.01*in4
          do i=1,4
             j=944+16*(k-1)+i+12
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          soaz(k)=0.01*in4
       enddo

       do i=1,4
          j=1840+i
          by4(i)=scan(j)
       enddo
       call by2int4(machine,by4,in4)
       hmsl=100.*in4

!-----------------------------------------------------------------------
!c  To decode earth observation.

       do k=1,56
          ij=1852+160*(k-1)
          ifqc(k)=0
          do kk=1,19
             do i=1,4
                j=ij+4*(kk-1)+i
                by4(i)=scan(j)
             enddo
             call by2int4(machine,by4,in4)
             brtmp(kk,k)=0.01*in4
             if( in4.lt.0 )ifqc(k)=1
          enddo
          do i=1,4
             j=ij+76+i
             by4(i)=scan(j)
          enddo
          call by2int4(machine,by4,in4)
          brtmp(20,k)=10.e10
!         brtmp(20,k)=0.0001*in4
!         radiate=0.0001*in4
!         brtmp(20,k)=20862.6/log(1.+36310465./radiate)
       enddo
       if( isqc.eq.0 )then
           ij=0
           do k=1,56
              ij=ij+ifqc(k)
           enddo
           if( ij.eq.56 )isqc=1
       endif

       return
       end

       subroutine dc_gen_inf(genbyte,machine,sateid,instrument,nline &
                            ,year,month,day,hour,minute,ier)
!***********************************************************************
! Subroutine/Function : dc_gen_inf
!
! Usage
!    call dc_gen_inf(genbyte,machine,sateid,instrument,nline
!   1               ,year,month,day,hour,minute,ier)
!
! Description      : To decode general information  of noaa-NN 1d data.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I    genbyte(96) byte        general information byte data array.
!    I     machine    integer     machine index,
!                                 =0, for HP, SGI machine.
!                                 .not. 0, for DEC, PC machine.
!    O     sateid     integer     satellite id.
!    O     instrument integer     instrument code.
!                                 =5, HIRS;  =6, MSU;  =10, AMSU-A;
!                                 =11, AMSU-B;  =12, MHS.
!    O     nline      integer     count of scan lines in this data set.
!    O     year       integer     start of data set year.
!    O     month      integer     start of data set month.
!    O     day        integer     start of data set day.
!    O     hour       integer     start of data set hour.
!    O     minute     integer     start of data set minute.
!    O     ier        integer     error message.
!                                 =0, success.
!                                 =1, failure.
!
! Modules Called : by2int4,  jumonday1
! 
! General Information Structure: 4*24 bytes
!  data    field  
!  type    description
!  ....    .................................................
!  C  3    1D dataset creation ID (NSS for NESDIS, etc)
!  C  1    filler
!  C  3    site of originating centre for 1B data
!          (MSC for CWB satellite center)
!  C  1    filler
!  I  4    level 1D format version number
!  I  4    level 1D format version year
!  I  4    level 1D format version day of year
!  I  4    count of header Records in this data set (=1)
!  I  4    satellite id (WMO code)
!  I  4    instrument code (5=HIRS; 6=MSU; 10=AMSU-A; 11=AMSU-B; 12=MHS)
!  I  4    10 * (nominal satellite altitude, km) (=8790)
!  I  4    nominal orbit period (seconds)
!  I  4    orbit number (at start of dataset)
!  I  4    start of data set year
!  I  4    start of data set day of year
!  I  4    start of data set UTC time of day in milliseconds
!  I  4    orbit number (at end of dataset)
!  I  4    end of data set year
!  I  4    end of data set day of year
!  I  4    end of data set UTC time of day in milliseconds
!  I  4    count of scan lines in this data set
!  I  4    count of missing scan lines
!  I  4    ATOVPP version number (values above 9000 indicate test vm)
!  I  4    instruments present (bit0=HIRS, bit1=MSU, bit2=AMUSU-A,
!          bit3=AMUSU-B, bit4=AVHRR)
!  I  4    version number of data set for antenna corrections
!          (=0 if data not corrected)
!  I  4    spare
!***********************************************************************
       implicit none
       byte genbyte(96),by4(4)
       integer machine,sateid,instrument,nline,ier
       integer year,month,day,hour,minute
       integer in4
       integer i,j,jday,mn

       ier=0

!-----------------------------------------------------------------------
!c  To decode satellite id

       do i=1,4
          j=i+24
          by4(i)=genbyte(j)
       enddo
       call by2int4(machine,by4,in4)
       sateid=in4

!-----------------------------------------------------------------------
!c  To decode general information.

       do i=1,4
          j=i+28
          by4(i)=genbyte(j)
       enddo
       call by2int4(machine,by4,in4)
       if( (in4.ne.5).and.(in4.ne.10).and.(in4.ne.11).and. &
           (in4.ne.12) )then
           ier=1
           print*,'Error in subroutine: dc_amsua.'
           print*,'This data header is not AMSU-A data header.'
           print*,'instrument code 5=HIRS, 10=AMSU-A, 11=AMSU-B, 12=MHS'
           print*,'instrument code: ',in4
           return
       endif
       instrument=in4

       do i=1,4
          j=i+44
          by4(i)=genbyte(j)
       enddo
       call by2int4(machine,by4,in4)
       year=in4     
       do i=1,4
          j=i+48
          by4(i)=genbyte(j)
       enddo
       call by2int4(machine,by4,in4)
       jday=in4     
       call jumonday1(year,jday,month,day) 
       do i=1,4
          j=i+52
          by4(i)=genbyte(j)
       enddo
       call by2int4(machine,by4,in4)
       mn=in4/60000
       hour=mn/60
       minute=mn-60*hour

       do i=1,4
          j=i+72
          by4(i)=genbyte(j)
       enddo
       call by2int4(machine,by4,in4)
       nline=in4     

       return
       end

       subroutine jumonday1(iyear,jday,month,iday)
!***********************************************************************
! Subroutine/Function : jumonday1
!
! Usage :
!    call jumonday1(iyear,jday,month,iday)
!
! Description      : To get month and day giving the day of year.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     iyear      integer     year.
!    I     jday       integer     the day of year.
!    O     month      integer     month.
!    O     iday       integer     day.
!
! Modules Called : none
!***********************************************************************

       implicit none
       integer iyear,jday,month,iday
       integer mon(12),mon1(12)
       data mon/31,59,90,120,151,181,212,243,273,304,334,365/
       data mon1/31,60,91,121,152,182,213,244,274,305,335,366/
       integer id,i1,i2,i3,i4,i

       id=0
       i1=iyear/4
       i4=iyear-4*i1
       if( i4.eq.0 )id=1
       i1=iyear/400
       i2=iyear-400*i1
       if( i2.eq.0 )id=1
       i1=iyear/100
       i3=iyear-100*i1
       if( (i3.eq.0).and.(i2.ne.0) )id=0
       if( id.eq.0 )then
           do i=1,12
              if( jday.le.mon(i) )then
                  month=i
                  if( month.eq.1 )then
                      iday=jday
                  else
                      iday=jday-mon(i-1)
                  endif
                  return
              endif
           enddo
       else
           do i=1,12
              if( jday.le.mon1(i) )then
                  month=i
                  if( month.eq.1 )then
                      iday=jday
                  else
                      iday=jday-mon1(i-1)
                  endif
                  return
              endif
           enddo
       endif
       return
       end

      subroutine by2int4(machine,by4,in4)
!***********************************************************************
! Subroutine/Function : by2int4
!
! Usage :
!    call by2int4(by4,in4)
!
! Description      : To transform 4-bytes to integer for HP, SGI machine.
!
! Arguments :
!  I/O/W   name,      type,       description
!    I     machine    integer     machine index,
!                                 =0, for HP, SGI machine.
!                                 .not. 0, for DEC, PC machine.
!    I     by4(4)     byte array  input 4 bytes.
!    O     in4        integer   output integer.
!
! Modules Called : none
!***********************************************************************

      integer in4,ib
      byte by4(4),by_tem(4)
      equivalence (ib,by_tem(1))
      if( machine.eq.0 )then
          by_tem(1)=by4(4)
          by_tem(2)=by4(3)
          by_tem(3)=by4(2)
          by_tem(4)=by4(1)
      else
          by_tem(1)=by4(1)
          by_tem(2)=by4(2)
          by_tem(3)=by4(3)
          by_tem(4)=by4(4)
      endif
      in4=ib
      return
      end
