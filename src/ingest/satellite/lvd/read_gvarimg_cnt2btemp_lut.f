       subroutine read_gvarimg_cnt2btemp_lut(csatid,
     &chtype,cnt2btemp,istatus)
c
       implicit none

       integer i,idummy
       integer n,nn,np
       integer icnt
       real*4  cnt2btemp(0:1023)
       real*4  rdummy
       integer istatus

       character*100 path
       character*255 filename
       character*6   csatid
       character*3   chtype

       istatus = 0

       nn=index(chtype,' ')-1
       if(nn.le.0)nn=3
       call get_directory('static',path,np)
       path=path(1:np)//'lvd/'//csatid//'lookup/'
       np=index(path,' ')-1
       filename=path(1:np)//chtype(1:nn)//'.lut'

       n=index(filename,' ')
       open(22,file=filename,status='old',form='formatted',
     &err=900)

       do i=0,1023
          read (22,*) idummy,rdummy,cnt2btemp(i)
c         read(22,*,err=901)icnt,rdummy,cnt2btemp(i),idummy
       enddo

       close(22)

       goto 1000

900    write(6,*)'Error opening file: ',filename(1:n)
       istatus = -1
       goto 1000

901    write(6,*)'Error reading file: ',filename(1:n)
       istatus = -1

1000   return
       end
