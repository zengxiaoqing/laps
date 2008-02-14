c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine covmat_humid(i4time) 

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c c History:
c Date         Name          Action
c --------     -----------   -------------------------------------------
c 03/02/2007   Ok-Yeon Kim   Created. 
c ----------------------------------------------------------------------


      implicit none

      include           'netcdf.inc'

cc-- variable for initialization --cc
   
      integer           max_files
      integer           nvars
      integer           nrec
      integer           ndims
      integer           nrecord
      integer           namelen
      integer           nav
      integer           ix,iy,iz
      integer           iblev
      real              dx,dy

      parameter         (max_files=5000)
      parameter         (nvars=50)         ! number of variables
      parameter         (nrec=1)           ! change this to generalize
      parameter         (ndims=7)          ! number of dimensions
      parameter         (nrecord=1)
      parameter         (namelen=132)
      parameter         (nav=1)
      parameter         (ix=93)
      parameter         (iy=65)
      parameter         (iz=19)
      parameter         (iblev=5)
      parameter         (dx=81.2705)
      parameter         (dy=81.2705)

      integer         i4time
      integer           num_of_ens       ! number of time lagged ensembles
      integer           iens
      integer           ilength
      integer           jlength
      integer           dlen,flen 
      integer           itimes(max_files)
      integer           bg_files
      integer           ihour
      integer           i,j,ij,jj,n,m,k,l 
      integer           ibkgd
      integer           ifcst_bkgd(3000)
      integer           i4timeinit(3000)
      integer           i4time_min_diff
      integer           valid_time_1,valid_time_2
      integer           indx_best_init,indx_best_fcst

      character*4       af 
      character*9       bkgd(3500),bkgd_ens(3500)
      character*4       fcst(3000,120),fcst_ens(3000,12)
      character*9       finit1,finit2      
      character*13      ens_names(max_files)
      character*200     iname(max_files)
      character*256     path_to_ens_bkg  ! path to time lagged ensembles
      character*256     names(max_files)
      character*256     ifnames_ens(max_files)
      character*200     fname

cc-- variable for readnet_nam --cc

c      integer           nvars,nrec,ndims,namelen,nav,nrecord
c      integer           ix,iy,iz,iblev
c      real              dx,dy
      integer         rcode              ! error code
      integer         recdim             ! recode dimension
      integer         ncid

      real            gh_sfc(ix,iy,nrec)
      real            gh(ix,iy,iz,nrec)

      real            rh_2mFH(ix,iy,nrec)
      real            rh(ix,iy,iz,nrec)
      real            rh_lbis(ix,iy,iblev,nrec)

      real            t_2mFH(ix,iy,nrec)
      real            t(ix,iy,iz,nrec)
      real            t_lbis(ix,iy,iblev,nrec)

      real            uw_10mFH(ix,iy,nrec)
      real            uw(ix,iy,iz,nrec)
      real            uw_lbis(ix,iy,iblev,nrec)

      real            vw(ix,iy,iz,nrec)
      real            vw_lbis(ix,iy,iblev,nrec)
      real            vw_10mFH(ix,iy,nrec)

      real            av(ix,iy,iz,nrec)
      real            pvv(ix,iy,iz,nrec)
      real            p_sfc(ix,iy,nrec)
      real            heli(ix,iy,nrec)
      real            cape_sfc(ix,iy,nrec)
      real            cape_lbis(ix,iy,nrec)
      real            cin_sfc(ix,iy,nrec)
      real            cin(ix,iy,nrec)
      real            bli_lbis(ix,iy,nrec)
      real            pli_lbis(ix,iy,nrec)
      real            pw(ix,iy,nrec)
      real            emspMSL(ix,iy,nrec)
      real            prMSL(ix,iy,nrec)
      real            cp_sfc(ix,iy,nrec)
      real            tp_sfc(ix,iy,nrec)
      integer         isoLevel(iz)
      integer         boundryLevel(iblev)
      real*8            valtime(nrec)
      real*8            reftime(nrec)
      character*1       origin(namelen)
      character*1       model(namelen)
      character*1       grid_type(namelen,nrec)
      character*1       j_dim(namelen,nrec)
      character*1       i_dim(namelen,nrec)
      integer*2         version
      integer*2         Ni(nav)
      integer*2         Nj(nav)
      real            La1(nav)
      real            La2(nav)
      real            Lo1(nav)
      real            Lo2(nav)
      real            Di(nav)
      real            Dj(nav)
      real            IntLat1(nav)
      real            IntLat2(nav)
      real            Lon0(nav)
      integer         start(ndims)            ! hyperslab starting index
      integer         count(ndims)            ! hyperslab count from start
      integer           vdims(ndims)            ! max # of var dims
      character*1024    strbuf                  ! string buffer for var
                                                !  and attr names
cc-- variable for covariance matrix --cc

      character*200     ofname(max_files)
      character*200     ofile
      character*256     indir
      character*9       a9
      character*7       c_ext
      character*6       a_ext
      character*5       b_ext
      character*4       q_ext

      integer           inum
      integer           len_dir,inlen
      integer         istatus
      real            rh_con(ix,iy,iz,20)
      real            t_con(ix,iy,iz,20)
      real            q_con(ix,iy,iz,20)
      real            q_inno(ix,iy,iz,20)
      real            q_mean(ix,iy,iz)
      real            q_mat(ix*iy*iz,20)
      real            q_cov(ix*iy*iz,ix*iy*iz) ! full covariance matrix
      real            qcovhor(9,9,ix,iy,iz)    ! horizontal covariance
      real            qcovert(iz,iz,ix,iy)     ! vertical covariance 
      real            aa,bb
 
      real              es(ix,iy,iz,20)
      real              e(ix,iy,iz,20)
      real              rho(ix,iy,iz,20)
      real              es0,lv,rv,t0,ratio

cc-- namelist --cc

      namelist /time_lagged_ens_nl/ path_to_ens_bkg,
     +                              num_of_ens,
     +                              ilength,
     +                              jlength


c--------------- End of Diclaration -----------------------------------


      write (6,*) '                                    '
      write (6,*) '++++++++++++++++++++++++++++++++++++'
      write (6,*) 'Start of Background Error Covariance'
      write (6,*) '++++++++++++++++++++++++++++++++++++'
      write (6,*) '                                    '

c
c Read analysis time
c

      write (6,*) 'Read analysis time:  ',i4time


c ----------------------------------------------------------------------
c Read time lagged ensemble membesr (background)
c ----------------------------------------------------------------------


      call get_directory('static',fname,dlen)
        open(22,file=fname(1:dlen)//'time_lagged_ens.nl'
     +         ,status='old')
        read(22,time_lagged_ens_nl)  
        close(22)

      call s_len(path_to_ens_bkg,dlen) 
        fname=path_to_ens_bkg(1:dlen)

      call get_file_times(fname,max_files,names,itimes,
     +                    bg_files,istatus)

c
c 
c

      do i=1,bg_files
         call s_len(names(i),j) 

        j=j-13
        if (j .ge. 0) then
           ens_names(i)=names(i)(j+1:j+13) 
        endif
      enddo

      
      ij=0
      ibkgd=1
      do n=1,bg_files
         if(n .eq. 1) then
            bkgd(ibkgd)=ens_names(n)(1:9)
            ifcst_bkgd(ibkgd)=1
            call i4time_fname_lp(bkgd(ibkgd),i4timeinit(ibkgd),
     +                           istatus)
         else
            finit1=ens_names(n-1)(1:9)
            finit2=ens_names(n)(1:9)
            if(finit1 .eq. finit2)then
                ifcst_bkgd(ibkgd)=ifcst_bkgd(ibkgd)+1
            else
                ibkgd=ibkgd+1
                bkgd(ibkgd)=ens_names(n)(1:9)
                ifcst_bkgd(ibkgd)=1
                call i4time_fname_lp(bkgd(ibkgd),i4timeinit(ibkgd),
     +                               istatus)
            endif
         endif       
            fcst(ibkgd,ifcst_bkgd(ibkgd))=ens_names(n)(10:13)
      enddo       

      do n=1,ibkgd ; do m=1,ifcst_bkgd(n)
         bkgd_ens(n)=bkgd(n) 
         fcst_ens(n,m)=fcst(n,m) 
      enddo ; enddo 

c
c
c

      n=1
      indx_best_init=0
      indx_best_fcst=0
      i4time_min_diff=100000
      do while(n.le.ibkgd)
         do jj=2,ifcst_bkgd(n)
            af=fcst(n,jj-1)
            read(af,'(i4)') ihour
            valid_time_1=i4timeinit(n)+ihour*3600
            af=fcst(n,jj)
            read(af,'(i4)') ihour
            valid_time_2=i4timeinit(n)+ihour*3600
            if(valid_time_1 .le. i4time.and.
     +         valid_time_2 .ge. i4time)then
                if(abs(i4timeinit(n)-i4time).lt.i4time_min_diff)then
                   i4time_min_diff=abs(i4timeinit(n)-i4time)
                   indx_best_init=n
                   indx_best_fcst=jj-1
                endif
            endif
          enddo
      n=n+1
      enddo

c ----------------------------------------------------------------------
c Accept time-phased ensemble members
c ----------------------------------------------------------------------
      

      do n=1,num_of_ens       ! except for analysis time for real time run

         ifnames_ens(n)=bkgd_ens(indx_best_init-n)
     +                //fcst_ens(indx_best_init,indx_best_fcst+n)
 
         write (6,*) 'Ensemble member(',n,')', TRIM(ifnames_ens(n))

      enddo 

c ----------------------------------------------------------------------
c Read NAM (NetCDF) file
c ----------------------------------------------------------------------

      do 300 n=1,num_of_ens 
         
         iname(n)=path_to_ens_bkg(1:dlen)//'/'//ifnames_ens(n)
         call s_len(iname(n),flen)

         write (6,*) 'member =',n,iname(n)(1:flen) 

         ncid=ncopn(iname(n),ncnowrit,rcode)

         call readnet_nam(ncid,ix,iy,iz,iblev,
     +                    gh_sfc,gh,rh_2mFH,rh,rh_lbis,t_2mFH,
     +                    t,t_lbis,uw_10mFH,uw,uw_lbis,vw,
     +                    vw_lbis,vw_10mFH,av,pvv,p_sfc,heli,
     +                    cape_sfc,cape_lbis,cin_sfc,cin,
     +                    bli_lbis,pli_lbis,pw,emspMSL,prMSL,
     +                    cp_sfc,tp_sfc,isoLevel,boundryLevel,
     +                    valtime,reftime,origin,model,grid_type,
     +                    j_dim,i_dim,Ni,Nj,La1,La2,Lo1,Lo2,
     +                    Di,Dj,IntLat1,IntLat2,Lon0,start,
     +                    count,vdims,strbuf,num_of_ens)

         do m=1,nrec      ! rh (relative humidity)
           do k=1,iz
             do j=1,iy
               do i=1,ix
                 rh_con(i,j,k,n)=rh(i,j,k,m)
              enddo
             enddo
           enddo
         enddo 

         do m=1,nrec      ! t (temperature)
           do k=1,iz
             do j=1,iy
               do i=1,ix
                 t_con(i,j,k,n)=t(i,j,k,m)
               enddo
             enddo
           enddo
         enddo

  300 enddo


c ----------------------------------------------------------------------
c Convert relative humidity to specific humidity
c ----------------------------------------------------------------------

      es0=6.11         ! constant for conversion from rh to sh
      lv=2.5*(10**6)
      rv=461.5
      t0=273.15
      ratio=0.62197 ! (rl/rw)

      do n=1,num_of_ens 
        do k=1,iz
          do j=1,iy
            do i=1,ix
              es(i,j,k,n)=es0*exp(lv/rv*(1/t0-1/t_con(i,j,k,n)))
              e(i,j,k,n)=(rh_con(i,j,k,n)/100)*es(i,j,k,n)
              rho(i,j,k,n)=ratio*(e(i,j,k,n)
     +                     /(isoLevel(k)+(e(i,j,k,n)*(ratio-1))))
              q_con(i,j,k,n)=rho(i,j,k,n)
            enddo
          enddo
        enddo
      enddo


c ----------------------------------------------------------------------
c Compute Ensemble Mean / Innovation
c ----------------------------------------------------------------------

      do k=1,iz ; do j=1,iy ; do i=1,ix
         q_mean(i,j,k)=0. 
      enddo ; enddo ; enddo 

      do n=1,num_of_ens
        do k=1,iz ; do j=1,iy ; do i=1,ix
           q_mean(i,j,k)=q_mean(i,j,k)+q_con(i,j,k,n)
        enddo ; enddo ; enddo
      enddo

      do k=1,iz ; do j=1,iy ; do i=1,ix
         q_mean(i,j,k)=q_mean(i,j,k)/float(num_of_ens)
      enddo ; enddo ; enddo 
         
      do n=1,num_of_ens
        do k=1,iz ; do j=1,iy ; do i=1,ix
           q_inno(i,j,k,n)=q_con(i,j,k,n)-q_mean(i,j,k)
        enddo ; enddo ; enddo
      enddo


      istatus=1
      call make_fnam_lp(i4time,a9,istatus)

      call get_directory('cov',indir,len_dir)
      call s_len(indir,inlen)

      b_ext='.q.em'
      ofname(1)=indir(1:inlen)//a9//b_ext
      call s_len(ofname(1),flen)

      write(ofile,'(a)') ofname(1)(1:flen)
      print*,ofile
      open(1,file=ofile,status='unknown',form='unformatted'
     +                 ,access='sequential')
      write(1) (((q_mean(i,j,k),i=1,ix),j=1,iy),k=1,iz)
      close(1)


      write (6,*)  '                                    '
      write (6,*)  'Compute ensemble mean and innovation'
      write (6,*)  '                                    '

      iens=num_of_ens-1

  
c ----------------------------------------------------------------------
c Full Error Covariance Matrix 
c ----------------------------------------------------------------------

      do n=1,num_of_ens ; inum=0
        do k=1,iz
          do j=1,iy
            do i=1,ix
              inum=inum+1
                q_mat(inum,n)=q_inno(i,j,k,n)
            enddo
          enddo
        enddo
      enddo
 
      do m=2,inum          ! inum = ix*iy*iz
         do n=1,m-1 ; aa = 0.
            do k=1,num_of_ens 
               aa = aa + (q_mat(m,k)*q_mat(n,k))
            enddo         ! k
               aa = aa / float(iens)
               q_cov(m,n) = aa
               q_cov(n,m) = aa
         enddo            ! n
      enddo               ! m


      do m=1,inum ; bb = 0.   ! inum = ix*iy*iz
         do k=1,num_of_ens
            bb = bb + (q_mat(m,k)*q_mat(m,k))
         enddo            ! k
            bb = bb / float(iens)
            q_cov(m,m) = bb
      enddo               ! m

      istatus=1
      call make_fnam_lp(i4time,a9,istatus)
      
      call get_directory('cov',indir,len_dir)
      call s_len(indir,inlen)

      c_ext='.q.bmat'
      ofname(1)=indir(1:inlen)//a9//c_ext
      call s_len(ofname(1),flen)

      write(ofile,'(a)') ofname(1)(1:flen)
      print*,ofile
      open(1,file=ofile,status='unknown',form='unformatted'
     +                 ,access='sequential',recl=ix*iy*iz)
      do l=1,inum
         write(1,rec=l) (q_cov(m,l),m=1,inum)
      enddo 
      close(1)


      write (6,*) '                                                 '
      write (6,*) 'Compute full error covariance matrix for humidity'
      write (6,*) '                                                 '


c ----------------------------------------------------------------------
c Localization of covariances from lagged ensemble forecasts 
c ----------------------------------------------------------------------
c
c Step 1 : Covariance (covhor) in the averaged Horizontal Length (ilength,jlength)
c

      call hor_cov(i4time,q_inno
     +                   ,num_of_ens
     +                   ,ix,iy,iz,isoLevel
     +                   ,dx,dy
     +                   ,ilength,jlength
     +                   ,qcovhor)

      istatus=1
      call make_fnam_lp(i4time,a9,istatus)

      call get_directory('cov',indir,len_dir)
      call s_len(indir,inlen)

      a_ext='.q.hor'
      ofname(1)=indir(1:inlen)//a9//a_ext
      call s_len(ofname(1),flen)

      write(ofile,'(a)') ofname(1)(1:flen)
c      write(ofile,'(a,a,a)') './',a9,'.horcov.bin'
      print*,ofile
      open(1,file=ofile,status='unknown',form='unformatted'
     +                 ,access='sequential')
      write(1) 
     +    (((((qcovhor(l,m,i,j,k),l=1,9),m=1,9),i=1,ix),j=1,iy),k=1,iz)
      close(1)


      write (6,*) '                                                '
      write (6,*) 'Covariance in the horizontal length for humidity'
      write (6,*) '                                                '


c
c Step 2 : Covariance (covert) in the vertical direction
c

      call vert_cov(i4time,q_inno
     +                    ,num_of_ens
     +                    ,ix,iy,iz,isoLevel
     +                    ,qcovert)

      istatus=1
      call make_fnam_lp(i4time,a9,istatus)

      call get_directory('cov',indir,len_dir)
      call s_len(indir,inlen)

      q_ext='.cov'
      ofname(1)=indir(1:inlen)//a9//q_ext
      call s_len(ofname(1),flen)

      write(ofile,'(a)') ofname(1)(1:flen)
c      write(ofile,'(a,a,a)') './',a9,'.vertcov.bin'
      print*,ofile
      open(1,file=ofile,status='unknown',form='unformatted'
     +                 ,access='sequential')
      write(1)
     +    ((((qcovert(k,l,i,j),k=1,iz),l=1,iz),i=1,ix),j=1,iy)
      close(1)


      write (6,*)  '                                    '
      write (6,*)  'Covariance in the vertical direction'
      write (6,*)  '                                    '


c ----------------------------------------------------------------------
c Call for double check ...
c ----------------------------------------------------------------------

 
      call recheck(i4time,qcovert,qcovhor,ix,iy,iz)    

      write (6,*)  '                                             '
      write (6,*)  'Done for covariance estimate and double check'
      write (6,*)  '                                             '
 


      return
      end

