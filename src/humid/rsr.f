       subroutine rsr (kan, rad, ii,jj)
       real rad (ii,jj)
       integer kan,ii,jj
       character*2 ckan



       write(ckan,2) kan
2      format (i2.2)


       open (1, file = '/data/peaks/laps/nest7grid/lapsdat/lsr/
     1lapssndr_rad_
     1'//ckan//'.dat', form = 'unformatted')
       read (1) rad

       close (1)

       return
       end

