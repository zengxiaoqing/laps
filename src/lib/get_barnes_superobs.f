
        subroutine get_barnes_superobs(ni,nj,nk,r_missing_data           ! I
     1                        ,n_var,i_var,max_obs,obs_barnes            ! I
     1                        ,i_ratio,i4time_sys                        ! I
     1                        ,superobs_barnes,ncnt_super)               ! O

        include 'barnesob.inc'
        type (barnesob) :: obs_barnes(max_obs)
        type (barnesob) :: superobs_barnes(max_obs)

        integer nsup(ni,nj,nk)
        real sumi(ni,nj,nk)
        real sumj(ni,nj,nk)
        real sumo(ni,nj,nk,n_var)
        real sumw(ni,nj,nk)
        real sumv(ni,nj,nk)
        real xnew(n_var)

	real residualu
	real sumsq_obs

        ncnt_super = 0
        residualu = 0.
        sumsq_obs = 0.
   
        nsup = 0
        sumi = 0.
        sumj = 0.
        sumo = 0.
        sumw = 0.
        sumv = 0.

        write(6,*)' get_barnes_superobs...'

        if(nk .gt. 1)then
            nprint = 20
        else
            nprint = 100
        endif

        r_ratio = i_ratio

        sumw1 = 0.

        do n = 1,max_obs    

            ri = obs_barnes(n)%i
            rj = obs_barnes(n)%j
            rk = obs_barnes(n)%k

!           i4time = obs_barnes(n)%i4time

            i_cell = min(nint((ri-1.)/r_ratio) * i_ratio + 1, ni)
            j_cell = min(nint((rj-1.)/r_ratio) * i_ratio + 1, nj)
            k_cell = rk

            nsup(i_cell,j_cell,k_cell) = nsup(i_cell,j_cell,k_cell) + 1 
            sumi(i_cell,j_cell,k_cell) = sumi(i_cell,j_cell,k_cell) + ri
            sumj(i_cell,j_cell,k_cell) = sumj(i_cell,j_cell,k_cell) + rj
            sumo(i_cell,j_cell,k_cell,:) = sumo(i_cell,j_cell,k_cell,:) 
     1                                 + obs_barnes(n)%value(:)

!           Note that time weight is being factored into the superob
            call get_time_wt(i4time_sys,obs_barnes(n)%i4time
     1                                 ,time_wt,istatus)

            sumw(i_cell,j_cell,k_cell) = sumw(i_cell,j_cell,k_cell) 
     1                                 + obs_barnes(n)%weight * time_wt
            sumv(i_cell,j_cell,k_cell) = sumv(i_cell,j_cell,k_cell) 
     1                                 + obs_barnes(n)%vert_rad_rat

            sumw1 = sumw1 + obs_barnes(n)%weight * time_wt

        enddo ! n

        sumw2 = 0.

        do k = 1,nk
        do j = 1,nj
        do i = 1,ni
            if(nsup(i,j,k) .gt. 0)then
                inew    = nint(sumi(i,j,k)   / float(nsup(i,j,k)))
                jnew    = nint(sumj(i,j,k)   / float(nsup(i,j,k)))
                xnew(:) =      sumo(i,j,k,:) / float(nsup(i,j,k))  
                vnew    =      sumv(i,j,k)   / float(nsup(i,j,k))  
                ncnt_super = ncnt_super + 1
                do l = 1,n_var
                    superobs_barnes(ncnt_super)%value(l) = xnew(l)
                enddo ! l
                superobs_barnes(ncnt_super)%i = inew
                superobs_barnes(ncnt_super)%j = jnew
                superobs_barnes(ncnt_super)%k = k
                superobs_barnes(ncnt_super)%weight = sumw(i,j,k)
                superobs_barnes(ncnt_super)%i4time = i4time_sys 
                superobs_barnes(ncnt_super)%vert_rad_rat = vnew 

                sumw2 = sumw2 + superobs_barnes(ncnt_super)%weight
            endif
        enddo ! i
        enddo ! j
        enddo ! k

        write(6,*)' i_ratio,max_obs,ncnt_super = ',
     1              i_ratio,max_obs,ncnt_super

        write(6,*)' sumw1,sumw2 = ',sumw1,sumw2

        return

        end

