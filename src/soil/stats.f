
      subroutine moment(data,n,
     &                  ave,adev,sdev,var,skew,curt,
     &                  istatus)
c
c
      dimension data(9)
      parameter (r_missing_data = +1e37)
      istatus=0
      if(n.le.1)then
            print *, 'n must be at least 2'
            istatus=1
            return
      endif
      s=0.
      do 11 j=1,n
        s=s+data(j)
11    continue
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      do 12 j=1,n
        s=data(j)-ave
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    continue
      adev=adev/n
      var=var/(n-1)
      if(var.ne.0.)then
        sdev=sqrt(var)
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.  
      else
        sdev=r_missing_data
        skew=r_missing_data
        curt=r_missing_data
c       print*,  'no skew or kurtosis when zero variance'
      endif
      return
      end

