
      subroutine skill_scores(contable,lun_out           ! I
     1                       ,frac_coverage              ! O (obs OR fcst)
     1                       ,frac_obs                   ! O 
     1                       ,frac_fcst                  ! O 
     1                       ,bias,ets)                  ! O

      integer,parameter :: k12 = selected_int_kind(12)

!     First index is observed, second index is forecast
!     0 is Yes, 1 is No
      integer (kind=k12) :: contable(0:1,0:1)

      integer (kind=k12) :: hits,misses,false_alarms,correct_negatives
     1                     ,total 
      
      hits              = contable(0,0)
      misses            = contable(0,1)
      false_alarms      = contable(1,0)
      correct_negatives = contable(1,1)

      rmiss = -999.

      if(minval(contable) .lt. 0)then
          write(6,*)' ERROR: skill_scores contable has negative values'
          frac_coverage = rmiss
          frac_obs = rmiss
          frac_fcst = rmiss
          bias = rmiss
          ets = rmiss
          return
      endif

      total = hits + misses + false_alarms + correct_negatives

      if(total .gt. 0)then
          frac_negatives = float(correct_negatives) / float(total)
          frac_coverage = 1.0 - frac_negatives

          frac_obs  = float(hits + misses)       / float(total)
          frac_fcst = float(hits + false_alarms) / float(total)

          accuracy = float(hits + correct_negatives) / float(total)
      else
          frac_negatives = rmiss
          frac_coverage = rmiss
          frac_obs = rmiss
          frac_fcst = rmiss
          accuracy = rmiss
          if(lun_out .gt. 0)then
              write(6,*)' WARNING: no data points in skill_scores'
          endif
      endif

      if(hits + misses .gt. 0)then
          bias = float(hits + false_alarms) / float(hits + misses)
          pod = float(hits) / float(hits + misses)
      else
          bias = rmiss
          pod = rmiss
      endif

      if(hits + false_alarms .gt. 0)then
          far = float(false_alarms) / 
     1          float(hits + false_alarms)      
      else
          far = rmiss
      endif

      if(total .gt. 0)then
          hits_random = float((hits + misses) * (hits + false_alarms)) 
     1            / float(total)

          denom = float(hits + misses + false_alarms) - hits_random
          if(denom .gt. 0.)then
              ets = (float(hits) - hits_random) / denom
          else
              ets = rmiss
          endif
      else
          hits_random = rmiss
          ets = rmiss
      endif

      if(lun_out .gt. 0)then
          write(lun_out,*)' Hits = ',hits
          write(lun_out,*)' Misses = ',misses
          write(lun_out,*)' False Alarms = ',false_alarms
          write(lun_out,*)' POD = ',pod
          write(lun_out,*)' FAR = ',far
          write(lun_out,*)' Correct Negatives = ',correct_negatives
          write(lun_out,*)' Hits Random = ',hits_random
          write(lun_out,*)' Total = ',total

          write(lun_out,*)' Accuracy = ',accuracy
          write(lun_out,*)' Bias = ',bias
          write(lun_out,*)' ETS = ',ets
      endif

      return
      end
