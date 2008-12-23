
      subroutine skill_scores(contable,lun_out)

!     First index is observed, second index is forecast
!     0 is Yes, 1 is No
      integer contable(0:1,0:1)

      integer hits,misses,false_alarms,correct_negatives,total
      
      hits              = contable(0,0)
      misses            = contable(0,1)
      false_alarms      = contable(1,0)
      correct_negatives = contable(1,1)

      total = hits + misses + false_alarms + correct_negatives

      accuracy = float(hits + correct_negatives) / float(total)
      bias = float(hits + false_alarms) / float(hits + misses)

      hits_random = ((hits + misses) * (hits + false_alarms)) / total
      ets = (hits - hits_random) / 
     1      (hits + misses + false_alarms - hits_random)

      write(lun_out,*)' Hits = ',hits
      write(lun_out,*)' Misses = ',misses
      write(lun_out,*)' False Alarms = ',false_alarms
      write(lun_out,*)' Correct Negatives = ',correct_negatives
      write(lun_out,*)' Total = ',total

      write(lun_out,*)' Accuracy = ',accuracy
      write(lun_out,*)' Bias = ',bias
      write(lun_out,*)' ETS = ',ets
   
      return
      end
