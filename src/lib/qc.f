

        function pct_rejected(n_good_qc,n_bad_qc)

        n_total = n_good_qc + n_bad_qc
 
        if(n_total .gt. 0)then
            pct_rejected = (float(n_bad_qc) / float(n_total)) * 100.       
        else
            pct_rejected = 0.
        endif

        return
        end

