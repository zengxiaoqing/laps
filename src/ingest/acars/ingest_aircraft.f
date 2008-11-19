
       program ingest_aircraft

!      Driver for aircraft ingest (PIN intermediate file)

!      Steve Albers      May-1999       Original Version

       write(6,*)
       write(6,*)' Call ingest_pireps'
       call ingest_pireps(istatus)

       write(6,*)
       write(6,*)' Call ingest_acars'
       call ingest_acars(istatus)

       write(6,*)
       write(6,*)' Call ingest_wisdom'
       call ingest_wisdom(istatus)

       end

