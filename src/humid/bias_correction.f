cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      function bias_correction (raw_data, satellite, sounder, 
     1     channel_index)

      implicit none


c     routine to fit GOES bias corrections into variational processing
c     note that optran is too warm.  to make data appear as optran values,
c     add the bias.
c     author birkenheuer     9/29/99

c     optran interface code
c     chief application GIMPAP/AFWA

      real bias_correction
      integer channel_index
      integer satellite         ! 8= goes 8,  10 = goes 10 etc
      integer sounder           ! 1= sounder, 0 = imager
      real bias_8 (22), bias_10(22)
      namelist /bias_coefficients_nl/ bias_8,bias_10
      save bias_8, bias_10, first_time
      real raw_data
      character*200 fname
      integer len
      integer first_time
      data first_time /0/



      if (first_time .eq. 0) then !get coef

c     read in namelist
         call get_directory('static',fname,len)
         open (23, file=fname(1:len)//'bias_coefficients.nl',
     1        status = 'old', err = 24)

         read(23,bias_coefficients_nl,end=24)
         close (23)

         first_time = 1         ! set to not first time
      endif


c     end read namelist

      bias_correction = raw_data

      if(satellite .eq. 8) then
         if(sounder .eq. 1) then ! correct raw_data to match Optran
            bias_correction = raw_data + bias_8(channel_index)
         endif
      endif

      
      if(satellite .eq.10) then
         if(sounder.eq.1) then
            bias_correction = raw_data + bias_10(channel_index)
         endif
      endif

      return

 24   write(6,*) 'WARNING no bias coefficients found'

      return

      end

