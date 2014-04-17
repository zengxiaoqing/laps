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
        subroutine climate_sm(v_lat,            !input
     1                julian_day,       !input
     1                standard_press,   !output
     1                tempertur_guess,  !ouput
     1                mixratio_guess,   !output
     1                istatus)          !ouput

c       $log: climate_sm.for,v $
c revision 1.1  1996/08/30  20:44:00  birk
c initial revision
c


c       dolan, m        20-aug-84       modified uw source code
c
c    title:             climate_sm
c
c    file:              profs4::user_devclimate.for
c
c    abstract:          this routine generates a climatological
c                       guess for satellite simulation/retrieval
c                       at 40 fixed pressure levels.  it is to
c                       be used as a guess profile into the
c                       vasar subsystem.
c
c    environment:       processor:      digital equipment corp.  vax 11/750
c                       o.s.:           vax/vms  version  3.6
c                       language:       vax-11 fortran (fortran 77 v5.0)
c
        implicit none

c       adapted to unix environment 10/14/94  (dan birkenheuer)
c


c-------------------------------------------------------------------------------
c    formal parameter declarations
c-------------------------------------------------------------------------------
       save 
       real  v_lat,
     1  standard_press(40),
     1  tempertur_guess(40),
     1  mixratio_guess(40)

        integer       julian_day,
     1          istatus


c-------------------------------------------------------------------------------
c    working variable declarations
c-------------------------------------------------------------------------------
        real  temp_tmp(40,2), !temporary temperature array
     1  temp_pres(40),  !temporary pressure array
     1  temp_mixrat(20,2),!temporary mixing ratio array
     1  temp_hgt(2),    !temporary height array
     1  fs2,            !in line function
     1  f,              !real variable for in line function definition
     1  sat_pt_sm,      !saturation function
     1  wt1,wt2,        !weights for interpolation
     1  wx,z10,ws,      !temporary mixing ratio variables
     1  alat,           !latitude
     1  relative_hum    !relative humidity
c
c.......temperatures are in degrees kelvin:
c
        integer       jan_tmp_15(40),         !15n latitude jan temp
     1  jan_tmp_30(40),         !30n latitude jan temp
     1  jan_tmp_45(40),         !45n latitude jan temp
     1  jan_tmp_60(40),         !60n latitude jan temp
     1  jan_tmp_75(40),         !75n latitude jan temp
     1  jul_tmp_15(40),         !15n latitude july temp
     1  jul_tmp_30(40),         !30n latitude july temp
     1  jul_tmp_45(40),         !45n latitude july temp
     1  jul_tmp_60(40),         !60n latitude july temp
     1  jul_tmp_75(40)         !75n latitude july temp
c
        integer
     1  jan_rel_hum_15(20),     !15n latitude jan rel hum
     1  jan_rel_hum_30(20),     !30n latitude jan rel hum
     1  jan_rel_hum_45(20),     !45n latitude jan rel hum
     1  jan_rel_hum_60(20),     !60n latitude jan rel hum
     1  jan_rel_hum_75(20),     !65n latitude jan rel hum
     1  jul_rel_hum_15(20),     !15n latitude july rel hum
     1  jul_rel_hum_30(20),     !30n latitude july rel hum
     1  jul_rel_hum_45(20),     !45n latitude july rel hum
     1  jul_rel_hum_60(20),     !60n latitude july rel hum
     1  jul_rel_hum_75(20)     !75n latitude july rel hum
c
        integer
     1          kelv_temp(40,5,2),      !3-dimen temp array in k
     1          iwvmr(20,5,2),          !3-dimen relative hum array
     1          ihite(5,2),             !height array( 5 levels )
     1          lr(15),
c     1          is2,                    !in line function not used
     1          jlat,
     1          jl,kk,i,j,k,jl1,jl2,
     1          it1,it2,iz1,iz2,iw1,iw2,kday,nmon,imon,
     1          lmon,nday,nyr,
     1          iz10


c-------------------------------------------------------------------------------
c    equivalence statements
c-------------------------------------------------------------------------------
        equivalence (jan_tmp_15(1),kelv_temp(1,1,1)),
     1      (jul_tmp_15(1),kelv_temp(1,1,2))
        equivalence (jan_tmp_30(1),kelv_temp(1,2,1)),
     1      (jul_tmp_30(1),kelv_temp(1,2,2))
        equivalence (jan_tmp_45(1),kelv_temp(1,3,1)),
     1      (jul_tmp_45(1),kelv_temp(1,3,2))
        equivalence (jan_tmp_60(1),kelv_temp(1,4,1)),
     1      (jul_tmp_60(1),kelv_temp(1,4,2))
        equivalence (jan_tmp_75(1),kelv_temp(1,5,1)),
     1      (jul_tmp_75(1),kelv_temp(1,5,2))
        equivalence (jan_rel_hum_15(1),iwvmr(1,1,1)),
     1      (jul_rel_hum_15(1),iwvmr(1,1,2))
        equivalence (jan_rel_hum_30(1),iwvmr(1,2,1)),
     1      (jul_rel_hum_30(1),iwvmr(1,2,2))
        equivalence (jan_rel_hum_45(1),iwvmr(1,3,1)),
     1      (jul_rel_hum_45(1),iwvmr(1,3,2))
        equivalence (jan_rel_hum_60(1),iwvmr(1,4,1)),
     1      (jul_rel_hum_60(1),iwvmr(1,4,2))
        equivalence (jan_rel_hum_75(1),iwvmr(1,5,1)),
     1      (jul_rel_hum_75(1),iwvmr(1,5,2))
c--------------------------------------------------------------------------------
c    data statements
c-------------------------------------------------------------------------------
        data temp_pres/
     1 .1,.2,.5,1.,1.5,2.,3.,4.,5.,7.,10.,15.,20.,25.,30.,50.,60.,
     2 70.,85.,100.,115.,135.,150.,200.,250.,300.,350.,400.,430.,475.,
     3 500.,570.,620.,670.,700.,780.,850.,920.,950.,1000./

c * 15n annual
        data jan_tmp_15/
     1 233,248,265,270,266,261,257,249,246,240,235,229,225,222,219,209,
     2 205,201,197,193,198,205,209,220,231,240,246,253,257,262,265,272,
     3 276,280,283,288,292,297,299,301/
c * 30n january
        data jan_tmp_30/
     1 235,250,264,269,263,257,250,245,242,236,230,225,222,219,216,209,
     2 207,204,203,205,207,210,212,216,225,233,240,246,250,254,257,264,
     3 268,272,275,280,283,285,287,288/
c * 45n january
        data jan_tmp_45/
     1 242,252,265,263,255,248,239,233,229,221,218,216,6*215,216,3*217,
     2 218,2*219,226,232,237,241,245,247,253,257,260,262,265,267,269,
     3 270,272/
c * 60n january
c
        data jan_tmp_60/
     1 250,254,260,248,241,236,230,225,221,218,216,213,2*211,212,214,
     2 2*215,216,7*217,223,229,232,237,239,246,250,252,253,256,258,259,
     3 2*258/
c * 75n january
        data jan_tmp_75/
     1 240,246,252,240,233,228,222,217,213,210,7*208,209,210,2*211,212,
     2 213,214,215,216,221,226,228,232,234,239,242,245,247,251,253,251,
     3 250,249/
c * 15n annual (repeated for symmetry)
        data jul_tmp_15/
     1 233,248,265,270,266,261,257,249,246,240,235,229,225,222,219,209,
     2 205,201,197,193,198,205,209,220,231,240,246,253,257,262,265,272,
     3 276,280,283,288,292,297,299,301/
c * 30n july
        data jul_tmp_30/
     1 231,248,266,272,268,263,255,250,246,241,234,229,226,223,220,214,
     2 211,209,207,204,203,204,209,221,232,241,248,255,259,264,267,273,
     3 277,280,283,288,292,297,300,303/
c * 45n july
        data jul_tmp_45/
     1 232,251,268,276,273,267,259,254,250,244,237,232,228,226,224,220,
     2 219,218,5*216,221,230,238,245,252,255,260,262,268,273,277,279,
     3 285,289,292,293,296/
c * 60n july
        data jul_tmp_60/
     1 226,250,273,277,275,273,265,259,253,245,238,234,231,229,227,
     2 10*225,230,238,244,248,253,256,263,266,270,271,276,280,284,285,
     3 288/
c * 75n july
        data jul_tmp_75/
     1 230,245,268,272,270,268,260,254,248,240,238,237,234,233,231,
     2 8*230,228,227,229,235,242,245,249,252,258,263,267,270,273,275,
     3 276,277,279/
c
c * 'iwvmr' contains relative humidity, organized like temperature

        data jan_rel_hum_15/
     1 3*15,2*20,25,4*30,5*35,5*75/
        data jan_rel_hum_30/
     1 3*15,2*20,6*30,2*35,40,45,50,60,70,75,80/
        data jan_rel_hum_45/
     1 3*15,20,30,2*35,2*40,45,2*50,2*55,60,65,2*70,73,77/
        data jan_rel_hum_60/
     1 5*15,20,40,45,50,25,20,50,60,2*65,2*70,2*75,80/
        data jan_rel_hum_75/
     1 4*15,20,37,40,42,2*45,2*50,52,55,20,2*60,65,70,80/
        data jul_rel_hum_15/
     1 3*15,2*20,25,4*30,5*35,5*75/
        data jul_rel_hum_30/
     1 3*15,20,2*30,2*35,2*40,2*45,50,55,2*60,2*65,75,80/
        data jul_rel_hum_45/
     1 3*15,20,25,6*30,35,40,45,50,55,60,65,70,75/
        data jul_rel_hum_60/
     1 3*15,20,30,35,40,43,47,50,52,57,60,20,65,3*70,72,75/
        data jul_rel_hum_75/
     1 4*15,20,30,35,37,40,45,47,50,55,20,40,65,20,75,80,85/
c
        data lr/40,37,35,31,28,26,25,24,23,20,18,16,15,13,11/
c
c * 1000 mb height
        data ihite/116,175,142,101,98, 116,119,116,84,101/
c------------------------------------------------------------------------------
c    in line function definitions
c------------------------------------------------------------------------------
        fs2(i)=.01*(i) ! replaced float with mixed mode 2/23/99 db
c        is2(f)=100.*f+.5  !not used
c...............................begin climate........................................

        alat=abs(v_lat)

c.......set pressure data values to variable name

        do i=1,40
            standard_press(i)=temp_pres(i)
        enddo

c.......check for too small of latitude

        if(alat.gt.15.) go to 100

        jl=1            !set <15 latitude to 15 lat coordinates
        kk=1
        go to 200

c.......check for too large of latitude

100     if(alat.lt.75.) go to 300

        jl=5            !set >75 latitude to 75 lat coordinates
        kk=2

c.......set temporary values for future computation

200     do k=1,kk
            do i=1,40
                temp_tmp(i,k)=kelv_temp(i,jl,k)
            enddo       !i loop

            do i=1,20
                temp_mixrat(i,k)=fs2(iwvmr(i,jl,k))
            enddo       !i loop

            temp_hgt(k)=ihite(jl,k)
        enddo           !k loop

        if(kk.ne.1) go to 400

        do i=1,40
            tempertur_guess(i)=temp_tmp(i,1)
        enddo   !i loop

        do i=1,20
c            mixratio_guess(i) = .001  !old statement
            mixratio_guess(i) = .001 + (i-1)*.000052631 !ramp up
            mixratio_guess(i+20)=temp_mixrat(i,1)
        enddo   !i loop

        z10=temp_hgt(1)
        go to 500

c.......determine month for proper indexing

300     jlat=nint(alat)
        jl1=jlat/15
        jl2=jl1+1
        wt1=(float(jl2*15)-alat)/15.    !set weights
        wt2=1.-wt1

        do k=1,2
            do  i=1,40
                it1=kelv_temp(i,jl1,k)
                it2=kelv_temp(i,jl2,k)
                temp_tmp(i,k)=wt1*float(it1)+wt2*float(it2)
            enddo       !i loop

            do i=1,20
                iw1=iwvmr(i,jl1,k)
                iw2=iwvmr(i,jl2,k)
                temp_mixrat(i,k)=wt1*fs2(iw1)+wt2*fs2(iw2)
            enddo       !i loop

            iz1=ihite(jl1,k)
            iz2=ihite(jl2,k)
            temp_hgt(k)=wt1*float(iz1)+wt2*float(iz2)
        enddo           !k loop

c.......determine calendar day for computation of data

400     kday=mod(julian_day,100000)

        call calday(kday,nyr,nmon,nday,lmon)

        if(v_lat.lt.0.) nmon=nmon+6
        if(nmon.gt.12) nmon=nmon-12
        imon=iabs(nmon-7)
        wt1=float(imon)/6.      !recalculate weights
        wt2=1.-wt1

c.......calculate guess temperatures and set to variable name

        do i=1,40
            tempertur_guess(i)=wt1*temp_tmp(i,1)+wt2*temp_tmp(i,2)
        enddo   !i loop

        do i=1,20
            j=i+20
            wx=wt1*temp_mixrat(i,1)+wt2*temp_mixrat(i,2)
c            mixratio_guess(i) = 0.001  !old statement
c            mixratio_guess(i) = .001 + (i-1)*.000052631 !ramp up
            mixratio_guess(i) = 0.0
            mixratio_guess(j)=wx
        enddo   !i loop

        z10=wt1*temp_hgt(1)+wt2*temp_hgt(2)

500     continue

        iz10=nint(z10)

c.......compute mixing ratios from relative humidities
c       and set to variable name

        do j=21,40
            ws=sat_pt_sm(standard_press(j),tempertur_guess(j))
            relative_hum=mixratio_guess(j)
            mixratio_guess(j)=relative_hum*ws
        enddo   !j loop

        istatus = 1
999     return
        end             !climate

        subroutine det_mix_rat_sm(press, temp, dewpt_dep,
     1           mix_ratio, nl)
c
c               07-jun-84       source code from wang
c                               (not on u. of w. source tape)
c
c               20-aug-84       modified

c               dewpt_dep       - dewpoint depression
c-------------------------------------------------------------------------------
c    variable declarations
c-------------------------------------------------------------------------------

        integer nl
        real  press(nl),
     1  temp(nl),
     1  dewpt_dep(nl),
     1  mix_ratio(nl),
     1  temp_dew                !dew pt temperature


        do 130 i = 1, nl
            temp_dew = temp(i) - dewpt_dep(i)
            if (temp_dew .gt. 253.0) go to 110
            es = vpice_sm(temp_dew)
            go to 120
110         es = satvap_sm(temp_dew)
120         mix_ratio(i) = 622.0 * es / press(i)
130     continue

        return
        end             !det_mix_rat_sm

      real function sat_pt_sm(press, temp)
c     
c     14-may-84       added from csu vas code
c     (not on u. of w. source tape)
      real  press(1),         !pressure
     1     temp(1),             !temperature
     1     depression(1),
     1     mix_ratio(1)         !mixing ratio
c-------------------------------------------------------------------------------
      
      depression(1) = 0.
      call det_mix_rat_sm(press,temp,depression, !calculate mixing ratio
     1     mix_ratio,1)
      sat_pt_sm=mix_ratio(1)
c     
c     14-may-84       return added
c     
      return
      end

        real function vpice_sm(temp)
c
c               07-jun-84       source code from wong
c                               (not on u. of w. source tape)
c
        dimension c(6)
c
        save
        real*8 c, t
c
        data c/ 0.7859063157d+00,  0.3579242320d-01,
     1      -0.1292820828d-03,  0.5937519208d-06,
     2       0.4482949133d-09,  0.2176664827d-10/
c
        t = temp - 273.16
        vplog = c(1) + t * (c(2) + t * (c(3) + t * (c(4) + t * (c(5) +
     1              t *  c(6)))))
        vpice_sm = 10.0 ** vplog
        return
        end

        function satvap_sm(temp)
c
c               07-jun-84       source code from wang
c                               (not on u. of w. source tape)
c
        dimension c(10)
c
        save
        real*8 c, s, t
c
        data c/ 0.9999968760d-00, -0.9082695004d-02,
     1       0.7873616869d-04, -0.6111795727d-06,
     2       0.4388418740d-08, -0.2988388486d-10,
     3       0.2187442495d-12, -0.1789232111d-14,
     4       0.1111201803d-16, -0.3099457145d-19/
        t = temp - 273.16
        s = 6.1078000000d+00 / (c(1) + t * (c(2) + t * (c(3) +
     1                              t * (c(4) + t * (c(5) +
     2                              t * (c(6) + t * (c(7) +
     3                              t * (c(8) + t * (c(9) +
     4                              t * c(10)))))))))) ** 8
        satvap_sm = s
        return
        end

        subroutine calday(nyrday,kyear,nmon,kalday,mon)

c   to obtain calendar date from 'nyrday' which is of the form 'yyddd'
c   kyear,nmon,kalday are binary integer year,month,day.  mon is
c   the month expressed in left-justified three-letter code.
        
        save
        dimension months(12)
        data months/4hjan ,4hfeb ,4hmar ,4hapr ,4hmay ,4hjun ,
     1  4hjul ,4haug ,4hsep ,4hoct ,4hnov ,4hdec /


        kyear = nyrday/1000
        julian_day = mod(nyrday,1000)
        jcode = 1115212
        if(mod(kyear,4) .eq. 0) jcode = 1115208
        do  m = 1,12
            month = m
            ndm = 31 - mod(jcode,4)
            julian_day = julian_day - ndm
            if(julian_day .le. 0) go to 120
            jcode = jcode/4
        enddo

 120    kalday = julian_day + ndm
        mon = months(month)
        nmon = month
        return
        end
