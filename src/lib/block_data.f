      block data

c     routine lapsgrid.f

        include 'lapsparms.cmn'
        data iflag_lapsparms_cmn /0/

c     routine laps_thermo.f

        COMMON/INDX/ P(70),T(70),TD(70),HT(70),PBECR(20,4),TDFCR(20,2)
     1              ,VEL(20),temdif(70),partem(70),pbe(70)
     #              ,DD85,FF85,DD50,FF50
        DATA PBE/70*0./


c     routine supmap.f

        COMMON/SUPMP1/PI,TOVPI,DTR,RTD,EPS,OV90,CON1,CON2,PART
        COMMON/SUPMP2/NPTS,MAXLAT,MINLAT,MAXLON,MINLON,PTS(200)
        COMMON/SUPMP3/POLONG,CONE,RLAT,RLON,JGR,ILF,SGN
        COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD
        COMMON/SUPMP5/PHIOC,SINO,COSO,SINR,COSR,IPROJ
        COMMON/SUPMP6/UMIN,UMAX,VMIN,VMAX,UEPS,VEPS
        COMMON/SUPMP7/PHIO,PHIA,IGRID,IDOT,ILTS
        COMMON/SUPMP8/U,V,U1,V1,U2,V2
        COMMON/SUPMPA/IIER
        COMMON/MAPCOL/MPCOL1,MPCOL2,MPCOL3,MPCOL4
        COMMON/MAPDAS/LDASH1,LDASH2,LDASH3,LDASH4
c        DATA                    !Default line intensities and dash patterns
c     1          MPCOL1,LDash1   /255,'1777'O/,  !Map lines
c     2          MPCOL2,LDash2   /128,'1756'O/,  !Grid lines
c     3          MPCOL3,LDash3   /192,'1777'O/,  !Limb lines
c     4          MPCOL4,LDash4   /255,'1777'O/   !Perimeter
        DATA                    !Default line intensities and dash patterns
     1          MPCOL1,LDash1   /255,1023/,  !Map lines
     2          MPCOL2,LDash2   /128,1006/,  !Grid lines
     3          MPCOL3,LDash3   /192,1023/,  !Limb lines
     4          MPCOL4,LDash4   /255,1023/   !Perimeter


c       COMMON/SUPMP1/PI,TOVPI,DTR,RTD,EPS,OV90,CON1,CON2,PART
c       COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD
        COMMON/SUPMP9/DS,DI,DSRDI

        DATA   CON1 / 1.00001/
        DATA   CON2 / 179.99999/
        DATA     DI / 16./
        DATA    DTR / 1.7453292519943E-2/
        DATA    EPS / 1.E-6/
        DATA   OV90 / 1.11111111111111E-2/
        DATA     PI / 3.1415926535898/
        DATA    RTD / 57.295779513082/
        DATA  TOVPI / 0.63661977236758/
        DATA   UOLD / 0.0 /
        DATA   VOLD / 0.0 /
        DATA PART   / 1.0 /          !SIZE OF PICTURE (90% OF SCREEN)
  
c     routine xsect.f
        logical l_convert
        common/lapsplot_omega/l_convert
        data l_convert /.true./

c     routine config_satellite_lvd (in file lapsgrid.f).
        include 'satellite_dims_lvd.inc'
        include 'satellite_common_lvd.inc'
        include 'sat_data_static_lvd.inc'
        data iflag_lvd_common /0/ 

c     routine config_ingest_rrv_parms
        include 'ingest_rrv_dims.inc'
        include 'ingest_rrv_common.inc'
        include 'ingest_rrv_constants.dat'
      end

