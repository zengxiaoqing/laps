      subroutine sat_sublatlon(ewc4,ewi4,nsc4,nsi4,f_time,imc,orbAt,
     &SatSubLAT,SatSubLON,istatus)
c
c
c
      Implicit None

      include 'instco.inc'

      Real*8        orbAt(336)
      Real          time_50,time50
      Real*8        t50_8
      Real*8        t
      Real*8        f_time
      Real*8        SatSubLAT,SatSubLON
      Real          pi,radtodeg

      Integer     ewc4,ewi4
      Integer     nsc4,nsi4
      Integer     INSTR
      Integer     time_spec(2)
      Integer     imc
      Integer     istatus

      istatus = -1
      INSTR=1          !1=Imager, 2=Sounder

      pi=3.141592653589793
      radtodeg=180.0/pi

      call bcd_to_int(orbAt(12),time_spec)
      time_50 = time50(time_spec)
      t50_8=time_50
      t = f_time /60. + 7305. * 24. * 60.

      call SETCON(INSTR,nsc4,nsi4,ewc4,ewi4)
      call LMODEL(t,t50_8,OrbAt,imc,SatSubLAT,SatSubLON)

      write(6,*)'  Sat Subpoint lat (deg) ',SatSubLAT*radtodeg
      write(6,*)'  Sat Subpoint lon (deg) ',SatSubLON*radtodeg

      istatus=1

      return
      end
