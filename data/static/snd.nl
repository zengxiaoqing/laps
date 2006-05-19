 &snd_nl
 path_to_raw_raob='/public/data/raob/netcdf',
 path_to_local_raob='/data/fxa/LDAD/rawinsonde/netCDF',
 path_to_raw_tower='/data/fxa/LDAD/mesonet/met-tower/netCDF',
 path_to_raw_drpsnd='/null/public/data/dropsonde/netcdf',
 path_to_raw_poessnd='/null/public/data/madis/point/POES/netcdf',
 path_to_raw_satsnd='/public/data/sat/nesdis/goes12/sounding/binary/',
 path_to_raw_radiometer='/public/data/madis/point/radiometer/netcdf',
 /
c SOUNDING ingest (ingest_sounding.exe)
c
c
c 'path_to_raw_raob'
c
c This works for NIMBUS, WFO, MADIS, CWB, and AFWA formats. If WFO or MADIS
c data are being used, the parameter 'c8_project_common' should be set to 'WFO'
c in the 'nest7grid.parms' namelist.
c
c
c 'path_to_raw_drpsnd'
c
c Note that dropsonde 'drpsnd' data works for AVAPS, CWB, and SND formats at 
c present. If SND format is the input, then 'c8_project' should be set to
c 'AIRDROP' in 'nest7grid.parms'. When we have the ASCII SND file in one 
c directory being the input, and SND as output in a second directory, the 
c software's main purpose is applying the time windowing to the observations.
c
c
c 'path_to_raw_poessnd'
c
c The 'poessnd' data is for MADIS POES data (if 'c8_project' is not AFWA)
c The default is for this to be "off" as this is still somewhat experimental.
