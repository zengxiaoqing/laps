 &snd_nl
 path_to_raw_raob='/public/data/raob/netcdf',
 path_to_local_raob='/data/fxa/LDAD/rawinsonde/netCDF',
 path_to_raw_tower='/public/data/tower/public/netcdf',
 path_to_raw_drpsnd='/public/data/dropsonde/netcdf',
 path_to_raw_poessnd='/public/data/madis/point/POES/netcdf',
 path_to_raw_satsnd='/public/data/sat/nesdis/goes12/sounding/binary/',
 path_to_raw_radiometer='/public/data/madis/point/radiometer/netcdf',
 l_fill_ht=.true.,
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
c Note that dropsonde 'drpsnd' data works for NIMBUS, AVAPS, CWB, and SND 
c formats at present. If SND format is the input, then 'c8_project' should 
c be set to 'AIRDROP' in 'nest7grid.parms'. When we have the ASCII SND file in 
c one directory being the input, and SND as output in a second directory, the 
c software's main purpose is applying the time windowing to the observations.
c
c
c 'path_to_raw_poessnd'
c
c The 'poessnd' data is for MADIS POES data (if 'c8_project' is not AFWA)
c
c
c 'path_to_raw_satsnd'
c
c For most values of 'c8_project' this refers to the path containing GOES 
c sounding data. If 'c8_project' is set to 'AFWA' the path refers to AFWA
c satellite sounding data instead.
c
c
c 'l_fill_ht'
c
c Flag to decide whether to calculate heights from the hypsometric equation if
c they aren't reported in the raw sounding data
