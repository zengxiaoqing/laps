 &pro_nl
 c8_blp_format='default',
 /

c PROFILER INGEST PARAMETERS
c
c 'c8_blp_format' - A value of 'default' means that we are using the 
c                  'c8_project' in 'nest7grid.parms' to specify the boundary
c                  layer profiler format. An override to this can be specified 
c                  as follows:
c
c                  'NIMBUS' denotes GSD NetCDF format and NIMBUS file timestamp
c                  'WFO' denotes MADIS Multi-Agency Profiler NetCDF format and
c                        WFO file timestamp
c                  'MADIS' denotes MADIS Multi-Agency Profiler NetCDF format and
c                        WFO file timestamp
c                  'RSA' is the RSA/LDAD format available on the RSA version
c                        of AWIPS
c