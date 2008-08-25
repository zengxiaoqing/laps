 &vrc_nl
 c_rawdata_path = '/public/data/radar/wsi/nowrad/netcdf/',
 aoml_path_in = '/nodata/home/fab/albers/aoml/sweep/',
 vrc_outdir = 'lapsprd',
 iwsimsng = 255,
 icheckint = 30,
 iwaittime = 180,
 iagethresh = 3601,
 /
c
c
c iagethresh - Only loop through the waiting if most recent unprocessed input
c              data is younger than this threshold age (seconds)
c
c vrc_outdir - set to 'lapsprd' if there's only one VRC output directory
c              set to 'rdr' if there are multiple VRC output directories
c
c Note: Program will not rewrite any output VRC files for times earlier than 
c       the latest one initially present.