 &cloud_drift_nl
 n_paths_drift=2,
 path_to_cloud_drift='/public/data/madis/point/HDW1h/netcdf',
                     '/public/data/madis/point/HDW/netcdf',
 cloud_drift_format='MADIS','MADIS',
 /

c Cloud Drift ingest (ingest_cloud_drift.exe)
c
c 'n_paths_drift' - Can be up to 10.
c
c 'path_to_cloud_drift' - Directory for each source of data
c
c 'cloud_drift_format' - Format of data contained in each directory
c                        Valid values are 'NESDIS', 'CWB_SATOB', 'CWB_HDSW', 
c                        'MADIS', and 'AFWA'
