 &cloud_drift_nl
 n_paths_drift=3,
 path_to_cloud_drift='/public/data/sat/nesdis/goes8/wind/cloud-drift/ascii',
                     '/public/data/sat/nesdis/goes9/wind/cloud-drift/ascii',
                     '/public/data/sat/nesdis/goes10/wind/cloud-drift/ascii',
 cloud_drift_format='NESDIS','NESDIS','NESDIS',
 /

c Cloud Drift ingest (ingest_cloud_drift.exe)
c
c 'n_paths_drift' - Can be up to 10.
c
c 'path_to_cloud_drift' - Directory for each source of data
c
c 'cloud_drift_format' - Format of data contained in each directory
c                        Valid values are 'NESDIS', 'CWB_SATOB', 'CWB_HDSW', 
c                        and 'AFWA'
