&lfmpost_nl

out_cdf=.true.,
out_grib=.false.,
out_v5d=.false.,

make_micro=.true.,
make_firewx=.true.,

verbose=.true.,
realtime=.true.,
write_to_lapsdir=.false.,
make_donefile=.false.,

precip_dt=10800,
c_m2z='rams',

n3d_pts_thr=99000000,
p3d_pts_thr=99000000,

&end

c c_m2z various options for reflectivity calculation 
c 'kessler', 'wrf', 'rams', synpolrad 'synp', future additions are 
c thompson 'thom' and others. 
c
c n3d_pts_thr: threshold number of input 3D grid points for reduced processing
c
c p3d_pts_thr: threshold number of output 3D grid points for reduced processing
