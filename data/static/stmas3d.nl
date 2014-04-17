&number_states
numstat=5
/

&stmas3d
ifbkgnd=0,
ifbound=1,
ifpcdnt=1,

penal0x=2.0,2.0,2.0,2.0,2.0,
penal0y=2.0,2.0,2.0,2.0,2.0,
penal0z=0.05,0.05,0.05,0.05,0.05,
penal0t=0.01,0.01,0.01,0.01,0.01,
pnlt0pu=0.1,
pnlt0pv=0.1,

numdims=4,
numgrid=0,0,0,0,
maxgrid=0,0,0,0,
fnstgrd=0,

w_cmpnnt=0,
u_cmpnnt=1,
v_cmpnnt=2,
pressure=3,
temprtur=4,
humidity=5,

cosstep=5,
midgrid=6,
finstep=5,
pnlt0hy=1.0,
taul_hy=0.5,
endhylv=0,
endgslv=0
/

c numstat:  total number of analysis state variables
c ifbkgnd:  option for use of background fields: 0 no; 1 yes
c ifbound:  option for apply a bound on analysis: 0 no; 1 yes
c ifpcdnt:  option for vertical coordinate: 1 isobaric
c
c penal0*:  smoothing parameter for all variables and dimensions
c pnlt0o*:  penalty parameters for geostrophic balance in x and y
c
c numdims:  dimensions of the analysis; design to run 2 or 3 D analysis
c numgrid:  numbers of gridpoints on the coarsest multigrid
c maxgrid:  numbers of gridpoints on the finest multigrid, 
c           if maxgrid(1) = 0, stmas uses a default multigrid setup
c fnstgrd:  total number of multigrid levels
c
c obradius:  observation influence radius if specified (no hardcoded in input_bg_obs.f90
c
c w_cmpnnt...: Indices of the number of analysis variables
c
c cosstep:  Number of minimization iterations requested
c midgrid:  a multigrid level where users can set different minimization iterations
c           any grid coarser than midgrid level uses cosstep iterations; otherwise, uses
c           finstep
c finstep:  number of minimization for multigrid levels finer than midgrid
c pnlt0hy:  starting weight for penalizing hydrostatic balance
c taul_hy:  ratio of reducing the hydrostatic weight between two adjacent multigrid levels
c endhylv:  the level after which no hydrostatic balance is applied
c endgslv:  the level after which no geostrophic balance is applied
c
