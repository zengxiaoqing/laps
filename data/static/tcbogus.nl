&tcbogus_nl
bogus=.false.
itc_no=1
dskm=15.
year=2002
month=07
day=09
hour=12.
cen_lat=24.2,20.,
cen_lon=120.2,125.,
alpha_i=0.85,0.8,
vmax=18.,60.,
rmax=50.,60.,
r_7deg=80.,300.,
speed=15.,20.,
direc=30.,45.,
&end

c tcbogus_nl variable descriptions:
c BOGUS =   The flag control whether do bogus or not.
c ITC_NO =  The number of the typhoon on the LAPS domain,now at most 3 typhoons.
c DSKM =    Background model  resolution (unit: km).
c           Following 4 parameters are date (utc)
c YEAR = 
c MONTH = 
c DAY =  
c HOUR =  
c CEN_LAT = Typhoon center for latitude
c CEN_LON = Typhoon center for longitude
c ALPHA_I = The exponential parameter; specifies the horizontal wind profile outside the RMW(Radius Maximum Wind) of the Rankine vortex.
c VMAX =    Typhoon maximum wind speed (unit:  m/s)
c RMAX =    The radius of typhoon maximum speed (unit: km)
c R_7DEG =  The radius of  15 m/s. (unit : km )
c SPEED =   The speed of typhoon  movement  (unit: km/hour )
c DIREC =   The direction of typhoon movement (unit : degree 0 ~ 359 )
c Guo-Ji  designed this bogus  program by Rankine vortex. By his experiences,
c ALPHA_I = 0.7, or 0.8, or 0.9 (small typhoon ---> big typhoon),RMAX=60,or 50
c (weak typhoon ---> strong typhoon)


