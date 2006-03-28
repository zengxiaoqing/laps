#############################################################################
#  map_utils.pm
#
#  This is a Perl module for converting (i,j) to (lat,lon) and vice versa
#  for any of 4 map projections:  Cylindrical Equidistant (lat-lon), 
#  Mercator, lambert conformal (secant and tangent), and polar 
#  stereographic.  
#
#  Usage
#
#    1.  First, you need to know the basic parameters of your map
#        projection and grid:
#          a.  TYPE:Projection type ("LC", "PS", "LL", or "ME" for
#              Lambert Conformal, Polar-Stereographic, Lat/lon cylindrical
#              equidistant, or mercator, respectively.
#          b.  KNOWN_LAT/KNOWN_LON:  You need to know the lat/lon (deg N/E)   
#              of one point on your grid.  Typically, this will be the center
#              or SW corner point, but the module allows the use of any
#              known point on the grid.
#          c.  KNOWN_RI/KNOWN_RJ:  You have to know the real (i,j) coordinate
#              of your known lat/lon.  The value ranges from (1:nx,1:ny), 
#              where nx and ny are the number of grid points in the E-W and
#              N-S direction, respectively.  For example, if KNOWN_LAT/
#              KNOWN_LON represent the SW corner, then KNOWN_RI/KNOWN_RJ
#              should be set to (1.0,1.0).  If you are providing the
#              center point, you should use ( 0.5*(nx+1), 0.5*(ny+1) ).
#              should be the SW corner of the grid
#          d.  DX:The grid spacing in meters at the true latitude (not used
#              for lat/lon grids)
#          e.  TRUELAT1: The latitude at which the grid spacing is true
#              (polar-stereo,
#              lambert conformal, and mercator grids only).  For lat/lon grids,
#              the first truelat value is set to the latitudinal grid spacing
#              in degrees
#          f.  TRUELAT2:The second true latitude if a lambert conformal
#          g.  STDLON:The standard longitude (LC and PS grids only)...
#              for lat/lon 
#              grids the standard longitude is set to the longitudinal
#              grid increment in degrees
#          h.  NX:The number of east-west points
#          i.  NY:The number of north south points
#            

#
#    2.  Once you know all of the above, you call the map_set subroutine,
#        passing in all of the above arguments (even if not used for your
#        particular projection type.  This routine returns a hash:
#
#       my %proj = &map_utils::map_set($TYPE,$KNOWNLAT,$KNOWNLON, 
#                   KNOWN_RI, KNOWN_RJ,$DX,$STDLON,
#                   $TRUELAT1,$TRUELAT2,$NX,$NY);
#   
#    3.  The %proj hash is then used as an input argument to the 
#        coordinate conversion routines.  To convert an (i,j) to
#        a lat/lon, you would call:
#
#        my ($lat,$lon) = &map_utils::ij_to_latlon($i,$j,%proj);
#
#        Or, to convert a lat/lon to an i/j:
#   
#        my ($i,$j) = &map_utils::latlon_to_ij($lat,$lon,%proj);
#
#  History:
#
#   3 Sep 2002:  Initial version checked into LAPS repository..B. Shaw
#
#
############################################################################

package map_utils;
require 5;
use strict;

# We need the trig module for things like pi, deg2rad, rad2deg,
# and other trig functions.

use Math::Trig;
use POSIX;
# Define earth constant

my $earth_radius_m = 6371200.0;

# Define hash referencing short character codes to their
# full description

my %proj_descr = ( LL => "CYLINDRICAL EQUIDISTANT",
                   ME => "MERCATOR",
                   LC => "LAMBERT CONFORMAL",
                   PS => "POLAR STEREOGRAPHIC",
                   RL => "ROTATED LAT-LON"
                 );

# Here is a hash cross-referencing the character string
# codes for each projection to the NCAR Graphics code.

my %proj_codes = ( LL => 0,
                   ME => 1,
                   LC => 3,
                   PS => 5,
                   RL => 3
                 );


#--------------------------------------------------------------------------
#  Sub for returning a hash to be populated
#--------------------------------------------------------------------------

sub map_hash {
  my %map_hash = ( proj => "XX",
                   latsw => -999.9,
                   lonsw => -999.9,
                   latne => -999.9,
                   lonne => -999.9,
                   latnw => -999.9,
                   lonnw => -999.9,
                   latse => -999.9,
                   lonse => -999.9,
                   latcen => -999.9,
                   loncen => -999.9,
                   dx => -999.9,
                   stdlon => -999.9, 
                   truelat1 => -999.9,
                   truelat2 => -999.9,
                   hemi => 0,
                   cone => -999.9,
                   polei => -999.9,
                   polej => -999.9,
                   rsw => -999.9,
                   rebydx => -999.9,
                   nx => -99,
                   ny => -99,
                   dellon => -999.9,
                 );
  return %map_hash;

} 
###########################################################################
#  Main map_set driver to set up map proj hash table
##########################################################################
sub map_set($$$$$$$$$$$) {

  # Get the input parameters
  my ($type, $knownlat, $knownlon, $kri, $krj, $dx, $stdlon, 
      $truelat1, $truelat2, $nx, $ny) = @_;

     $type = 'LC' if $type =~ /RL/i;

  my %proj = map_set_sub($type, $knownlat, $knownlon, $dx, $stdlon,
                         $truelat1,$truelat2,$nx,$ny);

  if (($kri ne 1.000) or ($krj ne 1.000)) {
    my $rswi = 2.0 - $kri;
    my $rswj = 2.0 - $krj; 
    my ($latsw,$lonsw) = ij_to_latlon($rswi,$rswj,%proj);
    
    # Call map_set again...
    %proj = map_set_sub($type, $latsw, $lonsw, $dx, $stdlon,
                         $truelat1,$truelat2,$nx,$ny);
  }
  return %proj;
}


# -------------------------------------------------------------------------
#  Sub for setting up hash table defining map structure
# -------------------------------------------------------------------------

sub map_set_sub($$$$$$$$$) {

  # Get the input parameters
  my ($type, $latsw, $lonsw, $dx, $stdlon, $truelat1, $truelat2,
      $nx, $ny) = @_;

  # Get the hash template
  my %proj = &map_hash();
 
  # Start populating the hash, doing some validity checks based 
  # on projection type

  $type = uc $type;

  # Check to make sure projection type is supported
  if ( ! exists $proj_descr{$type}) {
    print "map_set: Invalid projection type specified: $type\n";
    print "  Allowed types: \n";
    my @projection = each %proj_descr ;
    my $check = @projection;
    while ($check > 0) {
      print "  @projection\n";
      @projection = each %proj_descr ;
      $check = @projection;
    }
    die;
  }else{
   ${proj{proj}} = $type;
  }

  # Ensure latitude of southwest corner is between -90 and 90 degrees

  if (abs($latsw) > 90.0) {
    print "map_set: Invalid lat: $latsw\n";
    print "  -90.0 <= latsw <= 90.\n";
    die;
  }else{
    ${proj{latsw}} = $latsw;
  }

  # Ensure longitude of southwest corner is between -180 and 180 degrees
  if (abs($lonsw) > 180.0){
    print "map_set: Invalid SW corner lon: $lonsw\n";
    print "  -180 <= lonsw <= 180\n";
    die;
  }else{ 
   ${proj{lonsw}} = $lonsw;
  }

  # If not a latlon grid, ensure dx is positive.  Dx is not used
  # in the case of a latlon grid
  if ( ($type ne "LL") and ($dx <= 0)){
    print "map_set: Invalid dx: $dx\n:";
    print "  For projections other than lat/lon, dx must be set to\n";
    print "  a positive value in meters.\n";
    die;
  }else{
   ${proj{dx}} = $dx;
  }

  # Check stdlon, which should be set for all projection types
  # except mercator.  In the case of a latlon grid, this is set
  # to an increment
  if ($type ne "ME"){
    if (abs($stdlon) > 180.) {
      print "map_set: Invalid standard longitude: $stdlon\n";
      print "  For PS and LC : -180 <= stdlon <= 180\n";
      print "  For LL, this should be the delta-lon value\n";
      die;
      if ($type eq "LL"){
        if ($stdlon eq 0) {
          print "map_set:  Invalid delta-longitude for LL grid\n";
          die;
        }
      }
    }else{
      ${proj{stdlon}} = $stdlon;
    }
  }else{
    ${proj{stdlon}} = 0;
  } 
 
  # All projections use truelat1.  In the case of a Latlon grid, however,
  # this is really the delta-lat parameter.

  if (abs($truelat1) > 90) {
    print "map_set:  Invalid true latitude 1: $truelat1\n";
    die;
  }else{
    if (($type eq "LL") and ($truelat1 eq 0.)){
      print "map_set:  Delta-lat for LL grid must be non-zero!\n";
      die;
    }
  }
  ${proj{truelat1}} = $truelat1;

  # LC projection requires truelat2
  if ($type eq "LC"){
    if (abs($truelat2) > 90) {
      print "map_set:  Invalid true latitude 2: $truelat2\n";
      die;
     }else{
      ${proj{truelat2}} = $truelat2;
     }
  }

  # Check nx/ny
  if ($nx <= 1){
    print "map_set: Invalid nx: $nx\n";
    print "  Nx must be > 1\n";
    die;
  }else{
   ${proj{nx}} = $nx;
  }

  if ($ny <= 1){
    print "map_set: Invalid ny: $ny\n";
    print "  Ny must be > 1\n";
    die;
  }else{
   ${proj{ny}} = $ny;
  }

  # Fill in the rest of the proj hash
  
  # hemisphere parameter 
  if ($type ne "LL") {
    if ($truelat1 < 0.) {
      ${proj{hemi}} = -1.0;
    }else{
      ${proj{hemi}} = 1.0;
    }
    ${proj{rebydx}} = $earth_radius_m / $dx;
  } 

  # Case-dependent calls for final setup

  if ($type eq "PS") {   
    ${proj{cone}} = 1.0;
    my $reflon = ${proj{stdlon}} + 90.;
    my $scale_top = 1. + ${proj{hemi}} * sin(deg2rad(${proj{truelat1}}));
    my $ala1 = deg2rad(${proj{latsw}});
    ${proj{rsw}} = ${proj{rebydx}}*cos($ala1)*$scale_top/
                   (1.+${proj{hemi}}*sin($ala1)) ;

    # Find the pole point
    my $alo1 = deg2rad(${proj{lonsw}} - $reflon);
    ${proj{polei}} = 1. - ${proj{rsw}} * cos($alo1);
    ${proj{polej}} = 1. - ${proj{hemi}} * ${proj{rsw}} * sin($alo1);
    
  } elsif ($type eq "LC" ) {
    
    # Make sure truelat1 <= truelat2
    if ($truelat1 > $truelat2) {
       ${proj{truelat1}} = $truelat2;
       ${proj{truelat2}} = $truelat1;
       $truelat1 = ${proj{truelat1}};
       $truelat2 = ${proj{truelat2}};
    } 
    # Set cone factor 
    ${proj{cone}} = lc_cone(${proj{truelat1}},${proj{truelat2}});

    # Compute longitude differences to avoid the "cut" zone
    my $deltalon = ${proj{lonsw}} - ${proj{stdlon}};
    if ($deltalon > 180.) { $deltalon = $deltalon - 360. }
    if ($deltalon < -180.){ $deltalon = $deltalon + 360. }

    my $costl1 = cos(deg2rad($truelat1));

    # Radius to SW corner
    ${proj{rsw}} = ${proj{rebydx}} * $costl1/${proj{cone}} *
                   (tan(deg2rad(90.*${proj{hemi}}-${proj{latsw}})/2.) /
                    tan(deg2rad(90.*${proj{hemi}}-${proj{truelat1}})/2.) ) **
                    ${proj{cone}};
   
    # Find the pole point
    my $param1 = ${proj{cone}}*deg2rad($deltalon);
    ${proj{polei}} = 1. - ${proj{hemi}} * ${proj{rsw}} * sin($param1);
    ${proj{polej}} = 1. + ${proj{rsw}} * cos($param1); 
   
  }elsif($type eq "ME") {
    my $clain = cos(deg2rad(${proj{truelat1}}));
    ${proj{dellon}} = ${proj{dx}} / ( $earth_radius_m * $clain);
    ${proj{rsw}} =0.;
    if (${proj{latsw}} ne 0.) {
      ${proj{rsw}} = (log(tan(0.5*deg2rad(${proj{latsw}}+90.))))/${proj{dellon}};
    }

  }elsif($type eq "LL") {
    if (${proj{lonsw}} < 0.) { ${proj{lonsw}} = ${proj{lonsw}} + 360. }
  }else{
    print "Unknown projection: $type\n";
    return;
  } 

  # Call ij_to_latlon to fill in corners/center
  (${proj{latnw}},${proj{lonnw}}) = ij_to_latlon(1.,$ny,%proj);
  (${proj{latne}},${proj{lonne}}) = ij_to_latlon($nx,$ny,%proj);
  (${proj{latse}},${proj{lonse}}) = ij_to_latlon($nx,1.,%proj);
  (${proj{latcen}},${proj{loncen}}) = ij_to_latlon(($nx+1.)*.5,
                                                   ($ny+1.)*.5, %proj);
  return %proj;

}
########################################################################
sub latlon_to_ij(\$\$\%) {

  # Wrapper subroutine to compute the i/j point (1->nx,1->ny) from
  # a provided latitude and longitude for a given projection.  
  #
  # Arguments:
  #
  #     $lat = input latitude (float value in degrees N)
  #     $lon = input longitude (float value in degrees E)
  #     %proj = input map information hash as set up by map_set subroutine
  #
  #     Returns:  ($ri,$rj), real i/j values
  #
  #
  
  my ($ri,$rj);
  my ($lat,$lon,%proj) = @_;
  
  if ( ${proj{proj}} eq "PS" ) {
    ($ri,$rj) = llij_ps($lat,$lon,%proj);
  }elsif( ${proj{proj}} eq "LC") {
    ($ri,$rj) = llij_lc($lat,$lon,%proj);
  }elsif( ${proj{proj}} eq "LL") {
    ($ri,$rj) = llij_ll($lat,$lon,%proj);
  }elsif( ${proj{proj}} eq "ME") {
    ($ri,$rj) = llij_me($lat,$lon,%proj);
  }else{
    die "Unsupported projection type in latlon_to_ij: ${proj{proj}}\n";
  }

  return($ri,$rj);
}

#############################################################################
sub ij_to_latlon(\$\$\%) {

  # Wrapper subroutine to compute the latitude and longitude for 
  # a given map/grid projection and a given (i,j) coordinate.
  #
  # We define the coordinate to be (1,1) at the origin (SW corner) 
  #
  # Arguments:
  #
  #     $ri = input float i coordinate (E-W direction)
  #     $rj = input float j coordinate (N-S direction)
  #     %proj = input map information hash as set up by map_set subroutine
  #
  #     Returns:  ($lat,$lon), real lat/lon values in degrees N/E
  #

  my ($lat,$lon);
  my ($ri,$rj,%proj) = @_;

  if ( ${proj{proj}} eq "PS" ) {
    ($lat,$lon)=ijll_ps($ri,$rj,%proj);
  }elsif( ${proj{proj}} eq "LC") {
    ($lat,$lon)=ijll_lc($ri,$rj,%proj);
  }elsif( ${proj{proj}} eq "LL") {
    ($lat,$lon)=ijll_ll($ri,$rj,%proj); 
  }elsif( ${proj{proj}} eq "ME") {
    ($lat,$lon)=ijll_me($ri,$rj,%proj);
  }else{
    die "Unsupported projection type in latlon_to_ij: ${proj{proj}}\n";
  }

  return($lat,$lon);
}

############################################################################
sub llij_ps(\$\$\%proj) {

  # Subroutine to compute i/j from lat/lon for a polar-stereographic
  # grid projection.  Arguments are the same as for latlon_to_ij...in fact
  # ij_to_latlon calls this routine as necessary so the calling routine
  # has one interface to all required conversion routines.

  my ($ri,$rj);
  my ($lat,$lon,%proj) = @_;

  my $reflon = ${proj{stdlon}} + 90.;
  
  # Compute numerator term of map scale factor

  my $scale_top = 1. + ${proj{hemi}} * sin(deg2rad(${proj{truelat1}}));

  # Find radius to desired point
  my $ala = deg2rad($lat);
  my $rm = ${proj{rebydx}}*cos($ala)*$scale_top/(1.+${proj{hemi}}*sin($ala));
  my $alo = deg2rad($lon-$reflon);
  $ri = ${proj{polei}} + $rm * cos($alo);
  $rj = ${proj{polej}} + ${proj{hemi}} * $rm * sin($alo);
  
  return ($ri,$rj);

}

############################################################################
sub ijll_ps(\$\$\%proj) {

  # Subroutine to compute lat/lon from i/j for a polar-stereographic
  # grid projection.  Arguments are the same as for ij_to_latlon...in fact
  # ij_to_latlon calls this routine as necessary so the calling routine
  # has one interface to all required conversion routines.

  my ($lat,$lon);
  my ($ri,$rj,%proj) = @_;

  # Compute the reference longitude by rotating 90 degrees to the east
  # to find the longitude line parallel to the positive x-axis
  my $reflon = ${proj{stdlon}} + 90.;

  # Compute numerator term of map scale factor
  my $scale_top = 1. + ${proj{hemi}} * sin(deg2rad(${proj{truelat1}}));


  # Compute radius to point of interest
  my $xx = $ri - ${proj{polei}};
  my $yy = ($rj - ${proj{polej}}) * ${proj{hemi}};
  my $r2 = $xx**2 + $yy**2;

  # Now, the magic code
  if ($r2 eq 0) {
    $lat = ${proj{hemi}} * 90.;
    $lon = $reflon;
  }else{
    my $gi2 = (${proj{rebydx}} * $scale_top)**2.;
    $lat = rad2deg( ${proj{hemi}}*asin(($gi2-$r2)/($gi2+$r2)) );
    my $arccos = acos($xx/sqrt($r2));
    if ( $yy > 0) {
      $lon = $reflon + rad2deg($arccos);
    }else{
      $lon = $reflon - rad2deg($arccos);
    }
  }
 
  # Convert to a -180 -> 180 East convention

  if ($lon > 180.) {$lon = $lon - 360. }
  if ($lon < -180.) { $lon = $lon + 360. }

  return ($lat,$lon);
}
###############################################################################
sub llij_lc(\$\$\%) {
  #
  # Computes i/j from lat/lon for a Lambert Conformal Grid.
  #
 
  my ($ri,$rj);
  my ($lat,$lon,%proj) = @_;

  # Compute deltalon
  my $deltalon = $lon - ${proj{stdlon}};
  if ($deltalon > 180.) { $deltalon = $deltalon - 360. }
  if ($deltalon < -180.){ $deltalon = $deltalon + 360. }

  # Convert truelat1 to radians and get cosine for future use
  my $costl1 = cos(deg2rad(${proj{truelat1}}));

  # Compute the radius to the desired point
  my $rm = ${proj{rebydx}} * $costl1/${proj{cone}} *
                   (tan(deg2rad(90.*${proj{hemi}}-$lat)/2.) /
                    tan(deg2rad(90.*${proj{hemi}}-${proj{truelat1}})/2.) ) **
                    ${proj{cone}}; 

  my $param1 = ${proj{cone}} * deg2rad($deltalon);
  $ri = ${proj{polei}} + ${proj{hemi}} * $rm * sin($param1);
  $rj = ${proj{polej}} - $rm * cos($param1);

  # If in the southern hemisphere, we need to flip the i/j values such
  # that (1,1) is still the SW corner.

  if (${proj{hemi}} eq -1.) {
    $ri = 2. - $ri;
    $rj = 2. - $rj;
  }

  return ($ri,$rj);
}
#############################################################################
sub ijll_lc(\$\$\%) {

  # Subroutine to compute lat/lon from i/j for Lambert conformal maps.

  my ($lat,$lon);
  my ($ri,$rj,%proj) = @_;

  my $chi1 = deg2rad(90. - ${proj{hemi}}*${proj{truelat1}});
  my $chi2 = deg2rad(90. - ${proj{hemi}}*${proj{truelat2}});

  # Flip indices if we are in southern hemisphere
  if (${proj{hemi}} < 0.) {
    $ri = 2. - $ri;
    $rj = 2. - $rj;
  }

  # Compute square of radius to i/j
  my $xx = $ri - ${proj{polei}};
  my $yy = ${proj{polej}} - $rj;
  my $r2 = $xx**2. + $yy**2.;
  my $r = sqrt($r2)/${proj{rebydx}};

  # Convert to lat/lon
  if ($r2 eq 0.) {
    $lat = ${proj{hemi}} * 90.;
    $lon = ${proj{stdlon}};
  }else{
    $lon = ${proj{stdlon}} + rad2deg( 
           atan2(${proj{hemi}}*$xx,$yy)/${proj{cone}}); 
    $lon = fmod(($lon+360),360.);
    my $chi;
    if ($chi1 eq $chi2) { 
      $chi = 2.0 * atan(($r/tan($chi1))**(1./${proj{cone}}) *
                 tan($chi1*0.5));
    }else{
      $chi = 2.0 * atan(($r*${proj{cone}}/sin($chi1))**(1./${proj{cone}}) *
             tan($chi1*0.5)); 
    }
    $lat = (90.0 - rad2deg($chi)) * ${proj{hemi}};
  }
  if ($lon > 180.) { $lon = $lon - 360. }
  if ($lon < -180.){ $lon = $lon + 360. }
  return($lat,$lon);
} 
###############################################################################
sub llij_me(\$\$\%){

  # Subroutine that computes i/j from lat/lon for mercator maps

  my ($lat,$lon,%proj) = @_;
  my ($ri,$rj);
  my $deltalon = $lon - ${proj{lonsw}};
  if ($deltalon < -180.){$deltalon = $deltalon + 360. }
  if ($deltalon > 180.) {$deltalon = $deltalon - 360. }
  $ri = 1. + ($deltalon/rad2deg(${proj{dellon}}));
  $rj = 1. + (log(tan(0.5*deg2rad($lat+90.))))/${proj{dellon}} - ${proj{rsw}};
  return ($ri,$rj);
}
###############################################################################
sub ijll_me(\$\$\%){

  # Subroutine that computes lat/lon from i/j for mercator maps
  
  my ($ri,$rj,%proj) = @_;
  my ($lat,$lon);

  $lat = 2.0*rad2deg(atan(exp(${proj{dellon}}*(${proj{rsw}}+$rj-1.)))) -90.;
  $lon = ($ri-1.)*rad2deg(${proj{dellon}}) + ${proj{lonsw}};
  if ($lon > 180.) { $lon = $lon - 360. }
  if ($lon < -180.){ $lon = $lon + 360. }
  return ($lat,$lon);
}
##############################################################################
sub llij_ll(\$\$\%){

  # Subroutine that computes i/j from lat/lon for lat/lon grids

  my ($lat,$lon,%proj) = @_;
  my ($ri,$rj);
  my $latinc = ${proj{truelat1}};
  my $loninc = ${proj{stdlon}};

  # Compute deltalat/deltalon

  my $deltalat = $lat - ${proj{latsw}};

  # Account for possible issues around dateline
  my $lon360;
  if ($lon < 0.) {
    $lon360 = $lon + 360.;
  }else{
    $lon360 = $lon;
  }
  my $deltalon = $lon360 - ${proj{lonsw}};
  if ($deltalon < 0) { $deltalon = $deltalon + 360. }

  # Compute i/j
  $ri = $deltalon/$loninc + 1.;
  $rj = $deltalat/$latinc + 1.;

  return($ri,$rj);
}    
##############################################################################
sub ijll_ll(\$\$\%){

  # Subroutine to convert i/j to lat/lon for lat/lon grids

  my ($ri,$rj,%proj) = @_;
  my ($lat,$lon,$latinc,$loninc,$deltalat,$deltalon);

  $latinc = ${proj{truelat1}};
  $loninc = ${proj{stdlon}};

  $deltalat = ($rj-1.)*$latinc;
  $deltalon = ($ri-1.)*$loninc;
  $lat = ${proj{latsw}} + $deltalat;
  $lon = ${proj{lonsw}} + $deltalon;

  if ( (abs($lat) > 90.) or (abs($deltalon) > 360.) ) {
    # Off the earth!
    print "Warning:  $ri $rj is off this grid!\n";
    $lat = -999.;
    $lon = -999.;
  }else{
    $lon = $lon + 360.;
    $lon = fmod($lon,360.);
    if ($lon > 180.){ $lon = $lon - 360.};
  }
  return ($lat,$lon);
}
###############################################################################
sub lc_cone(\$\$) {

  # Computes the cone factor for a lambert conformal map
  my $cone;
  my ($truelat1,$truelat2) = @_;

  if ( ($truelat2 - $truelat1) > 0.01 ) {
    $cone = (log(cos(deg2rad($truelat1))) -
             log(cos(deg2rad($truelat2)))) / 
            (log(tan(deg2rad(90.-abs($truelat1))*0.5)) -
             log(tan(deg2rad(90.-abs($truelat2))*0.5)) );

  }else{
    $cone = sin(deg2rad(abs($truelat1)));
  }
  return $cone;
}
   

