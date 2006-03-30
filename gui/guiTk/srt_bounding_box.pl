
# ---------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# bounding_box.pl 
# 	Create a bounding box on a canvas that can be moved and 
#	resized allowing user to create a custom model domain.
#
# Author: Paula McCaslin   30 Jan 2003  Original Version
# --------------------------------------------------


#use warnings;
use strict 'subs';
use strict 'refs';
#use strict;
use vars qw($can $b_clear);

# ------------------------------------------
# creating_bbox
#
# ------------------------------------------
sub creating_bbox {
    use strict;

    my $mytitle="Domain Box Instructions";
    my $mymsg=
"To create a model domain use the domain bounding box to interactively select an area of interest using the following instructions. You'll also have the opportunity to fine-tune your domain after pressing 'Update Map'.\n
The following instructions ask you to press and hold down mouse button 1 while dragging the mouse to a new location.\n
o To create a domain box -- Press & hold the mouse while you drag the cursor. 
\t\t\tA domain box, with tabs, will be created.
o To resize a domain box -- Press & hold the mouse on any corner or side tag while you drag
\t\t\tthe mouse.
o To move a domain box -- Press & hold the mouse anywhere on the domain box while
\t\t\tyou drag the mouse.
\t\t\tOr, manually enter center lat/lon. This is the only place to modify these
\t\t\tparameters.\n
\t\t\tThe 'Clear' button will invoked.\n";
       info_dbox($mytitle,$mymsg);
}

# ------------------------------------------
# Set bindings to be able to draw a bounding box.
#
# 1=draw bbox
# 0=move/resize existing bbox
# ------------------------------------------
sub set_bindings_to_draw {
    my ($arg)=@ARG;

    
    # In case cursor is not refreshed (rare condition).
    $can->configure(-cursor => "crosshair");

    if ($arg eq 1) {
       # Create a bounding box, that can be resized and moved.
       $can->Tk::bind ("<Button-1>", [\&select_start, Ev('x'), Ev('y')]);
       $can->Tk::bind ("<B1-Motion>", [\&select_draw, Ev('x'), Ev('y')]);
       $can->Tk::bind ("<ButtonRelease-1>", [\&select_end, Ev('x'), Ev('y')]);

       # Disable the projection and center point widgets.
       set_snap_box_state(0);

    } elsif ($arg eq 0) {
       # Redefine bindings to display cursor output with lat,lon values.
       $can->Tk::bind ("<Button-1>", [\&can_xy_to_ll_proj, Ev('x'), Ev('y')]);
       $can->Tk::bind ("<B1-Motion>", [\&can_xy_to_ll_proj, Ev('x'), Ev('y')]);
       $can->Tk::bind ("<ButtonRelease-1>", \&delete_latlon_text);
       
       # Activate the projection and center point widgets.
       set_snap_box_state(1);
     
    } elsif ($arg eq 2) {
       # Create a nest box, that can be resized and moved.
       $can->Tk::bind ("<Button-1>", [\&select_start_nest, Ev('x'), Ev('y')]);
       $can->Tk::bind ("<B1-Motion>", [\&select_draw_nest, Ev('x'), Ev('y')]);
       $can->Tk::bind ("<ButtonRelease-1>", [\&select_end_nest, Ev('x'), Ev('y')]);

       $can->configure(-cursor => "sizing");
       $hint_msg="To DRAW A NEST, place the mouse cursor in the canvas.";
    }
}

# ----------------------------------------------------------------------
# USAGE:  select_start <x> <y>
#         select_draw <x> <y>
#         select_end <x> <y>
#
# Used to handle bounding box operations on the canvas.  "select_start"
# marks the starting point (one corner of the box).  "select_draw"
# draws the bounding box at each motion point. "select_end" removes the
# bindings, once you create a box it can further be resized and moved. 
#
#   $c->canvasx($x) coverts x space to canvas space.
# ----------------------------------------------------------------------
sub select_start {
    my ($c, $x, $y)=@ARG;
    $x0 = $c->canvasx($x);
    $y0 = $c->canvasy($y);

    # x,y values must not exceed +/-180, +/-90, respectively.
    $knownXmax=$img_xmax+($img_xmax-$x0);
    $knownYmax=$img_ymax+($img_ymax-$y0);
    $knownXmin=-$x0;
    $knownYmin=-$y0;
}

sub select_draw {
    my ($c, $x, $y)=@ARG;
    $x = $c->canvasx($x);
    $y = $c->canvasy($y);

    # x,y values must not exceed +/-180, +/-90, respectively.
    if ($x >= $knownXmax) { $x=$knownXmax; }
    if ($x <= $knownXmin) { $x=$knownXmin; }
    if ($y >= $knownYmax) { $y=$knownYmax; }
    if ($y <= $knownYmin) { $y=$knownYmin; }

    $can->delete("bbox");
    $can->createRectangle($x0, $y0, $x, $y, 
                          -outline => $colorW, 
                          -fill => $colorW, 
                          -stipple => 'transparent',
                          -tags => "bbox",
                          -width => 1);

    # Suggest projection type based on center point lat/lon.
    suggest_proj_type();
    create_tags($x0, $y0, $x, $y);
    $pix_per_gpoint=($x-$x0)/$nx_dim;
    bbox_corners_etc($x0, $y0, $x, $y);
}

sub select_end {
    # Redefine bindings to display cursor output.
    set_bindings_to_draw(0);

    if (&min_size_test()) {return;}

    # Update needs to be pressed, a parameter has changed.
    highlight_update_button(1);
}

# ------------------------------------------
# USAGE:  move_start <x> <y>
#         move_draw <x> <y>
#         draw_end <x> <y>
#
# Bounding box can be moved interactively. 
# ------------------------------------------

sub move_start {
    if (!$calc_boundingBox || !$calc_centerLatLon) { return; }

    my ($c, $x, $y)=@ARG;
    $x0 = $c->canvasx($x);
    $y0 = $c->canvasy($y);

    $can->configure(-cursor => "fleur");
    hide_tags();

    # x,y values must not exceed +/-180, +/-90, respectively.
    my ($l, $t, $r, $b) = $can->coords("bbox");
    $knownDistX=($r-$l)/2;
    $knownDistY=($b-$t)/2;
}

sub move_draw {
    if (!$calc_boundingBox || !$calc_centerLatLon) { return; }

    my ($c, $x, $y)=@ARG;
    $x = $c->canvasx($x);
    $y = $c->canvasy($y);

    # x,y values must not exceed +/-180, +/-90, respectively.
    my $stayput=0;
    my ($cnt_x, $cnt_y)=$can->coords("c_img");

    if ($cnt_x < 0) { $cnt_x=0; $stayput=1; }
    if ($cnt_y < 0) { $cnt_y=0; $stayput=1; }
    if ($cnt_x > $img_xmax) { $cnt_x=$img_xmax; $stayput=1; }
    if ($cnt_y > $img_ymax) { $cnt_y=$img_ymax; $stayput=1; }
    if ($stayput) {
        $l=$cnt_x+$knownDistX; 
        $r=$cnt_x-$knownDistX; 
        $t=$cnt_y+$knownDistY; 
        $b=$cnt_y-$knownDistY;
        $can->coords("bbox", $l, $t, $r, $b);
        $can->coords("c_img", $cnt_x, $cnt_y);
        return;
    }

    # Get bounding box coords.
    my ($l, $t, $r, $b) = $can->coords("bbox");
    my $dist_x=$x-$x0;
    my $dist_y=$y-$y0;

    # Add new distances.
    $l=$l+$dist_x; 
    $r=$r+$dist_x; 
    $t=$t+$dist_y; 
    $b=$b+$dist_y;

    # Set bounding box coords.
    $can->coords("bbox", $l, $t, $r, $b);
    centerpoint_draw($l, $t, $r, $b);
    $x0=$x;
    $y0=$y;

    # Calculate UL & LR corners (and grid spacing).
    bbox_corners_etc($l,$t,$r,$b);
}

sub draw_end {

    # Force the user to select (or re-select) the projection type 
    # only AFTER the user is completely satisfied with the domain.
    if ($globalMap == 1) {$proj_label="None";}
   
    # If bind ButtonRelease contains bbox or c_img, then return because
    # user is just dragging and ButtonReleasing the cursor while outputing
    # lat,lon values.
    if ($globalMap == 0 && ($_[1] eq 'bbox' || $_[1] eq 'c_img') ) {return;} 

    $can->configure(-cursor => "crosshair");
    update_tags();

    # If domain_mode =1 and a resize, then check min size.
    if ($domain_mode == 1 && $_[1] ne 'bbox' && $_[1] ne 'c_img' ) { 
      if (min_size_test()) {return;} 
    }

    #if ($domain_mode == 1) { center_bbox(); }
    if ($calc_boundingBox) { highlight_update_button(1); }
}

# ------------------------------------------
# USAGE:  resize_start <x> <y>
#         resize_draw <x> <y>
#         draw_end <x> <y>
#
# Rubber box can be resized interactively. 
# ------------------------------------------

sub resize_start {
    if (!$calc_boundingBox) { return; }

    my ($c, $x, $y, $flag)=@ARG;
    $x0 = $c->canvasx($x);
    $y0 = $c->canvasy($y);

    if($flag eq "l_edge") {
       $can->configure(-cursor => 'left_side');
    }elsif($flag eq "r_edge") {
       $can->configure(-cursor => 'right_side');
    }elsif($flag eq "t_edge") {
       $can->configure(-cursor => 'top_side');
    }elsif($flag eq "b_edge") {
       $can->configure(-cursor => 'bottom_side');
    }elsif($flag eq "l_top") {
       $can->configure(-cursor => 'top_left_corner');
    }elsif($flag eq "r_top") {
       $can->configure(-cursor => 'top_right_corner');
    }elsif($flag eq "l_bot") {
       $can->configure(-cursor => 'bottom_left_corner');
    }elsif($flag eq "r_bot") {
       $can->configure(-cursor => 'bottom_right_corner');
    }

    hide_tags();
}

sub resize_draw {
    if (!$calc_boundingBox) { return; }

    my ($c, $x, $y, $flag)=@ARG;
    $x = $c->canvasx($x);
    $y = $c->canvasy($y);

    # x,y values must not exceed +/-180, +/-90, respectively.
    my $stayput=0;
    my ($cnt_x, $cnt_y)=$can->coords("c_img");
    if ($cnt_x < 0) { $cnt_x=0; $stayput=1; }
    if ($cnt_y < 0) { $cnt_y=0; $stayput=1; }
    if ($cnt_x > $img_xmax) { $cnt_x=$img_xmax; $stayput=1; }
    if ($cnt_y > $img_ymax) { $cnt_y=$img_ymax; $stayput=1; }
    if ($stayput) {

        my ($l, $t, $r, $b) = $can->coords("bbox");
        $knownDistX=($r-$l)/2;
        $knownDistY=($b-$t)/2;

        $l=$cnt_x+$knownDistX; 
        $r=$cnt_x-$knownDistX; 
        $t=$cnt_y+$knownDistY; 
        $b=$cnt_y-$knownDistY;
        $can->coords("bbox", $l, $t, $r, $b);
        $can->coords("c_img", $cnt_x, $cnt_y);

        return;
    }

    my $dist_x=$x-$x0;
    my $dist_y=$y-$y0;

    if ($calc_centerLatLon) { 
      # Move one side of box and thus move the center point.
      $opp_dist_x=$opp_dist_y=0;
    } else {
      # Move opposite sides of box and thus do Not move 
      # or calculate the center point.
      $opp_dist_x=-$dist_x;
      $opp_dist_y=-$dist_y;
    }

    # Get bounding box coords.
    my ($l, $t, $r, $b) = $can->coords("bbox");

    # Set bounding box coords.
    if($flag eq "l_edge") {
       $l+=$dist_x; 
       $r+=$opp_dist_x;

    } elsif($flag eq "r_edge") {
       $l+=$opp_dist_x; 
       $r+=$dist_x;

    }elsif($flag eq "t_edge") {
       $t+=$dist_y; 
       $b+=$opp_dist_y;

    }elsif($flag eq "b_edge") {
       $t+=$opp_dist_y; 
       $b+=$dist_y;

    }elsif($flag eq "l_top") {
       $l+=$dist_x; 
       $r+=$opp_dist_x;
       $t+=$dist_y; 
       $b+=$opp_dist_y;

    }elsif($flag eq "r_top") {
       $l+=$opp_dist_x; 
       $r+=$dist_x;
       $t+=$dist_y; 
       $b+=$opp_dist_y;

    }elsif($flag eq "r_bot") {
       $l+=$opp_dist_x; 
       $r+=$dist_x;
       $t+=$opp_dist_y; 
       $b+=$dist_y;

    }elsif($flag eq "l_bot") {
       $l+=$dist_x; 
       $r+=$opp_dist_x;
       $t+=$opp_dist_y; 
       $b+=$dist_y;
    }
    $can->coords("bbox", $l, $t, $r, $b);


    $x0=$x;
    $y0=$y;

    if ($calc_centerLatLon) { 
        centerpoint_draw($l, $t, $r, $b);
    }

    bbox_corners_etc($l,$t,$r,$b);
}

# ------------------------------------------
# USAGE:  min_size_test <x> <y>
#
# The left and right size of the bbox cannot be equal or less than 4 pixels.
# The top and bottom size of the bbox cannot be equal or less than 4 pixels.
# ------------------------------------------
sub min_size_test {
    use strict;

    my ($l, $t, $r, $b) = $can->coords("bbox");
    my $mymsg;
    my $mytitle="Domain Box Info";
    if (  $l == "" || $t == "" || $r == "" || $b == "" || # Creating domain test.
          ($r-$l) <  4 || ($b-$t) <  4 ) {                # Resizing domain test.
       # Domain box is null.

       $mymsg="The domain bounding box that you tried to create is too small to continue. The following instructions ask you to press and hold down mouse button 1 while dragging the mouse to a new location.\n\no To create a domain box -- Press & hold the mouse while you drag the cursor. \n\t\t\tA domain box, with tabs, will be created.\no To resize a domain box -- Press & hold the mouse on any corner or side tag while you drag\n\t\t\tthe mouse.\no To move a domain box -- Press & hold the mouse anywhere on the domain box while\n\t\t\tyou drag the mouse.\n\n\t\t\tThe 'Clear' button will invoked.\n";

       info_dbox($mytitle,$mymsg);
       render_map_clear();
       # Error.
       return(1);

    } elsif ( ($r-$l) < 10 || ($b-$t) < 10 ) {           # Resizing domain test.
       # Domain box is too small.
       
       $mymsg="The domain box is not large enough to continue. Resize domain box or press 'Clear'.";
       info_dbox($mytitle,$mymsg);
       # Error.
       return(1);
    } 

    # Success.
    return(0);
}

# ------------------------------------------
# Continually update the center point as bounding box is 
# interactively moved.
# ------------------------------------------
sub centerpoint_draw {
    my ($l, $t, $r, $b)=@_;

    # Update Cursor.
    my $cent_x=(($r-$l) / 2) +$l;
    my $cent_y=(($b-$t) / 2) +$t;
    $can->coords("c_img", $cent_x, $cent_y);
 
    if ($globalMap) {
        cylindrical_xy2ll($cent_x, $cent_y);
    }
}
 
# -----------------------------------------------------
# Calculate lat&lon values for display entry boxes based on x&y.
#
# -----------------------------------------------------
sub cylindrical_xy2ll {
    my ($cent_x, $cent_y)=@ARG;
  
    # With x,y known calculate lat,lon.
    $grid_cen_lon_cmn=x_to_cyl_lon($cent_x);
    $grid_cen_lat_cmn=y_to_cyl_lat($cent_y);

    # In order to have show std lon, true lats with  
    # globalMap these need to be set here.
    set_true_lats();
}

# -----------------------------------------------------
# wrap the Cylindrical x&y values for user entered lat&lon values
#
# -----------------------------------------------------
sub wrap_cylindrical_ll2xy {
    my ($null,$zoom_val)=@ARG;
 
    if ($grid_cen_lon_cmn eq $grid_cen_lon_orig && 
        $grid_cen_lat_cmn eq $grid_cen_lat_orig) { return;}
    if (!($mw->focusCurrent eq $clonEnt ||
          $mw->focusCurrent eq $clatEnt)) { return; }
 
    cylindrical_ll2xy(0, $zoom_val);

    $grid_cen_lon_orig=$grid_cen_lon_cmn;
    $grid_cen_lat_orig=$grid_cen_lat_cmn;
}

# -----------------------------------------------------
# Calculate x&y values for user entered lat&lon values
#
# Then call place_bbox to move, zoom and possibly
# center the bounding box.
# -----------------------------------------------------
sub cylindrical_ll2xy {
    my ($null,$zoom_val)=@ARG;


    # Update needs to be pressed, a parameter has changed.
    highlight_update_button(1);

    # With lat,lon known calculate x,y.
    my $cent_x=($grid_cen_lon_cmn + 180.) / (360 / $img_xmax);
    #my $cent_x=($grid_cen_lon_cmn + 360.) / (720 / $img_xmax);
    my $cent_y=$img_ymax - (($grid_cen_lat_cmn + 90.) / (180 / $img_ymax));

    # Check x,y after calculations.
    if ($cent_x >= $img_xmax) { $cent_x=$img_xmax; }
    if ($cent_y >= $img_ymax) { $cent_y=$img_ymax; }
    if ($cent_x <= 0) { $cent_x=0; }
    if ($cent_y <= 0) { $cent_y=0; }

    place_bbox($cent_x,$cent_y,$zoom_val); 
  
    # If center lat changes, then update initial values of truelat1, truelat2.
    set_true_lats(); 
}

# -----------------------------------------------------
# Calculate lat&lon values for display entry boxes based on x&y.
#
# -----------------------------------------------------
sub x_cylindrical_xy2ll {
    my ($cent_x, $cent_y)=@ARG;
  
    # With x,y known calculate lat,lon.
    $stdlon  =x_to_cyl_lon($cent_x);
    $truelat1=y_to_cyl_lat($cent_y);

    # In order to have show std lon, true lats with  
    # globalMap these need to be set here.
    set_center_latlon();
}

# -----------------------------------------------------
# wrap the Cylindrical x&y values for user entered lat&lon values
#
# -----------------------------------------------------
sub x_wrap_cylindrical_ll2xy {
    my ($null,$zoom_val)=@ARG;
 
    if ($stdlon   eq $stdlon_orig && 
        $truelat1 eq $truelat1_orig) { return;}
    if (!($mw->focusCurrent eq $lonEnt ||
          $mw->focusCurrent eq $latEnt)) { return; }
 
    x_cylindrical_ll2xy(0, $zoom_val);

    $stdlon_orig=$stdlon;
    $truelat1_orig=$truelat1;
}

# -----------------------------------------------------
# Calculate x&y values for user entered lat&lon values
#
# Then call place_bbox to move, zoom and possibly
# center the bounding box.
# -----------------------------------------------------
sub x_cylindrical_ll2xy {
    my ($null,$zoom_val)=@ARG;


    # Update needs to be pressed, a parameter has changed.
    highlight_update_button(1);

    # With lat,lon known calculate x,y.
    my $cent_x=($stdlon + 180.) / (360 / $img_xmax);
    #my $cent_x=($grid_cen_lon_cmn + 360.) / (720 / $img_xmax);
    my $cent_y=$img_ymax - (($truelat1 + 90.) / (180 / $img_ymax));

    # Check x,y after calculations.
    if ($cent_x >= $img_xmax) { $cent_x=$img_xmax; }
    if ($cent_y >= $img_ymax) { $cent_y=$img_ymax; }
    if ($cent_x <= 0) { $cent_x=0; }
    if ($cent_y <= 0) { $cent_y=0; }

    place_bbox($cent_x,$cent_y,$zoom_val); 
  
    # If center lat changes, then update initial values of truelat1, truelat2.
    set_center_latlon();
}

# ------------------------------------
# place_bbox
#
# Place the boundingbox with new information.
# If zoom is set the bounding box will be scaled.
# ------------------------------------
sub place_bbox {
    #if (!$calc_boundingBox) { return; }
    my ($cnt_x,$cnt_y,$zoom_val)=@ARG;

    # Get bounding box coords.
    my ($l, $t, $r, $b) = $can->coords("bbox");
    my $bbox_width =$r-$l;
    my $bbox_height=$b-$t;
    my $half_w=(($bbox_width ) / 2) * $zoom_val;
    my $half_h=(($bbox_height) / 2) * $zoom_val;

    # Set bounding box, scale to zoom_val.
    $can->coords("bbox", $cnt_x-$half_w, $cnt_y-$half_h,
                         $cnt_x+$half_w, $cnt_y+$half_h);
    update_tags();
    center_bbox(); 

    # Calculate UL & LR corners (and grid spacing, if =1)
    # Recalculate l,t,r,b.
    ($l, $t, $r, $b) = $can->coords("bbox");
    bbox_corners_etc($l,$t,$r,$b);
}

# ------------------------------------
# center_bbox
#
# Center the boundingbox with new information.
# ------------------------------------
sub center_bbox {

    my ($cnt_x, $cnt_y)=$can->coords("c_img");
    $xMoveto=($cnt_x - ($can->width / 2)) / $img_xmax;
    $yMoveto=($cnt_y - ($can->height/ 2)) / $img_ymax;
    if ($xMoveto <= 0) {$xMoveto=0;}
    if ($yMoveto <= 0) {$yMoveto=0;}

    # Adjust the scrollbar sliders in order to view centerpoint icon.
    $can->xviewMoveto($xMoveto);
    $can->yviewMoveto($yMoveto);
}

# ----------------------------------------------------------------------
# USAGE:  create_tags
#
# Adds active tags to the corners and sides of the bounding box,
# allowing the user to move and resize the bounding box.
# ----------------------------------------------------------------------
sub create_tags {
  my ($l, $t, $r, $b)=@_;

  $can->delete("rtags");

  # Get bounding box coords.
  my $cent_x=(($r-$l) / 2) +$l;
  my $cent_y=(($b-$t) / 2) +$t;

  if ($globalMap) {

      # Update cylindrical centerpoint (get lat,lon from x,y).
      if ($perform_calc) {cylindrical_xy2ll($cent_x, $cent_y);}

  } else {

      # Update non-cylindrical centerpoint (get x,y from lat,lon).
      ($cent_x, $cent_y)=ll_proj_to_xy($grid_cen_lat_cmn, $grid_cen_lon_cmn);

      # Get bounding box coords.
      my $bbox_width =$r-$l;
      my $bbox_height=$b-$t;
      my $half_w=($bbox_width / 2);
      my $half_h=($bbox_height/ 2);

      # Set bounding box, scale to zoom_val.
      $l=$cent_x-$half_w; 
      $r=$cent_x+$half_w;
      $t=$cent_y-$half_h;
      $b=$cent_y+$half_h;
      $can->delete("bbox");
      $can->createRectangle($l, $t, $r, $b,
                          -outline => $colorW, 
                          -fill => $colorW, 
                          -stipple => 'transparent', # Tk800.020 (vs $transparent)
                          -tags => "bbox",
                          -width => 1);
  }

  # Swap vars, if user drags mouse backwards and up (west and north).
 if(1) {
  if ($r < $l) {
     my $tmp=$l;
     $l=$r;
     $r=$tmp;
  }
  if ($b < $t) {
     my $tmp=$t;
     $t=$b;
     $b=$tmp;
  }
 }

  my $sz=6;
  $can->createImage($cent_x, $cent_y,
                                  -image => $centerPoint,
                                  -tags => ["c_img", "rtags"]);

  $can->createRectangle( ($l+0), ($t+0), ($l+$sz), ($t+$sz), 
                                  -outline => $colorW,
                                  -fill => $colorW,
                                  -tags => ["l_top", "rtags"]);
  $can->createRectangle( ($r-$sz), ($t+$sz), ($r+0), ($t+0), 
                                  -outline => $colorW,
                                  -fill => $colorW,
                                  -tags => ["r_top", "rtags"]);
  $can->createRectangle( ($l+0), ($b-$sz), ($l+$sz), ($b+0), 
                                  -outline => $colorW,
                                  -fill => $colorW,
                                  -tags => ["l_bot", "rtags"]);
  $can->createRectangle( ($r-$sz), ($b-$sz), ($r+0), ($b+0), 
                                  -outline => $colorW,
                                  -fill => $colorW,
                                  -tags => ["r_bot", "rtags"]);
  $cent_x=$cent_x-3;
  $cent_y=$cent_y-3;

  $can->createRectangle( $cent_x, ($t+0), ($cent_x+$sz), ($t+$sz), 
                                  -outline => $colorW,
                                  -fill => $colorW,
                                  -tags => ["t_edge", "rtags"]);
  $can->createRectangle( $cent_x, ($b+0), ($cent_x+$sz), ($b-$sz), 
                                  -outline => $colorW,
                                  -fill => $colorW,
                                  -tags => ["b_edge", "rtags"]);
  $can->createRectangle( ($l-0), $cent_y, ($l+$sz), ($cent_y+$sz), 
                                  -outline => $colorW,
                                  -fill => $colorW,
                                  -tags => ["l_edge", "rtags"]);
  $can->createRectangle( ($r-0), $cent_y, ($r-$sz), ($cent_y+$sz), 
                                  -outline => $colorW,
                                  -fill => $colorW,
                                  -tags => ["r_edge", "rtags"]);
  
  @bbox_widget_list=("l_top", "r_top", "l_bot", "r_bot",
                     "t_edge", "b_edge", "l_edge", "r_edge");

#   The Tk::Widget inherits from Tk::Canvas. 
#   There are differences with ->bind vs ->Tk::bind; bind is either 
#   available as a Tk method and as a Tk::Canvas method. If you want 
#   to apply the Tk method, you have to use $canvas->Tk::bind, otherwise 
#   use just $canvas->bind. The same applies for other methods (e.g. the 
#   lower() and raise() methods).  - Slaven Rezic 

  # Bind transparent bbox and center icon to move the bounding box.
  foreach my $ttag ("bbox", "c_img") {
     $can->bind ($ttag, "<Button-1>", [\&move_start, Ev('x'), Ev('y')]);
     $can->bind ($ttag, "<B1-Motion>", [\&move_draw, Ev('x'), Ev('y')]);
     $can->bind ($ttag, "<ButtonRelease-1>", [\&draw_end, $ttag]);
  }

  # Bind corner and edge tags to resize the bounding box.
  foreach $ttag (@bbox_widget_list) {
     $can->bind ($ttag, "<Button-1>", [\&resize_start, Ev('x'), Ev('y'), $ttag ]);
     $can->bind ($ttag, "<B1-Motion>", [\&resize_draw, Ev('x'), Ev('y'), $ttag]);
     $can->bind ($ttag, "<ButtonRelease-1>", [\&draw_end, $ttag]);
  }

}

# ----------------------------------------------------------------------
# USAGE:  hide_tags
#
# Temporarily hides tags on the bounding box.
# but not the centerpoint icon $can->coords("c_img", 0, 0);
# ----------------------------------------------------------------------
sub hide_tags {
    foreach my $ttag (@bbox_widget_list) {
         $can->coords($ttag, -1, -1, -1, -1);
    }
}

# ----------------------------------------------------------------------
# USAGE:  update_tags
#
# Updates the locations of the active tags on the bounding box.
# ----------------------------------------------------------------------
sub update_tags {

  if (!$calc_boundingBox) { return; }

  # Suggest projection type based on new center point lat/lon.
  suggest_proj_type();

  # Get bounding box coords.
  my ($l, $t, $r, $b) = $can->coords("bbox");
  #my $cent_x=(($r-$l) / 2) +$l;
  #my $cent_y=(($b-$t) / 2) +$t;
  my $cent_x=(($r+$l) / 2);
  my $cent_y=(($b+$t) / 2);
  my $sz=6;

  $can->coords("c_img",  $cent_x, $cent_y);
  $can->coords("l_top",  $l,       $t,      ($l+$sz), ($t+$sz));
  $can->coords("r_top", ($r-$sz),  $t,       $r,      ($t+$sz));
  $can->coords("r_bot", ($r-$sz), ($b-$sz),  $r,       $b);
  $can->coords("l_bot",  $l,      ($b-$sz), ($l+$sz),  $b);

  $cent_x=$cent_x-3;
  $cent_y=$cent_y-3;

  $can->coords("t_edge", $cent_x, $t,      ($cent_x+$sz), ($t+$sz));
  $can->coords("l_edge", $l,      $cent_y, ($l+$sz),      ($cent_y+$sz));
  $can->coords("r_edge", $r,      $cent_y, ($r-$sz),      ($cent_y+$sz));
  $can->coords("b_edge", $cent_x, $b,      ($cent_x+$sz), ($b-$sz));
}

# -----------------------------------------
# bbox_corners_etc
#
# As the user interactively defines the bounding box, 
# calculate the values for
#
#  -corner lat,lon values calculated from x,y
#  -grid spacing
#  -horizontal NX NY grid values.
# -----------------------------------------

sub bbox_corners_etc {

    my ($l,$t,$r,$b)=@_;


    # Calculate bbox width and height.
    my $bbox_height=$b-$t;
    my $bbox_width =$r-$l;

    # With x,y known calculate lat,lon.
    if ($globalMap) {
        
       # Cylindircal.
       $lonSW=x_to_cyl_lon($l);
       $lonNE=x_to_cyl_lon($r);
       $latNE=y_to_cyl_lat($t);
       $latSW=y_to_cyl_lat($b);

    } else {
      
       # Non-cylindircal.
       ($latSW,$lonSW)=xy_to_ll_proj($l, $b);
       ($latNE,$lonNE)=xy_to_ll_proj($r, $t);
    }

    # Use the LAPS tools great_circle_distance to get domain distances.
    # Use domain distances to calc grid spacing.
    if ($calc_gridSpacing eq 1) {

       # Convert degrees to radians.
       my $rad_l = $deg2rad * $lonSW;
       my $rad_b = $deg2rad * $latSW;
       my $rad_r = $deg2rad * $lonNE;
       my $rad_t = $deg2rad * $latNE;
       my $rad_x = $deg2rad * $grid_cen_lon_cmn;
       my $rad_y = $deg2rad * $grid_cen_lat_cmn;

       # Treat PS carefully, the great_circle_distance 
       # is invalid at polar caps, +/-90. 
       if ($polar_cap) {
          $rad_x = $rad_r;
          if ($grid_cen_lat_cmn >= 0) {
             $rad_y = $rad_b;
          } else {
             $rad_y = $rad_t;
          }
       }

       $xdist=(&laps_tools::great_circle_distance
              ($rad_y,$rad_l,$rad_y,$rad_r)/1000);
       $ydist=(&laps_tools::great_circle_distance
              ($rad_t,$rad_x,$rad_b,$rad_x)/1000);

       if ($domain_mode > 1) { 
         $grid_spacing_dkm=sprintf ("%.1f", ($xdist/ ($nx_dim -1) ) );
       }
       $grid_spacing_km =sprintf ("%.1f", ($xdist/ ($nx_dim -1) ) );
       if ($grid_spacing_km > 0) { $ny_dim=int($ydist/$grid_spacing_km)+1; }


    } elsif ($grid_val_restrict == 1) {
       # If $grid_val_restrict eq edit bbox, then updating nx_dim & ny_dim is a must.

       $nx_dim=int($bbox_width /$pix_per_gpoint + 0.5);
       $ny_dim=int($bbox_height/$pix_per_gpoint + 0.5);
       # Set nx_ddim too.
       set_dnx();

       # Highlight background if values change.
       if ($nx_ddim != $nx_orig) { $nx_ent->configure(-bg => $colorY2); }
       if ($ny_dim  != $ny_orig) { $ny_ent->configure(-bg => $colorY2); }
    } 

    # Store current grid spacing var for use again.
    $dx_orig=$grid_spacing_dkm;
    $grid_pts=$nx_ddim*$ny_dim;

    display_ll_ur($lonSW,$lonNE,$latNE,$latSW);   
}

# _______________ NESTING ROUTINES _____________________________
#

# ----------------------------------------------------------------------
# USAGE:  select_start_nest <x> <y>
#         select_draw_nest <x> <y>
#         select_end_nest <x> <y>
#
# Used to handle bounding box operations on the canvas.  "select_start"
# marks the starting point (one corner of the box).  "select_draw"
# draws the bounding box at each motion point. "select_end" removes the
# bindings, once you create a box it can further be resized and moved. 
#
#   $c->canvasx($x) coverts x space to canvas space.
# ----------------------------------------------------------------------
sub select_start_nest {
    my ($c, $x, $y)=@ARG;
    $x0 = $c->canvasx($x);
    $y0 = $c->canvasy($y);
}

sub select_draw_nest {
    my ($c, $x, $y)=@ARG;
    $x = $c->canvasx($x);
    $y = $c->canvasy($y);

    $can->delete("nestbox"); 
    $can->createRectangle($x0, $y0, $x, $y, 
                          -outline => $colorY2, 
                          -tags => "nestbox",
                          -width => 1);

    bbox_corners_nest($x0, $y0, $x, $y, 0);
}

sub select_end_nest {
    # Redefine bindings to display cursor output.
    set_bindings_to_draw(0);

    my ($l, $t, $r, $b) = $can->coords("nestbox");
    if ( bbox_corners_nest($l, $t, $r, $b, 1) ) {
	delete_nest();
        $hint_msg=
"*** New nest CANNOT BE completely outside of it's parent domain! ***";
#"$hint_msg\t*** New nest $nest_id_label CANNOT BE completely outside of it's parent domain! ***";
	return;
    }
    ($l, $t, $r, $b) = $can->coords("nestbox");
    create_tags_nest($l, $t, $r, $b);

    fill_nest_table_entry();

    if (min_size_test()) {return;}
    $can->configure(-cursor => "crosshair");

    # Display all nests.
    if ($num_nests > 1) { draw_all_nests(); }

}

sub render_nest_label {

    my ($l, $t, $r, $b)=@_;
    $can->createText($r, $t,
                          -fill => $colorY2, 
                          -text => $nest_id_label,
                          -anchor => 'ne',
                          -font => 'helvetica 11 bold', 
                          -tags => "nest_lab");
}
sub remove_nest_label {
    if ($can->gettags('nest_lab')) { $can->delete("nest_lab"); }
}

# ------------------------------------------
# USAGE:  move_start_nest <x> <y>
#         move_draw_nest <x> <y>
#         draw_end_nest <x> <y>
#
# Bounding box can be moved interactively. 
# ------------------------------------------

sub move_start_nest {
    my ($c, $x, $y)=@ARG;
    $x0 = $c->canvasx($x);
    $y0 = $c->canvasy($y);

    hide_tags_nest();
    $can->configure(-cursor => "fleur");
}

sub move_draw_nest {
    my ($c, $x, $y)=@ARG;
    $x = $c->canvasx($x);
    $y = $c->canvasy($y);

    # Get bounding box coords.
    my ($l, $t, $r, $b) = $can->coords("nestbox");
    my $dist_x=$x-$x0;
    my $dist_y=$y-$y0;

    # Add new distances.
    $l=$l+$dist_x; 
    $r=$r+$dist_x; 
    $t=$t+$dist_y; 
    $b=$b+$dist_y;


    # If nest border values exceed the bbox border values, return.
    my ($stop)=stop_drawing_nest($l,$t,$r,$b);
    if ($stop) { return;}

    # Set nest box coords.
    $can->coords("nestbox", $l, $t, $r, $b);

    $x0=$x;
    $y0=$y;
}

sub draw_end_nest {
    $can->configure(-cursor => "crosshair");

    # Calculate UL & LR corners (and set nestbox grid).
    my ($l, $t, $r, $b) = $can->coords("nestbox");
    bbox_corners_nest($l,$t,$r,$b, 1);
    update_tags_nest();
    render_nest_label($l, $t, $r, $b);
  
    # Display all nests.
    my $HOLD=$nest_index;
    if ($num_nests > 1) { draw_all_nests(); }
    set_nest_index($HOLD+1);
    #assign_nl_vars_nest(); # Need to preserve these new values.
}

# ------------------------------------------
# stop_drawing_nest 
#
# The nestbox border values must not exceed the bbox border values. 
# ------------------------------------------
sub stop_drawing_nest {
  my ($l, $t, $r, $b)=@_;

    # Get the nestbox corners in lat,lon (from x,y).
    ($latSW,$lonSW)=xy_to_ll_proj($l, $b);
    ($latNE,$lonNE)=xy_to_ll_proj($r, $t);

    # Get the nestbox corners in i,j (from lat,lon).
    #ptm my ($lli,$llj) = &map_utils::latlon_to_ij($latSW,$lonSW,%parent_hash);
    #ptm my ($uri,$urj) = &map_utils::latlon_to_ij($latNE,$lonNE,%parent_hash);
    my $parent_i=$nl_var{PARENT_ID}[$nest_index];
    my ($lli,$llj) = &map_utils::latlon_to_ij($latSW,$lonSW,%{ $ArrayOfHashes[$parent_i]} );
    my ($uri,$urj) = &map_utils::latlon_to_ij($latNE,$lonNE,%{ $ArrayOfHashes[$parent_i]} );

    # Get nearest integer.
    $lli=int($lli+.5);
    $llj=int($llj+.5);
    $uri=int($uri+.5);
    $urj=int($urj+.5);

    #test_nest_borders();

    my $error=0;
    # grid_buffer=4.
    if ($lli < 1+$grid_buffer || 
        $llj < 1+$grid_buffer || 
        $uri > $nx_nest-$grid_buffer || 
        $urj > $ny_nest-$grid_buffer){ 
       # The nest border values exceed the bbox border values.
       $can->configure(-cursor => "X_cursor");
       
       # Error.
       $error=1;
    } elsif ($can->cget(-cursor) eq "X_cursor") {
       $can->configure(-cursor => "fleur");
    } else {
       # Keep drawing.
    }

    return($error);
}

# -----------------------------------
# test_nest_borders
#
# The nest border values exceed the bbox border values
# with grid_buffer=4.
# -----------------------------------
sub test_nest_borders {

    my $exceeds=0;
    if ($nest_lli < 1+$grid_buffer || 
        $nest_llj < 1+$grid_buffer ||
        $nest_uri > $nx_nest-$grid_buffer || 
        $nest_urj > $ny_nest-$grid_buffer){ 
        $exceeds=1;
        if ($nest_lli < 1+$grid_buffer){ $nest_lli = 1+$grid_buffer;}
        if ($nest_llj < 1+$grid_buffer){ $nest_llj = 1+$grid_buffer;}
        if ($nest_uri > $nx_nest-$grid_buffer){ $nest_uri = $nx_nest-$grid_buffer;}
        if ($nest_urj > $ny_nest-$grid_buffer){ $nest_urj = $ny_nest-$grid_buffer;}
    }
    return($exceeds);
}
 
# ------------------------------------------
# USAGE:  resize_start_nest <x> <y>
#         resize_draw_nest <x> <y>
#         draw_end_nest <x> <y>
#
# Rubber box can be resized interactively. 
# ------------------------------------------
sub resize_start_nest {

    my ($c, $x, $y, $flag)=@ARG;
    $x0 = $c->canvasx($x);
    $y0 = $c->canvasy($y);

    if($flag eq "l_edge_nest") {
       $can->configure(-cursor => 'left_side');
    }elsif($flag eq "r_edge_nest") {
       $can->configure(-cursor => 'right_side');
    }elsif($flag eq "t_edge_nest") {
       $can->configure(-cursor => 'top_side');
    }elsif($flag eq "b_edge_nest") {
       $can->configure(-cursor => 'bottom_side');
    }elsif($flag eq "l_top_nest") {
       $can->configure(-cursor => 'top_left_corner');
    }elsif($flag eq "r_top_nest") {
       $can->configure(-cursor => 'top_right_corner');
    }elsif($flag eq "l_bot_nest") {
       $can->configure(-cursor => 'bottom_left_corner');
    }elsif($flag eq "r_bot_nest") {
       $can->configure(-cursor => 'bottom_right_corner');
    }

    hide_tags_nest();
}

sub resize_draw_nest {

    my ($c, $x, $y, $flag)=@ARG;
    $x = $c->canvasx($x);
    $y = $c->canvasy($y);

    # Get cursor distance travled.
    my $dist_x=$x-$x0;
    my $dist_y=$y-$y0;

    # Get bounding box coords.
    my ($l, $t, $r, $b) = $can->coords("nestbox");

    # Set bounding box coords.
    if($flag eq "l_edge_nest") {
       $l+=$dist_x;

    } elsif($flag eq "r_edge_nest") {
       $r+=$dist_x;

    }elsif($flag eq "t_edge_nest") {
       $t+=$dist_y; 

    }elsif($flag eq "b_edge_nest") {
       $b+=$dist_y;

    }elsif($flag eq "l_top_nest") {
       $t+=$dist_y; 
       $l+=$dist_x;

    }elsif($flag eq "r_top_nest") {
       $t+=$dist_y; 
       $r+=$dist_x;

    }elsif($flag eq "r_bot_nest") {
       $b+=$dist_y;
       $r+=$dist_x;

    }elsif($flag eq "l_bot_nest") {
       $b+=$dist_y;
       $l+=$dist_x;
    }

    # If nest border values exceed the bbox border values, return.
    my ($stop)=stop_drawing_nest($l,$t,$r,$b);
    if ($stop) { return;}
   
    # Draw nest box.
    $can->coords("nestbox", $l, $t, $r, $b);

    $x0=$x;
    $y0=$y;

    # Comment out, its too much processing.
    # Get corner and center information.
    # bbox_corners_nest($l,$t,$r,$b, 0);

}

# ---------------------------------------
# remove_nestbox 
# ---------------------------------------
sub remove_nestbox {
   if ($can->gettags('nestbox')) { $can->delete("nestbox"); }
   if ($can->gettags('rtags_nest')) { $can->delete("rtags_nest"); }
}

# ----------------------------------------------------------------------
# USAGE:  create_tags_nest
#
# Adds active tags to the corners and sides of the bounding box,
# allowing the user to move and resize the bounding box.
# ----------------------------------------------------------------------
sub create_tags_nest {
  my ($l, $t, $r, $b)=@_;

  # Check for valid nest, else delete.
  if ($l==$r && $t==$b) { delete_nest(); return }

  render_nest_label($l, $t, $r, $b); 

  # Get bounding box coords.
  my $cent_x=(($r-$l) / 2) +$l;
  my $cent_y=(($b-$t) / 2) +$t;

  # Get bounding box coords.
  my $bbox_width =$r-$l;
  my $bbox_height=$b-$t;
  my $half_w=($bbox_width / 2);
  my $half_h=($bbox_height/ 2);

  # Set bounding box, scale to zoom_val.
  $l=$cent_x-$half_w; 
  $r=$cent_x+$half_w; 
  $t=$cent_y-$half_h;
  $b=$cent_y+$half_h;

  # Swap vars, if user drags mouse backwards and up (west and north).
  if($projection_type eq 'LC' or $projection_type eq 'RL') {
     if ($r < $l) {
        my $tmp=$l;
        $l=$r;
        $r=$tmp;
     }
     if ($b < $t) {
        my $tmp=$t;
        $t=$b;
        $b=$tmp;
     }
  }

  # Delete old tags.
  remove_nestbox();
  $can->createRectangle($l, $t, $r, $b,
                        -outline => $colorY2, 
                        -fill => $colorY2, 
                        -stipple => 'transparent',
                        -tags => "nestbox",
                        -width => 1);

  my $sz=6;
  $can->createImage($cent_x, $cent_y,
                                  -image => $centerPoint_dotY,
                                  -tags => ["c_img_nest", "rtags_nest"]);
  $can->createRectangle( ($l+0), ($t+0), ($l+$sz), ($t+$sz), 
                                  -outline => $colorY2,
                                  -fill => $colorY2,
                                  -tags => ["l_top_nest", "rtags_nest"]);
  $can->createRectangle( ($r-$sz), ($t+$sz), ($r+0), ($t+0), 
                                  -outline => $colorY2,
                                  -fill => $colorY2,
                                  -tags => ["r_top_nest", "rtags_nest"]);
  $can->createRectangle( ($l+0), ($b-$sz), ($l+$sz), ($b+0), 
                                  -outline => $colorY2,
                                  -fill => $colorY2,
                                  -tags => ["l_bot_nest", "rtags_nest"]);
  $can->createRectangle( ($r-$sz), ($b-$sz), ($r+0), ($b+0), 
                                  -outline => $colorY2,
                                  -fill => $colorY2,
                                  -tags => ["r_bot_nest", "rtags_nest"]);
  $cent_x=$cent_x-3;
  $cent_y=$cent_y-3;

  $can->createRectangle( $cent_x, ($t+0), ($cent_x+$sz), ($t+$sz), 
                                  -outline => $colorY2,
                                  -fill => $colorY2,
                                  -tags => ["t_edge_nest", "rtags_nest"]);
  $can->createRectangle( $cent_x, ($b+0), ($cent_x+$sz), ($b-$sz), 
                                  -outline => $colorY2,
                                  -fill => $colorY2,
                                  -tags => ["b_edge_nest", "rtags_nest"]);
  $can->createRectangle( ($l-0), $cent_y, ($l+$sz), ($cent_y+$sz), 
                                  -outline => $colorY2,
                                  -fill => $colorY2,
                                  -tags => ["l_edge_nest", "rtags_nest"]);
  $can->createRectangle( ($r-0), $cent_y, ($r-$sz), ($cent_y+$sz), 
                                  -outline => $colorY2,
                                  -fill => $colorY2,
                                  -tags => ["r_edge_nest", "rtags_nest"]);

  @nestbox_widget_list=("l_top_nest", "r_top_nest", "l_bot_nest", "r_bot_nest",
                      "t_edge_nest", "b_edge_nest", "l_edge_nest", "r_edge_nest");

  # Bind transparent bbox and center icon to move the bounding box.
  foreach my $ttag ("nestbox", "c_img_nest") {
     $can->bind ($ttag, "<Button-1>", [\&move_start_nest, Ev('x'), Ev('y')]);
     $can->bind ($ttag, "<B1-Motion>", [\&move_draw_nest, Ev('x'), Ev('y')]);
     $can->bind ($ttag, "<ButtonRelease-1>", [\&draw_end_nest, $ttag]);
  }

  # Bind corner and edge tags to resize the bounding box.
  foreach $ttag (@nestbox_widget_list) {
     $can->bind ($ttag, "<Button-1>", [\&resize_start_nest, Ev('x'), Ev('y'), $ttag ]);
     $can->bind ($ttag, "<B1-Motion>", [\&resize_draw_nest, Ev('x'), Ev('y'), $ttag]);
     $can->bind ($ttag, "<ButtonRelease-1>", [\&draw_end_nest, $ttag]);
  }

}

# ----------------------------------------------------------------------
# USAGE:  hide_tags_nest
#
# Temporarily hides tags on the nest box.
# ----------------------------------------------------------------------
sub hide_tags_nest {
    foreach my $ttag (@nestbox_widget_list) {
         $can->coords($ttag, -1, -1, -1, -1);
    }
    $can->coords("c_img_nest", -1, -1);
    remove_nest_label();
}

# ----------------------------------------------------------------------
# USAGE:  update_tags_nest
#
# Updates the locations of the active tags on the bounding box.
# ----------------------------------------------------------------------
sub update_tags_nest {

  # Get bounding box coords.
  my ($l, $t, $r, $b) = $can->coords("nestbox");
  my $cent_x=(($r+$l) / 2);
  my $cent_y=(($b+$t) / 2);
  my $sz=6;

  $can->coords("c_img_nest", $cent_x, $cent_y);

  $can->coords("l_top_nest",  $l,       $t,      ($l+$sz), ($t+$sz));
  $can->coords("r_top_nest", ($r-$sz),  $t,       $r,      ($t+$sz));
  $can->coords("r_bot_nest", ($r-$sz), ($b-$sz),  $r,       $b);
  $can->coords("l_bot_nest",  $l,      ($b-$sz), ($l+$sz),  $b);

  $cent_x=$cent_x-3;
  $cent_y=$cent_y-3;

  $can->coords("t_edge_nest", $cent_x, $t,      ($cent_x+$sz), ($t+$sz));
  $can->coords("l_edge_nest", $l,      $cent_y, ($l+$sz),      ($cent_y+$sz));
  $can->coords("r_edge_nest", $r,      $cent_y, ($r-$sz),      ($cent_y+$sz));
  $can->coords("b_edge_nest", $cent_x, $b,      ($cent_x+$sz), ($b-$sz));
}

# -----------------------------------------
# bbox_corners_nest
#
# As the user interactively defines the bounding box, 
# calculate the values for
#  -corner lat,lon values calculated from x,y
# -----------------------------------------
sub bbox_corners_nest {

    my ($l,$t,$r,$b,$switch)=@_;

    #if($projection_type eq 'LC') {
    if(0) {
       if ($r > $l) {
          my $tmp=$l;
          $l=$r;
          $r=$tmp;
       }
       if ($b > $t) {
          my $tmp=$t;
          $t=$b;
          $b=$tmp;
       }
    }
   
    # Get the nestbox corners in lat,lon (from x,y).
    ($latSW,$lonSW)=xy_to_ll_proj($l, $b);
    ($latNE,$lonNE)=xy_to_ll_proj($r, $t);

    # Get the nestbox corners in i,j (from lat,lon).
    my ($lli,$llj) = &map_utils::latlon_to_ij($latSW,$lonSW,%parent_hash);
    my ($uri,$urj) = &map_utils::latlon_to_ij($latNE,$lonNE,%parent_hash);
   
    # Get nearest integer.
    $lli=int($lli+.5);
    $llj=int($llj+.5);
    $uri=int($uri+.5);
    $urj=int($urj+.5);

    # grid_buffer=4.
    if ($switch) {
        if ( ($lli < 1+$grid_buffer && $uri < 1+$grid_buffer) ||
             ($llj < 1+$grid_buffer && $urj < 1+$grid_buffer) ||
             ($lli > $nx_nest-$grid_buffer && $uri > $nx_nest-$grid_buffer) ||
             ($llj > $ny_nest-$grid_buffer && $urj > $ny_nest-$grid_buffer) ) {
		return(1);
        }

        if ($lli < 1+$grid_buffer){$lli=1+$grid_buffer; $switch=2;} 
        if ($llj < 1+$grid_buffer){$llj=1+$grid_buffer; $switch=2;}
        if ($uri > $nx_nest-$grid_buffer){$uri=$nx_nest-$grid_buffer; $switch=2;}
        if ($urj > $ny_nest-$grid_buffer){$urj=$ny_nest-$grid_buffer; $switch=2;}
    }

    # Only when releasing B1, set the nestbox corners exactly to the MOAD grid.
    if ($switch) {

       # Get the nestbox corners in lat,lon (from i,j).
       ($latSW,$lonSW) = &map_utils::ij_to_latlon($lli,$llj,%parent_hash);
       ($latNE,$lonNE) = &map_utils::ij_to_latlon($uri,$urj,%parent_hash);
   
       # Get the nestbox corners in x,y (from new lat,lon).
       ($r, $t)=ll_proj_to_xy($latNE, $lonNE);
       ($l, $b)=ll_proj_to_xy($latSW, $lonSW);
   
       # Account for offset that centered the bbox in the canvas. 
       $l= $l+$b_xoff;
       $r= $r+$b_xoff;
       $t= $t+$b_yoff;
       $b= $b+$b_yoff;
   
       # Set the nestbox corners exactly to the MOAD grid.
       $can->coords("nestbox", $l, $b, $r, $t);
    }

    if ($switch == 2) {
       $can->configure(-cursor => "X_cursor");
    } elsif ($can->cget(-cursor) eq "X_cursor") {
       $can->configure(-cursor => "crosshair");
    }

    $grid_pts=${parent_hash{nx}} * ${parent_hash{ny}};
    display_ll_ur($lonSW,$lonNE,$latNE,$latSW);   

    # Calculate the centerPoint lat,lon and i,j.
    my $cent_x=(($r-$l) / 2) +$l;
    my $cent_y=(($b-$t) / 2) +$t;
    my ($c_lat,$c_lon)=xy_to_ll_proj($cent_x, $cent_y);
    my ($ii,$jj) = &map_utils::latlon_to_ij($c_lat,$c_lon,%parent_hash);
    $ii=int($ii+.5);
    $jj=int($jj+.5);
    $c_lat=sprintf ("%.2f",$c_lat);
    $c_lon=sprintf ("%.2f",$c_lon);

    $nest_l[$nest_index]=$l;
    $nest_t[$nest_index]=$t;
    $nest_r[$nest_index]=$r;
    $nest_b[$nest_index]=$b;

    $nest_lli=$lli;
    $nest_llj=$llj;
    $nest_uri=$uri;
    $nest_urj=$urj;
    sync_nl_vars_nest();

    $hint_msg="$hint_msg\t Centerpoint for $nest_id_label: $c_lat, $c_lon";

    # Success.
    return(0);
}
 
# -----------------------------------------
# delete_latlon_text_nest
#
# Clear cursor values.
# -----------------------------------------

sub delete_latlon_text_nest { 
     $can->delete("c_line_nest"); 
}

# -----------------------------------------
# render_bbox_label
#
# -----------------------------------------
sub render_bbox_label {

    ($l, $t, $r, $b) = $can->coords("bbox");
    $can->createText($r, $t,
                          -fill => $colorW, 
                          -text => "d01",
                          -anchor => 'ne',
                          -font => 'helvetica 11 bold', 
                          -tags => "bbox_lab");

}

sub remove_bbox_label {
   if ($can->gettags('bbox_lab')) { $can->delete("bbox_lab"); }
}

# -----------------------------------------
# render_parent_bbox
#
# Recreate a bounding box for the parent domain.
# -----------------------------------------
sub render_parent_bbox {

    my $i=$parent_of_nest-1;
    if($i == 0){
       
       ($l, $t, $r, $b) = $can->coords("bbox");

    } else {
       if (!defined $nest_l[$i] || $nest_l[$i]=="") { return;}

       $l=$nest_l[$i];
       $t=$nest_t[$i];
       $r=$nest_r[$i];
       $b=$nest_b[$i];
    }

    remove_parentbox();
    $can->createRectangle($l, $t, $r, $b,
                          -outline => $colorY2, 
                          -dash => '---',
                          -tags => "parentbox",
                          -width => 2);
    $can->createText($r, $t,
                          -fill => $colorY2, 
                          -text => $parent_of_nest_label[$nest_id_num],
                          -anchor => 'ne',
                          -font => 'helvetica 11 bold', 
                          -tags => "parentbox");

}

sub remove_parentbox {
   if ($can->gettags('parentbox')) { $can->delete("parentbox"); }
}

# -----------------------------------------
# convert_nest_vars_to_nestbox 
# 
#
# -----------------------------------------
sub convert_nest_vars_to_nestbox {
    my ($switch)=@_;
    my ($ll,$tt,$rr,$bb);

  #  if ($nest_l[$nest_index]=="") { 
       if (!defined %parent_hash) { return; } # No MOAD.
       # Get the nestbox corners in lat,lon (from i,j).
       ($latSW,$lonSW)=&map_utils::ij_to_latlon($nest_lli,$nest_llj,%parent_hash);
       ($latNE,$lonNE)=&map_utils::ij_to_latlon($nest_uri,$nest_urj,%parent_hash);
      
       # Get the nestbox corners in x,y (from new lat,lon).
       ($r, $t)=ll_proj_to_xy($latNE, $lonNE);
       ($l, $b)=ll_proj_to_xy($latSW, $lonSW);
      
       # Account for offset that centered the bbox in the canvas. 
       $ll=$nest_l[$nest_index]= $l+$b_xoff;
       $rr=$nest_r[$nest_index]= $r+$b_xoff;
       $tt=$nest_t[$nest_index]= $t+$b_yoff;
       $bb=$nest_b[$nest_index]= $b+$b_yoff;
  #  }

    if ($switch == 1) {
       $can->coords($nestBox[$nest_index], $ll, $tt, $rr, $bb);
       $can->coords($nestLab[$nest_index], $tt, $rr);
    } else {
       render_old_nestbbox();
    }

}

# -----------------------------------------
# render_old_nestbox
#
# Recreate nestbox.
# -----------------------------------------
sub render_old_nestbbox {

    if (!defined $nest_l[$nest_index]) { return;}
    if ($nest_l[$nest_index]=="") { return;}

    $l=$nest_l[$nest_index];
    $t=$nest_t[$nest_index];
    $r=$nest_r[$nest_index];
    $b=$nest_b[$nest_index];

    $can->createRectangle($l, $t, $r, $b,
                          -outline => $colorW, 
                          -tags => "oldnestbox",
                          -width => 1);

    $nest_id_label_old="d0$nest_id_num"; 
    $can->createText($r, $t,
                          -fill => $colorW, 
                          -text => $nest_id_label_old,
                          -anchor => 'ne',
                          -font => 'helvetica 11 bold', 
                          -tags => "oldnestbox");
}

sub remove_oldnestbox {
   if ($can->gettags('oldnestbox')) { $can->delete("oldnestbox"); }
   remove_nest_label();
}

sub remove_nest_graphics {
      remove_nestbox();
      remove_parentbox();
      remove_oldnestbox();
      remove_bbox_label();
      remove_nest_label();
      $mw->update;
      $mw->idletasks;
} 


# -----------------------------------------
# adjust_nestbox 
#
# Recalculate the size of the bounding box if
# the grid spacing, or the n,y dims are adjusted.
# then display new bbox. 
# -----------------------------------------
 
sub adjust_nestbox {
   
   # If 'Erase' was pressed, ignore sub adjust_nestbox,
   # as there is no nestbox to adjust.
   if ($nest_lli eq "") { return; } 


   # Test to see if nest is a parent. 
   # If so, do not proceed without confirmation.
   for ($ix=2; $ix <= $num_nests; $ix++) {
      if ($nest_id_num == $nl_var{PARENT_ID}[$ix]) {
         # This nest domain is a pArEnT.

         if(ask_delete_dependants() eq 'Cancel'){ 
            # Restore previous values and return.
            assign_nl_vars_nest();
            return;
         }
      }
   }

   # Work with current nest. 
   map_utils_setup_nest($nl_var{PARENT_ID}[$nest_index]);

   # Adjust nest.
   $nest_l[$nest_index]="";
   $nest_r[$nest_index]="";
   $nest_t[$nest_index]="";
   $nest_b[$nest_index]="";
   convert_nest_vars_to_nestbox(1);

   if (test_nest_borders()) { 
         print "PROBLEM: NEST domain exdeeds boundry Parent domain.\n";
         $nest_l[$nest_index]="";
         $nest_r[$nest_index]="";
         $nest_t[$nest_index]="";
         $nest_b[$nest_index]="";
         convert_nest_vars_to_nestbox(1);
   }
   create_tags_nest($nest_l[$nest_index], $nest_t[$nest_index], 
                    $nest_r[$nest_index], $nest_b[$nest_index]);
 
   sync_nl_vars_nest();
   assign_nl_vars_nest();

   # Display all nests.
   my $HOLD=$nest_index;
   if ($num_nests > 1) { reload_nests(); }
   $nest_index=$HOLD;
   $nest_id_num=$nest_index+1;

}

# _______________ NESTING ROUTINES _____________________________
#

# -----------------------------------------
# display_ll_ur 
#
# Trim the floating point value to 2 significant digits.
# -----------------------------------------
sub display_ll_ur {
    my ($L, $R, $T, $B)=@_;

    # Format to 2 sig digits.
    for $arrg ($L, $R, $T, $B) {
      $arrg=sprintf ("%.2f", $arrg);
    }
    
    # Make sure range is correct.
    if ($T >  90) { $T= 180-$T; }
    if ($T < -90) { $T=-180-$T; }

    # If domain's UR lat is < LL lat, then swap the value pairs.
 if(1) {
    if ($T < $B ) { 
      my $tmp=$T; 
      $T=$B;
      $B=$tmp; 

      $tmp=$R; 
      $R=$L;
      $L=$tmp; 
    }
 }

    # Display info in suggestion box. 
    $hint_msg="Corners LL: $B, $L, UR: $T, $R\t\tNo. Grid Points: $grid_pts";
}

# -----------------------------------------
# delete_latlon_text
#
# Clear cursor values.
# -----------------------------------------

sub delete_latlon_text { $can->delete("c_line"); }

# -----------------------------------------
# can_xy_to_ll_proj 
#
# Drag cursor to display lat/lon values.
# 
# Get the projection lat,lon pair from the canvas's cursor x,y pair.
# Make an external call to 
# `$GUI_EXE/pwrap_ll_xy_convert.exe $type $arg1 $arg2 $setsupData`
#      if $type=0, then the conversion is xy_to_ll and args are x,y.
#      if $type=1, then the conversion is ll_to_xy and args are lat,lon.
# 
# The following variables are sourced from file 
#  /tmp/vector_instructions.tk:
# $scaleby_x, $scaleby_y, $scale_xMin, $scale_yMax.
#
# And, $c_xmax is the greater of ($c_width,$c_height).
#
#    $x= ($x/$c_xmax/$scaleby_x) + $scale_xMin;
#    $y=(($y/$c_ymax/$scaleby_y) - $scale_yMax) * -1;
#
# -----------------------------------------

sub can_xy_to_ll_proj {
    if ($globalMap) { return; }

    my ($c, $c_x, $c_y)=@ARG;
    $c_x=$c->canvasx($c_x); 
    $c_y=$c->canvasy($c_y); 
    my $x=$c_x-$b_xoff;
    my $y=$c_y-$b_yoff;

    # Calculate.
    $x=$x/$c_xmax;
    $y=$y/$c_ymax;

    if( ($projection_type eq 'LC' or $projection_type eq 'RL') 
         && ${parent_hash{latcen}}<0) {
	$x=$x_frac-$x;
	$y=$y_frac-$y;
    }
    $x= ($x/$scaleby_x) + $scale_xMin;
    $y=(($y/$scaleby_y) - $scale_yMax) * -1;

    # Call exe.
    my ($p_lat,$p_lon)=split (/,/, `$pwrap_ll_xy_convert_exe 0 $x $y $setsupData`);
    $p_lat=sprintf ("%.2f",$p_lat);
    $p_lon=sprintf ("%.2f",$p_lon);
    delete_latlon_text();

    my ($ii,$jj) = &map_utils::latlon_to_ij($p_lat,$p_lon,%{ $ArrayOfHashes[0]} );
    $ii=int($ii+.5);
    $jj=int($jj+.5);
    if ($ii<=0 || $ii>$nl_var{XDIM}[0] ||
        $jj<=0 || $jj>$nl_var{YDIM}[0]) {$ii=$jj="-";}
    my $my_anchor='w';
    if ($c_x > $c_xmax*.80 ) { $my_anchor='e'; }
    $can->createText($c_x, $c_y, -fill => 'white', 
            #-text => "    $p_lat,$p_lon\n     i:$ii,j:$jj\n x:$c_x,y:$c_y",
             -text => "    $p_lat,$p_lon\n     i:$ii,j:$jj",
             -font => $legend_font,
             -anchor => $my_anchor,
             -tags => "c_line");
}

# -----------------------------------------
# xy_to_ll_proj 
#
# Get the projection lat,lon pair from the x,y pair.
# Make an external call to 
# `$GUI_EXE/pwrap_ll_xy_convert.exe $type $arg1 $arg2 $setsupData`
#      if $type=0, then the conversion is xy_to_ll and args are x,y.
#      if $type=1, then the conversion is ll_to_xy and args are lat,lon.
#
# Conversion code without dealing with LC.
#$x= ($x/$c_xmax/$scaleby_x) + $scale_xMin;
#$y=(($y/$c_ymax/$scaleby_y) - $scale_yMax) * -1;
# -----------------------------------------

sub xy_to_ll_proj {
    if ($globalMap) { return; }


    my ($c_x, $c_y)=@ARG;
    my $x=$c_x-$b_xoff;
    my $y=$c_y-$b_yoff;

    # Calculate.
    if ($c_xmax==0 ) { die ($hint_msg="divide by zero.\n",watch_cursor(0),reset_vars()); exit;}
    $x=$x/$c_xmax;
    $y=$y/$c_ymax;

    if( ($projection_type eq 'LC' or $projection_type eq 'RL') 
         && ${parent_hash{latcen}}<0) {
	$x=$x_frac-$x;
	$y=$y_frac-$y;
    }
    $x= ($x/$scaleby_x) + $scale_xMin;
    $y=(($y/$scaleby_y) - $scale_yMax) * -1;

    # Call exe.
    my ($p_lat,$p_lon)=split (/,/, `$pwrap_ll_xy_convert_exe 0 $x $y $setsupData`);
    return($p_lat,$p_lon);
}

# -----------------------------------------
# ll_proj_to_xy  
#
#
# Get the projection x,y pair from the lat,lon pair.
# Make an external call to 
# `$pwrap_ll_xy_convert_exe $type $arg1 $arg2 $setsupData`
#      if $type=0, then the conversion is xy_to_ll and args are x,y.
#      if $type=1, then the conversion is ll_to_xy and args are lat,lon.
#
# Conversion code without dealing with LC.
# $c_x=int(0.001+ ($x - $scale_xMin)*$c_xmax*$scaleby_x);
# $c_y=int(0.001+ ($scale_yMax - $y)*$c_ymax*$scaleby_y);
# -----------------------------------------
sub ll_proj_to_xy {
    my ($p_lat,$p_lon)=@ARG;

    # Call exe.
    my ($x,$y)=split (/,/, `$pwrap_ll_xy_convert_exe 1 $p_lat $p_lon $setsupData`);

    # Calculate.
    $c_x=($x - $scale_xMin)*$scaleby_x;
    $c_y=($scale_yMax - $y)*$scaleby_y;

    if( ($projection_type eq 'LC' or $projection_type eq 'RL') 
         && ${parent_hash{latcen}}<0) {
	$c_x=$x_frac-$c_x;
	$c_y=$y_frac-$c_y;
    }
    $c_x=int(0.001+ ($c_x*$c_xmax));
    $c_y=int(0.001+ ($c_y*$c_ymax));

    return($c_x,$c_y);
}

# ------------------- 
#
# ------------------- 
sub grid_spacing_adjusted {

  if ($dx_orig eq $grid_spacing_dkm) { return;}
  if ($mw->focusCurrent ne $dx_ent) { return;}

  if ($grid_spacing_dkm <= 0) { $grid_spacing_dkm=0.1; }
  $pix_per_gpoint=$pix_per_gpoint / $dx_orig * $grid_spacing_dkm;

  adjust_bbox();
  $dx_orig =$grid_spacing_dkm;

  if($nmm) { display_degrees(); }
}

sub display_degrees {

  # Calc and display grid distance in degrees.
  $nx_degrees=$grid_spacing_km * 0.00920;
  $nx_degrees= sprintf ("%.3f",$nx_degrees);
  $ny_degrees=$grid_spacing_km * 0.009120;
  $ny_degrees= sprintf ("%.3f",$ny_degrees);

  # Update grid distance label.
  $gpoint_dist_lab->configure(-text => 
               "Dist between Grid Points(km): (or $nx_degrees deg)");
  $hint_msg="$hint_msg\t\t In radian degrees NX: $nx_degrees, NY: $ny_degrees";
}

# ------------------- 
#
# ------------------- 
sub dims_adjusted {

  if ($nx_orig eq $nx_ddim && $ny_orig eq $ny_dim) { return;}

  if (!($mw->focusCurrent eq $nx_ent ||
        $mw->focusCurrent eq $ny_ent)) { return; }

  adjust_bbox(); 
  $nx_orig=$nx_ddim;
  $ny_orig=$ny_dim;
}
 
# -----------------------------------------
# adjust_bbox 
#
# Recalculate the size of the bounding box if
# the grid spacing, or the n,y dims are adjusted.
# then display new bbox. 
# -----------------------------------------
 
sub adjust_bbox {

    set_nx();
    set_dx();

    # Update needs to be pressed, a parameter has changed.
    highlight_update_button(1);

    # Recalc x,y pixel dist, either dims or pix_per_gpoint changed.
    my $xdist=($nx_dim * $pix_per_gpoint) / 2;
    my $ydist=($ny_dim * $pix_per_gpoint) / 2;

    # Get center point x,y.
    my ($cent_x, $cent_y)=ll_proj_to_xy($grid_cen_lat_cmn, $grid_cen_lon_cmn);

    # Calculate change in the size of the bounding box.
    $cent_x= $cent_x+$b_xoff;
    $cent_y= $cent_y+$b_yoff;
    my $l=$cent_x - $xdist;
    my $r=$cent_x + $xdist;
    my $t=$cent_y - $ydist;
    my $b=$cent_y + $ydist;

    # Set bounding box coords.
    $can->coords("bbox", $l, $t, $r, $b);

    # Update tags & Calculate UL & LR corners (but not grid spacing).
    # Tags are disabled regardless.
    $calc_gridSpacing=0;
    bbox_corners_etc($l,$t,$r,$b);
}

# -----------------------------------------
# render_old_bbox
#
# Recreate bounding box and tags.
# -----------------------------------------
sub render_old_bbox {

    $can->createRectangle($bbox_L, $bbox_T, $bbox_R, $bbox_B, 
                          -outline => $colorW, 
                          -fill => $colorW, 
                          -stipple => 'transparent',
                          -tags => "bbox",
                          -width => 1);
}


# -----------------------------------------
# make_bbox_from_nlist 
#
# Make a bounding box from info in a namelist
# calculate the size of the bounding box based
# namelist values then display it. 
#
# Called from sub load_domain found in srt_user_interface.pl.
# -----------------------------------------

sub make_bbox_from_nlist {

   # Set map_utils hash.
   if (map_utils_setup_d01()) {
       # Cannot work with namelist file values.
       fail_dbox("Bad $namelist_arg", 
       "Cannot load $ROOT_TEMPLATES/$domain_select. \nFile possibly contains inconsistent domain values.\nCannot continue!"); 
       $domain_select="";
       $domain_mb->cget(-menu)->invoke($domain_mode-1);
       return;
   }

   # Redefine bindings to display cursor output.
   set_bindings_to_draw(0);

   # Enable the Update button.
   set_snap_box_state(1);
   render_new_map();

   # Disable widgets.
   set_latlon_wdgt_state();
}

# -------------------------------------------------------------
#
# Get cylindircal longitude value from x, knowing the image size.
# -------------------------------------------------------------
sub x_to_cyl_lon {
    my ($x)=@ARG;
    my $lon=($x * (360 / $img_xmax) - 180);
    #my $lon=($x * (720 / $img_xmax) - 360);
    $lon=sprintf ("%.2f", $lon);
    
    return $lon;
}
# -------------------------------------------------------------
#
# Get cylindircal latitude value from y, knowing the image size.
# -------------------------------------------------------------
sub y_to_cyl_lat {
    my ($y)=@ARG;
    my $lat=(($img_ymax - $y) * (180 / $img_ymax) - 90);
    $lat=sprintf ("%.2f", $lat);
    
    return $lat;
}

### Return 1 to the calling use statement ###
1;

