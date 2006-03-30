
#------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# srt_horiz_editor.pl 
# 	Creates the Horizontal Grid Editor panel, its widgets and tools
#       and logic for switching between domains.
#
# Author: Paula McCaslin   30 Jan 2003  Original Version
# ---------------------------------------------------------------------------


#use warnings;
#use strict;
use strict 'subs';
use strict 'refs';
use subs qw/native_optionmenu/;

# ----------------------------------
# Grid Editor Frame (h_frame)
#
# ----------------------------------
sub create_horizGrid_panel {

  $h_frame=$panel2->Frame()
              ->pack(-expand => 1, -fill => 'both',
                     -padx => 5, -pady => 3);

  #--------- Create Toolbar --------
  $he_toolbar=$h_frame->Frame()
             ->pack(-expand => 0, -fill => 'x');

  $b_zoomIn=$he_toolbar->Button(-image => $zoom_in,
                -relief => 'flat',
                -highlightthick => 0,
                -takefocus => 0,
                -command => [\&zoom_in_out, 2])
       ->pack(-side => 'left');

  $b_zoomOut=$he_toolbar->Button(-image => $zoom_out,
                -relief => 'flat',
                -highlightthick => 0,
                -takefocus => 0, 
                -command => [\&zoom_in_out, -2])
       ->pack(-side => 'left', -padx => 1); 

  $b_info=$he_toolbar->Button(-image => $info_im,
                -relief => 'flat',
                -highlightthick => 0,
                -takefocus => 0, 
                -command => [\&creating_bbox])
       ->pack(-side => 'left', -anchor => 'n'); 

  # Add balloon message to inform user
  $balloon->attach($b_zoomIn, 
     -msg => "Zoom in"); 
  $balloon->attach($b_zoomOut, 
     -msg => "Zoom out"); 
  $balloon->attach($b_info,
     -msg => "Instructions on how to create \na domain bounding box.");


  #------- map display window ---------
  $he_map=$h_frame->Frame(-relief => 'sunken',
                          -bd => 2)
                  ->pack(-expand => 1, -fill => 'both', 
                         -side => 'left', -anchor => 'nw');

  #panel grip used to resize 
  $h_frame->Adjuster(-widget => $he_map, -side => 'left')
                   ->pack(-fill => 'both', -side => 'left');

  #map widgets window
  $he_sel=$h_frame->Scrolled('Pane', 
                             (-relief => 'sunken', -bd => 2,
                             -width => 418, 
                             -sticky => 'nwes'),
                             -scrollbars => 'os')
                  ->pack(-side => 'left', -anchor => 'nw',
                         -expand => 1, -fill => 'both');
  $he_sel->Subwidget('xscrollbar')->configure(-width => 10);
  #$he_sel->packPropagate(1);

  #------- Create Canvas for Map -------------

  $scroll_can = $he_map->Scrolled('Canvas', 
                                       (-cursor => 'crosshair', 
                                        -bg => $bg_canvas, 
                                        -width => 340), 
                                  -scrollbars => 'osoe');
  $scroll_can->Subwidget('xscrollbar')->configure(-width => 10);
  $scroll_can->Subwidget('yscrollbar')->configure(-width => 10);

  # NOTE: pack canvas in a separate line to be able to expand & fill!
  $scroll_can->pack(-expand => 1, -fill => 'both',
                    -padx => 0, -pady => 0);
  $can=$scroll_can->Subwidget('canvas');
  $c_photo=$can->Photo;

  #$dmap_image=$ce_dmap_image;
  #$c_photo=$dmap_image;

  create_hori_widgets();
}

sub create_hori_widgets {

  #---- Create NoteBook & Panels ----
  $h_nb = $he_sel->NoteBook(-backpagecolor => $bg_canvas,
                           )
         ->pack(-expand => 1, -fill => 'both', 
                -padx => 0, -pady => 0);

  my @panel_tag= ('',
               'Domain Grid',
               'Nest Grid',
               );


  my $h_ival=1;
  $h_panel1=$h_nb->add($h_ival, -label => $panel_tag[$h_ival],
                                -raisecmd => [\&present_nest_widgets, 0]
                      );
  $h_ival=2;
  $h_panel2=$h_nb->add($h_ival);
  $h_nb->raise($h_ival);
  $h_nb->pageconfigure($h_ival, -label => $panel_tag[$h_ival],
                                -raisecmd => [\&present_nest_widgets, 1],
                      );

  # Set to 'null' since nests do not apply to LAPS.
  if ($model_name eq "LAPS" ) { $h_nb->pageconfigure(2, -label => "              "); } 

  $nest_panel=0;
  $h_nb->raise(1);

  #---

  $balloon->attach($h_nb,
              -balloonmsg => \$msg,
              -motioncommand => sub {
                  my($h_nb,$x,$y)=@ARG;
                  $x-=$h_nb->rootx; # Convert screen coords to widget coords.
                  $y-=$h_nb->rooty;
                  my $tab = $h_nb->identify($x,$y);
                  if (defined $tab && $tab == 2) {
                      #$num="tpanel$tab";
                      $msg="Create and edit nests\n for your domain.";
                      #if ($model_name eq "LAPS" ) { 
                      #$msg="Create and edit nests\n for your domain. \n\nDoes NOT apply to LAPS."; }

                      0; # show balloon
                  } else {     
                      1; # cancel balloon
                  }
              });

  #---

  create_parent_ui();
  create_nest_ui(); 

}

sub create_parent_ui {

  #--------- Horizontal Information -------
  $h_sel=$h_panel1->Frame(-relief => 'flat')
                ->pack(-anchor => 'w', -expand => 1, -fill => 'both'); 


  #--------- BCD Information -------------
  $bcd=$h_sel->Frame(-relief => 'groove', -bd => 2)
             ->pack(-anchor => 'w',
                    -padx => 5, -pady => 5,
                    -ipadx => 10, -ipady => 7);
  
    #--------- Grid Information -------
    $bcd->gridRowconfigure(0, -minsize => 3); 

    $bcd->gridColumnconfigure( 0, -minsize => 5); 
    $bcd->gridColumnconfigure( 1, -minsize => 185); 
    $bcd->gridColumnconfigure( 5, -minsize => 2); 

    #--------- BCD Label -------------
    $bcd->Label(-text => "Map Files")
        ->place( -x => 5, -y => -7);

    #----- Activate menulist as it is created.
   if (0) {
    $global_mb=global_menu_list();
    $global_choice=$dmap_image;
    $global_choice='Global Cylindrical';

    # Add balloon messages to inform.
    $balloon->attach($global_mb, -msg => "Choice of Global, Norther Hem, \nor Southern Hem.");
    $roww=2; 
   }

    $roww=1; 
    $bcd_mb=bcd_menu_list();
    $bcd_choice=$bcdFileDefault;
    update_bcd_list();

    # Add balloon messages to inform.
    $balloon->attach($bcd_mb, -msg => "List of cartographic map files.");


if(0){
    $roww=3; 
    my $clr_lab=$bcd->Label(-text => "Edit Color for:",
                         -fg => $colorN) 
                 ->grid( -row => $roww, 
                         -column => 1, 
                         -sticky => 'w');

    my @map_line_labs = 
               ("Lat, Lon Reference Lines", # Lat, Lon reference lines (white)
                "Continent & State Lines",  # Continents, State line (coral)
                "User Selected Maps",       # BCD lines (khaki)
    );
    my @map_lines_arg = ("latlon_lines",   # Lat, Lon reference lines (white)
                         "map_lines",      # Continents, State line (coral)
                         "bcd_lines",      # BCD lines (khaki)
    );

    $linetype_mb=$bcd->Menubutton(-indicator => 0,
                             -tearoff => 0,
                             -relief => 'raised',
                             -bd => 2,
                             -text => "Line Type",
                             -width => 23,
                             -justify => 'right',
                             )
                   ->grid( -row => $roww,
                           -column => 3,
                           -columnspan => 3,
                         );

    # Add balloon messages to inform.
    $balloon->attach($linetype_mb, -statusmsg => "Select color for map lines.");

    my $i=0;
    until ($i>2) {
        $linetype_mb->command(-label => $map_line_labs[$i],
                              -command => [\&selectColor, $map_lines_arg[$i]]
                          );
        $i++;
    }
}

  
  #--------- Projection Information -------
  $proj=$h_sel->Frame(-relief => 'groove', -bd => 2) 
              ->pack(-anchor => 'w',
                     -padx => 5, -pady => 5,
                     -ipadx => 10, -ipady => 7);

  
    #--------- Grid Information -------
    $proj->gridRowconfigure(0, -minsize => 3); 
    $proj->gridRowconfigure(6, -minsize => 3); 

    $proj->gridColumnconfigure( 0, -minsize => 5); 
    $proj->gridColumnconfigure( 1, -minsize => 210); 
    $proj->gridColumnconfigure( 2, -minsize => 10); 
    $proj->gridColumnconfigure( 3, -minsize => 59); 
    $proj->gridColumnconfigure( 6, -minsize => 3); 
  

    my $rrow=2;
    #--------- Projection Selector ----------
    $maproj_mb_lab=$proj->Label( -text =>'Map Projection:', -fg => $colorN) 
                      ->grid( -row => $rrow, 
                              -column => 1, 
                              -columnspan => 2, 
                              -sticky => 'w');

    $maproj_mb=$proj->Menubutton( -indicator => 1,
                             -tearoff => 0,
                             -relief => 'raised',
                             -bd => 2,
                             -width => 15,
                             -justify => 'right',
                             -activebackground => $update_color,
                             -textvariable => \$proj_label)
                   ->grid( -row => $rrow,
                           -column => 3,
                           -columnspan => 3,
                           -sticky => 'w');

    # Add balloon messages to inform.
    $balloon->attach($maproj_mb, -statusmsg => "Select projection type.");

    $i=0;
    until ($i>$max_proj_choices) {
        $maproj_mb->command( -label => $proj_choices[$i][0],
                           -command => [\&set_proj_vars, $i] );
        $i++;
    }


    #---- Standard Lon Widget -------
    $proj->Label( -text => "Standard Longitude:", -fg => $colorN)
        ->grid( -row => ++$rrow,
                -column => 1,
                -sticky => 'w');

    $lonEnt=$proj->Entry( -width => 8,
                          -justify => 'right',
                          -textvariable => \$stdlon,
                          -fg => $colorG,
                          )
        ->grid( -row => $rrow,
                -column => 4,
                -sticky => 'e');

    my $lon=$proj->Frame( -highlightthickness => 0) 
              ->grid( -row => $rrow, 
                      -column => 5, 
                      -sticky => 'w');
  
      $lon_up=$lon->Button( -bitmap => $small_up, 
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'stdlon', 0.1, $lonEnt]
                  ) ->pack();
      $lon_dw=$lon->Button( -bitmap => $small_down,
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'stdlon', -.1, $lonEnt]
                  ) ->pack();

    #---- Standard Lat Lat2 Widgets -------
  
    $proj->Label( -text => "True Latitude 1:", -fg => $colorN)
        ->grid( -row => ++$rrow,
                -column => 1,
                -sticky => 'w');

    $latEnt=$proj->Entry( -width => 8,
                          -justify => 'right',
                          -textvariable => \$truelat1)
        ->grid( -row => $rrow,
                -column => 4,
                -sticky => 'e');

    my $lat=$proj->Frame( -highlightthickness => 0) 
              ->grid( -row => $rrow, 
                      -column => 5, 
                      -sticky => 'w');
  
      $lat_up=$lat->Button( -bitmap => $small_up, 
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'truelat1', 0.1,
                                $latEnt]) ->pack();

      $lat_dw=$lat->Button( -bitmap => $small_down,
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'truelat1', -.1,
                                $latEnt]) ->pack();
    #---
  
    $proj->Label( -text => "True Latitude 2:", -fg => $colorN)
         ->grid( -row => ++$rrow,
                 -column => 1,
                 -sticky => 'w');

    $lat2Ent=$proj->Entry( -width => 8,
                           -justify => 'right',
                           -textvariable => \$truelat2)
         ->grid( -row => $rrow,
                 -column => 4,
                 -sticky => 'e');

    my $lat2=$proj->Frame( -highlightthickness => 0) 
                       ->grid( -row => $rrow, 
                               -column => 5, 
                               -sticky => 'w');
  
      $lat2_up=$lat2->Button(-bitmap => $small_up, 
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'truelat2', 0.1, $lat2Ent]) ->pack();
      $lat2_dw=$lat2->Button( -bitmap => $small_down,
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'truelat2', -.1, $lat2Ent]) ->pack();

    # Add blank line.
    ++$rrow;
  
    #---
    $proj->Label(-text => "Centerpoint Longitude:", -fg => $colorN) 
         ->grid( -row => ++$rrow,
                 -column => 1,
                 -sticky => 'w');

    $clonEnt=$proj->Entry(-width => 8,
                          -justify => 'right',
                          -textvariable => \$grid_cen_lon_cmn)
         ->grid( -row => $rrow, 
                 -column => 4, 
                 -sticky => 'e');
   
    $clon=$proj->Frame( -highlightthickness => 0) 
               ->grid( -row => $rrow, 
                       -column => 5,
                       -sticky => 'w');
  
      $clon_up=$clon->Button(-bitmap => $small_up, 
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'grid_cen_lon_cmn', 0.1,
                                 $clonEnt]) ->pack();
  
      $clon_dw=$clon->Button(-bitmap => $small_down,
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'grid_cen_lon_cmn', -.1,
                                 $clonEnt]) ->pack();
    #---
 
    $proj->Label( -text => "Centerpoint Latitude:", -fg => $colorN)
         ->grid( -row => ++$rrow,
                 -column => 1,
                 -sticky => 'w');

    $clatEnt=$proj->Entry(-width => 8,
                          -justify => 'right',
                          #-validate => 'key',
                          #-validatecommand => \&valid_numeric,
                          #-invalidcommand => sub {$mw->bell},
                          -textvariable => \$grid_cen_lat_cmn)
         ->grid( -row => $rrow,
                 -column => 4,
                 -sticky => 'e');
    my $clat=$proj->Frame( -highlightthickness => 0) 
               ->grid( -row => $rrow, 
                       -column => 5, 
                       -sticky => 'w');
  
      $clat_up=$clat->Button(-bitmap => $small_up, 
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'grid_cen_lat_cmn', 0.1,
                                  $clatEnt]) ->pack();
    
      $clat_dw=$clat->Button(-bitmap => $small_down,
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'grid_cen_lat_cmn', -.1,
                                  $clatEnt]) ->pack();
    

    #--------- Projection Label -------------
    $proj->Label(-text => "Projection (degrees)")
         ->place(-x => 5, -y => -7);
  
     
  
  #--------- Horiz Grid Information -------------
  $h_grid=$h_sel->Frame(-relief => 'groove', -bd => 2)
           ->pack(-anchor => 'w',
                  -padx => 5, -pady => 5,
                  -ipadx => 10, -ipady => 7,
                  -after => $proj);
  @pInfo_h_grid=$h_grid->packInfo();
  
    $h_grid->gridRowconfigure(0, -minsize =>  3); 

    $h_grid->gridColumnconfigure( 0, -minsize => 5); 
    $h_grid->gridColumnconfigure( 1, -minsize => 224); 
    $h_grid->gridColumnconfigure( 3, -minsize => 67); 
    $h_grid->gridColumnconfigure( 6, -minsize => 9); 


    $rrow=0;
    #--------- Horizontal Dims X, Y, Resolution Widgets --------
  
    $h_grid->Label( -text => "Horizontal Dimension X:", -fg => $colorN) 
                   ->grid( -row => ++$rrow, 
                           -column => 1, 
                           -sticky => 'w');
  
    $nx_ent=$h_grid->Entry( -width => 8, 
                            -justify => 'right',
                            -textvariable => \$nx_ddim) 
                   ->grid( -row => $rrow, 
                           -column => 4, 
                           -sticky => 'e');

   my $nx=$h_grid->Frame( -highlightthickness => 0) 
                     ->grid( -row => $rrow, 
                             -column => 5, 
                             -sticky => 'w');

      $nx_up=$nx->Button(-bitmap => $small_up, 
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'nx_ddim', 1, $nx_ent]) ->pack();
      $nx_dw=$nx->Button(-bitmap => $small_down,
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'nx_ddim',-1, $nx_ent]) ->pack();
    #---
    
    $h_grid->Label( -text => "Horizontal Dimension Y:", -fg => $colorN) 
                   ->grid( -row => ++$rrow, 
                           -column => 1, 
                           -sticky => 'w');
  
    $ny_ent=$h_grid->Entry( -width => 8,
                          -justify => 'right',
                          -textvariable => \$ny_dim) 
                 ->grid( -row => $rrow, 
                         -column => 4, 
                         -sticky => 'e');

    my $ny=$h_grid->Frame( -highlightthickness => 0) 
                     ->grid( -row => $rrow, 
                             -column => 5, 
                             -sticky => 'w');

      $ny_up=$ny->Button(-bitmap => $small_up, 
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'ny_dim', 1, $ny_ent]) ->pack();
      $ny_dw=$ny->Button(-bitmap => $small_down,
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'ny_dim',-1, $ny_ent]) ->pack();
    #---

    $gpoint_dist_lab=$h_grid->Label(-text => "Distance between Grid Points (km):", 
                           -fg => $colorN) 
                   ->grid(-row => ++$rrow, 
                          -column => 1, 
                          -columnspan => 3,
                          -sticky => 'w');

    $dx_ent=$h_grid->Entry(-width => 8, 
                           -justify => 'right',
                           -textvariable => \$grid_spacing_dkm) 
                   ->grid(-row => $rrow, 
                          -column => 4, 
                          -sticky => 'e');
    my $dx=$h_grid->Frame(-highlightthickness => 0) 
               ->grid(-row => $rrow, 
                      -column => 5, 
                      -sticky => 'w');

      $dx_up=$dx->Button(-bitmap => $small_up, 
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'grid_spacing_dkm', 1,
                                $dx_ent]) ->pack();
      $dx_dw=$dx->Button(-bitmap => $small_down,
                  -highlightthickness => 0, 
                  -takefocus => 0,
                  -command => [ \&increment_value, 'grid_spacing_dkm', -1,
                                $dx_ent]) ->pack();
  
    #--------- Grid Label -------------------
    $h_grid->Label(-text => "Grid")
           ->place(-x => 5, -y => -7);
    $h_grid->packForget();

    #------Hide Grid Widgets ----------------
    $h_grid_hide=$h_sel->Frame(-relief => 'flat', -bd => 2,
                               -height => 123,
                               -width => 380)
                       ->pack(@pInfo_h_grid, -fill => 'y');
  
    # Add blank line.
    ++$rrow;
  
  #--------- Mode Information -------------
  $mode=$h_sel->Frame(-relief => 'groove', -bd => 2)
              ->pack(-anchor => 'w', 
                     -padx => 5, -pady => 5,
                     -ipadx => 10, -ipady => 5,
                                      );
  @pInfo_mode=$mode->packInfo();
  
    #----- Selectors --------

    # Buttons to restrict centerpoint, etc.
    #--------- Grid Information -------
    $mode->gridRowconfigure(0, -minsize => 5); 
    $mode->gridColumnconfigure( 1, -minsize => 17); 
    $mode->gridColumnconfigure( 3, -minsize => 17); 
    $mode->gridColumnconfigure( 5, -minsize => 17); 

    $grid_val_restrict=1;
    $rb_edit_bbox=$mode->Radiobutton(
                    -text => 'Domain Bounding Box',
                    -value => 1, 
                    -variable => \$grid_val_restrict,
                    -fg => $colorN,
                    -activeforeground => $colorN,
                    #-command => [\&set_gridSpacing, 1],
                    -command => [\&restrict_grid_var_calc],
                    -indicatoron => 1)->grid(-row => 1, -column => 2);
  
    $rb_edit_parms=$mode->Radiobutton(
                    -text => 'Grid & Projection Values',
                    -value => 0,
                    -variable => \$grid_val_restrict,
                    -fg => $colorN,
                    -activeforeground => $colorN,
                    -command => [\&restrict_grid_var_calc],
                    -indicatoron => 1)->grid(-row => 1, -column => 4);

    # Add balloon messages to inform
    $balloon->attach($rb_edit_bbox,
       -msg => "Change domain bounding box size\nby interactively using the cursor.");
    $balloon->attach($rb_edit_parms,
       -msg => "Change domain bounding box size\nby editing 'Grid' values.");


    #--------- Mode Label -------------
    $mode->Label(-text => "Fine scale editing mode")
         ->place( -x => 5, -y => -7);


  #------Action Buttons Info ----------------
  $action=$h_sel->Frame(-relief => 'groove', -bd => 2)
                ->pack(-anchor => 'w',
                       -padx => 5, -pady => 5,
                       -ipadx => 10, -ipady => 7,
                       );
  
    $action->gridRowconfigure(0, -minsize => 5); 
    $action->gridColumnconfigure(0, -minsize => 10); 
    $action->gridColumnconfigure(2, -minsize => 11); 
    $action->gridColumnconfigure(4, -minsize => 11); 
    $action->gridColumnconfigure(6, -minsize => 11); 
    $action->gridColumnconfigure(8, -minsize => 11); 
  
    #--------- Button Choices ---------------
    $b_clear=$action->Button(-text => "Clear", -width => 6,
                              -command => [\&render_map_clear])
                     ->grid(-row => 1, -column => 1); 

    $b_startover=$action->Button(-text => "Start Over", -width => 9,
                              -state => 'disabled',
                              -command => [\&render_old_map])
                     ->grid(-row => 1, -column => 3);
  
    $b_resetvars=$action->Button(-text => "Reset Values", -width => 10,
                              -command => [\&reinstate_gen_map_vars])
                     ->grid(-row => 1, -column => 5);
  
    $b_update=$action->Button(-text => "Update Map", -width => 9,
                              -command => [\&wrap_render_new_map])
                     ->grid(-row => 1, -column => 7);
  
  
    # Add balloon messages to inform.
    $balloon->attach($b_clear,
       -msg => "Clear everything.");
    $balloon->attach($b_startover,
       -msg => "Completely undo all 'Update Map' command(s)\nbut leave Domain bounding box.");
    $balloon->attach($b_resetvars,
       -msg => "Reset current Domain bounding box and Grid\nvalue changes since the 'Update Map'\nwas last pressed.");
    $balloon->attach($b_update,
       -msg => "Render new map in projection space\n defined by the domain bounding box.");


    #------ MOAD Actions Label ------------------
    $action->Label(-text => "Actions")
           ->place( -x => 5, -y => -7);
  
    # Need these lists in order to enable/disable widgets.
    @hori_wdgt_list  =($dx_ent,  $dx_up,    $dx_dw, 
                       $nx_ent,  $nx_up,    $nx_dw, 
                       $ny_ent,  $ny_up,    $ny_dw,
                       $latEnt,  $lat_up,   $lat_dw,
                       $clonEnt, $clon_up,  $clon_dw,
                       $clatEnt, $clat_up,  $clat_dw);

    @lon_wdgt_list   =($lonEnt,  $lon_up,   $lon_dw);  # For ME.
    @lat_wdgt_list   =($latEnt,  $lat_up,   $lat_dw);
    @lat2_wdgt_list  =($lat2Ent, $lat2_up,  $lat2_dw); # For ME and PS. 

    @cenlat_wdgt_list=($clonEnt, $clon_up,  $clon_dw);
    @cenlon_wdgt_list=($clatEnt, $clat_up,  $clat_dw);

    @entry_wdgt_list =($lonEnt, $latEnt, $lat2Ent, $clonEnt, $clatEnt);
    @entry_wdgt_list2=($dx_ent, $nx_ent, $ny_ent);

    bind_latlon_entries(1);

}  

# _______________ BEGIN NESTING ROUTINES _____________________________

sub create_nest_ui {

  #--------- Horizontal Information -------
  my $h_sel_nest=$h_panel2->Frame(-relief => 'flat')
                ->pack(-anchor => 'w', -expand => 1, -fill => 'both'); 

  #--------- Parm Information ---------
  $nest_fr=$h_sel_nest->Frame(-relief => 'groove', -bd => 2) 
              ->pack(-anchor => 'w',
                     -padx => 5, -pady => 5,
                     -ipadx => 10, -ipady => 10);

    #--------- Grid Information -------
    $nest_fr->gridRowconfigure(0, -minsize => 3); 
    $nest_fr->gridRowconfigure(5, -minsize => 10); 
    $nest_fr->gridColumnconfigure( 0, -minsize => 2); 
    $nest_fr->gridColumnconfigure( 1, -minsize => 232); 

    my $rrow=2;
    #--------- Domain Selector ----------
    $nestid_mb_lab=$nest_fr->Label( -text =>'Domain ID:', -fg => $colorN) 
                      ->grid( -row => $rrow, 
                              -column => 1, 
                              -columnspan => 2, 
                              -sticky => 'w');

    $nestid_mb=$nest_fr->Menubutton( -indicator => 1,
                             -tearoff => 0,
                             -relief => 'raised',
                             -bd => 2,
                             -width => 4,
                             -activebackground => $update_color,
                             -anchor => 'e',
                             -textvariable => \$nest_id_label)
                   ->grid( -row => $rrow,
                           -column => 6,
                           -columnspan => 2,
                           -sticky => 'w');

    $i=2;
    until ($i> $max_nests) {
        $nestid_mb->command(-label => $i,
                            -foreground => $colorG,
                            -command => [\&set_nest_index, $i] );
        $i++;
    }

    #--------- Nest Parent ID Widgets -------
    $nest_parent_text="Parent ID:";
    $nest_fr->Label( -text => $nest_parent_text, -fg => $colorN)
         ->grid( -row => ++$rrow,
                 -column => 1,
                 -columnspan => 3,
                 -sticky => 'w');

    $nestparent_mb=$nest_fr->Menubutton( -indicator => 1,
                             -tearoff => 0,
                             -relief => 'raised',
                             -bd => 2,
                             -width => 4,
                             -activebackground => $update_color,
                             -anchor => 'e',
                             -textvariable => \$parent_of_nest)
                   ->grid( -row => $rrow,
                           -column => 6,
                           -columnspan => 2,
                           -sticky => 'w');
    update_nest_parent_mb();
    
    #--------- Nest Grid Spacing Widgets -------
    
    $nest_spacing_label="Grid Spacing Ratio to Parent:";
    $nest_fr->Label( -textvariable => \$nest_spacing_label, -fg => $colorN)
         ->grid( -row => ++$rrow,
                 -column => 1,
                 -columnspan => 4,
                 -sticky => 'w');

    $nest_ratio="3";
    $nestratio_mb=$nest_fr->Menubutton( -indicator => 1,
                             -tearoff => 0,
                             -relief => 'raised',
                             -bd => 2,
                             -width => 4,
                             -activebackground => $update_color,
                             -anchor => 'e',
                             -textvariable => \$nest_ratio)
                   ->grid( -row => $rrow,
                           -column => 6,
                           -columnspan => 2,
                           -sticky => 'w');
    $i=1;
    $nest_ratio_max=5;
    until ($i> $nest_ratio_max) {
        $nestratio_mb->command(-label => $i,
                               -command => [\&wrap_set_nest_ratio, $i] );
        $i++;
    }

    #--------- Nest LL I,J Widgets -------
    $rrow++;
    $nest_fr->Label( -text => "Lower Left I,J:", -fg => $colorN)
         ->grid( -row => ++$rrow,
                 -column => 1,
                 -sticky => 'w');

    $lli_ent=$nest_fr->Entry(-width => 6,
                          -justify => 'right',
                          -textvariable => \$nest_lli)
         ->grid( -row => $rrow,
                 -column => 4,
                 -sticky => 'e');
    my $lli_nest=$nest_fr->Frame( -highlightthickness => 0) 
               ->grid( -row => $rrow, 
                       -column => 5, 
                       -sticky => 'w');
  
      $lli_up=$lli_nest->Button(-bitmap => $small_up, 
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'nest_lli', 1,
                                  $lli_ent]) ->pack();
    
      $lli_dw=$lli_nest->Button(-bitmap => $small_down,
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'nest_lli', -1,
                                  $lli_ent]) ->pack();
    
    #---
 
#    $nest_fr->Label( -text => "Lower Left J:", -fg => $colorN)
#         ->grid( -row => ++$rrow,
#                 -column => 1,
#                 -sticky => 'w');

    $llj_ent=$nest_fr->Entry(-width => 6,
                          -justify => 'right',
                          -textvariable => \$nest_llj)
         ->grid( -row => $rrow, 
                 -column => 6, 
                 -sticky => 'e');
   
    my $llj_nest=$nest_fr->Frame( -highlightthickness => 0) 
               ->grid( -row => $rrow, 
                       -column => 7,
                       -sticky => 'w');
  
      $llj_up=$llj_nest->Button(-bitmap => $small_up, 
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'nest_llj', 1,
                                 $llj_ent]) ->pack();
  
      $llj_dw=$llj_nest->Button(-bitmap => $small_down,
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'nest_llj', -1,
                                 $llj_ent]) ->pack();


    #--------- Nest UR I,J Widgets -------
    $nest_uri_label="Upper Right I, J:";
    $nest_fr->Label( -textvariable => \$nest_uri_label, -fg => $colorN)
         ->grid( -row => ++$rrow,
                 -column => 1,
                 -sticky => 'w');

    $uri_ent=$nest_fr->Entry(-width => 6,
                          -justify => 'right',
                          -textvariable => \$nest_uri)
         ->grid( -row => $rrow,
                 -column => 4,
                 -sticky => 'e');
    my $uri_nest=$nest_fr->Frame( -highlightthickness => 0) 
               ->grid( -row => $rrow, 
                       -column => 5, 
                       -sticky => 'w');
  
      $uri_up=$uri_nest->Button(-bitmap => $small_up, 
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'nest_uri', 1,
                                  $uri_ent]) ->pack();
    
      $uri_dw=$uri_nest->Button(-bitmap => $small_down,
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'nest_uri', -1,
                                  $uri_ent]) ->pack();
    #---
#    $nest_urj_label="Upper Right J:";
#    $nest_fr->Label( -textvariable => \$nest_urj_label, -fg => $colorN)
#         ->grid( -row => ++$rrow,
#                 -column => 1,
#                 -sticky => 'w');

    $urj_ent=$nest_fr->Entry(-width => 6,
                          -justify => 'right',
                          -textvariable => \$nest_urj)
         ->grid( -row => $rrow, 
                 -column => 6, 
                 -sticky => 'e');
   
    my $urj_nest=$nest_fr->Frame( -highlightthickness => 0) 
               ->grid( -row => $rrow, 
                       -column => 7,
                       -sticky => 'w');
  
      $urj_up=$urj_nest->Button(-bitmap => $small_up, 
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'nest_urj', 1,
                                 $urj_ent]) ->pack();
  
      $urj_dw=$urj_nest->Button(-bitmap => $small_down,
                    -highlightthickness => 0, 
                    -takefocus => 0,
                    -command => [\&increment_value, 'nest_urj', -1,
                                 $urj_ent]) ->pack();

    @nest_widgets_list=($lli_ent, $lli_up, $lli_dw,
                        $llj_ent, $llj_up, $llj_dw,
                        $uri_ent, $uri_up, $uri_dw,
                        $urj_ent, $urj_up, $urj_dw);

    #--------- Projection Label -------------
    $nest_fr->Label(-text => "Nest Parameters")
         ->place(-x => 5, -y => -7);

    # Bind key events to update the nestbox and its attributes.
    for my $event ("<Return>", "<Double-Return>") { #, "<Leave>") {
       $lli_ent->Tk::bind  ($event, [\&adjust_nestbox]);
       $llj_ent->Tk::bind  ($event, [\&adjust_nestbox]);
       $uri_ent->Tk::bind  ($event, [\&adjust_nestbox]);
       $urj_ent->Tk::bind  ($event, [\&adjust_nestbox]);
    }
    for my $event ("<ButtonRelease>") {
       $lli_up->Tk::bind   ($event, [\&adjust_nestbox]);
       $lli_dw->Tk::bind   ($event, [\&adjust_nestbox]);
       $llj_up->Tk::bind   ($event, [\&adjust_nestbox]);
       $llj_dw->Tk::bind   ($event, [\&adjust_nestbox]);
       $uri_up->Tk::bind   ($event, [\&adjust_nestbox]);
       $uri_dw->Tk::bind   ($event, [\&adjust_nestbox]);
       $urj_up->Tk::bind   ($event, [\&adjust_nestbox]);
       $urj_dw->Tk::bind   ($event, [\&adjust_nestbox]);
    }

   #-----------------------------------------#

  #--------- Summary Information ---------
  $n_summary_fr=$h_sel_nest->Frame(-relief => 'groove', -bd => 2) 
              ->pack(-anchor => 'w',
                     -padx => 5, -pady => 5,
                     -ipadx => 10, -ipady => 10);

    #------ Nest Summary Table ---------------

    $nest_text=$n_summary_fr->Scrolled('Text', 
                                       (-height => 8,
                                        -width => 49,
                                        -padx => 1,
                                        -wrap => 'none',
                                        -highlightthickness => 0, 
                                        -spacing1=> 2,
                                        -exportselection => 1,
                                        -fg => $colorN,
                                        -relief => 'flat'),
                             -takefocus => '0',
                             -scrollbars => 'osoe')
               #->grid( -row => 1, -column => 1, -sticky => 's');
                ->pack(-ipadx => 5, -pady => 0, -side => 'bottom');
   $nest_text->Subwidget('xscrollbar')->configure(-width => 10);
   $nest_text->Subwidget('yscrollbar')->configure(-width => 10);
   $nest_text->tagConfigure('color', -foreground => $normal_color);

    #------ Summary Actions Label ------------------
    $n_summary_fr->Label(-text => "Summary of Domains")
           ->place( -x => 5, -y => -7);

   #-----------------------------------------#

   #------Action Buttons Info ----------------
   my $n_action=$h_sel_nest->Frame(-relief => 'groove', -bd => 2)
                ->pack(-anchor => 'w',
                       -padx => 5, -pady => 5,
                       -ipadx => 10, -ipady => 7,
                       );
  
    $n_action->gridRowconfigure(0, -minsize => 7); 
    $n_action->gridRowconfigure(2, -minsize => 7); 
    $n_action->gridColumnconfigure(0, -minsize => 8); 
    $n_action->gridColumnconfigure(2, -minsize => 8); 
    $n_action->gridColumnconfigure(4, -minsize => 8); 
    $n_action->gridColumnconfigure(6, -minsize => 8); 
  
    #--------- Button Choices ---------------
    $n_restore=$n_action->Button(-text => "Restore All Nests", -width => 16,
                             -command => [\&reinstate_gen_nest_vars])
                     ->grid(-row => 1, -column => 1); 

    $n_del=$n_action->Button(-text => "Delete Nest", -width => 12,
                             -command => [\&delete_nest])
                     ->grid(-row => 1, -column => 3); 

    $n_erase=$n_action->Button(-text => "Erase Box", -width => 12,
                             -command => [\&erase_nestbox])
                     ->grid(-row => 1, -column => 5); 

    set_button_state(0,$n_restore);
    # Add balloon messages to inform.
    $balloon->attach($n_restore,
       -msg => "Restore all nest previously written to disk.");
    $balloon->attach($n_del,
       -msg => "Delete all information about this nest\n(including ALL dependants).");
    $balloon->attach($n_erase,
       -msg => "Erases the nest bounding box and\nallows you to draw another one\n(leaving Parent ID and Ratio to Parent unchanged).");
    @nest_actions=($n_del, $n_erase);

    #------ Nests Actions Label ------------------
    $n_action->Label(-text => "Actions")
           ->place( -x => 5, -y => -7);
}

# -----------------------------------
# set_nest_index
#
# Set map nest menubar widget to the domain nest ID number.
# -----------------------------------
sub set_nest_index { 
    
    $nestid_mb->configure(-bg => $bg_color);
    set_button_state(0,$n_restore);

    # Set index.
    ($nest_id_num)=@_; 
    $nest_id_label="d0$nest_id_num"; 
    $nest_index=$nest_id_num-1;
   
    # If user requested additional nest w/o filling with values.  
    sync_nl_maxs_nest();

  if($nest_panel) {

    #$nestBox[$nest_index]="nestBox";
    $nestBox[$nest_index]=$nest_id_label;
    $nestLab[$nest_index]=$nest_id_label;

    # Update menus and labels.
    $nest_parent_text="Parent ID:";
    update_nest_parent_mb();
    configure_domain_id_mb();

    if ($nest_id_num == 1) { $nest_id_label=""; return; } # Reset widgets, only.

    # Hide tags.
    hide_tags_nest();

    # Update values.
    assign_nl_vars_nest();
    sync_nl_vars_nest();


    # Set nest: new, or existing.
    if ($parent_of_nest eq "") {

       $nestparent_mb->configure(-bg => $colorY2);

       set_button_list_state(0,@nest_widgets_list); 
       remove_parentbox();

       # --- New nest ---
       # Allow a user to draw nest.
       set_bindings_to_draw(2);
       set_button_list_state(1,@nest_actions);

       # For nest_id_num=2, the only parent choice =1.
       if ($nest_id_num == 2) { 
          set_parent_of_nest(1); 
          $hint_msg="To DRAW A NEST, place the mouse cursor in the canvas.";
          return;
       } 

       # For all other nest_id_nums, ask user.
       my $my_guess=$nest_id_num-1;
       set_parent_of_nest($my_guess);
       $hint_msg="To DRAW A NEST, place the mouse cursor in the canvas.";

    } else { #if ($nest_l[$nest_index] > 0) {

       # Update values.
       #assign_nl_vars_nest();
       #sync_nl_vars_nest();

       # Reset parent hash.
       (%parent_hash)=(%{ $ArrayOfHashes[$nest_index]} );


       # --- Existing nest(s) ---
       # Color active subnest yellow. 
       create_tags_nest($nest_l[$nest_index], $nest_t[$nest_index], 
                        $nest_r[$nest_index], $nest_b[$nest_index]);

       # Color parent nest yellow. 
       render_parent_bbox();
       update_parent_id_label();
       set_button_state(1,$n_erase);

       # Test to see if nest is a parent. 
       # If so, do not proceed without confirmation.
       for ($ix=2; $ix <= $num_nests; $ix++) {
          if ($nest_id_num == $nl_var{PARENT_ID}[$ix]) {
             # This nest domain is a pArEnT.
             set_button_state(0,$n_erase);
             hide_tags_nest();
          }
       }
 
       # Allow a user edit nest.
       set_bindings_to_draw(0);
    }

  }
}

# -----------------------------------
# wrap_draw_all_nests
#
# Called from present_nest_widgets.
# -----------------------------------
sub wrap_draw_all_nests {

     draw_all_nests();
     update_nest_id_mb();
     set_nest_index($num_nests);
}

# -----------------------------------
# draw_all_nests
#
# Called from either:
# wrap_draw_all_nests (new) AND
# delete_dependants (modified).
# -----------------------------------
sub draw_all_nests { 
     my($my_end)=@_;
     watch_cursor(1);
     remove_oldnestbox();

     for (2 .. $num_nests) { 
       # Existing nest.
       $nest_id_num=$_; 
       $nest_index=$nest_id_num-1;
       assign_nl_vars_nest();

       # Set parent hash for this nest using map_utils.
       if (!defined $nest_spacing[$parent_of_nest]) { 
          set_nest_ratio($nl_var{RATIO_TO_PARENT}[$nest_index]); 
       }
       map_utils_setup_nest($parent_of_nest); 
       convert_nest_vars_to_nestbox(0);
       fill_nest_table_entry();

     }
     #if ($can->gettags('oldnestbox')) { $can->lower("oldnestbox"); }
     watch_cursor(0);
}

# -----------------------------------
# activate_nest_widgets 
#
# 
# -----------------------------------
sub activate_nest_widgets {

     # Activate these widgets for user.
     set_button_list_state(1,@nest_widgets_list); 
     update_nest_id_mb();
}

# -----------------------------------
# update_nest_id_mb 
#
# Set the menubar colors equal to the number of nests.
# -----------------------------------
sub update_nest_id_mb {

    if ($nest_index < 1) {
       $nestid_mb->entryconfigure(0, -foreground => 'blue');
       return;
    }

    # Color list. 
    for ($ix=2; $ix <= $max_nests; $ix++) {
        if($ix > $num_nests){
           $my_color='blue';
        } else {
           $my_color=$normal_color;
        }

        if($ix > $num_nests+1 || 
           $nl_var{DOMAIN_ORIGIN_LLI}[$ix-2] eq ""){
           $my_state='disabled';
        } else {
           $my_state='normal';
        }

        $nestid_mb->entryconfigure($ix-2, -foreground => $my_color,
                                          -state => $my_state);
    }

}

# -----------------------------------
# update_nest_parent_mb 
#
# Set the menubar choices equal to the number of nests-1.
# -----------------------------------
sub update_nest_parent_mb {

  foreach ($nestparent_mb) {
    # Clear list. 
    if (defined $nestparent_mb && $nestparent_mb->cget(-menu)) {
       $nestparent_mb->cget(-menu)->delete(0, 'end'); }

    # Build list. 
    for ($i=1; $i < $nest_id_num; $i++) {
       $nestparent_mb->command(-label => $i,
                   -command => [\&set_parent_of_nest, $i] ); }
  }

  update_force_loc_mb();
}

# -----------------------------------
# set_parent_of_nest
#
# Set map nest menubar widget to the domain nest ID number.
# -----------------------------------
sub set_parent_of_nest { 

    # User must have an null sub nest before this can be set.
    if (($parent_of_nest != @_) && 
         $nest_panel && 
         $nest_lli ne "") { 
       my $disp = info_dbox("Confirm Choice", 
"The selection of the Parent ID changes the projections.
You must press 'Erase Box' before this can be set.");
       return;   
    } 

    ($parent_of_nest)=@_; 
    $nl_var{PARENT_ID}[$nest_index]=$parent_of_nest;
    $parent_of_nest_label[$nest_id_num]="d0$parent_of_nest"; 
    $nestparent_mb->configure(-bg => $bg_color);

    # Set parent hash for this nest using map_utils.
    # --------------------------------
    map_utils_setup_nest($parent_of_nest); 

    # Helpful parent domain info, after map_utils is called.
    update_parent_id_label();

    # Assign a nest ratio value -- as is is a function of the parent id.
    if ($nest_ratio eq "") { 
       set_nest_ratio(3);

    } else {
       # Assume its ratio hasn't changed, but the value 
       # for grid spacing could, so update it regardless.
       set_nest_ratio($nest_ratio); 
    }

    # Increase num_nests, if user sets nest_parent.
    sync_nl_maxs_nest();
 
    # Nest borders.
    render_parent_bbox();

     
} 

sub update_parent_id_label {

    # Helpful info about parent domain, filled in once map_utils is called.
    my $parent_index=$parent_of_nest-1;
    $nest_parent_text=
"Parent ID (nx:$nest_nx[$parent_index], ny:$nest_ny[$parent_index], dx:$nest_spacing[$parent_index] km):";      
}

# -----------------------------------
# wrap_set_nest_ratio
#
# Assign a value to both nest_ratio AND display is resulting nest_spacing value.
# -----------------------------------
sub wrap_set_nest_ratio { 

    # User must have an null sub nest before this can be set.
    if ($nest_index < $num_nests-1 && 
       ($nl_var{RATIO_TO_PARENT}[$nest_index] != @_) && 
        $nest_panel && 
        $nest_lli ne "") { 
       my $disp = info_dbox("Confirm Choice", 
"The selection of the Ratio to Parent changes the projections.
You must press 'Erase Box' before this can be set.");
       return;   
    } 
    set_nest_ratio(@_); 
}

# -----------------------------------
# set_nest_ratio
#
# Assign a value to both nest_ratio AND display is resulting nest_spacing value.
# -----------------------------------
sub set_nest_ratio { 
    ($nest_ratio)=@_; 
    if ($nest_ratio) {
       $nestratio_mb->configure(-bg => $bg_color);

       my $parent_index=$parent_of_nest-1;
       my $my_nest_spacing=sprintf ("%.3f", $nest_spacing[$parent_index]/$nest_ratio);
       $nest_spacing_label=
          "Grid Spacing Ratio to Parent ($my_nest_spacing km):";

       # Set global vars.
       $nl_var{RATIO_TO_PARENT}[$nest_index]=$nest_ratio;
       $nest_spacing[$nest_index]=$my_nest_spacing;


    } else {
       $nest_spacing_label="Grid Spacing Ratio to Parent:";
    }
}

# -------------------------------------------------
# map_utils_setup_nest
#
# See "$ROOT_INSTALL/etc/map_utils.pm";
#
# This is a Perl module for converting (i,j) to (lat,lon) and vice versa
# for any of 4 map projections:  Cylindrical Equidistant (lat-lon),
# Mercator, lambert conformal (secant and tangent), and polar
# stereographic.
#
# -------------------------------------------------
sub map_utils_setup_nest {
   my ($n_parent_id)=@_;
   my $parent_index=$n_parent_id-1;
   my $d_x;

      my($known_lat,$known_lon) = 
           &map_utils::ij_to_latlon(
                $nl_var{DOMAIN_ORIGIN_LLI}[$parent_index],
                $nl_var{DOMAIN_ORIGIN_LLJ}[$parent_index],
                %{$ArrayOfHashes[$parent_index]} );
      my($known_ri,$known_rj)=(1,1);

      $nx_nest=(($nl_var{DOMAIN_ORIGIN_URI}[$parent_index] -
                 $nl_var{DOMAIN_ORIGIN_LLI}[$parent_index]) *
                 $nl_var{RATIO_TO_PARENT}[$parent_index]) + 1;
   
      $ny_nest=(($nl_var{DOMAIN_ORIGIN_URJ}[$parent_index] -
                 $nl_var{DOMAIN_ORIGIN_LLJ}[$parent_index]) *
                 $nl_var{RATIO_TO_PARENT}[$parent_index]) + 1;
      if ($parent_index == 0) { 
	  $nest_spacing[0]=$grid_spacing_km; 
          if ($nest_ratio) {
	     $nest_spacing[1]=sprintf ("%.3f", $grid_spacing_km/$nest_ratio);
          }
      }
      $d_x=$nest_spacing[$parent_index]*1000.;

      my $my_catch=Tk::catch {%{ $ArrayOfHashes[$nest_index] } =
           &map_utils::map_set($projection_type,
                               $known_lat,$known_lon,
                               $known_ri, $known_rj, $d_x, 
                               $stdlon,
                               $truelat1,$truelat2,
                               $nx_nest,$ny_nest);
      };

   if ($my_catch eq "") { 
     warn("There is a major problem with the namelist file: $template_nl.\n"), return(1);
   }

   # Recalculate the corners base on gridded map.
   my ($latSW,$lonSW) = 
      &map_utils::ij_to_latlon(1,1,%{ $ArrayOfHashes[$nest_index]} );
   my ($latNE,$lonNE) = 
      &map_utils::ij_to_latlon($nx_nest,$ny_nest,%{ $ArrayOfHashes[$nest_index]} );
   display_ll_ur($lonSW,$lonNE,$latNE,$latSW);   
   $hint_msg="$hint_msg";  # todo: finish this.


   # Reset parent hash.
   (%parent_hash)=(%{ $ArrayOfHashes[$nest_index]} );

   $nest_nx[$parent_index]=${parent_hash{nx}};
   $nest_ny[$parent_index]=${parent_hash{ny}};
 
   $known_ri = ($nest_nx[$parent_index]+1.)*0.5;   # Exact center of grid
   $known_rj = ($nest_ny[$parent_index]+1.)*0.5;   # Exact center of grid
   ($nest_cenlat[$parent_index],$nest_cenlon[$parent_index]) = 
      &map_utils::ij_to_latlon($known_ri,$known_rj,%parent_hash);

   # Success.
   return(0);
}

# ---------------------------------------
# reset_nest_values
# 
# Delete all nest values.
# ---------------------------------------
sub reset_nest_values {

   # --- Reset Nest Interface and vars, if displayed.
   $nest_index=1;
   $num_nests=1; 
   $nest_id_label=$nest_id_num=""; 
   $nest_parent_text="Parent ID:";
   set_button_list_state(0,@nest_actions); 

   # highlight brite yellow, so user will select a value.
   $nestid_mb->configure(-bg => $colorY2);
   
   # Reset current Nesting variables.
   for $key (@nesting_vars) {
      $nl_var_max{$key}=1; #There is always at least 1, the moad.
      for $i (1 .. $max_nests) { 
        $nl_var{$key}[$i]=""; }
   }
 
   for $i (1 .. $max_nests) { 
      if (defined $nest_l[$i]) { 
         $nest_l[$i]=$nest_t[$i]=$nest_r[$i]=$nest_b[$i]=""; } }
 
   assign_nl_vars_nest();
   remove_nest_graphics();
   empty_nest_table();

}

# ---------------------------------------
# erase_nestbox
# 
# ---------------------------------------
sub erase_nestbox {

  # Clear values.
  $nest_lli=$nest_llj=$nest_uri=$nest_urj="";  
  sync_nl_vars_nest;

  if ($nest_index <= 0) { return; }

  # Erase all nestboxes.
  remove_nestbox();
  remove_oldnestbox();
  set_button_list_state(0,@nest_widgets_list); 
  
  # Redraw nests except for the one just erased.
  if ($nest_index > 1) { 
     reload_nests(1);
  }

  # Get ready to draw again.
  set_bindings_to_draw(2);

  # Recover, if desired.
  if ($nl_var{NUM_DOMAINS}[0] > 1) { 
     set_button_state(1,$n_restore);
  }

}

# -----------------------------------
# reinstate_gen_nest_vars
#
# Redraw nests including the one(s) just erased or edited,
# if they have been 'saved' to disk.
# -----------------------------------
sub reinstate_gen_nest_vars{

     set_button_state(0,$n_restore);
     if ($nl_var{NUM_DOMAINS}[0] > 1) { 
        $domain_mode=2;
        render_old_map();
        present_nest_widgets(1);
     }
}

# ---------------------------------------
# reload_nests 
# 
# Redraw all sub nests.
# No need to recalculate map_utils_setup_nest vars. 
# ---------------------------------------
sub reload_nests {
  my ($arg)=@_;
  

  watch_cursor(1);
  remove_oldnestbox();

  my $max_nests=$num_nests-$arg;
  for (2 .. $max_nests) { 
     $nest_id_num=$_; 
     $nest_index=$nest_id_num-1;
     render_old_nestbbox();
  }
  $nest_id_num++;
  $nest_index++;
  if ($can->gettags('parentbox')) { $can->raise("parentbox"); }

  watch_cursor(0);

}

# ---------------------------------------
# delete_nest
# 
# Permentaly remove nest and sub nests.
# ---------------------------------------
sub delete_nest {

   remove_nest_label();
   remove_nestbox();
   $nest_index--;
   delete_dependants(); # And, in this case the current nest too.

   if ($nest_index >= 1) { 
      assign_nl_vars_nest();

   } else {
      reset_nest_values();
      # Redefine bindings to display cursor output.
      set_bindings_to_draw(0);
      set_button_list_state(0,@nest_widgets_list); 
   }
    # Display all nests.
    if ($num_nests > 1) { draw_all_nests(); }

    # Nest borders.
    render_parent_bbox();

    # Recover, if desired.
    if ($nl_var{NUM_DOMAINS}[0] > 1) { 
       set_button_state(1,$n_restore);
    }

}

# ---------------------------------------
# ask_delete_dependants
# 
# Delete nest dependant then continue.
# ---------------------------------------
sub ask_delete_dependants {
    
    # Make sure that sub nests are all valid (check for null data).
    reduce_nl_maxs_nest();

    my $nest_range;
    if($num_nests <= 2) { 
       $nest_range="d02";
    } elsif (!$nest_panel) {
       $nest_range="d02 - d0$num_nests";
    } else {
       $nest_range=$nest_index+2;
       $nest_range="d0$nest_range - d0$num_nests";
    }

    my $ans;
    if ($num_nests > 1) { 
       $ans=okcancel_dbox("Delete Dependant Nests", 
"You are editing a parent. If you continue ALL 
dependant nests ($nest_range) will be deleted.");

       if($ans eq 'Ok'){ 
          if($nest_panel) {
             # Delete some nests.
             $nest_index++;
             delete_nest(); 
             remove_nest_graphics();

          } else {
             # Delete all nest values.
             reset_nest_values(); 

             # Redefine bindings to display cursor output.
             set_bindings_to_draw(0);
          }
       }

    }

    return($ans);
}

# ---------------------------------------
# delete_dependants
# 
# Permentaly remove dependant.
# Called from delete_nest and adjust_nestbox.
# ---------------------------------------
sub delete_dependants {

   my $my_begin=$nest_index;
   for($nest_index=$num_nests; 
       $nest_index > $my_begin; $nest_index--) {

      # Reset current Nesting variables.
      $i=$nest_index;
      for $key (@nesting_vars) {
         $nl_var_max{$key}=$nl_var_max{$key}-1;
         $nl_var{$key}[$i]=""; 
      }

      if (defined $nest_l[$i]) { 
         $nest_l[$i]=$nest_t[$i]=$nest_r[$i]=$nest_b[$i]=""; 
      }

      reduce_nl_maxs_nest();
   }
   
   if ($nest_index >= 1) {
      set_nest_index($num_nests);
   }
   update_nest_id_mb();

}

# ---------------------------------------
# fill_nest_table_entry 
#
# Summarize a list of nests and their properties in a formatted table.
# ---------------------------------------
sub fill_nest_table_entry {

    my $i=$nest_index;
    my $parent_i=$nl_var{PARENT_ID}[$nest_index];

    $nest_nx[$i]=(($nl_var{DOMAIN_ORIGIN_URI}[$i]-$nl_var{DOMAIN_ORIGIN_LLI}[$i])
       * $nl_var{RATIO_TO_PARENT}[$i]) +1;
    $nest_ny[$i]=(($nl_var{DOMAIN_ORIGIN_URJ}[$i]-$nl_var{DOMAIN_ORIGIN_LLJ}[$i])
       * $nl_var{RATIO_TO_PARENT}[$i]) +1;
  
    # Calculate center lat/lon for current nest, too.
    # Information is not known before now.
    $known_ri=((($nl_var{DOMAIN_ORIGIN_URI}[$i] -
                 $nl_var{DOMAIN_ORIGIN_LLI}[$i]) +1.)*0.5 +
                 $nl_var{DOMAIN_ORIGIN_LLI}[$i]); 
    $known_rj=((($nl_var{DOMAIN_ORIGIN_URJ}[$i] -
                 $nl_var{DOMAIN_ORIGIN_LLJ}[$i]) +1.)*0.5 +
                 $nl_var{DOMAIN_ORIGIN_LLJ}[$i]);
    ($nest_cenlat[$i],$nest_cenlon[$i]) = 
       &map_utils::ij_to_latlon($known_ri,$known_rj,%{ $ArrayOfHashes[$parent_i] } );


    # Clear current selection, if any with (S)electionClear;
    $nest_text->SelectionClear();


    # Fill nest table information.
    # ---------------------------------

    # Allow update to sigma text widget.
    $nest_text->configure(-state => 'normal'); 
    $nest_text->delete("1.0", "end");
    $nest_text->insert("end", 
       "        NX, NY\tSPACE\tPA RA\tLL IJ\tUR IJ\tPTS\tCENT LAT LON\n");
    $nest_text->tagAdd('color', "1.0","2.0");

    # @nest_disp_ary{$key}="1"; # Element [0] is always =1.
    # Format nest table.
    for $count (0 .. $num_nests-1) {
       $i=$count;
       $dm=$count+1;
       $dm="d0$dm";
       $nest_text->insert("end", sprintf (
          "%s: %s, %s\t%s\t",
          $dm,$nest_nx[$i],$nest_ny[$i],$nest_spacing[$i]) );
       if($count==0) {
         $nest_text->insert("end", sprintf (
          "%s   %s\t%s, %s\t%s, %s\t",
          1,1,1,1,
          $nl_var{XDIM}[$count],
          $nl_var{YDIM}[$count],
          ) );
       } else {
         $nest_text->insert("end", sprintf (
          "%s   %s\t%s, %s\t%s, %s\t",
          $nl_var{PARENT_ID}[$count],
          $nl_var{RATIO_TO_PARENT}[$count],
          $nl_var{DOMAIN_ORIGIN_LLI}[$count],
          $nl_var{DOMAIN_ORIGIN_LLJ}[$count],
          $nl_var{DOMAIN_ORIGIN_URI}[$count],
          $nl_var{DOMAIN_ORIGIN_URJ}[$count]) );
       }
       $nest_text->insert("end", sprintf (
          "%s\t%.2f,%.2f\t%s\n",
          $nest_nx[$i]*$nest_ny[$i],
          $nest_cenlat[$i],$nest_cenlon[$i], $dm
          ));
    }

    #$nest_text->see('1.0');

    # Return sigma text widget to orig state.
    $nest_text->configure(-state => 'disabled');

}

# ---------------------------------
# write_nest_table 
#
# Write nest information to a file.
# ---------------------------------

sub write_nest_table {

    # Write nest information to a file.
    # ---------------------------------
    my $nest_info="$ROOT_TEMPLATES/$domain_select/nest_info.txt"; 
    open(NIT,">$nest_info") or 
          warn("Difficulty writing $nest_info.\n"), return(1);

    print NIT ("        NX, NY\tSPACE\tPA RA\tLL IJ\tUR IJ\t\tPTS\tCENT LAT LON\n");

    # Format nest table.
    for $count (0 .. $num_nests-1) {
       $i=$count;
       $dm=$count+1;
       $dm="d0$dm";
       printf NIT ("%s: %s, %s\t%s\t",
          $dm,$nest_nx[$i],$nest_ny[$i],$nest_spacing[$i]);
       if($count==0) {
         printf NIT ("%s   %s\t%s, %s\t%s, %s\t\t",
          1,1,1,1,
          $nl_var{XDIM}[$count], $nl_var{YDIM}[$count]);
       } else {
         printf NIT ("%s   %s\t%s, %s\t%s, %s\t\t",
          $nl_var{PARENT_ID}[$count],
          $nl_var{RATIO_TO_PARENT}[$count],
          $nl_var{DOMAIN_ORIGIN_LLI}[$count],
          $nl_var{DOMAIN_ORIGIN_LLJ}[$count],
          $nl_var{DOMAIN_ORIGIN_URI}[$count],
          $nl_var{DOMAIN_ORIGIN_URJ}[$count]);
       }
       printf NIT ("%s\t%.2f,%.2f\t%s\n",
          $nest_nx[$i]*$nest_ny[$i],
          $nest_cenlat[$i],$nest_cenlon[$i], $dm);
    }
    close(NIT);

}

# ---------------------------------------
# empty_nest_table 
#
# Empty formatted table.
# ---------------------------------------
sub empty_nest_table {

    # Clear current selection, if any with (S)electionClear;
    $nest_text->SelectionClear();
    $nest_text->configure(-state => 'normal'); 
    $nest_text->delete("1.0", "end");
    $nest_text->insert("end", 
       "        NX, NY\tSPACE\tPA RA\tLL IJ\tUR IJ\tPTS\n");
    $nest_text->tagAdd('color', "1.0","2.0");
    $nest_text->insert("end", "No sub nests exist.");
    $nest_text->configure(-state => 'disabled');

}

# -------------------------------------------
# access_to_nest_tabs
#
# The nest tabs need to be managed. Sometimes allowing
# access other times not.
# -------------------------------------------
sub access_to_nest_tabs {

   if ($globalMap or $projection_type eq 'RL' or $model_name eq "LAPS" ) { 
     $h_nb->pageconfigure(2, -state => 'disabled');
   } else {
     $h_nb->pageconfigure(2, -state => 'normal');
   }

}

# _______________ END NESTING ROUTINES _____________________________

# -----------------------------------
# bind_latlon_entries_init {
#
# -----------------------------------
sub bind_latlon_entries_init {
  my ($arg)=@ARG;
  
  if ($arg) {

    # Bind key events to update the bbox and its attributes.
    for my $event ("<Return>") {
       $lonEnt->Tk::bind ($event, [\&x_wrap_cylindrical_ll2xy, 1]);
       $latEnt->Tk::bind ($event, [\&x_wrap_cylindrical_ll2xy, 1]);
    }
    for my $event ("<Double-Return>") {
       $lonEnt->Tk::bind ($event, [\&x_wrap_render_new_map]);
       $latEnt->Tk::bind ($event, [\&x_wrap_render_new_map]);
    }
    for my $event ("<Leave>") {
       $lonEnt->Tk::bind ($event, [\&x_wrap_cylindrical_ll2xy, 1]);
       $latEnt->Tk::bind ($event, [\&x_wrap_cylindrical_ll2xy, 1]);
    }
    for my $event ("<ButtonRelease>") {
       $lon_up->Tk::bind ($event, [\&x_wrap_cylindrical_ll2xy, 1]);
       $lon_dw->Tk::bind ($event, [\&x_wrap_cylindrical_ll2xy, 1]);
       $lat_up->Tk::bind ($event, [\&x_wrap_cylindrical_ll2xy, 1]);
       $lat_dw->Tk::bind ($event, [\&x_wrap_cylindrical_ll2xy, 1]);
    }

  } else {

    #'Un-bind'.
    for my $event ("<Return>", "<Leave>") {
       $lonEnt->Tk::bind ($event, sub{ Tk::NoOp });
       $latEnt->Tk::bind ($event, sub{ Tk::NoOp });
    }
    for my $event ("<ButtonRelease>") {
       $lon_up->Tk::bind ($event, sub{ Tk::NoOp });
       $lon_dw->Tk::bind ($event, sub{ Tk::NoOp });
       $lat_up->Tk::bind ($event, sub{ Tk::NoOp });
       $lat_dw->Tk::bind ($event, sub{ Tk::NoOp });
    }
  }


}

# -----------------------------------
# bind_latlon_entries {
#
# -----------------------------------
sub bind_latlon_entries {
  my ($arg)=@ARG;
  
  if ($arg) {

    # Bind key events to update the bbox and its attributes.
    for my $event ("<Return>") {
       $clonEnt->Tk::bind ($event, [\&wrap_cylindrical_ll2xy, 1]);
       $clatEnt->Tk::bind ($event, [\&wrap_cylindrical_ll2xy, 1]);
       $dx_ent->Tk::bind  ($event, [\&grid_spacing_adjusted]);
       $nx_ent->Tk::bind  ($event, [\&dims_adjusted]);
       $ny_ent->Tk::bind  ($event, [\&dims_adjusted]);
       $latEnt->Tk::bind  ($event, [\&truelat_adjusted]);
       $lat2Ent->Tk::bind ($event, [\&truelat_adjusted])
    }
    for my $event ("<Double-Return>") {
       $clonEnt->Tk::bind ($event, [\&wrap_render_new_map]);
       $clatEnt->Tk::bind ($event, [\&wrap_render_new_map]);
       $dx_ent->Tk::bind  ($event, [\&wrap_render_new_map]);
       $nx_ent->Tk::bind  ($event, [\&wrap_render_new_map]);
       $ny_ent->Tk::bind  ($event, [\&wrap_render_new_map]);
       $latEnt->Tk::bind  ($event, [\&wrap_render_new_map]);
       $lat2Ent->Tk::bind ($event, [\&wrap_render_new_map]);
       $lonEnt->Tk::bind  ($event, [\&wrap_render_new_map]);
    }
    for my $event ("<Leave>") {
       $clonEnt->Tk::bind ($event, [\&wrap_cylindrical_ll2xy, 1]);
       $clatEnt->Tk::bind ($event, [\&wrap_cylindrical_ll2xy, 1]);
       $dx_ent->Tk::bind  ($event, [\&grid_spacing_adjusted]);
       $nx_ent->Tk::bind  ($event, [\&dims_adjusted]);
       $ny_ent->Tk::bind  ($event, [\&dims_adjusted]);
       $latEnt->Tk::bind  ($event, [\&truelat_adjusted]);
       $lat2Ent->Tk::bind ($event, [\&truelat_adjusted]);
    }
    for my $event ("<ButtonRelease>") {
       $clon_up->Tk::bind ($event, [\&wrap_cylindrical_ll2xy, 1]);
       $clon_dw->Tk::bind ($event, [\&wrap_cylindrical_ll2xy, 1]);
       $clat_up->Tk::bind ($event, [\&wrap_cylindrical_ll2xy, 1]);
       $clat_dw->Tk::bind ($event, [\&wrap_cylindrical_ll2xy, 1]);
       $dx_up->Tk::bind   ($event, [\&grid_spacing_adjusted]);
       $dx_dw->Tk::bind   ($event, [\&grid_spacing_adjusted]);
       $nx_up->Tk::bind   ($event, [\&adjust_bbox]);
       $nx_dw->Tk::bind   ($event, [\&adjust_bbox]);
       $ny_up->Tk::bind   ($event, [\&adjust_bbox]);
       $ny_dw->Tk::bind   ($event, [\&adjust_bbox]);
       $lat_up->Tk::bind  ($event, [\&truelat_adjusted]);
       $lat_dw->Tk::bind  ($event, [\&truelat_adjusted]);
       $lat2_up->Tk::bind ($event, [\&truelat_adjusted]);
       $lat2_dw->Tk::bind ($event, [\&truelat_adjusted]);
    }

  } else {

    #'Un-bind'.
    for my $event ("<Return>", "<Leave>") {
       #$clonEnt->Tk::bind ($event, [""]);
       $clonEnt->Tk::bind ($event, sub{ Tk::NoOp });
       $clatEnt->Tk::bind ($event, sub{ Tk::NoOp });
    }
    for my $event ("<ButtonRelease>") {
       $clon_up->Tk::bind ($event, sub{ Tk::NoOp });
       $clon_dw->Tk::bind ($event, sub{ Tk::NoOp });
       $clat_up->Tk::bind ($event, sub{ Tk::NoOp });
       $clat_dw->Tk::bind ($event, sub{ Tk::NoOp });
    }
  }


}

# -----------------------------------
# reinstate_gen_map_vars
#
# If the variables to generate a map change, 
# then let user reset them if desired.
# -----------------------------------
sub reinstate_gen_map_vars{

     # Reset all the app variables from model namelist. 
     assign_nl_vars(1);
     # (Set_nx and set_dx are called from adjust_bbox just below here.)
  
     # Reset display variables.
     $pix_per_gpoint=$pix_per_gpoint / $dx_orig * $grid_spacing_dkm;
     adjust_bbox();
     update_tags();

     # Un-highlight update button.
     highlight_update_button(0);

     # Enable tabs.
     if ($domain_mode == 2) {access_to_tabs(1);}
}

# -----------------------------------
# nyi
# -----------------------------------
sub global_menu_list {
    use vars qw($global_choice @global_image);

    @global_image = ('Global Cylindrical',
                     'North Hem Polar Stereo',
                     'South Hem Polar Stereo');
    my $rrow=1;
    my $nom_lab=$bcd->Label(-text => "Source image:", 
                            -fg => $colorN) 
                    ->grid( -row => $rrow, 
                            -column => 1, 
                            -sticky => 'w');

    my $nom = native_optionmenu(
        $bcd,
        \$global_choice,
        [sub {
                  switch_worlds();
             }, 
        ],
        @global_image,
    );
    $nom->grid(-row => $rrow, 
               -column => 3, 
               -columnspan => 3,
               -sticky => 'w');

    my $menu = $nom->cget(-menu);
    for my $i (0 .. scalar(@global_image)-1 ) {
        $menu->entryconfigure($i, -columnbreak => 1) unless $i % 7;
    }


    return($nom);
}

sub bcd_menu_list {
    use vars qw($bcd_choice);

    my $rrow=2;
    my $nom_lab=$bcd->Label(-text => "User selected Map File:", 
                            -fg => $colorN) 
                    ->grid( -row => $rrow, 
                            -column => 1, 
                            -sticky => 'w');

    my $nom = native_optionmenu(
        $bcd,
        \$bcd_choice,
        [sub {
               update_bcd_list();
           }, 
        ],
        @bcd_found,
    );
    $nom->grid(-row => $rrow, 
               -column => 3, 
               -columnspan => 3,
               -sticky => 'w');

    my $menu = $nom->cget(-menu);
    for my $i (0 .. scalar(@bcd_found)-1 ) {
        $menu->entryconfigure($i, -columnbreak => 1) unless $i % 7;
    }


    return($nom);
}

sub native_optionmenu {
    my($parent, $varref, $command, @optionvals)=@ARG;

    $$varref = $optionvals[0];

    my $mb = $parent->Menubutton(
        -textvariable       => $varref,
        -indicatoron        => 1,
        -width              => 20,
        -relief             => 'raised',
        -borderwidth        => 2,
        -highlightthickness => 2,
        -anchor             => 'c',
        -direction          => 'flush',
        -justify            => 'right',
    );
    my $menu = $mb->Menu(-tearoff => 0);
    $mb->configure(-menu => $menu);

    my $callback = ref($command) =~ /CODE/ ? [$command] : $command;

    foreach (@optionvals) {
	$menu->radiobutton(
            -label     => $ARG,
            -variable  => $varref,
            -command   => [@$callback, $ARG],
            -selectcolor => 'red',
            -indicatoron => 0,
        );
    }

   $mb;

} # end native_optionmenu

# ---------------------------------------
# update_bcd_list 
#
# Once the user selects a Map File, 1) expand variable by adding
# the full path name and 2) deselect all but the current choice.
# ---------------------------------------
sub update_bcd_list {

   my $tmp=join ('', $GUI_MAP, '/', $bcd_choice, '.bcd');
   $bcdFileName=eval('$tmp'); 
   $bcdFileName=$tmp; 

   for my $i (0 .. scalar(@bcd_found)-1 ) {
      $bcd_mb->entryconfigure($i, -bg => $bg_color);
   }
   $bcd_mb->entryconfigure($bcd_choice, -bg => $update_color);

   # Update needs to be pressed, a parameter has changed.
   if(Exists($b_update) ) {highlight_update_button(1);}
}
 

# ---------------------------------------
# adjust_scrollregion 
# 
# Set the scrollbars length and width to be the size of the
# current image and the scrollbar x & y sliders accordingly.
# ---------------------------------------
sub adjust_scrollregion {
    $can->configure(-scrollregion => [ 0,0,$ARG[0],$ARG[1] ] );
  # print "sub adjust_scrollregion \n";
  # $can->configure(-scrollregion => [($dmap_image->width/2),0, $ce_dmap_image->width, $ce_dmap_image->height] );
}

# -----------------------------------
# present_grid_editor 
#
# Select Create new domain, load existing domain, or update domain
#     1) Clear  -- set grid editor background to a global map image, or
#     2) Reload -- set grid editor to empty canvas in order to draw a map.
#     0) Update -- set grid editor to empty canvas in order to update a map.
# -----------------------------------
sub present_grid_editor {
  my ($arg)=@ARG;

  # Clean the canvas.
  $can->delete("all"); 

  # Un-highlight update button.
  highlight_update_button(0);


  if ($arg == 1) { 

    # Reset map and vars as 'Clear' has been pressed.
    # ------------------------
   
    # New Domain.
    $map_setup_init=1;
    $firstTime=1;
    $globalMap=1;
    set_bindings_to_draw(1);
    #bind_latlon_entries(1);
    bind_latlon_entries_init(1);

    #$c_photo->read($dmap_image, -shrink);
    #$can->delete('canv_map_im');
    $c_photo->blank;
    $c_photo->copy($dmap_image, -shrink);
    $can->createImage(0, 0, 
                      -anchor => 'nw',
                      -image => $c_photo,
                      -tags => "canv_map_im");

    # Reset zoom factor
    $zoom_factor=0;

    # Enable zoom mechanism.
    set_zoom_icon_state(1);
    $img_xmax=$c_photo->width();
    $img_ymax=$c_photo->height();
    adjust_scrollregion($img_xmax,$img_ymax);

    # Default canvas scrollbar locations.
    $can->xviewMoveto(0.004);
    $can->yviewMoveto(0.005);

    # Calculate the center lat,lon point
    $calc_centerLatLon=1;
    $grid_val_restrict=1;

    # Activate gridSpacing.
    set_gridSpacing(1);

    # Hide grid editor widgets (though these are already 
    # hidden when application is initially started).
    present_grid_editor_widgets();

  } elsif ($arg == 2) { 

    # Reload bbox because either 'Start Over' has been pressed.
    # or, the user has loaded an existing domain.
    # ------------------------

    # Existing Domain.
    $map_setup_init=1;
    $firstTime=1;
    $globalMap=1;
    set_bindings_to_draw(0);
    bind_latlon_entries_init(0);
    bind_latlon_entries(1);

    $can->createImage(0, 0, 
                      -anchor => 'nw',
                      -image => $c_photo,
                      -tags => "canv_map_im");

    # Enable zoom mechanism.
    set_zoom_icon_state(2);
    $img_xmax=$c_photo->width();
    $img_ymax=$c_photo->height();
    adjust_scrollregion($img_xmax,$img_ymax);


    # Activate gridSpacing.
    set_gridSpacing(1);

    # Update scrollbar sliders in order to view the centerpoint icon
    if ($z_factor < 0) { $z_factor=abs (1/$z_factor);}
    cylindrical_ll2xy(0, $z_factor);

    # Calculate the center lat,lon point
    $calc_centerLatLon=1;
    #set_button_list_state(1,@cenlon_wdgt_list);
    #set_button_list_state(1,@cenlat_wdgt_list);
    $grid_val_restrict=1;

    # Hide grid editor widgets.
    present_grid_editor_widgets();


  } elsif ($arg == 0) { 

    # Create map and set vars as 'Update' has been pressed.
    # ------------------------

    if ($prev_editor_arg=$arg) { return;}

    # Existing Domain.
    $globalMap=0;
    $can->createText(50, 50, 
                -anchor => 'w',
                -text => "...creating and loading domain file.",
                -fill => $colorW,
                -tags => "canv_msg", 
                (defined $error_font ? (-font => $error_font) : ())
                );

    # Disable zoom mechanism.
    set_zoom_icon_state(0);
    $can->configure(-scrollregion => [ 0,0,0,1]);

    # Take focus away from Centerpoint lat,lon.
    # (If user were to press 'Return' boundingbox would update.)
    switch_focus_main();

    # Do NOT Calculate the center lat,lon point.
    $calc_centerLatLon=0;
    $grid_val_restrict=1;

    # Disable gridSpacing.
    set_gridSpacing(0);

    # Show grid editor widgets.
    present_grid_editor_widgets();

    #'Un-bind' action associated with entries.
    bind_latlon_entries_init(0);
    bind_latlon_entries(0);
  }

  $prev_editor_arg=$arg;
  set_latlon_wdgt_state();
  access_to_nest_tabs();

}

# -------
# present_grid_editor_widgets
#
# Show and hide the grid editor depending on editing mode.
# -------
sub present_grid_editor_widgets {
 
  if ($globalMap == 1) {
     # Hide grid editor widgets.
     $mode->packForget();
     $h_grid_hide->pack(@pInfo_h_grid, -after => $proj, -fill => 'y');
     $h_grid->packForget();
     #set_button_state(0,$linetype_mb);

  } else {
     # Show grid editor widgets.
     $mode->pack(@pInfo_mode, -before => $action);
     $h_grid->pack(@pInfo_h_grid, -after => $proj);
     $h_grid_hide->packForget();
     #set_button_state(1,$linetype_mb);
  }

}
# -----------------------------------
# zoom_in_out
#
# NOTE: 
# -zoom interpolates and adds pixels
# -subsample subtracts pixels without interpolation
# -----------------------------------
sub zoom_in_out {
    my ($incre)=@ARG;

    watch_cursor(1);
    
    # Clean the canvas slate every time.
    $can->delete('canv_map_im');
    $c_photo->blank;

    # Maintain these values.
    $zoom_factor=$zoom_factor+($incre); 
    my $scale=abs($zoom_factor);

    # Activate zoom icons.
    set_zoom_icon_state(1);

    # Size the image.
    if ($zoom_factor > 0) {
       # Increase map background image (by replicating pixels).
       $c_photo->copy($dmap_image, -zoom => $scale, $scale , -shrink);
       ##$c_photo->copy($dmap_image_double, -shrink);

       # Limit maximum size (by disabling zoomIn button).
       if ($zoom_factor >= 4) { 
           $b_zoomIn->configure(-state => 'disabled'); 
           ##$c_photo->copy($dmap_image_doubledouble, -shrink);
       }

    } elsif ($zoom_factor < 0) {
       # Decrease map background image (by sub-sampling pixels).
       $c_photo->copy($dmap_image, -subsample => $scale, $scale, -shrink);
       ##$c_photo->copy($dmap_image_half, -shrink);

       # Limit minimum size (by disabling zoomOut button).
       if ($zoom_factor <= -2) { $b_zoomOut->configure(-state => 'disabled'); }

    } elsif ($zoom_factor == 0) {
       # Load default map background.
       $c_photo->copy($dmap_image, -shrink);
    }

    # Update image.
    $can->createImage(0, 0, 
                      -anchor => 'nw',
                      -image => $c_photo,
                      -tags => "canv_map_im");
    $can->lower('canv_map_im');
    $img_xmax=$c_photo->width();
    $img_ymax=$c_photo->height();
    adjust_scrollregion($img_xmax,$img_ymax);

    if (!($latSW eq "" || $lonSW  eq "" || $latNE  eq "" || $lonNE eq "" )){ 
      # Update location of center point icon accordingly.
      if ($incre < 0) { $incre=abs (1/$incre);}
      cylindrical_ll2xy(0, $incre);
    }
    
    watch_cursor(0);
}


# ---------------------------------------
# set_zoom_icon_state
#
# Set the zoom icons either disabled or normal
# based on present_grid_editor values and when 
# user selects Clear or Start Over.
# ---------------------------------------
sub set_zoom_icon_state {
    my ($state)=@ARG;

    if ($zoomForgotten) {
       # Restore buttons that apply.
       $b_info->packForget();
       $b_zoomIn->pack(-side => 'left');
       $b_zoomOut->pack(-side => 'left', -padx => 1);
       $b_info->pack(-side => 'left', -anchor => 'n');
    }

    if ($state == 0) {
       # Set buttons to disable state but save previous state.
       $b_zoomInState=$b_zoomIn->cget(-state);
       $b_zoomOutState=$b_zoomOut->cget(-state);

       # Remove buttons that do not apply.
       $zoomForgotten=1;
       $b_zoomIn->packForget();
       $b_zoomOut->packForget();

    } elsif ($state == 1) {
       # Set buttons to normal state.
       $zoomForgotten=0;
       $b_zoomIn->configure( -state => 'normal');
       $b_zoomOut->configure(-state => 'normal');

       $b_zoomInState='normal';
       $b_zoomOutState='normal';

    } else {
       # Set buttons to their previous state.
       $zoomForgotten=0;
       $b_zoomIn->configure( -state => $b_zoomInState);
       $b_zoomOut->configure(-state => $b_zoomOutState);
    }
}

sub set_gridSpacing {
  my ($arg)=@ARG;

  $calc_gridSpacing=$arg;
  restrict_grid_var_calc();
}

# ----------------------------------------------------------
# restrict_grid_var_calc: 
#
# Grid choices, one of the two need to be held constant
# allowing the user to do fine-tune adjustments with the 
# bbox OR the entry widgets.
# ----------------------------------------------------------
sub restrict_grid_var_calc {
  use vars qw($grid_val_restrict);
   
  if ($grid_val_restrict == 1) { # Edit using the bbox.

      # Activate bounding box. 
      $calc_boundingBox=1;
      if (!$globalMap) { update_tags(); }

      # Disable entry box vars. 
      foreach $wdgt (@hori_wdgt_list, @lon_wdgt_list) {
         $wdgt->configure(-state => 'disabled', -fg => $disabled_color); }

      # Switch focus to the mainwindow, away from any Entrybox.
      switch_focus_main();
       

  } elsif ($grid_val_restrict == 0) { # Edit by typing in parameters values.

      # Disable bounding box. 
      $calc_boundingBox=0;
      hide_tags();
      hide_tags_nest();

      # Activate (refresh) entry box vars.
      foreach $wdgt (@hori_wdgt_list) {
         $wdgt->configure(-state => 'normal', -fg => $normal_color); }


      # Sync bbox determined values of nx,xy and _orig values.
      $nx_orig=$nx_dim;
      $ny_orig=$ny_dim;

      # Allow user to drag & display lat/lon values.
      set_bindings_to_draw(0);
  }

  present_grid_editor_widgets();

  # Set entries as editable, or not.
  set_lat2_lon_wdgt_state();
 
}

# ---------------------------------------
# present_nest_widgets
#
# Edit nests by showing/hiding the nest ui.
# ---------------------------------------
sub present_nest_widgets {
   ($nest_panel)=@_;

   if ($nest_panel) {
      # Show 'Nest Domain' panel.

      # If user changes map, press 'Update Map'.
      if($b_update_highlight) {
            ## The map needs to be updated.
            my $ans=okcancel_dbox("Press 'Update Map'", 
"The domain map parameters have changed. Press 'Update Map' to continue.
By pressing 'Ok' the 'Update Map' button will be invoked automatically.");

            # Update the map, after calling sync_nl_vars.
            if($ans eq 'Ok'){ 
		wrap_render_new_map(); 
            } else {
		$nest_panel=0;
		$h_nb->raise(1);
		return;
            }
      }

      # Disable bounding box and add a label for it.
      $calc_boundingBox=0;
      hide_tags();
      render_bbox_label();

      # Display all nests.
      if ($num_nests > 1) { 
         wrap_draw_all_nests(); 
         set_button_list_state(1,@nest_actions);
      } else {
         set_button_list_state(0,@nest_widgets_list); 
      }

  } else {
      # Show 'MOAD Domain' panel.
 
      if ($domain_select eq "") { return; }

      # Reset fine scale editing mode.
      restrict_grid_var_calc();

      # Remove yellow colored objects.
      remove_nestbox();
      remove_parentbox();
      remove_nest_label();
      if ($num_nests > 1) { 
         #sync_nl_vars_nest();
         reload_nests(0); 
      }
  }

  # Switch focus to the mainwindow, away from any Entrybox.
  switch_focus_main();
}

# ---------------------------------------
# set_lat2_lon_wdgt_state
# ---------------------------------------
sub set_lat2_lon_wdgt_state{
  use vars qw($grid_val_restrict);

  # Disable entry box vars. 
  set_button_list_state(0,@lat2_wdgt_list);
  if ( ($projection_type eq 'LC' or $projection_type eq 'RL') && $grid_val_restrict == 0) {
     # Enable entry box var for Lambert, only. 
     set_button_list_state(1,@lat2_wdgt_list);
  }

  if ($grid_val_restrict == 0) {
     if ($projection_type eq 'ME') {
        # Disable entry box var for Mercator, only. 
        # Changing this value has no real effect on ME projections,
        # so don't allow the user to edit it.
        set_button_list_state(0,@lon_wdgt_list);
     } else {
        # Enable entry box vars. 
        set_button_list_state(1,@lon_wdgt_list);
     }
  }

}

# ---------------------------------------
# set_latlon_wdgt_state
#
# The behavior of the domain box when moving or resizing depends on the map 
# projection.  
#  
# For a Mercator projection, the box can be moved to any location or resized 
# in any dimension.  
#
# For a Lambert projection, the box cannot be moved and can only be resized 
# around the central point.  This means that dragging one side of the domain 
# box to resize it results in the opposite side moving in the opposing direction 
# at the same time.  
#
# For a Polar Stereographic projection, the box can be moved up and down along 
# the central longitude between the equator and pole but cannot be moved left 
# or right.  When resizing a domain box on a Polar Stereographic projection, 
# the northern and southern sides move independently but the eastern and western 
# sides move in concert to keep the central longitude constant. 
# ---------------------------------------
sub set_latlon_wdgt_state{

return;

  if (!$globalMap) {

    if (0) {
     if ($projection_type eq 'PS') {
        set_button_list_state(0,@cenlat_wdgt_list);
        set_button_list_state(1,@cenlon_wdgt_list);
     } elsif ($projection_type eq 'LC' or $projection_type eq 'RL') {
        set_button_list_state(0,@cenlat_wdgt_list);
        set_button_list_state(0,@cenlon_wdgt_list);
     } elsif ($projection_type eq 'ME') {
        set_button_list_state(1,@cenlat_wdgt_list);
        set_button_list_state(1,@cenlon_wdgt_list);
     }
    }

    # Until code is adapted for discussion above ...turn off widgets.
    set_button_list_state(0,@cenlat_wdgt_list);
    set_button_list_state(0,@cenlon_wdgt_list);

  }

}


# ---------------------------------------
# render_map_clear
# 
# Display global map and set the scrollbar sliders 
# to 0 when user presses "Clear".
# ---------------------------------------
sub render_map_clear {

   # Set grid editor to a globalMap and 
   # suppress extra grid editing widgets.
   if($domain_mode == 1 && $globalMap) {
      $can->delete("bbox", "rtags"); 
      set_bindings_to_draw(1);
   } else {
      switch_domains($domain_mode);
   }

   # Activate gridSpacing.
   set_gridSpacing(1);

   reset_vars();
}

# ---------------------------------------
# reset_vars
# 
# ---------------------------------------
sub reset_vars {

   # Reset vars.
   $startover_vars_init=0;
   $map_setup_init=0;
   $nx_dim=$ny_dim=$nx_orig=$ny_orig=$NXmax;
   $grid_spacing_dkm=$dx_orig=0;
   $grid_spacing_m_cmn=0;
   $polar_cap=0;
   $prev_projection_type="";
   $reset_projection_type="";
   #$d_num="d02";
  
   watch_cursor(0);

   # These are set to null after a clear. An error dialog
   # box will appear if these values are not set.
   $latSW=$lonSW=$latNE=$lonNE="";

   # Set Projection Type: None and Map file: US_Counties
   $proj_label="None";
   $bcd_choice=$bcd_prev=$bcdFileDefault;
   update_bcd_list();
   $maproj_mb_lab->configure(-text => "Map Projection:"); 

   # Un-highlight buttons.
   highlight_update_button(0);
   set_button_highlight(0,$maproj_mb);

   # Don't let user "Start Over" a cleared canvas.
   set_button_state(0,$b_startover);

   # Set the vertical sigma scheme to current.
   $scheme_choice=$vscheme[0];

   # Reset HGrid panel to MOAD Domain, not Nest Doamin.
   reset_nest_values();
   if($nest_panel) {$h_nb->raise(1);}
}

# -----------------------------------
# render_old_map
#
# Display global map when user presses "Start Over".
# -----------------------------------
sub render_old_map {

 if ($domain_mode == 1) { 
   watch_cursor(1);

   # Set grid editor to a globalMap and 
   # suppress extra grid editing widgets.
   present_grid_editor(2);

   # Reset vars
   $startover_vars_init=0;
   $map_setup_init=0;
   $nx_dim=$ny_dim=$nx_ddim=$NXmax;
   $proj_label="None";

   # Recreate bounding box and tags
   render_old_bbox();

   # ...and tags AND lat,lon.
   $perform_calc=1;
   my ($LL, $TT, $RR, $BB) = $can->coords("bbox");
   create_tags($LL, $TT, $RR, $BB);

   # Set truelat1 and truelat2 values.
   set_true_lats();

   bbox_corners_etc($LL,$TT,$RR,$BB);
   $polar_cap=0;

   # No restrictions.
   restrict_grid_var_calc();
   suggest_proj_type();
   $can->itemconfigure('c_img', -image => $centerPoint);
   center_bbox();
  
   # Highlight update button.
   highlight_update_button(1);

   # Activate cen lat/lon buttons.
   #set_button_list_state(1,@cenlon_wdgt_list);
   #set_button_list_state(1,@cenlat_wdgt_list);
   set_button_list_state(1,@lon_wdgt_list);
   set_button_list_state(1,@lat_wdgt_list);

   #'Bind' action associated with entries.
   bind_latlon_entries_init(1);
   #bind_latlon_entries(1);

   # Disable tabs.
   access_to_tabs(0);

   watch_cursor(0);

 } elsif ($domain_mode == 2) { 

   # Read template namelist and 'merge' with install namelist.
   $panel_index=1; # Must be = 1.
   load_namelist();
   $panel_index=2;
 
   # Make a bounding box from namelist (without dynamic user input).
   make_bbox_from_nlist();

   # Reset fine scale editing mode.
   $grid_val_restrict=1;
   restrict_grid_var_calc();

   # Enable tabs.
   access_to_tabs(1); 
 }
}

# -----------------------------------
# check_domain_args
#
# Check the user input values, esp for domain related values.
# -----------------------------------
sub check_domain_args {

    # Error_status=0, no errors.
    my $error_status=0;

    # Projection label needs to be set.
    # --------------------------------
    if($proj_label eq "None") {
        # map projection 'None' selected

        $error_status=1;
        my $msg="The map projection type needs to be selected.\nWould you like the suggested type: $proj_choices[$sug_proj_idx][0]?";
        my $disp = yesno_dbox("Value Not Selected", $msg, "Yes"); 
        if($disp eq 'Yes') {
          set_proj_vars($sug_proj_idx);
          $error_status=0;
        } else {
          $error_status=1;
          return ($error_status);   
        }

    }

    # Calculate grid value to test it against threshold value.
    # --------------------------------
    my $grid_size=int(($grid_spacing_km)*sqrt($nx_dim*$nx_dim+$ny_dim*$ny_dim)); 

    if($bcd_choice eq $bcdUSCounties || $bcd_choice eq 'US_RiverDrainBas'  ) {

        if($grid_cen_lat_cmn >=  50 || $grid_cen_lat_cmn <= 20 ||
            $grid_cen_lon_cmn >= -50) {
           # Cent lat,lon is not over the US.
           my $disp = yesno_dbox("Confirm Choice", "The selection of the $bcd_choice map background file \nis not applicable. Would you like to switch to $bcdFileDefault.", "Yes");
           if($disp eq 'Yes') { 
              $bcd_choice=$bcdFileDefault;
              update_bcd_list();
           } else {
              $error_status=1;
              return ($error_status);   
           }
        }

        my  $g_threshold=5000;  # km
        if($grid_size >= $g_threshold) {
           # Warn user that domain is very large to display US county info.
           my $disp = okcancel_dbox("Confirm Choice", "The selection of the $bcd_choice map background file and\nthe large domain size will take a little extra time to process.");
           if($disp ne 'Ok') { 
               $error_status=1;
               return ($error_status);   
           }
        }
    }

    # Get bounding box coords.
    # --------------------------------
    if ($domain_mode == 1) { ($l, $t, $r, $b) = $can->coords("bbox"); }

    # Recalculate PS with polar_cap flag set, if appropriate.
    # --------------------------------
    if ($domain_mode == 1 && $projection_type eq 'PS') { 
       $polar_cap=0;
       if ( ($latNE >= 89 || $latSW <= -89) && abs ($grid_cen_lat_cmn) > 66.5) {
          $polar_cap=1;
          bbox_corners_etc($l,$t,$r,$b);
          set_proj_vars(0);
       }
    }

 
    # Lambert cannot accept a 0 latitude value.
    # --------------------------------
    if( ($projection_type eq 'LC' or $projection_type eq 'RL') && $grid_cen_lat_cmn <3 && $grid_cen_lat_cmn >-3) {

        $error_status=1;
        my $msg="The center latitude for the Lambert projection cannot be so close to zero.\nMake the latitude value at least 3 degrees North or South.";
        my $disp = info_dbox("Center Latitude Error", $msg); 
    }

    # Test estimated grid_spacing vals.
    # --------------------------------
    if ($domain_mode == 1) { bbox_corners_etc($l,$t,$r,$b); }


    # Mercator must have a 0 true latitude2 value, esp for supmap.
    # --------------------------------
    if ($projection_type eq 'ME') {
       $truelat2=0;
    }

    return ($error_status);   
}


# -----------------------------------
# suggest_proj_type
#
# Estimate the best guess for projection information
# based on users Center Point selection.
# -----------------------------------
sub suggest_proj_type {

    if($domain_mode != 1) { return; };
 
    # Make suggestion.
    # ------------------- 
    if (abs ($grid_cen_lat_cmn) <= 23.5) {
       # If selected center pt is near equator choose: Mercator 
       $sug_proj_idx=2; 

    } elsif (abs ($grid_cen_lat_cmn) > 70) {
       # If selected center pt poleward of 70 choose: Polar stereographic 
       $sug_proj_idx=0; 
 
    } else {
       # Recommend: Lambert (but any will work)
       $sug_proj_idx=1; 
    }
    $maproj_mb_lab->configure(-text => 
                "Map Projection: (suggest $proj_choices[$sug_proj_idx][2])"); 

   
    # Enforce suggestion.
    # ------------------- 
    my $ii=0;
    until ($ii>2) {
        $maproj_mb->entryconfigure($ii, -state => 'normal');
        $ii++;
    }

    if (abs ($grid_cen_lat_cmn) >= 60) {
       # Disable ME
       $maproj_mb->entryconfigure(2, -state => 'disabled');
       # Choose None.
       if ($projection_type eq 'ME') { $proj_label="None"; }
    }

    if ($latNE >= 89 or $latSW <= -89 or
        abs ($grid_cen_lat_cmn) >= 90) {
       # Disable ME & LC.
       $maproj_mb->entryconfigure(1, -state => 'disabled');
       $maproj_mb->entryconfigure(2, -state => 'disabled');
       # Choose PS.
       set_proj_vars(0);
    }

}

# -------------------------------------------
# wrap_render_new_map
# 
# -------------------------------------------
sub wrap_render_new_map {

   # For bind...
   #cylindrical_ll2xy(0, 1);


   # If nest is a parent do not proceed without confirmation,
   # as changing the parent effects the depenpants.
   if (!$globalMap) {
     if ($b_update_highlight && 
       $num_nests > 1 &&
       ask_delete_dependants() eq 'Cancel'){ 
         # Restore previous values and return.
         reinstate_gen_map_vars();
         return(1);
     }
   } else {
     # Delete all nest values.
     reset_nest_values(); 
   }
   
  # sync_nl_vars();
   if (check_gen_map_vars(0) || $globalMap) {
       # Disable 'Next' button, render map, then enable 'Next' button.
       set_button_state(0,$cntl_next);
 
       render_new_map();
       set_button_state(1,$cntl_next);
       return(0);

   } else {
       highlight_update_button(0);
       return(1);
   }

   # Again, suggest projection type.
   suggest_proj_type();
}

# -------------------------------------------
# render_new_map
#
# Graphically display user selected domain 
# by drawing map vector lines on a canvas 
# when user presses "Update Map".
# -------------------------------------------
sub render_new_map {

   if ($domain_mode == 1 && $mw->cget(-cursor) eq 'watch') {
      # Busy processing something. 
      return;
   } 

   watch_cursor(1);

   # Make sure user input values are within the bounds
   # --------------------------------
   if (check_domain_args() ) {
       watch_cursor(0);
       return(1); 
   }

   # Store zoom factor and bounding box location for use in 
   # subroutine render_old_map -- if user chooses to "Start Over".
   # --------------------------------
   unless ($startover_vars_init) { 
       $startover_vars_init=1;
       $z_factor=$zoom_factor;

       # The bbox needs to be saved.
       my @box_sides=$can->coords("bbox");
       ($bbox_L, $bbox_T, $bbox_R, $bbox_B) = @box_sides;
 
       # The bbox needs to be scaled to normal size then saved.
       if ($zoom_factor > 0) { foreach (@box_sides) { $_=$_/$zoom_factor; }
       } elsif ($zoom_factor < 0) { foreach (@box_sides) { $_=$_*$zoom_factor; }
       }
       ($orig_bbox_L, $orig_bbox_T, $orig_bbox_R, $orig_bbox_B) = @box_sides;
   }

   watch_cursor(2);

   # Create blank canvas window.
   # --------------------------------
   if ($globalMap) { 
       # Create blank canvas.
       present_grid_editor(0); 
   } else {
       # Clear canvas contents.
       #$can->delete("all"); 
   }


   # Set bbox to fit within the canvas.
   # Maintain bbox's original ratio.
   # --------------------------------
   $mw->update();
   $mw->idletasks();
   my $c_width =($can->width() - 10); 
   my $c_height=($can->height() - 10);
   if ($c_width <= 1 ) {
     $c_width =375 - 15;
     $c_height=500 - 15;
   }

   if($c_width < $c_height ) {
      $c_xmax=$c_ymax=$c_width - 0;
   } else {
      $c_xmax=$c_ymax=$c_height - 0; 
   }

   # Set map_utils hash.
   # --------------------------------
   map_utils_setup_d01(); 

   # Create arg list.
   # --------------------------------
   set_rot_lat_args();
   $current_args = "$ncarg_idx $lat2_arg $stdlon $rotation_arg $latNE $lonNE $latSW $lonSW $bcdFileName $mapVectorFile $setsupData";

   # Create new map background file /tmp/vector_instructions.tk.
   # --------------------------------

   # Remove previous map background file
   #if (-e $mapVectorFile) { @result = `rm -rf $mapVectorFile`; }
   if(-e $mapVectorFile){unlink $mapVectorFile or 
          warn("Can't remove $mapVectorFile: $!");}

   # Create map background file
   print $current_args if $Debug;
   watch_cursor(2);

   my $my_catch=Tk::catch { `$gen_map_exe $current_args`; };
   if ($my_catch eq "") { 
      # There is a major problem with the executable file.
      fail_dbox("Executable Problem", 
"Error with map background executable command: 
$gen_map_exe 
\t$current_args 
Possible 'Exec format error'.
\n\tClear will be invoked!"); 
      run_sys::run_sys("$gen_map_exe $current_args",1);
      render_map_clear();
   }

   # Read map background file /tmp/vector_instructions.tk which contain the 
   #        sub vector_instructions AND 
   #        sub ll_instructions.
   # --------------------------------
   if(-e $mapVectorFile && !-z _) {
       chmod 0666, $mapVectorFile;
       chmod 0666, $setsupData;

       do $mapVectorFile;
       #$can->delete("canv_msg"); 
       $can->delete("all"); 

       # Draw lat,lon location vectors.
       ll_instructions($can,$c_xmax,$c_ymax);

       # Draw map vectors.
       vector_instructions($can,$c_xmax,$c_ymax);

   }

   watch_cursor(3);

   # Clean up stuff.
   # --------------------------------
   my ($LL, $TT, $RR, $BB) = $can->coords("bbox");
   create_tags($LL, $TT, $RR, $BB);

# ---

   # Create new box in canvas which is 50 pixels bigger than bbox.
   ($LL, $TT, $RR, $BB) = $can->coords("box_line");
   my $b_xcen_new=(($RR+$LL)/2);
   my $b_ycen_new=(($BB+$TT)/2);

   # Calculate offset.
   my $b_xcen =($c_width/2);
   my $b_ycen =($c_height/2);
   $b_xoff=($b_xcen - $b_xcen_new);
   $b_yoff=($b_ycen - $b_ycen_new);

   # Center bbox in canvas.
   $can->move("all", $b_xoff, $b_yoff);

   # Create tags but Do Not re-compute center Lat Lon
   update_tags();
   ($LL, $TT, $RR, $BB) = $can->coords("bbox");
   $pix_per_gpoint=($RR-$LL)/$nx_dim;

   # Show bbox.
   $can->itemconfigure('c_img', -image => $centerPoint_dot);
   $can->delete('fake_area');

# ---

   if (!$grid_val_restrict) { 
       # Keep 'bbox' yellow after its been redrawn.
       restrict_grid_var_calc(); 
   } elsif ($firstTime && $domain_mode == 1) {
       # Update the dims.
       bbox_corners_etc($LL,$TT,$RR,$BB);

       # Set map_utils hash, again.
       map_utils_setup_d01(); 
   } 

   watch_cursor(2);

   # Store current_args; if only the window size changes 
   # then it is not necessary to recalculate the map.
   $previous_args = $current_args;
   $firstTime=0; 

   watch_cursor(0);

   # Allow user to select "Start Over".
   my $orig_bbox="$ROOT_TEMPLATES/$domain_select/.orig_bbox.txt"; 
   if ($domain_mode == 1 || $orig_bbox_exists) {set_button_state(1,$b_startover);}

 
   # Un-highlight the update button.
   highlight_update_button(0);

   # Draw nests.
   if ($num_nests > 1) { draw_all_nests(); }

   # Define reset_projection_type.
   $reset_projection_type=$projection_type;

   # Success.
   return(0);
}

# -------------------------------------------------
# map_utils_setup_d01 
#
# See "$ROOT_INSTALL/etc/map_utils.pm";
#
# This is a Perl module for converting (i,j) to (lat,lon) and vice versa
# for any of 4 map projections:  Cylindrical Equidistant (lat-lon),
# Mercator, lambert conformal (secant and tangent), and polar
# stereographic.
#
# -------------------------------------------------
sub map_utils_setup_d01 {

   $known_ri = ($nx_dim+1.)*0.5;   # Exact center of grid
   $known_rj = ($ny_dim+1.)*0.5;   # Exact center of grid
   my $d_x=$grid_spacing_km * 1000.;

if ($Debug) {
print "
			       $projection_type,
                               $grid_cen_lat_cmn,$grid_cen_lon_cmn,
                               $known_ri, $known_rj, $d_x, 
                               $stdlon,
                               $truelat1,$truelat2,
                               $nx_dim,$ny_dim\n\n
                               $nx_ddim,$grid_spacing_dkm\n\n";
}
   
   @ArrayOfHashes= ();
   my $my_catch=Tk::catch {%{$ArrayOfHashes[0]} = 
	   &map_utils::map_set($projection_type,
                               $grid_cen_lat_cmn,$grid_cen_lon_cmn,
                               $known_ri, $known_rj, $d_x, 
                               $stdlon,
                               $truelat1,$truelat2,
                               $nx_dim,$ny_dim);
   };

  # Make sure user input values are within the bounds
  if ($my_catch eq "") { 
     warn("There is a major problem with the namelist file: $template_nl.\n"), 
     return(1);
   }

   # Fill hash for nesting subroutines.
   (%parent_hash)=(%{ $ArrayOfHashes[0]} );

   # Recalculate the corners base on gridded map.
   ($latSW,$lonSW) = &map_utils::ij_to_latlon(1,1,%{ $ArrayOfHashes[0]} );
   ($latNE,$lonNE) = &map_utils::ij_to_latlon($nx_dim,$ny_dim,%{ $ArrayOfHashes[0]} );
   display_ll_ur($lonSW,$lonNE,$latNE,$latSW);   

   # Success.
   return(0);
}

# -------------------------------------------------
# switch_worlds 
#
# Not yet implemented.
# -------------------------------------------------
sub switch_worlds {

    if($global_choice eq $global_image[0]) {
         $dmap_image=$cy_dmap_image; 
     } elsif($global_choice eq $global_image[1]) {
         $dmap_image=$nh_dmap_image; 
     } elsif($global_choice eq $global_image[2]) {
         $dmap_image=$sh_dmap_image; 
     }
     present_grid_editor(1);
     set_button_state(0,$maproj_mb);
     $ncarg_idx=$proj_choices[$sug_proj_idx][3];
     set_proj_widget("$proj_choices[$sug_proj_idx][1]");
}

#_____________________________________________________
#
# ------------------------------------
# highlight_update_button
#
# Because a parameter has changed 'Update' needs 
# to be pressed, so set the color of the 'Update' 
# button to 'update' or 'normal'.
#
# NOTE: subroutine is called each time a mapping
# parameter is changed. This happens a lot.
# ------------------------------------
sub highlight_update_button {
    my ($value)=@ARG;

    $b_update_highlight=$value;

    set_button_highlight($value,$b_update);
    if ($value == 1 && $proj_label eq "None") { 
       set_button_highlight($value,$maproj_mb); }


    # Highlight the entrybox of the values that have changed.
    if ($value == 1) { 
       if ($panel_index == 2 && !$globalMap) { 
          my $who=$mw->focusCurrent; 

          # Disable tabs.
          access_to_tabs(0);

          # Color signals that display does NOT
          # update until 'Update' is pressed.
          #foreach (@entry_wdgt_list) { 
          foreach ($lonEnt, $latEnt, $lat2Ent,
                   @entry_wdgt_list2) { 
             if ($who eq $_) { $who->configure(-bg => $colorY2); } }

          if ($who eq $clatEnt) { $who->configure(-bg => $colorY2); }
          if ($who eq $clonEnt) { $who->configure(-bg => $colorY2); }

       }
    } else {

       # Un-highlight all entryboxes.
       foreach (@entry_wdgt_list, @entry_wdgt_list2) { 
         set_button_highlight($value,$_); }
    }

    # Change reset button too, unless globalMap.
    unless ($globalMap) {
     if ($projection_type eq $reset_projection_type) { 
        set_button_state($value,$b_resetvars);
     }
    }
}

# ------------------------------------
# set_snap_box_state
#
# Set the bcd background file, projection type 
# selectors AND center point widgets
# to 'disable' or 'normal'.
# ------------------------------------
sub set_snap_box_state {
    my ($state)=@ARG;

    set_button_highlight($state,$b_update);
    set_button_state($state,$b_update);

    #set_button_list_state($state,@cenlon_wdgt_list);
    #set_button_list_state($state,@cenlat_wdgt_list);
    set_button_state($state,$maproj_mb);
    set_button_state($state,$bcd_mb);

    set_button_list_state($state,@lon_wdgt_list);
    set_button_list_state($state,@lat_wdgt_list);
}

#_____________________________________________________
#
# Confirm projection values. 
#_____________________________________________________

# -----------------------------------
# set_proj_vars
#
# Set domain projection variables.
# -----------------------------------
sub set_proj_vars {
   ($p_index)=@ARG;

   # Reset map projection menubar widget.
   set_button_highlight(0,$maproj_mb);

   # Set proj label & proj index; Set the 
   # proj type for the map_utils plotting routine.
   $proj_label=$proj_choices[$p_index][0];
   $ncarg_idx=$proj_choices[$p_index][3];
   $prev_projection_type=$projection_type; # old
   $reset_projection_type=$projection_type; # old
   $projection_type=$proj_choices[$p_index][4]; # new

   # Set truelat1 and truelat2 values (i.e. $proj_choices[$p_index][1 or 2]).
   if ($panel_index!=1) {
     set_true_lats();
   }

   # Set nx_ddim.
   if ($panel_index!=1) {
     set_dnx();
     set_ddx();
   }

   if (!$globalMap) { 
      adjust_bbox();
      $dx_orig = $grid_spacing_dkm;
      $nx_orig = $nx_ddim;
   }

   access_to_nest_tabs();
   highlight_update_button(1);

   # Disable "Reset Values" button.
   #set_button_state(0,$b_resetvars);

   # Set entries as editable, or not.
   set_latlon_wdgt_state();
   set_lat2_lon_wdgt_state();
}

# R.Rozumalski- Added Rotated Lat-Lon Proj for use with NMM WRF. 
# When comparing the LC to the RL domains, the area covered by
# the LC domain, using some value of NX and DX that we'll call
# NX_LC and DX_LC, should be almost identical to that of the
# RL domain using value of NX_RL and DX_RL, where:
# NX_RL = (NX_LC + 1) / 2, and
# DX_RL  = SQRT ( 2* DX_LC**2 )
# NX_RL and DX_RL are the values displayed in the window 
# AND output to the NMM wrfsi.nl file. 

sub set_dx {
   $grid_spacing_km  = ( $projection_type eq 'RL' ) 
      ?  sqrt ( 0.5 * $grid_spacing_dkm ** 2 ) : $grid_spacing_dkm;
   $grid_spacing_km = sprintf ("%.1f",$grid_spacing_km);
}

sub set_ddx {
   if ($projection_type eq 'RL') {
      $grid_spacing_dkm=sqrt ( 2 * $grid_spacing_km ** 2 );

      # Highlight nx_ddim entrybox.
      $dx_ent->configure(-bg => $colorY2);

   } elsif ($prev_projection_type eq 'RL') {
      $prev_projection_type='';
      $grid_spacing_dkm = $grid_spacing_km;

      # Highlight nx_ddim entrybox.
      $dx_ent->configure(-bg => $colorY2);

   } else {
      $grid_spacing_dkm = $grid_spacing_km;
   }
   $grid_spacing_dkm = sprintf ("%.1f",$grid_spacing_dkm);
}

sub set_nx {
   $nx_dim = ( $projection_type eq 'RL') 
      ?  $nx_ddim * 2 - 1 : $nx_ddim;
}

sub set_dnx {
   if ($projection_type eq 'RL') {
      $nx_ddim = int( ($nx_dim +1) / 2 ); 

      # Highlight nx_ddim entrybox.
      $nx_ent->configure(-bg => $colorY2);

   } elsif ($prev_projection_type eq 'RL') {
      $prev_projection_type='';
      $nx_ddim = $nx_ddim * 2 -1;

      # Highlight nx_ddim entrybox.
      $nx_ent->configure(-bg => $colorY2);

   } else {
      $nx_ddim = $nx_dim;
   }
}

# -----------------------------------
# set_center_latlon
#
# Set the domain projection variables:
# center lat/lon to true lat/lon.
# 
# -----------------------------------
sub set_center_latlon {

     if ($globalMap) {
        $grid_cen_lat_cmn=$truelat2=$truelat1;
        $grid_cen_lon_cmn=$stdlon;
     }
}

# -----------------------------------
# set_true_lats
#
# Set the domain projection variables:
# truelat1, truelat2 to the center lat
# when user selects a projection type 
# and when domain (bbox) moves.
# 
# Also, set true lon, aka stdlon.
# -----------------------------------
sub set_true_lats {

     if ($globalMap) {
        # First guess values for std lats/lon.
        $truelat2=$truelat1=$grid_cen_lat_cmn;
        $stdlon=$grid_cen_lon_cmn;
     }
     
     if ($projection_type eq 'PS') {
        # PolarStereo Latitude 2 must be +/- 90.0";
        if (!$polar_cap) {
           if($grid_cen_lat_cmn >= 0) {
             $truelat2=90;
           } else {
             $truelat2=-90;
           }

           if (!$globalMap) { $lat2Ent->configure(-bg => $colorY2); }

        } else {
           if ($grid_cen_lat_cmn >= 0) {
              #$truelat1=$latSW if defined $latSW;;
              $truelat1=60;
           } else {
              #$truelat1=$latNE if defined $latNE;;
              $truelat1=-60;
           }

           if (!$globalMap) { $latEnt->configure(-bg => $colorY2); }
        }


     } elsif($projection_type eq 'ME') {
        # Mercator Latitude 2 is ignored.
        $truelat2=0;
        if (!$globalMap) { $lat2Ent->configure(-bg => $colorY2); }

     } else { 
        # $projection_type eq 'LC' or 'RL'.
        if ($truelat2 eq 0 or $truelat2 eq 90 or $truelat2 eq -90) {
           $my_val=$nl_var{MOAD_STAND_LATS}[1];
           if($model_name eq "LAPS"){$my_val=$nl_var{STANDARD_LATITUDE2}[0];}
           my $disp = yesno_dbox("Reset True Latitude 2", 
"True Latitude 2 has changed. Would you like to use the value $my_val", "Yes");
           if($disp eq 'Yes') { $truelat2=$my_val; }
           check_truelat_val();
           $lat2Ent->configure(-bg => $colorY2); 
        }
     }
}

# -------------------------------------------
# set_rot_lat_args
#
# Set rotation_arg and lat_arg for call to setsup.
# -------------------------------------------
sub set_rot_lat_args {

   if ($projection_type eq 'LC' or $projection_type eq 'RL') {

      $rotation_arg=$truelat2; 
      $lat2_arg=$truelat1;
   } else {

      $rotation_arg=0.;
      $lat2_arg=$truelat2;
   }
}

# -------------------------------------------
# truelat_adjusted
#
# The user adjusted values, check them and 
# highlight that update needs to be pressed.
# -------------------------------------------
sub truelat_adjusted {

  if ($t1_orig eq $truelat1 && $t2_orig eq $truelat2) { return;}
  if (!($mw->focusCurrent eq $latEnt ||
        $mw->focusCurrent eq $lat2Ent)) { return; }
 
  highlight_update_button(1);
  check_truelat_val();
}

# -------------------------------------------
# check_truelat_val 
#
# True latitude values need to be within range
# when the user types in values.
# -------------------------------------------
sub check_truelat_val {

  # Lambert 
  my $my_lmbrt=0; 
  if($projection_type eq 'LC' or $projection_type eq 'RL' ) { $my_lmbrt=1; }

  my $truelat_val="True Latitude Value Error";
  my $msg_hgh="The True Latitude absolute value cannot be greater than or equal to 90.";
  my $msg_low="The True Latitude absolute value cannot be less than or equal to 0.";
  my $msg_low2="The True Latitude absolute value cannot be less than 0.";
  my $msg_low3="The True Latitude absolute value cannot be greater than 0.";
 
  if($grid_cen_lat_cmn >= 0) {

    if ($truelat1 >= 90) {
      info_dbox($truelat_val, $msg_hgh);
      $truelat1=89.9; 
      if($projection_type eq 'ME') {$truelat1=88;} 
    } elsif ($truelat2 >= 90 && $my_lmbrt) {
      info_dbox($truelat_val, $msg_hgh);
      $truelat2=89.9; 
    } elsif ($truelat1 < 0) {
      info_dbox($truelat_val, $msg_low2);
      $truelat1= 0;
    } elsif ($truelat2 <= 0 && $my_lmbrt) {
      info_dbox($truelat_val, $msg_low);
      $truelat2= 0.1; 
    }

  } else {

    if ($truelat1 <= -90) {
      info_dbox($truelat_val, $msg_hgh);
      $truelat1=-89.9; 
      if($projection_type eq 'ME') {$truelat1=-88;} 
    } elsif ($truelat2 <= -90 && $my_lmbrt) {
      info_dbox($truelat_val, $msg_hgh);
      $truelat2=-89.9; 
    } elsif ($truelat1 > 0) {
      info_dbox($truelat_val, $msg_low3);
      $truelat1= 0;
    } elsif ($truelat2 >= 0 && $my_lmbrt) {
      info_dbox($truelat_val, $msg_low);
      $truelat2= -0.1; 
    }
  }

  $t1_orig=$truelat1;
  $t2_orig=$truelat2;

}

### Return 1 to the calling use statement ###
1;
