# ---------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# srt_vert_editor.pl 
# 	Creates the Vertical Grid Editor panel, its widgets and tools.
#
# Author: Paula McCaslin   30 Jan 2003  Original Version
# ---------------------------------------------------------------------------

#use warnings;
#use strict;
use strict 'subs';
use strict 'refs';
use vars qw(@vscheme $restrict_sigma $model_name);


# ----------------------------------
# set_vert_vars
#
# Vertical editor vars.
# ----------------------------------

sub set_vert_vars {
  package v;

  use vars qw($scheme_choice
  $logp $sig_edit_mode $pbot_mb $temp_k $gamma $v_constants 
  @default_sig_levels);


  $scheme_choice=$vscheme[0]; #'Current namelist'
  $logp=1;
  $sig_edit_mode=0;
  $pbot_mb=1013;
  if ($system_tool::model_name eq "LAPS") { $pbot_mb=101300; }
  $temp_k=288.0;
  $gamma=0.0065; #constant
  my $R=287;  #dry air constant
  my $g=9.81; #gravity constant
  $v_constants=($R*$gamma/$g); #0.190163

  @default_sig_levels= (1.000, 0.993, 0.980, 0.966, 0.950, 0.933,
              0.913, 0.892, 0.869, 0.844, 0.816, 0.786,
              0.753, 0.718, 0.680, 0.639, 0.596, 0.550,
              0.501, 0.451, 0.398, 0.345, 0.290, 0.236,
              0.188, 0.145, 0.108, 0.075, 0.046, 0.021,
              0.000);

  if ($my_nmm) {
  @default_sig_levels= 
        (1.00000 , 0.99505 , 0.98888 , 0.98153 , 0.97178 ,
         0.95971 , 0.94419 , 0.92420 , 0.90110 , 0.87509 ,
         0.84423 , 0.80683 , 0.76266 , 0.71174 , 0.65178 ,
         0.59570 , 0.54331 , 0.49439 , 0.44877 , 0.40624 ,
         0.36665 , 0.32982 , 0.29560 , 0.26384 , 0.23439 ,
         0.20712 , 0.18189 , 0.15859 , 0.13710 , 0.11729 ,
         0.09908 , 0.07991 , 0.06258 , 0.04595 , 0.02937 ,
         0.01482 , 0.00211 , 0.00000);
  }

  if ($system_tool::model_name eq "LAPS") {
  @default_sig_levels= 
        (110000., 105000., 100000., 95000., 90000., 85000.,
          80000., 75000., 70000., 65000., 60000., 55000.,
          50000., 45000., 40000., 35000., 30000., 25000.,
          20000., 15000., 10000., );
  }


}

# ----------------------------------
# create_vertGrid_panel 
#
# Vert Grid Editor Frame (v_frame).
# ----------------------------------

sub create_vertGrid_panel {
  use strict;
  use vars qw($panel3 $ve_sel $bg_canvas $vcan);

  # Set vars.
  set_vert_vars();

  # Set log pressure of top and bottom mb levels.
  set_logp_top_bot(); 

  my $v_frame=$panel3->Frame()
              ->pack(-expand => 1, -fill => 'both',
                     -padx => 5, -pady => 5);

  my $ve_toolbar=$v_frame->Frame(-height => 26)
                         ->pack(-fill => 'x');

  #------- vert display window ---------
  my $ve_map=$v_frame->Frame(-relief => 'sunken', -bd => 2)
                  ->pack(-side => 'left', -anchor => 'nw',
                         -expand => 0, -fill => 'y');

#  #panel grip used to resize 
#  $v_frame->Adjuster(-widget => $ve_map, -side => 'left')
#                   ->pack(-side => 'left', -fill => 'both');

  #panel grip used to resize 
  $v_frame->Frame(-width => 8)
                   ->pack(-side => 'left', -fill => 'both');

  #vert widgets window
  $ve_sel=$v_frame->Frame(-relief => 'sunken', -bd => 2)
                  ->pack(-side => 'left', -anchor => 'nw',
                         -expand => 1, -fill => 'both');


  #------- Create Canvas for Vert Levels -------------
  my $scroll_vcan = $ve_map->Scrolled('Canvas', 
                         (-cursor => 'crosshair white', 
                         -bg => $bg_canvas, 
                         -width => 340, 
                         -relief => 'groove', -bd => '2'),
                -scrollbars => 'osoe');
  $scroll_vcan->Subwidget('xscrollbar')->configure(-width => 10);
  $scroll_vcan->Subwidget('yscrollbar')->configure(-width => 10);


  # NOTE: pack canvas in a separate line to be able to expand & fill!
  $scroll_vcan->pack(-expand => 1, -fill => 'both',
                    -padx => 5, -pady => 5);

  $vcan=$scroll_vcan->Subwidget('canvas');

  create_vert_wdgets();
}


sub create_vert_wdgets {
  use strict;
  use vars qw($ve_sel
              $balloon $bold_font $thin_font $small_up $small_down
              $logp $pbot_mb $ptop_mb $sig_edit_mode
              $temp_k $ve_mb
             );
  use vars qw(
             $gsig
             $nz_ent $nz_lab $nz_label $nz_vert
             $pb_ent $pt_ent $temp_ent
             $roww
             $scaled_bot $scaled_top
             $sig_bold_font
             $sig_text $sig_view
             $view_edit_lab
             $v_grid
             @logp_list
             @nz_list
             );
  use vars qw($typein_lab);
  use vars qw($colorN);
  use vars qw($sigma_msg);
  use vars qw($logp_button);

   # --- Generate Levels ---
   $gsig=$ve_sel->Frame(-relief => 'groove', -bd => 2);
   if (!$restrict_sigma){
         $gsig->pack(-anchor => 'w', -fill => 'none',
                  -padx => 4, -pady => 10);
   }
  
   $gsig->gridRowconfigure( 0, -minsize => 20); 
   $gsig->gridRowconfigure( 3, -minsize =>  8); 
   $gsig->gridRowconfigure( 5, -minsize => 12); 
   $gsig->gridColumnconfigure( 0, -minsize => 10); 
   $gsig->gridColumnconfigure( 2, -minsize => 25); 
   $gsig->gridColumnconfigure( 4, -minsize => 10); 

   $roww=1;
   $ve_mb=create_vert_scheme_list();

   # Add balloon messages to inform.
   $sigma_msg="Reload namelist\nINSTALLROOT/templates/your_domain/$system_tool::namelist_arg.";
   $balloon->attach($ve_mb, -msg => \$sigma_msg);

   $roww++;
    #$nz_label="Number of Levels:"; This label is set elsewhere within code.
    $nz_lab=$gsig->Label(-text => $nz_label, -fg => $colorN) 
                   ->grid(-row => $roww, 
                          -column => 1, -columnspan => 2,
                          -sticky => 'w');

    my $nz_frame=$gsig->Frame(-highlightthickness => 0) 
                       ->grid(-row => $roww++, 
                              -column => 3, 
                              -sticky => 'e',
                              -ipadx => 2);

    $nz_ent=$nz_frame->Entry(-width => 8, 
                           -justify => 'right',
                           -textvariable => \$nz_vert) 
                   ->pack(-side => 'left',
                          -anchor => 'c');

  
    my $nz=$nz_frame->Frame(-highlightthickness => 0) 
                       ->pack();

      my $nz_up=$nz->Button(-bitmap => $small_up, 
                        -highlightthickness => 0, 
                        -takefocus => 0,
                        -command => [\&increment_value, 'nz_vert', 1, $nz_ent])
                    ->pack();
      my $nz_dw=$nz->Button(-bitmap => $small_down,
                        -highlightthickness => 0, 
                        -takefocus => 0,
                        -command => [\&increment_value, 'nz_vert', -1, $nz_ent])
                    ->pack();

  
    # Widgets are disabled when vert scheme is Default or Current.
    @nz_list =($nz_ent, $nz_up, $nz_dw);


   #--------- Generate Label -------------
   my $my_values_lab2='Generate Full-Sigma (WRF Eta) Values';
   if ($system_tool::model_name eq "LAPS") { $my_values_lab2='Generate Pressure Values'; }
   $gsig->Label(-text => $my_values_lab2)
          ->place( -x => 5, -y => -7);


   # --- Edit Levels ---
   my $esig=$ve_sel->Frame(-relief => 'groove', -bd => 2);
   if (!$restrict_sigma){
            $esig->pack(-anchor => 'w', -fill => 'none',
                  -padx => 4, -pady => 0, -ipadx => 31);
   }
  
    $esig->gridRowconfigure(0, -minsize => 20); 
    $esig->gridRowconfigure(4, -minsize =>  4); 
    $esig->gridRowconfigure(6, -minsize => 12); 

    $esig->gridColumnconfigure( 0, -minsize => 1); 

   $roww=1; 
   $typein_lab=$esig->Label(-text => 
                "Type in levels (separate values with a blank space):",
               -fg => $colorN, 
               -justify => 'left')
       ->grid( -row => $roww++, 
               -column => 1, -columnspan => 1,
               -sticky => 'w');

   $sig_text=$esig->Scrolled('Text', 
                                       (-height => 6,
                                        -width => 44,
                                        -padx => 10,
                                        -wrap => 'word',
                                        -highlightthickness => 0, 
                                        -spacing1=> 2,
                                        -exportselection => 1,
                                        -font => $thin_font,
                                        -relief => 'groove'),
                             -takefocus => '0',
                             -scrollbars => 'oe')
              ->grid(-row => $roww++,
                     -column => 1, -columnspan => 1,
                     -sticky => 'e');
  $sig_text->Subwidget('xscrollbar')->configure(-width => 10);
  $sig_text->Subwidget('yscrollbar')->configure(-width => 10);



  #@sig_font=$sig_text->configure(-font);
  #$sig_bold_font="$sig_font[3] bold";
  $sig_bold_font="$bold_font";
  
  $roww++;
  $sig_view=$esig->Button(-text => 'View Levels', 
                                 -width => 14, -justify => 'right',
                                 -state => 'disabled',
                                 -command => [\&typed_in_levels])
          ->grid(-row => $roww++, 
                 -column => 1, 
                 -sticky => 'e',
                 -ipady => 2);

    # Add balloon messages to inform.
    $balloon->attach($sig_view,
       -msg => "View values listed in text \nwindow after they are filtered.");

   #--------- Edit Levels Label -------------
   my $my_values_lab3='View/Edit (WRF Eta) Levels';
   if ($system_tool::model_name eq "LAPS") { $my_values_lab3='View/Edit Pressure Levels';}
   $view_edit_lab=$esig->Label(-text => $my_values_lab3)
        ->place( -x => 5, -y => -7);


  #--------- Vert Grid Information -------------
  $v_grid=$ve_sel->Frame(-relief => 'groove', -bd => 2)
             ->pack(-anchor => 'w',
                    -padx => 5, -pady => 10);
  
    $v_grid->gridRowconfigure( 0, -minsize => 20); 
    $v_grid->gridRowconfigure( 5, -minsize =>  8); 
    $v_grid->gridRowconfigure( 8, -minsize => 12); 

    $v_grid->gridColumnconfigure( 0, -minsize => 10); 
    $v_grid->gridColumnconfigure( 2, -minsize => 95); 
    $v_grid->gridColumnconfigure( 3, -minsize => 90); 
    $v_grid->gridColumnconfigure( 5, -minsize => 12); 
    
  #--------- Vertical Levels Widgets ----------
    $roww=1;


    $roww++;
    my $pt_lab=$v_grid->Label(-text => "Pressure at top of model (mb):", 
                              -fg => $colorN);
    if (!$restrict_sigma){
            $pt_lab->grid( -row => $roww, 
                           -column => 1, 
                           -sticky => 'w');
    }

    $pt_ent=$v_grid->Entry(-width => 8, 
                           -justify => 'right',
                           -textvariable => \$ptop_mb);
    if (!$restrict_sigma){
            $pt_ent->grid( -row => $roww,
                           -column => 3, -columnspan => 1,
                           -sticky => 'e');
    }

    my $pt=$v_grid->Frame(-highlightthickness => 0);
      if (!$restrict_sigma){
                $pt->grid( -row => $roww++, 
                           -column => 4, 
                           -sticky => 'w');
      }

      my $incre=1;
      if ($system_tool::model_name eq "LAPS") { $incre=1000;}
      my $pt_up=$pt->Button(-bitmap => $small_up, 
                      -highlightthickness => 0, 
                      -takefocus => 0,
                      -command => [\&increment_value, 'ptop_mb', $incre, $pt_ent])
                    ->pack();
      $incre=-1;
      if ($system_tool::model_name eq "LAPS") { $incre=-1000;}
      my $pt_dw=$pt->Button(-bitmap => $small_down,
                      -highlightthickness => 0, 
                      -takefocus => 0,
                      -command => [\&increment_value, 'ptop_mb', $incre, $pt_ent])
                    ->pack();

    my $pb_lab=$v_grid->Label(-text => "*Representative MSL surface pressure (mb):", 
                              -fg => $colorN) 
                   ->grid(-row => $roww, 
                          -column => 1, -columnspan => 2, 
                          -sticky => 'w');

    $pb_ent=$v_grid->Entry(-width => 8, 
                           -justify => 'right',
                           -textvariable => \$v::pbot_mb) 
                   ->grid(-row => $roww,
                          -column => 3, -columnspan => 1,
                          -sticky => 'e');

    my $pb=$v_grid->Frame(-highlightthickness => 0) 
                       ->grid(-row => $roww++, 
                              -column => 4, 
                              -sticky => 'w');

      $incre=10;
      if ($system_tool::model_name eq "LAPS") { $incre=1000;}
      my $pb_up=$pb->Button(-bitmap => $small_up, 
                     -highlightthickness => 0, 
                     -takefocus => 0,
                     -command => [\&increment_value, 'v::pbot_mb', $incre, $pb_ent])
                    ->pack();
      $incre=-10;
      if ($system_tool::model_name eq "LAPS") { $incre=-1000;}
      my $pb_dw=$pb->Button(-bitmap => $small_down,
                     -highlightthickness => 0, 
                     -takefocus => 0,
                     -command => [\&increment_value, 'v::pbot_mb', $incre, $pb_ent])
                    ->pack();

    my $temp_lab=$v_grid->Label(
                          -text => "*Representative surface temperature (K):", 
                          -fg => $colorN) 
                   ->grid(-row => $roww, 
                          -column => 1, -columnspan => 2, 
                          -sticky => 'w');

    $temp_ent=$v_grid->Entry(-width => 8, 
                          -justify => 'right',
                          -textvariable => \$v::temp_k) 
                   ->grid(-row => $roww,
                          -column => 3, -columnspan => 1,
                          -sticky => 'e');

    #$main::mw->traceVariable(\$temp_k, 'w' => [\&char_checker], 'temp_k');


    my $temp=$v_grid->Frame(-highlightthickness => 0) 
                    ->grid( -row => $roww++, 
                            -column => 4, 
                            -sticky => 'w');

      my $temp_up=$temp->Button(-bitmap => $small_up, 
                        -highlightthickness => 0, 
                        -takefocus => 0,
                        -command => [\&increment_value, 'v::temp_k', 1, $temp_ent])
                      ->pack();
      my $temp_dw=$temp->Button(-bitmap => $small_down,
                        -highlightthickness => 0, 
                        -takefocus => 0,
                        -command => [\&increment_value, 'v::temp_k', -1, $temp_ent])
                      ->pack();

    $roww++;
    if ($system_tool::model_name ne "LAPS") { 

      my $logp_lab=$v_grid->Label(-text => "Display levels in:", 
                                -fg => $colorN) 
                   ->grid(-row => $roww, 
                          -column => 1, 
                          -sticky => 'w');


      $logp_button=$v_grid->Checkbutton(  -text => 'Log Pressure',
                           -justify => 'right',
                           -selectcolor => 'yellow',
                           -variable => \$v::logp, 
                           -command => [\&set_sig_disp_mode],
                          );

       $logp_button->grid(-row => $roww++,
                          -column => 3, -columnspan => 2, 
                          -sticky => 'e');
    }
 
    # Widgets are disabled when display mode is not in log pressure.
    @logp_list=($pt_ent,  $pt_up,  $pt_dw, 
                $pb_ent,  $pb_up,  $pb_dw, 
                $temp_ent, $temp_up, $temp_dw);

    #--------- Vert Label -------------
    $v_grid->Label(-text => "Vertical Parameters")
           ->place( -x => 5, -y => -7);

 

  for my $event ("<ButtonRelease>") {
     $nz_up->Tk::bind ($event, [\&nz_vert_adjusted]);
     $nz_dw->Tk::bind ($event, [\&nz_vert_adjusted]);
     $pt_up->Tk::bind ($event, [\&pressure_adjusted]);
     $pt_dw->Tk::bind ($event, [\&pressure_adjusted]);
     $pb_up->Tk::bind ($event, [\&pressure_adjusted]);
     $pb_dw->Tk::bind ($event, [\&pressure_adjusted]);
     $temp_up->Tk::bind ($event, [\&temps_adjusted]);
     $temp_dw->Tk::bind ($event, [\&temps_adjusted]);
  } 
  for my $event ("<Return>","<Leave>") {
     $nz_ent->Tk::bind ($event, [\&nz_vert_adjusted]);
     $pt_ent->Tk::bind ($event, [\&pressure_adjusted]);
     $pb_ent->Tk::bind ($event, [\&pressure_adjusted]);
     $temp_ent->Tk::bind ($event, [\&temps_adjusted]);
  }

  # Create bindings in order to 
  # drag cursor to display sigma/pressure values.
  $vcan->Tk::bind ("<Button-1>", [\&vcan_xy_to_hgt, Ev('x'), Ev('y')]);
  $vcan->Tk::bind ("<B1-Motion>", [\&vcan_xy_to_hgt, Ev('x'), Ev('y')]);
  $vcan->Tk::bind ("<ButtonRelease-1>", [\&delete_vgrid_text]);
 
#  # Original setup.
#  set_vert_scheme_current();

}

# ---------------------------------------
# set_logp_top_bot 
#
# Set log pressure of top and bottom mb levels.
# ---------------------------------------
sub set_logp_top_bot {

  use vars qw($ptop_mb);
  if (defined $ptop_mb) {
     $scaled_top = log($ptop_mb);
  }
  if (defined $v::pbot_mb) {
     $scaled_bot = log($v::pbot_mb);
  }

}

# ---------------------------------------
# set_vert_scheme_current
#
# Original setup, new level values.
# ---------------------------------------

sub set_vert_scheme_current {

  # Load and display current levels.
  reset_vert_wdgets();
  $scheme_choice=$vscheme[0]; #'Current namelist'
  select_vert_scheme($vscheme[0]); #'Current namelist'
}

# ---------------------------------------
# reset_vert_wdgets
#
# ---------------------------------------

sub reset_vert_wdgets {

  # Display current scheme label.
  $v::pbot_mb=1013;
  if ($system_tool::model_name eq "LAPS") { $v::pbot_mb=101300; }
  mark_vlevel_axis();

  # Set logp display mode.
  set_sig_disp_mode();

  # Set edit mode to no-edit.
  set_sig_edit_mode();

  # Reset entry values.
  set_vert_vars();
}

sub create_vert_scheme_list {
    use vars qw($scheme_choice);

    @vscheme = ('Reread model namelist', 
                'Load original values', 
                'Browse...',
                'Edit levels',
                'Calc Linear in sigma', 
                'Calc Square Root in sigma', 
                'Calc 1/3 Linear 2/3 Sqrt in sigma', 
                );

    if ($model_name eq "LAPS") {
    @vscheme = ('Reload model namelist', 
                'Load original values', 
                'Browse...',
                'Edit levels',
                );
    }

    my $nom_lab=$gsig->Label(-text => "Choose what you want to do:", 
                               -fg => $colorN) 
                       ->grid( -row => $roww, 
                               -column => 1, 
                               -sticky => 'w');

    my $nom = native_optionmenu(
        $gsig,
        \$scheme_choice,
        [sub {
                my $menu = $Tk::event->W;
                select_vert_scheme($scheme_choice);
             }, 
        ],
        @vscheme,
    );
    $nom->grid(-row => $roww++, 
               -column => 3, 
               -sticky => 'e');

    my $menu = $nom->cget(-menu);
    for my $i (0 .. scalar(@vscheme) ) {
        $menu->entryconfigure($i);
    }
    return($nom);
}

# -------------------------------------------
# Subroutines to display vertical grid lines.
#
# -------------------------------------------


# -------------------------------------------
# set_sig_disp_mode
#
# Display levels in sigma or pressure values. 
# Set vertical pressure and temp widgets to 
# normal or disabled.
# -------------------------------------------

sub set_sig_disp_mode {

    switch_focus_main();

    foreach my $widgt (@logp_list) {
     if ($v::logp == 0) {
       $widgt->configure( -state => 'disabled', -fg => $disabled_color);
     } else {
       $widgt->configure( -state => 'normal', -fg => $normal_color);
     }
    }

    activate_vert_scheme();
    mark_vlevel_axis();
}

# -------------------------------------------
# activate_vert_scheme
#
# Get the vertical scheme type, then 
# call sub select_vert_scheme
# -------------------------------------------

sub activate_vert_scheme {

    if ($v::sig_edit_mode) { 
       typed_in_levels();
    } else { 
       # Load and display current levels.
       my $vert_label=$ve_mb->cget(-text);
       if ($vert_label eq $vscheme[2]) { 
          typed_in_levels();
          return; 
       }
       select_vert_scheme($vert_label);
    }
}

# -------------------------------------------
# select_vert_scheme
#
# Process user request for vertical scheme.
# Note: subroutine draw_vert_lines is called from 
# sub get_vert_lines.
# -------------------------------------------

sub select_vert_scheme {
    my ($scheme)=@ARG;

    set_button_list_state(1,@nz_list);

    if($v::sig_edit_mode == 1) {
       $v::sig_edit_mode=0; 
       set_sig_edit_mode();
    } else {
       check_vpressures();
       check_vtemps();
    }
    check_nz_vert();

    if ($scheme eq $vscheme[0]) {
       #'Reload Domain values'
       $sigma_msg=
"Reload namelist 
INSTALLROOT/templates/$domain_select/$system_tool::namelist_arg";
       set_button_list_state(0,@nz_list);
       @sig_levels=@nl_sig_levels;
       check_sigmas();
       get_vert_lines();

    } elsif ($scheme eq $vscheme[1]) {
       #'Load original values'
       $sigma_msg=
"Load original default values from 
SOURCE_ROOT/data/static/$system_tool::namelist_arg";
       set_button_list_state(0,@nz_list);
       @sig_levels=@v::default_sig_levels;
       check_sigmas();
       get_vert_lines();

    } elsif ($scheme eq $vscheme[4]) {
       #'Calc Linear in sigma'
       $sigma_msg=
"Automatic calculation of a Linear 
distribution of values for sigma";
       calc_linear_levels();
       check_sigmas();
       get_vert_lines();

    } elsif ($scheme eq $vscheme[5]) {
       #'Calc Square Root in sigma'
       $sigma_msg=
"Automatic calculation of Square Root
distribution of values for sigma";
       calc_sqrt_levels();
       check_sigmas();
       get_vert_lines();

    } elsif ($scheme eq $vscheme[6]) {
       #'Calc Linear Sqrt in sigma'
       #'Calc 1/3 Linear 2/3 Sqrt', 
       $sigma_msg=
"Automatic calculation of a distribution
of sigma values where the top 1/3 are
Linear and bottom 2/3 are Square Root";
       calc_lsqrt_levels();
       check_sigmas();
       get_vert_lines();

    } elsif ($scheme eq $vscheme[2]) {
       #'Browse for file with level data',
       $sigma_msg="Browse for file on your system\ncontaining vertical level data";
       set_button_list_state(0,@nz_list);
       browse_for_sigma_file();

    } elsif ($scheme eq $vscheme[3]) {
       #'Edit levels'
       $sigma_msg="Edit vertical levels using\n the text editor below";
       $v::sig_edit_mode=1; 
       set_sig_edit_mode();

    }

}

# --- Set vertical levels editing widgets to normal or disabled. --- 

sub set_sig_edit_mode {
    use vars qw($act_nz_vert $typein_lab);

    if ($v::sig_edit_mode) { 
       $view_edit_lab->configure(-text =>'Edit Levels');
       $typein_lab->configure(-text => 
                "Type in levels (separate values with a blank space):");
       $sig_textState='normal';
       $sig_text->configure(-state => $sig_textState, -relief => 'sunken');
       set_button_state(1,$sig_view);
       set_button_list_state(0,@nz_list);

       # Adjust the nz label; &set_nz_label().
       $nz_vert=$act_nz_vert;
       $nz_lab->configure(-text => "Number of Levels:");

    } else {
       $view_edit_lab->configure(-text =>'View Levels');
       $typein_lab->configure(-text => "");
       $sig_textState='disabled';
       $sig_text->configure(-state => $sig_textState, -relief => 'groove');
       set_button_state(0,$sig_view);
       set_button_list_state(1,@nz_list);
    }

    # Clear current selection, if any with (S)electionClear;
    $sig_text->SelectionClear();
}

# --- User entered new values, check validity and update display. --- 

sub pressure_adjusted {
    if ($ptop_orig eq $ptop_mb &&
        $pbot_orig eq $v::pbot_mb) { return;}
    if (!($mw->focusCurrent eq $pt_ent ||
          $mw->focusCurrent eq $pb_ent)) { return;}

    check_vpressures();
    mark_vlevel_axis();
    activate_vert_scheme();
}

sub temps_adjusted {
    if ($tempk_orig eq $v::temp_k) { return;}
    if ($mw->focusCurrent ne $temp_ent) { return;}

    check_vtemps();
    activate_vert_scheme();
}

sub nz_vert_adjusted { 
    if ($nz_vert_orig eq $nz_vert) { return;}
    if ($mw->focusCurrent ne $nz_ent) { return;}

    check_nz_vert();
    activate_vert_scheme();
}

# --- Make sure values entered by user are valid. --- 

sub check_vpressures {

    my $my_ckeck=0;
    if ($model_name ne "LAPS") {
       if ($ptop_mb <   10) {$ptop_mb=  10; $my_ckeck=1;}
       if ($ptop_mb >  250) {$ptop_mb= 250; $my_ckeck=1;}
       if ($v::pbot_mb <  750) {$v::pbot_mb= 750; $my_ckeck=1;}
       if ($v::pbot_mb > 1100) {$v::pbot_mb=1100; $my_ckeck=1;}
    }

    if ($ptop_mb > $v::pbot_mb) {
        info_dbox("Pressure Value Error", "Pressure top cannot be greater than pressure bottom.");
        $ptop_mb=100;
        $v::pbot_mb=1013;
        if ($model_name eq "LAPS") { 
           $ptop_mb=10000;
           $v::pbot_mb=101300; }
    }
    $pbot_orig=$v::pbot_mb;
    $ptop_orig=$ptop_mb;
    if ($my_ckeck = 1) { $hint_msg="Pressure value was outside of acceptable range."; }
}

sub check_vtemps {

    my $my_ckeck=0;
    if ($v::temp_k < 200) {$v::temp_k=200; $my_check=1;}
    if ($v::temp_k > 320) {$v::temp_k=320; $my_check=1;}
    $tempk_orig=$v::temp_k;
    if ($my_ckeck = 1) { $hint_msg="Temperature value was outside of acceptable range."; }
}

sub check_nz_vert {
    
    my $my_ckeck=0;
    if ($nz_vert <  15) {$nz_vert= 15; $my_check=1; $range="less than 15";}
    if ($nz_vert > 100) {$nz_vert=100; $my_check=1; $range="greater than 100";}
    $nz_vert_orig=$nz_vert;
    if ($my_ckeck = 1) { $hint_msg="Number of vertical values were $range."; }
}

# -------------------------------------------
# check_sigmas
#
# Check level values. The list must include 0 
# and 1. Delete values less than 0 or greater 
# than 1 or duplicate values.
# -------------------------------------------

sub check_sigmas {


    # Levels list must include 0 and 1.
    if ($model_name eq "WRF") { push(@sig_levels, (0, 1) ); }

    # Sort list into descending order.
    @sig_levels=sort reverse_numerically @sig_levels;

    # Clean up list.
    my @clean_levels=();
    @clean_levels=@sig_levels;
    my $i=0;
    my $j=0;

    for my $ith (0 .. scalar(@sig_levels)-1) {

       if ( 
       ($model_name eq "WRF"  && ($sig_levels[$ith] > 1 || $sig_levels[$ith] < 0)) ||
       ($model_name eq "LAPS" && ($sig_levels[$ith] > 110500 || $sig_levels[$ith]+0.1 < $ptop_mb))
         ) {

          # Delete values less than MIN level or greater than MAX level.
          splice(@clean_levels, $i-$j, 1);
          $j++;

       } elsif ($sig_levels[$ith] == $sig_levels[$ith-1] ) {

          # Delete duplicate values.
          splice(@clean_levels, $i-$j, 1);
          $j++;
       }

       $i++;

    } 

    @sig_levels=@clean_levels;
    set_nz_label();
}

# -------------------------------------------
# set_nz_label
#
# Update vert levels label for Linear Sqrt, only.
# (For more info see subroutine calc_lsqrt_levels).
# -------------------------------------------

sub set_nz_label {

    # Determine value for nz_vert.
    $act_nz_vert=scalar(@sig_levels);


    # The nz_vert label is not accurate for 'Calc 1/3 Linear 2/3 Sqrt in sigma', 
    # do not update nz_vert var, however update nz label.
    if ( ($ve_mb->cget(-text) eq $vscheme[44]) && !($v::sig_edit_mode) ) { 
      $nz_label="Number of Levels (actual $act_nz_vert):";
    } else {
      $nz_vert=$act_nz_vert;
      $nz_label="Number of Levels:";
    } 
    $nz_lab->configure(-text => $nz_label);

    if ($nz_vert <  15) {
        info_dbox("Number of Levels Error", "The Number of levels should not be less than 15. 
There were only $nz_vert found.");
       
    }

}

# -------------------------------------------
# reverse_numerically
#
# Sort all values from 1 to 0 in 
# descending order.
# -------------------------------------------

sub reverse_numerically {$b <=> $a}

sub calc_linear_levels {

    for(my $i = 0; $i <= $nz_vert-1; $i++ ) {
        my $sigma  = ($i / ($nz_vert-1));
        push(@sig_levels, "$sigma ");
    }
}

sub calc_sqrt_levels {

    for(my $i = 0; $i <= $nz_vert-1; $i++ ) {
        my $sigma  = sqrt($i / ($nz_vert-1) );
        push(@sig_levels, "$sigma ");
    }
}

# -------------------------------------------
# calc_lsqrt_levels
#
# Combination of 1/3 linear and 2/3 square root method.
#
# Since these two methods overlap each other, the 
# value for number of elements in the array is usually 
# larger than nz_vert requested by the user.
# -------------------------------------------

sub calc_lsqrt_levels {

    for (my $i = 0; $i <= $nz_vert-1; $i++ ) {
        my $sigma = ($i / ($nz_vert-1));
        my $sigma2 = sqrt($i / ($nz_vert-1) );

        my $plevel  = $sigma * ($v::pbot_mb - $ptop_mb) + $ptop_mb;
        my $plevel2 = $sigma2 * ($v::pbot_mb - $ptop_mb) + $ptop_mb;

        if($plevel  < 333.33) { push(@sig_levels, "$sigma "); } 
        if($plevel2 > 333.33) { push(@sig_levels, "$sigma2 "); }
    }
}

# -------------------------------------------
# typed_in_values
#
# Get values from text edit box.
# -------------------------------------------

sub typed_in_levels {

    my $sinput= $sig_text->get("1.0", "end");
    chomp($sinput);
    @sig_levels=split / /, $sinput;

    # Load and display current levels.
    check_sigmas();
    get_vert_lines();
}

# ----------------------------------------
# get_vert_lines 
#
# Convert levels to log pressure if logp eq 1.
# Update name list variables.
# Determine min and max values.
# Display vertical lines.
# ----------------------------------------
sub get_vert_lines {

    #if (!defined $v::ptop_mb) { return; }

    @hgt_levels=();
    @vert_levels=();

    # Allow update to text widget.
    $sig_text->configure(-state => 'normal'); 
    $sig_text->delete("1.0", "end");

    $paula=0; 
    foreach my $arg (@sig_levels) {
       # Format levels for display text.
       $fmt="%0.0f, ";
       if ($model_name eq "WRF") {$fmt="%.3f, ";}
       $sig_text->insert("end", sprintf ($fmt, $arg) );

       # Reformat level for output array.
       $fmt="%0.0f";
       if ($model_name eq "WRF") {$fmt="%.3f";}
       $arg=sprintf ($fmt, $arg);  

       # Create the vertical lines.
       if ($model_name eq "LAPS") {
           #$ptop_mb=10000;
           $plevel = $arg;
           my $new_arg = sprintf ("%.3f", 
             ( (log($plevel)-$scaled_top)/($scaled_bot-$scaled_top)) );

           # Compute height values.
           my $Z=sprintf ("%.0f", ($v::temp_k/$v::gamma)
                         *(1.0-($plevel/$v::pbot_mb)**($v::v_constants)) 
                 );
           push(@hgt_levels, "$Z ");
           push(@vert_levels, "$new_arg ");

           $paula++;

       } elsif ($v::logp eq 1) {

           $plevel = $arg * ($v::pbot_mb - $ptop_mb) + $ptop_mb;
           my $new_arg = sprintf ("%.3f", 
             ( (log($plevel)-$scaled_top)/($scaled_bot-$scaled_top)) );

           # Compute height values.
           my $Z=sprintf ("%.0f", ($v::temp_k/$v::gamma)
                         *(1.0-($plevel/$v::pbot_mb)**($v::v_constants)) 
                 );
           push(@hgt_levels, "$Z ");
           push(@vert_levels, "$new_arg ");
       } else {

           push(@vert_levels, "$arg ");
       }
    }

    $sig_text->see('end');

    # Return text widget to orig state.
    $sig_text->configure(-state => $sig_textState);

    # Update name list variables.
    my $keey='LEVELS';
    if ($model_name eq "LAPS") { $keey='PRESSURES'; }
    $nl_var_max{$keey} = scalar(@sig_levels);  # num of values
    $nl_var{$keey} = [@sig_levels];            # array of values

    # Determine min and max values.
    min_max_sigma_heights();
 
    # Display vertical lines.
    draw_vert_lines();
}

# -------------------------------------
# min_max_sigma_heights 
#
# Calculate height distance extremes (meters).
# Save the corresponding sigma min and max values, in
# order to display indicator arrows showing min, max. 
# -------------------------------------
sub min_max_sigma_heights {


    if ($v::logp ne 1) { 
       # Don't calc height extremes unless log press is set.
       $hint_msg="";
       return; 
    }

    # Calc height extremes.
    $min_dx_hgt=99999.9;
    $max_dx_hgt=0.0;
    $zlev_max=scalar(@hgt_levels);
    if ($zlev_max < 3) {
        $min_dx_idx= "1.0";
        $min_dx_idx2="1.0";
        $max_dx_idx= "1.0";
        $max_dx_idx2="1.0";
        return;
    };

    # Word length =7 is a function of these chars: "0.9500,^", 
    my $word_len=7;

    my $idx;
    for $idx (0 .. ($zlev_max-2)) {
       my $dx_hgt=$hgt_levels[$idx+1]-$hgt_levels[$idx]; 

       if ($min_dx_hgt > $dx_hgt) { 
           $min_dx_hgt=$dx_hgt;
           $min_dx_sig=$vert_levels[$idx];
           $min_dx_idx=$idx*$word_len;
           $min_dx_idx2=$min_dx_idx + $word_len-1;
       } elsif ($max_dx_hgt < $dx_hgt) { 
           $max_dx_hgt=$dx_hgt;
           $max_dx_sig=$vert_levels[$idx];
           $max_dx_idx=$idx*$word_len;
           $max_dx_idx2=$max_dx_idx + $word_len-1;
       }
    }

    # Set text widget indices (e.g. 1.020 to 1.026) in order to
    # indicate the min, max thicknesses. Once the indices are
    # calculated they are used to chage the font/color of the
    # value.
    $min_dx_idx= "1.$min_dx_idx";
    $min_dx_idx2="1.$min_dx_idx2";
    $max_dx_idx= "1.$max_dx_idx";
    $max_dx_idx2="1.$max_dx_idx2";


    $min_dx_hgt=sprintf ("%.0f",$min_dx_hgt);
    $max_dx_hgt=sprintf ("%.0f",$max_dx_hgt);

    if($panel_index == 3 and $model_name ne "LAPS") { vert_hint_msg(); };
    return;
}

sub vert_hint_msg {
    use vars qw/$hint_msg $min_dx_hgt $min_dx_hgt/;

    $hint_msg="Sigma height distance extremes (m):   min=$min_dx_hgt, max=$max_dx_hgt";
    return;
}


# --------------------------------------------
# draw_vert_lines 
#
# 
# Draw vertical lines from level values. 
# --------------------------------------------

sub draw_vert_lines {

    my $bd=30;
    my $cx1=105;
    my $cx2=$cx1+120;
    my $cx3=300;
    my $cy;
  
    my $y_height=360;
  
    # Delete previous lines.
    $vcan->delete('level_lines', 'level_arrows');

    # Allows line colors to taper off into background.
    @vert_line_color=($bg_canvas, 'gray50', 'gray60', 'gray70', 'gray80');

    # Create the vertical lines.
    foreach my $color (@vert_line_color) {
      my $i=0;
      foreach my $lev (@vert_levels) {
        $cy=$lev*$y_height+$bd;
        $vcan->createLine($cx1, $cy, $cx2, $cy,
              -fill => $color, -width => 1.0, 
              -tags => 'level_lines');

        # Create height labels: This loops several times 
        # (because many colors are used), but we need label 
        # created only once, then display every 4th value.
        if ($color eq $bg_canvas && $i%4 == 0) {
          #if ($hgt_levels[$i] >=0) { 
           #Suppress values below 0m.
           $vcan->createText($cx3, $cy,
                 -fill => 'white', -text => $hgt_levels[$i], 
                 -anchor => 'e',   -tags => 'level_lines');
          #}
        }
        $i++;

      }
      $cx1+=3;
      $cx2-=3;
    }

    # Put a label near the highest vertical line.
    $vcan->createText($cx3, $cy, 
                 -fill => 'white', #-text => $hgt_levels[$zlev_max-1], 
                 -anchor => 'e',   -tags => 'level_lines');

    if ($v::logp eq 1) {

      if ($model_name eq "LAPS") { @sig_levels=(); return; }

      # Indicate sigma height distance extremes text values.
      $sig_text->tagConfigure('bold', -font => $sig_bold_font,
                                      -foreground => 'navy');
      $sig_text->tagAdd('bold', $min_dx_idx, $min_dx_idx2);
      $sig_text->tagAdd('bold', $max_dx_idx, $max_dx_idx2);

      # Indicate sigma height distance extremes line values.

      # Add dashes and arrows to demark lines; largest and smallest dx
      # -dash '1 8' translates to one blank pixel for every 8 pixels.

      # Max thickness level   
      $cy=$max_dx_sig*$y_height+$bd;
      $vcan->createLine($cx1, $cy, $cx2, $cy,
                -fill => $bg_canvas, -dash => '1 8', -width => 1.0, 
                -tags => 'level_lines');
      $vcan->createLine($cx1-25, $cy-3, $cx1-15, $cy-3, 
                -arrow => 'last', -arrowshape => [5, 5, 5],
                -fill => 'white', -width => 1.0, 
                -tags => 'level_arrows');
      # Min thickness level   
      $cy=$min_dx_sig*$y_height+$bd;
      $vcan->createLine($cx1, $cy, $cx2, $cy,
                -fill => $bg_canvas, -dash => '1 2', -width => 1.0, 
                -tags => 'level_lines');
      # Add a blank segment in front of arrow
      $vcan->createLine($cx1-28, $cy, $cx1-25, $cy,
                -fill => $bg_canvas, -width => 1.0, 
                -tags => 'level_arrows');
      $vcan->createLine($cx1-25, $cy, $cx1-15, $cy, 
                -arrow => 'last', -arrowshape => [4, 4, 4],
                -fill => 'white', -width => 1.0, 
                -tags => 'level_arrows');
    }

    # Clear array.
    @sig_levels=();
}

# --------------------------------------------
# mark_vlevel_axis 
#
# Mark vertical level axis.
# --------------------------------------------

sub mark_vlevel_axis {

    # Reset log pressure of top and bottom mb levels.
    set_logp_top_bot(); 

    # Set pixel reference points.
    my $bd=30;
    my $cx1=115;
    my $cx2=$cx1+100;
    my $cx3=300;
  
    my $y_height=360;
    my $cy1=$bd;
    my $cy2=$bd+$y_height;
    use vars qw($colorW);
  
    # Create axis box only once.
    if (!$levels_box) {
       $levels_box=$vcan->createRectangle($cx1-30, $cy2, $cx2+30, $cy1, 
             -outline => $colorW, -width => 1, -tags => 'levels_box');
    }

    # Previous labels need to be deleted.
    $vcan->delete('level_labels');
  
    # Axis labels are outside box, so reset cx1.
    $cx1=$cx1-40;


    my $i;
    if ($v::logp eq 1) {
       # Create pressure labels from e.g. "1100 to 50 mb".

       my $press_units='(mb)          Log p';
       if ($model_name eq "LAPS") { $press_units='(pa)          Log p'; }
       $vcan->createText($cx1+15, $cy2+40, -fill => $colorW, 
                            -text => $press_units,
                            -tags => 'level_labels');
       $vcan->createText($cx3, $cy2+40, -fill => $colorW, 
                            -text => 'Height AGL (m)', -anchor => 'e',
                            -tags => 'level_labels');
       my $botp;
       if($v::pbot_mb >= 1100){
         $botp=1100;
       } elsif($v::pbot_mb > 1000){
         $botp=1000;
       } elsif($v::pbot_mb > 900){
         $botp=900;
       } elsif($v::pbot_mb > 800){
         $botp=800;
       } else{
         $botp=700;
       }

       my $my_interval=100;
       #if ($model_name eq "LAPS") {$botp=$v::pbot_mb; $my_interval=10000;}
       if ($model_name eq "LAPS") {$botp=110000; $my_interval=10000;}

       for($i = $botp; $i >= $my_interval; $i -= $my_interval ) {
          my $scaled_logp=(log($i)-$scaled_top)/($scaled_bot-$scaled_top);
          $cy=$scaled_logp*$y_height+$bd;
          $vcan->createText($cx1, $cy, -fill => $colorW, -text => $i, 
                            -anchor => 'e',
                            -tags => 'level_labels');
       }

       # Create axis box only once.
       if ($model_name eq "LAPS") {
         if (!$levels_box2) {
          $levels_box2=$vcan->createRectangle($cx1+10, $cy2+20, $cx2+30, $cy1, 
                -outline => $colorW, -width => 1, -tags => 'levels_box2');
         }
       }

       # At top of domain add label, e.g. "50 mb".
       $cy=$bd;
       $vcan->createText($cx1, $cy, -fill => $colorW, -text => int($ptop_mb), 
                            -anchor => 'e',
                            -tags => 'level_labels');

       if ($model_name eq "LAPS") { 
          ## At bottom of domain add label, e.g. "950 mb".
          #$vcan->createText($cx1-30, $cy2, -fill => $colorW, -text => $v::pbot_mb, 
          #                  -anchor => 'e',
          #                  -tags => 'level_labels');

          # At bottom of domain add label, e.g. "0 m".
          $vcan->createText($cx3, $cy2, -fill => $colorW, -text => "0m ", 
                            -anchor => 'e',
                            -tags => 'level_labels');
       } else {
          # At bottom of domain add label, e.g. "950 mb".
          $vcan->createText($cx1, $cy2+10, -fill => $colorW, -text => $v::pbot_mb, 
                            -anchor => 'e',
                            -tags => 'level_labels');
       }

    } else {
       # Create sigma labels from 0 to 1.
       $vcan->createText($cx1, $cy2+40, -fill => $colorW, 
                            -text => 'Sigma (WRF Eta) Values', -anchor => 'w',
                            -tags => 'level_labels');

       for(my $i = 0.0; $i <= 1.0; $i += .2 ) {
          $cy=$i*$y_height+$bd;
          my $ii=sprintf ("%.1f",$i);
          $vcan->createText($cx1, $cy, -fill => $colorW, -text => $ii, 
                            -anchor => 'e',
                            -tags => 'level_labels');
       }
    }

}

# --------------------------------------------
# vcan_xy_to_hgt 
#
# Drag cursor to display sigma/pressure values.
# --------------------------------------------

sub vcan_xy_to_hgt {
    my ($c, $c_x, $c_y)=@ARG;
    my $x=$c->canvasx($c_x);
    my $y=$c->canvasy($c_y);
  
    my $bd=30;
    my $cx1=85;
    my $cx2=$cx1+160;

    my $y_height=360;
    my $cy1=$bd;
    my $cy2=$bd+$y_height;
 
    # Clear cursor value.
    delete_vgrid_text();

    if ($x < $cx1 || $x > $cx2) {return;}
    if ($y < $cy1 || $y > $cy2) {return;}
 
    $cy=($y-$bd)/$y_height;

    # Format cursor value.
    if ($v::logp eq 1) {
       $cy=exp(($cy*($scaled_bot-$scaled_top))+$scaled_top);
       $cy=sprintf ("%.0f",$cy);
    } else {
       $cy=sprintf ("%.3f",$cy);
    } 
    if ($model_name eq "LAPS") { $cy=(int($cy/100)*100); }
    $vcan->createText($c_x, $c_y, -fill => 'white', 
                                 -text => "   $cy",
                                 -anchor => 'w',
                                 -font => $legend_font,
                                 -tags => "sigma_val");
}

# --------------------------------------------
# delete_vgrid_text;  Clear cursor value.
# --------------------------------------------
sub delete_vgrid_text { $vcan->delete("sigma_val"); }

# --------------------------------------------
# browse_for_sigma_file {
#
# Browse directories for a file with level data.
# Values above 1 and below 0 are filtered out
# as are words.
# --------------------------------------------
sub browse_for_sigma_file {

    my $file = $mw->getOpenFile();
    if ($file eq "") {return;}

    # Read and parse vert levels from a user file.
    open(NL,"$file") or die "Can't open file: $file\n";
    while (<NL>) {
	chomp($_);
	if (m/,/) {
	   @lines=split /,/, $_;
	} else {
	   @lines=split;
	}
	push(@sig_levels, @lines);

    } 
    close(NL);
    
    # Load and display current vert levels.
    check_sigmas();
    get_vert_lines();
}

### Return 1 to the calling use statement ###
1;
