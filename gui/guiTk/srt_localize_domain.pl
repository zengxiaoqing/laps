# ---------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# srt_localize_domain.pl 
#	Creates the Localize Domain panel, its widgets and tools 
#	and activates model localization.
#
# Author: Paula McCaslin   30 Jan 2003  Original Version
# ---------------------------------------------------------------------------


#use warnings;
#use strict;
use strict 'subs';
use strict 'refs';

use vars qw(%nl_var %nl_var_orig %nl_var_max
            @nl_var @nl_var_array $nl_var_max 
            @num_nl_entries $num_nl_sections @nl_section
            $nl_sig_levels $geog_path $ROOT_GEOG);


# ----------------------------------
# create_locParms_panel 
#
# writeNamelist Frame (w_frame)
# ----------------------------------

sub create_locParms_panel {
use vars qw($normal_color $disabled_color 
           $update_color $bg_color $colorN);

  my $w_pane=$panel4->Frame() 
                    ->pack(-expand => 1, -fill => 'both',
                           -padx => 10, -pady => 15);

  #--------- Namelist Information -------
  $w_sel=$w_pane->Frame() ->pack(-anchor => 'nw');

  $w_sel->gridRowconfigure(0, -minsize => 10); 
  $w_sel->gridRowconfigure(1, -minsize =>  5); 

  my ($iwidth, $justify, $col, $sec_lab, $ffont); 
  my ($i, $j, $k);

 if ($model_name eq 'WRF') {
  for $i (3,4) { 
    if ($i eq 3) {$iwidth = 10; $justify='right'; $col=1;  #hgridspec
                  $ffont = $bold_font;
                  $sec_lab= "Horizontal Grid Spec"; }
    if ($i eq 4) {$iwidth = 35; $justify='left';  $col=2;  #sfcfiles
                  $ffont = $thin_font;
                  $sec_lab= "Static Geographical Data Files"; }

    #--------- Namelist Information -------
    my $wentry=$w_sel->LabFrame(-label => $sec_lab, -labelside => 'acrosstop') 
                  ->grid( -row => 1,
                          -column => $col,
                          -sticky => 'nsew', # fill area
                          -padx => 4);
  
    #--------- Localization Selector ----------
    my $w_key=$wentry->Frame() ->pack(-padx => 12, -pady => 5);

    #--------- Grid Information -------
    $w_key->gridColumnconfigure( 2, -minsize => 2); 

    my ($state,$state_color);
    if ($i eq 3) {
      #hgridspec
      $state_color=$disabled_color;
      $state='disabled'; 
      $relief='sunken';
    } else {
      #sfcfiles
      $state_color=$normal_color;
      $state='normal';
      $relief='sunken';
    }

    my $ro=0;
    if ($i eq 4) {

        $w_key->Label(-text => "GEOG_DATAROOT",
		      -fg => $colorN) 
	      ->grid( -row => $ro,
		      -column => 1,
		      -sticky => 'w');
        
        my $w_path_ent=$w_key->Entry(
                      -width => $iwidth,
		      -justify => $justify,
                      -relief => 'flat',
                      -state => 'disabled',
                      -selectbackground => 'gray80',
		      -textvariable => \$ROOT_GEOG) 
	      ->grid( -row => $ro, 
		      -column => 3, 
		      -columnspan => 2, 
		      -sticky => 'w');

        # Add balloon messages to inform user.
        $balloon->attach($w_path_ent, -msg => 
"To change this path go to 'Path Preferences' 
under the 'Edit' pull-down menu."); 

        $ro++; 
        $w_key->Label()->grid();
        $ro++; 
     }

    for $j (1 .. $num_nl_entries[$i]) {
          
       #Only use 0th element for now as there are no nests.
       #my $max=$nl_var_max{$key}; # all elements
       my $max=1;                   # 0th element

       my $key=$nl_var_array[$i][$j]; 
       for ($k=0; $k < $max; $k++) {
 
          $w_key->Label(-text => $key, 
                        -fg => $colorN) 
                ->grid( -row => $ro,
                        -column => 1,
		        -columnspan => 2, 
                        -sticky => 'w');

          if ($key eq 'TOPTWVL_PARM_WRF' || $key eq 'SILAVWT_PARM_WRF') {
              $state_color=$normal_color;
              $state='normal';
              $relief='sunken';
          }

          #$w_tag=$w_key->Entry(-width => $iwidth,
          if ($i eq 3) {
             if ($key =~ /PARENT|ORIGIN/) {
                $justify='left';} else {$justify='right';
             }
          }
          my $x_view=$w_key->Entry(-width => $iwidth,
                        -justify => $justify,
                        -state => $state,
                        -relief => $relief,
                        -fg => $state_color,
                        -font => $ffont,
                        -textvariable => \$nl_var{$key}[$k]) 
                ->grid( -row => $ro++, 
                        -column => 4, 
                        -sticky => 'e');
          # Make entry contents viewable.
          if ($i eq 4) { $x_view->xview((length $geog_path) - 25); }
          if ($i eq 6) { $x_view->xview((length $si_path) - 25); }


       }
       #if ($i eq 3) { push @l_tag_list, $w_tag; }
    }


    if ($i eq 4) {
        $w_key->Label(-text => "")->grid(-row => $ro++);
        $ro++;
        $geog_path_button=$w_key->Button(-text => "Update Path", 
                       -command => [\&set_path_to_geog])
	      ->grid( -row => $ro,
		      -column => 2,
		      -columnspan => 3,
 		      -sticky => 's');
 
        $geog_path_button->Tk::bind ("<ButtonRelease-2>", [\&restore_path_to_geog]);
        # Add balloon messages to inform user.
        $gmsg=
"Replace the string: '$geog_path'
with: '$ROOT_GEOG'.
Use MouseButton1 to make this change,
use MouseButton2 to undo this change."; 
        $balloon->attach($geog_path_button, -msg => \$gmsg);

    } 

  }

  # Change the Entry box state to 'normal' for the following variable:
  #  TOPTWVL_PARM_WRF & SILAVWT_PARM_WRF
  #my $nn = scalar(@l_tag_list);
  #foreach my $entri ($w_tag_list[$nn-1], $w_tag_list[$nn-2] ) {
  #      $entri->configure(-state => 'normal', -fg => $normal_color); }


 } elsif ($model_name eq 'LAPS') {

    my $ro=1;
    $i=1;
    if ($i eq 1) {$iwidth = 10; $justify='right'; $col=1;  #hgridspec
                  $ffont = $thin_font;
                  $sec_lab= "Grid and Namelist Specifications"; }

    #--------- Namelist Information -------
    my $wentry=$w_sel->LabFrame(-label => $sec_lab, -labelside => 'acrosstop') 
                  ->grid( -row => 1,
                          -column => $col,
                          -columnspan => 3,
                          -sticky => 'nsew', # fill area
                          -padx => 4);
  
    #--------- Localization Selector ----------
    my $w_key=$wentry->Frame() ->pack();

    #--------- Grid Information -------
    $w_key->gridRowconfigure( 0, -minsize => 15); 
    $w_key->gridRowconfigure( 1, -minsize =>  8); 
    $w_key->gridRowconfigure( 3, -minsize =>  8); 
    $w_key->gridRowconfigure( 9, -minsize => 35); 
    $w_key->gridColumnconfigure( 0, -minsize => 50); 
    $w_key->gridColumnconfigure( 1, -minsize => 270); 
    $w_key->gridColumnconfigure( 5, -minsize => 50); 
    #$w_key->gridColumnconfigure( 0, -minsize => 15); 
    #$w_key->gridColumnconfigure( 2, -minsize => 15); 
    #$w_key->gridColumnconfigure( 5, -minsize => 40); 
    #$w_key->gridColumnconfigure(11, -minsize => 15); 
    #$w_key->gridColumnconfigure(17, -minsize => 15); 

    my ($state,$state_color);
    if ($i eq 1) {
      #hgridspec
      $state_color=$disabled_color;
      #$state='disabled'; 
      $state='normal';
      $relief='sunken';
    } else {
      #sfcfiles
      $state_color=$normal_color;
      $state='normal';
      $relief='sunken';
    }


    $col_1=1; 
    $col_2=4; 

    if (0) {
    for $j (3 .. 13) {
          my $key=$nl_var_array[$i][$j]; 
          my $max=1;                   # 0th element
 
          $w_key->Label(-text => $key, 
                        -fg => $colorN) 
                ->grid( -row => $ro,
                        -column => $col_1,
		        -columnspan => 2, 
                        -sticky => 'w');

          my $x_view=$w_key->Entry(-width => $iwidth,
                        -justify => $justify,
                        -state => $state,
                        -relief => $relief,
                        -fg => $state_color,
                        -font => $ffont,
                        -textvariable => \$nl_var{$key}[$k]) 
                ->grid( -row => $ro++, 
                        -column => $col_2, 
                        -sticky => 'e');
    }
    }

  $roww=1; 
  $nl_typein_lab=$w_key->Label(-text => "",
               -fg => $colorN, 
               -justify => 'left')
       ->grid( -row => $roww++, 
               -column => 1, -columnspan => 4,
               -sticky => 'w');


  $nl_text=$w_key->Scrolled('TextUndo',  # subclass of Text widget.
                             ( -height => 20,
                               -width => 80,
                               -padx => 10,
                               -wrap => 'none',
                               -highlightthickness => 0, 
                               -spacing1=> 2,
                               -exportselection => 1,
                               -font => $thin_font,
                               -relief => 'groove'),
                             -takefocus => '0',
                             -scrollbars => 'se')
                  ->grid( -row => $roww++, 
                        -column => $col_1, 
		        -columnspan => 4, 
                        -sticky => 'e');
  $nl_text->configure(-state => 'normal', -relief => 'sunken', -bd => 2); 

  $nl_text->Subwidget('xscrollbar')->configure(-width => 10);
  $nl_text->Subwidget('yscrollbar')->configure(-width => 10);
#  $nl_text->bindtags([$nl_text, ref($nl_text), $nl_text->toplevel, 'all']);
#  $nl_text->bind('<p>' => sub{make($nl_text,$wlist,split //,$nl_text->get('0.0','end'));Tk->break});

  #@btags = $nl_text->bindtags;
  #$nl_text->bindtags( [ @btags[1,0,2,3] ] );
 
  $roww++;
  $roww++;
  $geog_path_button=$w_key->Button(-text => 'Update Geog Path', 
                                 -width => 16, -justify => 'right',
                                 -command => [\&set_path_to_geog])
          ->grid(-row => $roww, 
                 -column => 2, 
                 -sticky => 'e',
                 -ipadx => 10,
                 -ipady => 2);
  $geog_path_button->Tk::bind ("<ButtonRelease-2>", [\&restore_path_to_geog]);

  $w_key->Button(-text => 'Undo', 
                               -width => 8, -justify => 'right',
                               -command => sub{ $nl_text->undo; })
          ->grid(-row => $roww, 
                 -column => 3, 
                 -sticky => 'e',
                 -ipadx => 0,
                 -ipady => 1);
  $w_key->Button(-text => 'Redo', 
                               -width => 8, -justify => 'right',
                               -command => sub{ $nl_text->redo; })
          ->grid(-row => $roww, 
                 -column => 4, 
                 -sticky => 'e',
                 -ipadx => 0,
                 -ipady => 1);

    # Add balloon messages to inform.
        $gmsg=
"To replace the geography dir path string
use MouseButton1; to undo this change 
use MouseButton2."; 
    $balloon->attach($geog_path_button, -msg => \$gmsg);

    if (0) {
# ----  GEOG_DATAROOT interface.
        $iwidth = 35; 
        $justify='left';
        $ro=1; $col_1+=6; $col_2+=6; 
        $w_key->Label(-text => "GEOG_DATAROOT",
		      -fg => $colorN) 
	      ->grid( -row => $ro,
		      -column => $col_1,
		      -sticky => 'w');
        
        my $w_path_ent=$w_key->Entry(
                      -width => $iwidth,
		      -justify => $justify,
                      -state => 'disabled',
                      -relief => 'flat',
                      -selectbackground => 'gray80',
                      -font => $bold_font,
		      -textvariable => \$ROOT_GEOG) 
	      ->grid( -row => $ro, 
		      -column => $col_2, 
		      -columnspan => 2, 
		      -sticky => 'w');

        # Add balloon messages to inform user.
        $balloon->attach($w_path_ent, -msg => 
"To change this path go to 'Path Preferences' 
under the 'Edit' pull-down menu."); 
    }

# ----  

    $w_key->Label(-text => "")->grid(-row => $ro++);
    $ro++; 
    $state_color=$normal_color;
    $state='normal';
    if (0) {
    for $j (23 .. 30) {
          my $key=$nl_var_array[$i][$j]; 
          my $max=1;                   # 0th element
 
          $w_key->Label(-text => $key, 
                        -fg => $colorN) 
                ->grid( -row => $ro,
                        -column => $col_1,
		        -columnspan => 2, 
                        -sticky => 'w');

          my $x_view=$w_key->Entry(-width => $iwidth,
                        -justify => $justify,
                        -state => $state,
                        -relief => $relief,
                        -fg => $state_color,
                        -font => $ffont,
                        -textvariable => \$nl_var{$key}[$k]) 
                ->grid( -row => $ro++, 
                        -column => $col_2, 
                        -sticky => 'e');
    }
# ----  
        $ro++;
        $geog_path_button=$w_key->Button(-text => "Update Path", 
                       -command => [\&set_path_to_geog])
	      ->grid( -row => $ro,
		      -column => $col_2,
		      -columnspan => 3,
 		      -sticky => 's');
 
        $geog_path_button->Tk::bind ("<ButtonRelease-2>", [\&restore_path_to_geog]);
        # Add balloon messages to inform user.
        $gmsg=
"Replace the string: '$geog_path'
with: '$ROOT_GEOG'.
Use MouseButton1 to make this change,
use MouseButton2 to undo this change."; 
        $balloon->attach($geog_path_button, -msg => \$gmsg);
    }
# ----  


 } #WRF or LAPS

}

# ----------------------------------
# nl_text_fill 
#
# Fill text window editor with namelist
# ----------------------------------

sub nl_text_fill {

  #add colors
  $nl_text->tagConfigure( 'tag0',   -foreground => 'navy'); 
  $nl_text->tagConfigure( 'tag1',   -foreground => 'gray45', -font => $bold_font );
  $nl_text->tagConfigure( 'tag2',   -foreground => 'black', -font => $bold_font );


  # Allow input to the Text Window.
  $nl_typein_lab->configure(-text => "");
  $nl_text->configure(-state => 'normal'); 
  $nl_text->delete("1.0", "end");

  # Input data to the Text Window.
  $k=0;
  for $i (1 .. $num_nl_sections) {
    $nl_text->insert("end", sprintf ("&$nl_section[$i]\n"), 'tag0'); #section header
    for $j (1 .. $num_nl_entries[$i]) {
     $key=$nl_var_array[$i][$j]; 

     # Dont allow edits to grid and projection parms.
     if($key =~ m/C6|STAND|GRID|_L_CMN|NK_LAPS/ ) { 
        $nl_text->insert("end", sprintf ("$key = "), 'tag0') ;
        $nl_text->insert("end", sprintf ("$nl_var{$key}[$k]\n"), 'tag1' ) ;
     } else {
        $nl_text->insert("end", sprintf ("$key = "), 'tag0') ;
        $nl_text->insert("end", sprintf ("$nl_var{$key}[$k]\n"), 'tag2' ) ;
     }
    }
    $nl_text->insert("end", sprintf ("/\n"), 'tag0'); #section close
  }
      
  # Don't let user edit a red-colored parameter.
  $nl_text->tagBind ('tag0', "<Any-Enter>" => sub{ shift->configure(-cursor=>"X_cursor red"); });
  $nl_text->tagBind ('tag0', "<Any-Leave>" => sub{ shift->configure(-cursor=>"left_ptr"); });
  $nl_text->tagBind ('tag0', "<Key>", \&nl_text_null);
  $nl_text->tagBind ('tag1', "<Any-Enter>" => sub{ shift->configure(-cursor=>"X_cursor red"); });
  $nl_text->tagBind ('tag1', "<Any-Leave>" => sub{ shift->configure(-cursor=>"left_ptr"); });
  $nl_text->tagBind ('tag1', "<Key>", \&nl_text_null);
}

# ----------------------------------
# nl_text_null
#
# Fill text window editor with namelist
# ----------------------------------

sub nl_text_null {

    shift->markSet('insert','end');
}

# ----------------------------------
# nl_text_readit 
#
# Write namelist to TEMPLATE and MOAD_DATAROOT
# ----------------------------------

sub nl_text_readit {

      $install_nl=0; # This is an existing namelist, so don't update any $nl_var_orig{}[].
      @nl_lines = split /\n/, $nl_text->get("1.0", "end");
      get_namelist_array(@nl_lines);
      nl_text_fill();
}

# ----------------------------------
# nl_text_save 
#
# Write namelist to TEMPLATE and MOAD_DATAROOT
# ----------------------------------

sub nl_text_save {

  $nl_typein_lab->configure(-text => "");
  $nl_text->configure(-state => 'disabled', -relief => 'groove', -bd => 2, -bg => $bg_color); 
  $nl_save_but->configure(-state => 'disabled');
  $nl_reset_but->configure(-state => 'disabled');

  my $rsp=yesno_dbox("Accept Namelist Changes", 
        "Changes to this namelist are not checked. 
Do you want to continue?
(Selecting 'No' will reset all changes to their previous value.)", "No"); 
  if ($rsp eq "No"){ 
     nl_text_fill();
     $hint_msg="All namelist values are reset.";
     return(1); 
  }

  # Write template namelist (and dataroot namelist, if it exists).
  write_textWindow_namelist();

  # Read template namelist and 'merge' with install namelist.
  my $panel_index_save=$panel_index;
  $panel_index=1; # Must be = 1.
  load_namelist();
  $panel_index=$panel_index_save;

  highlight_update_button(0);

}

# ----------------------------------
# create_initControls_panel 
#
# initControls Frame (e_pane)
# ----------------------------------

sub create_initControls_panel {
use vars qw($d_panel1 
           $normal_color $disabled_color 
           $update_color $bg_color $colorN);

  my $e_pane=$d_panel1->Frame()
                    ->pack(-expand => 1, -fill => 'both',
                           -padx => 10, -pady => 15);

  #--------- Namelist Information -------
  $e_sel=$e_pane->Frame() ->pack(-anchor => 'nw');

  $e_sel->gridRowconfigure(0, -minsize => 10); 
  $e_sel->gridRowconfigure(1, -minsize =>  5); 

  my ($iwidth, $justify, $row, $col, $sec_lab, $ffont); 
  my ($i, $j, $k);

  foreach $i (5,6) {
    if ($i eq 2) {$iwidth = 5; $justify='right';  $row=1, $col=1;  #hgridspec
                  $ffont = $bold_font;
                  $sec_lab= "File Time Specifications"; }
    if ($i eq 5) {$iwidth = 10; $justify='right'; $row=1, $col=2;  #sfcfiles
                  $ffont = $bold_font;
                  $sec_lab= "Interpolate Controls"; }
    if ($i eq 6) {$iwidth = 45; $justify='left';  $row=1, $col=3;  #sfcfiles
                  $ffont = $thin_font;
                  $sec_lab= "Standard Initialization Paths"; }

    #--------- Namelist Information -------
    my $eentry=$e_sel->LabFrame(-label => $sec_lab, -labelside => 'acrosstop') 
                  ->grid( -row => $row,
                          -column => $col,
                          -sticky => 'nsew', # fill area
                          -padx => 4);
  
    #--------- Localization Selector ----------
    my $e_key=$eentry->Frame() ->pack(-padx => 12, -pady => 5);

    #--------- Grid Information -------
    $e_key->gridColumnconfigure(2, -minsize => 2); 
    if($i eq 6) {$e_key->gridRowconfigure(0, -minsize => 67);}

    my ($state,$state_color);
    if ($i ne 3) {
      $state_color=$normal_color;
      $state='normal';
      $relief='sunken';
    }

    my $ro=1;
    if ($i eq 6) {
        $e_key->Label(-text => "EXT_DATAROOT",
		      -fg => $colorN) 
	      ->grid( -row => $ro,
		      -column => 1,
		      -sticky => 'w');
        
        my $e_path_ent=$e_key->Entry(
                      -width => $iwidth, 
		      -justify => $justify,
                      -relief => 'flat',
                      -state => 'disabled',
                      -selectbackground => 'gray80',
		      -textvariable => \$ROOT_EXT) 
	      ->grid( -row => $ro, 
		      -column => 3, 
		      -columnspan => 2, 
 	              -sticky => 'w');

        # Add balloon messages to inform user.
        $balloon->attach($e_path_ent, -msg => 
"To change this path go to 'Path Preferences' 
under the 'Edit' pull-down menu."); 

        $ro++; 
        $e_key->Label()->grid();
        $ro++; 
    }

    for $j (1 .. $num_nl_entries[$i]) {
          
       my $max=1;                   # 0th element

       my $key=$nl_var_array[$i][$j]; 
       if ($key eq 'LEVELS' || $key eq 'PTOP_PA') { 
       } else {
       for ($k=0; $k < $max; $k++) {
 
 

          $e_key->Label(-text => $key, 
                        -fg => $colorN) 
                ->grid( -row => $ro,
                        -column => 1,
		        -columnspan => 2, 
                        -sticky => 'w');

          if ($key eq 'TOPTWVL_PARM_WRF' || $key eq 'SILAVWT_PARM_WRF') {
              $state_color=$normal_color;
              $state='normal';
              $relief='sunken';
          }

          my $x_view=$e_key->Entry(-width => $iwidth,
                        -justify => $justify,
                        -state => $state,
                        -relief => $relief,
                        -fg => $state_color,
                        -font => $ffont,
                        -textvariable => \$nl_var{$key}[$k]) 
                ->grid( -row => $ro++, 
                        -column => 4, 
                        -sticky => 'e');
          # Make entry contents viewable.
          if ($i eq 6) { $x_view->xview((length $si_path) - 15); }
       }
       }
    }

    if ($i eq 6) {
        $e_key->Label(-text => "")->grid(-row => $ro++);
        $ro++;
        $si_path_button=$e_key->Button(-text => "Update SI Path", 
                       -command => [\&set_path_to_si])
	      ->grid( -row => $ro,
		      -column => 2,
		      -columnspan => 3,
		      -sticky => 's');

        $si_path_button->Tk::bind ("<ButtonRelease-2>", [\&restore_path_to_si]);
        # Add balloon messages to inform user.
        $balloon->attach($si_path_button, -msg => \$si_msg);
    } 

  }   

}

# --------------------------------------
# set_path_to_si 
#
# Replace directory path name from surface files 
# with new one.
# --------------------------------------
sub set_path_to_si {
     my $keyy;

     # Strip off pathname to si files (section=6). 
     for $j (1 .. $num_nl_entries[6]) {
        $keyy=$nl_var_array[6][$j]; 
        # Store orig.
        if ($nl_sivar{$keyy}[0] ne $nl_var{$keyy}[0]){
            $nl_sivar{$keyy}[0]=$nl_var{$keyy}[0];}
        $nl_var{$keyy}[0]="'$ROOT_EXT/extprd'";
     }
  
     $si_path=$ROOT_EXT;

     # Un-highlight.
     set_button_highlight(0,$si_path_button); 
}

# --------------------------------------
# restore_path_to_si 
#
# Strip off current directory path name from surface files 
# replace with new one.
# --------------------------------------
sub restore_path_to_si {

     my $keyy;
     # Strip off pathname to si files (section=6). 
     for $j (1 .. $num_nl_entries[6]) {
        $keyy=$nl_var_array[6][$j]; 
        $nl_var{$keyy}[0]=$nl_sivar{$keyy}[0];
     }

     $ARG=$nl_var{LBCPATH}[0];
     s/\/[^\/]*$//;                # Strip off directory path.
     s/'//g;                       # Strip off quote(s).
     $si_path=$ARG;
}

# --------------------------------------
# set_path_to_geog 
#
# Strip off current directory path name from surface files 
# replace with new one.
# --------------------------------------
sub set_path_to_geog {


     $nl_text->see("32.0"); # View last row.
     $nl_text->Tk::Text::FindAll("-exact", "-nocase", "$geog_path"); 
     if ($geog_path eq $ROOT_GEOG) {
        info_dbox("Geographical Directory Error", 
          "Path 1:$geog_path is equal to \nPath 2:$ROOT_GEOG\n
To define 'Path 1 go to 'Path Preferences' 
under 'Edit' pull-down menu.");
        return(1);
     }

     my $rsp=yesno_dbox("Change Geog Dir Path", 
        "Change $geog_path to \n$ROOT_GEOG?", "Yes"); 
     if ($rsp eq "No"){ return(1); }

     my $keyy;
     # Strip off pathname to surface files (section=4). 
     if ($model_name ne "LAPS") {
       for $j (1 .. $num_nl_entries[4]) {
        $keyy=$nl_var_array[4][$j]; 
        $nl_var{$keyy}[0] =~ s/'$geog_path\//'$ROOT_GEOG\//;

       }
     } else {
       for $j (1 .. 50) {
        $keyy=$nl_var_array[1][$j]; 
        $nl_var{$keyy}[0] =~ s/'$geog_path\//'$ROOT_GEOG\//;

       }
     }

     set_button_highlight(0,$geog_path_button);
     $geog_path_hold=$geog_path; 
     $geog_path=$ROOT_GEOG; 
     
     # Success.
     nl_text_fill();
     $nl_text->see("32.0"); # View last row.
     return(0);
}

# --------------------------------------
# restore_path_to_geog 
#
# Strip off current directory path name from surface files 
# replace with new one.
# --------------------------------------
sub restore_path_to_geog {

     if ($geog_path_hold eq "") {
        $hint_msg="There is nothing to undo...";
        return(1);
     }
     $nl_text->see("32.0"); # View last row.
     $nl_text->Tk::Text::FindAll("-exact", "-nocase", "$geog_path"); 
     if ($geog_path eq $geog_path_hold) {
        info_dbox("Geographical Directory Error", 
           "$geog_path is equal to \n$geog_path_hold");
        set_button_highlight(0,$geog_path_button);
        return(1);
     }

     my $rsp=yesno_dbox("Change Geog Dir Path", 
        "Change $geog_path to \n$geog_path_hold?", "Yes"); 
     if ($rsp eq "No"){ return(1); }

     my $keyy;
     # Strip off pathname to surface files (section=4). 
     if ($model_name ne "LAPS") {
        for $j (1 .. $num_nl_entries[4]) {
           $keyy=$nl_var_array[4][$j]; 
           $nl_var{$keyy}[0] =~ s/'$geog_path\//'$geog_path_hold\//;

        }
     } else {
       for $j (1 .. 50) {
        $keyy=$nl_var_array[1][$j]; 
        $nl_var{$keyy}[0] =~ s/'$geog_path\//'$geog_path_hold\//;

       }
     }

     $geog_path=$geog_path_hold; 
     
     # Success.
     nl_text_fill();
     $nl_text->see("32.0"); # View last row.
     return(0);
}

# --------------------------------------
# check_geog_path 
#
# Check the validity of the directory path to sfcfiles (i.e. geog).
# --------------------------------------
sub check_geog_path {
     # Strip off quote(s).

     if ($geog_path eq "" || $geog_path=~ m/ /) {
        my $rsp=info_dbox("Static Geographical Directory Error", 
           "Null character was found. You must enter an absolute \ndirectory path to model geography data."); 
       return(1);
     }
     if ( !-d "$geog_path") {
        my $rsp=yesno_dbox("Static Geographical Directory Error", 
           "Path to $geog_path does not exist. \nAccept path >$geog_path< anyway?", "No"); 
        if ($rsp eq "No"){ 
            $geog_path=$ROOT_GEOG;
            return(1);
        } 
     }
     return(0);
}

# --------------------------------------
# sync_geog_path 
#
# Sync the the directory path to sfcfiles (i.e. geog) 
# if the entered value differs from orig value.
# --------------------------------------

sub sync_geog_path {

      if ($ROOT_GEOG ne $geog_path) { 
         #set_button_highlight(1,$geog_path_button);
         $msg="The GEOG_DATAROOT differs from the path\nto geography files listed above. By pressing\n'Yes' the 'Update Path' button will be invoked.";
         my $ans=yesnocancel_dbox("Press 'Update Path'",$msg);

         if($ans eq 'Yes'){ 
           # Press button for user.
           #if ( !$geog_path_button->invoke() ) { return(1); };
           $mw->update;
           $mw->idletasks;
           sleep(1);
           write_namelist();
           $cntl_next->invoke();
           return(1);
         } elsif($ans eq 'No'){ 
           #set_button_highlight(0,$geog_path_button); 
           return(0);
         } else { 
           return(1);
         }
      };

      # Success.
      return(0);
} 


# ----------------------------------
# wrap_set_geog_dataroot
#
# Browse for a new EXT_DATAROOT. If necessary files exist, 
# then great. If not, then copy them from the ROOT_INSTALL.
# ----------------------------------
sub wrap_set_geog_dataroot {

    $ROOT_GEOG=browse4dir($ROOT_GEOG);
    set_geog_dataroot();
}


# ----------------------------------
# set_geog_dataroot
#
# ----------------------------------
sub set_geog_dataroot {

    my $err_stat=0;
    if (!-d "$ROOT_GEOG") {
       info_dbox("Directory Error", "The GEOG_DATAROOT entered does not exist.
It will be replaced with $GEOG_DATAROOT.");
       $ROOT_GEOG="$GEOG_DATAROOT"; 

       $err_stat=1;
    }

    if ((-e "$ROOT_GEOG/topo_30s") && (-e "$ROOT_GEOG/islope") &&
        (-e "$ROOT_GEOG/landuse_30s") && (-e "$ROOT_GEOG/greenfrac") ){
       # Necessary files exist.

    } else {

       # Necessary files do not exist.
       my $rsp=info_dbox("Directory Error", 
"Files are missing from
$ROOT_GEOG.\n
It is suggested that you ftp the necessary files to this location,
or browse for the file location, if the files exist on your system."); 
       #$err_stat=1;
    }

    if ($geog_path ne $ROOT_GEOG) {
      set_button_highlight(1,$geog_path_button); }
    return($err_stat);
}


#___________________________________________
#

# ----------------------------------
# create_localizeDomain_panel 
#
# localizeDomain Frame (l_frame)
# ----------------------------------

sub create_localizeDomain_panel {
use vars qw($panel5 
           $normal_color $disabled_color 
           $update_color $bg_color $colorN);

use vars qw($l_frame);

  $l_frame=$panel5->Frame() 
                  ->pack(-expand => 1, -fill => 'both',
                         -padx => 20,  -pady => 15);

  #--------- Execute Command -------
  my $l_execute=$l_frame->Frame(-relief => 'groove', -bd => 2) 
                 ->pack(-anchor => 'nw', 
                        -pady => 15, -ipadx => 60, -ipady => 60);

    $l_execute->gridRowconfigure(3, -minsize => 40); 
    $l_execute->gridRowconfigure(4, -minsize => 20); 

    $loc_frame=$l_execute->Frame()
                         ->grid(-row => 2, -column => 1,
                                -columnspan => 1, -sticky => 'nw');
  
    $execCmd = $loc_frame->ExecuteCommand(-command => 'hostname')
                         ->pack(-fill => 'both');

    #--------- Force Domain Localization Widgets -------
    $l_force=$l_execute->Frame();

    $forceLocalization=1;
    $forceDomLoc_mb=$l_force->Menubutton( -indicator => 1,
                             -tearoff => 0,
                             -relief => 'raised',
                             -bd => 2,
                             -width => 4,
                             -activebackground => $update_color,
                             -anchor => 'e',
                             -textvariable => \$forceLocalization)
             ->pack(-fill => 'both', -side => 'right', -padx => 10);

    $l_force->Label( -text => "Force domain localization of Domain ID:", -fg => $colorN)
             ->pack(-fill => 'both', -side => 'right');

    update_force_loc_mb();

    $saveExecCmd = $loc_frame->Button(
                           -width => 14,
                           -command => [\&save_ExecCmd],
                           -text => 'Save Command')
                          ->pack(-fill => 'none', -anchor => 'e', -padx => 5);
            #->pack(-fill => 'both', -side => 'right', -padx => 5);

    $refreshExecCmd = $loc_frame->Button(
                           -width => 14,
                           -command => [\&localization_command],
                           -text => 'Refresh');
   #$refreshExecCmd->pack(-fill => 'none', -anchor => 'e', -padx => 5, -pady => 5);


    $balloon->attach($execCmd->Subwidget('doit'), -msg => "Press to run window_domain_rt.pl
Press everytime the namelist projection parameters change."); 
    $balloon->attach($saveExecCmd, -msg => "Press to save window_domain_rt.pl to a file.");
   #$balloon->attach($refreshExecCmd, -msg => "Press to create window_domain_rt.pl\nwith command line arguments.");


    #--------- Execute Command Label --------
    $l_execute->Label(-text => "Localize - using window_domain_rt.pl")
        ->place( -x => 5, -y => -7);
}

# -----------------------------------
# present_fwidget
#
# Set the menubar choices equal to the number of nests-1.
# -----------------------------------
sub present_fwidget {
    
    if ($show_forceWidget) {
	# Show fwidget
	$l_force->grid(-row => 4, -column => 1, 
 		       -columnspan => 2, -sticky => 'we');
        update_force_loc_mb();
    } else {
	# Hide fwidget
	$l_force->gridForget();
    }

}

# -----------------------------------
# update_force_loc_mb
#
# Set the menubar choices equal to the number of nests-1.
# -----------------------------------

sub update_force_loc_mb {

    # Reset this to 1.
    $forceLocalization=1;

    # Clear list. 
    if (defined $forceDomLoc_mb && $forceDomLoc_mb->cget(-menu)) {
       $forceDomLoc_mb->cget(-menu)->delete(0, 'end'); }

    # Build list. 
    if ($nest_id_num > 1) {
    $forceDomLoc_mb->command(-label => "None", 
			     -command => [\&force_localization, "None"] );
    }
    for ($i=1; $i < $nest_id_num+1; $i++) {
       $forceDomLoc_mb->command(-label => $i, 
				-command => [\&force_localization, $i] );
    }
}

sub force_localization {
   $forceLocalization=@_[0];
   if (@_[0] eq "None") {
	# Unsetenv var.
	$ENV{FORCE_LOCALIZATION}="";
   } else {
	# Setenv var.
	$ENV{FORCE_LOCALIZATION}=@_[0];
   }
   print " oooo FORCE_LOCALIZATION $ENV{FORCE_LOCALIZATION} oooo \n";
   update_force_loc_mb();
}

# -----------------------------------------
# localization_command 
#
# Confirm localization command. It is created only
# after the domain name has been selected.
# -----------------------------------------
sub localization_command {

use vars qw($sys_cmd_update);

    my $d_flag="$dataroot_select/$domain_select";
    my $c_flag="";
    if (!-e $d_flag) { $c_flag="-c"; }
    my $cmd0="$ROOT_INSTALL/etc/window_domain_rt.pl"; 
    my $cmd1=" -w $window_domain_arg2";
    my $cmd2=" -s $ROOT_SOURCE";
    my $cmd3=" -i $ROOT_INSTALL";
    my $cmd4=" -d $d_flag";
    my $cmd5=" -t $ROOT_TEMPLATES/$domain_select";
    my $cmd6=" $c_flag";
    my $cmd7=" $geog_path";

    # Create localization command.
    $sys_cmd="$cmd0 $cmd1 $cmd2 $cmd3 $cmd4 $cmd5 $cmd6";

    # Display localization command.
my $localiza_cmd="
$cmd0
\tModel:\t\t   $cmd1
\t$root_env_label[0]:\t   $cmd2
\t$root_env_label[1]:\t   $cmd3
\t$root_env_label[2]:\t   $cmd4
\tDomain:\t\t   $cmd5
\tConfigure:\t   $cmd6\n
Your Geography directory path is: $cmd7\n\n";

    # Enter localization command.
    $execCmd->configure(-display_command => $localiza_cmd);
    $execCmd->configure(-command => $sys_cmd);
    localization_button_highlight();

}


# ----------------------------------
# localization_button_highlight 
#
# Highlight "Run Localization" button"
# ----------------------------------
sub localization_button_highlight {
    $execCmd->Subwidget('doit')->configure(-background => $colorY2); #'#ffffaa'); # yellow
}


# ----------------------------------
# localization_began 
#
# ----------------------------------
sub localization_began {

    set_button_state(0,$cntl_back);
    set_button_state(0,$cntl_next);

    switch_focus_main();
    $hint_msg="NOTE: Localization can take SEVERAL minutes, depending on the size of the domain...";
}

# ----------------------------------
# localization_done 
#
# Check for successful domain localization
# then change buttons accordingly.
# ----------------------------------
sub localization_done {
  
    $static_select=
      "$dataroot_select/$domain_select/static/static.$window_domain_arg";

    # Check for successful domain localization.
    if(-e $static_select) {
       # Valid and localized domain - so add a checkmark.
       add_checkmark();

       # Enable "Interpolate Data" button.
       interp_data_but_state(1); 

       # Allow user to create graphics.
       reset_domain_graphics_vars();

       # For NMM link static.wrfsi.rotlat to static.wrfsi.d01 
       if ($projection_type eq 'RL') {
          link( "$dataroot_select/$domain_select/static/static.$window_domain_arg",
            "$dataroot_select/$domain_select/static/static.wrfsi.d01");
       }

    }

    set_button_state(1,$cntl_back);
    set_button_state(1,$cntl_next);
    $cntl_next->focus;
}
# ----------------------------------
# save_ExecCmd
#
# Write the command line to a file in $ROOT_INSTALL/user_scripts.
# ----------------------------------

sub save_ExecCmd { # write_grib_prep_cmd 


  # Create directory, if necessary.
  my $user_scripts="$ROOT_INSTALL/user_scripts";
  my $user_scripts_file="$user_scripts/window_domain_rt_$domain_select.sh";

  if (!-d $user_scripts) {
      mkdir $user_scripts, 0777
      or (&info_dbox("Make Directory Error", "Cannot make $user_scripts.
Change write permissions on INSTALLROOT (chmod), then press 'Write File' again.")), 
         return(1);
  }

  # Write command.
  open(LCMD,">$user_scripts_file") or 
    &fail_dbox("Failure","Can't open file: $user_scripts_file."), return(1);
  print LCMD "# The following command runs window_domain_rt.pl to localize domain.\n";
  print LCMD "$sys_cmd\n";
  close(LCMD);
  chmod 0775, $user_scripts_file;

  # Let user know file was written.
  &info_dbox("Wrote File", "Wrote file 'window_domain_rt_$domain_select.sh' located in directory\n$user_scripts containing the command:\n\n$sys_cmd"); 

  return(0);
}

#___________________________________________
#

# ----------------------------------
# create_domainGraphics_panel 
#
# If NCARG && NCL_COMMAND are not set (and valid), then
# for now do NOT create or display geog graphics interface.
# ----------------------------------
sub create_domainGraphics_panel {

  my $dg_frame=$panel6->Frame()
                      ->pack(-expand => 1, -fill => 'both',
                             -padx => 20,  -pady => 15);

  #--------- Localization Output -------
  my $dg_output=$dg_frame->Frame(-relief => 'groove', -bd => 2) 
                 ->pack(-anchor => 'nw',
                        -pady => 15, -ipady => 70);

    $dg_output->gridRowconfigure(0, -minsize =>   0); 
    $dg_output->gridRowconfigure(2, -minsize =>   0); 
    $dg_output->gridRowconfigure(4, -minsize =>   2); 
    $dg_output->gridRowconfigure(6, -minsize =>  10); 
    $dg_output->gridRowconfigure(8, -minsize => 170); 
    $dg_output->gridColumnconfigure(0, -minsize => 100); 
    $dg_output->gridColumnconfigure(2, -minsize =>  30); 
    $dg_output->gridColumnconfigure(4, -minsize =>  30); 
    $dg_output->gridColumnconfigure(6, -minsize => 130); 


  @display_choice = (
            'All',
            'Terrain elevation',
            'Land use',
            'Top Layer of soil',
            'Bottom Layer of soil',
            'Annual min greenness fraction',
            'Annual max greenness fraction',
            'Annual Soil Temp',
            'Terrain slope index',
            'Max snow albedo',
            'Max snow albedo interpolated',
            'Land water mask',
  );

  @render_choice = qw(
            meta
            avc
            use
            stl
            sbl
            gnn
            gnx
            tmp
            slp
            alb
            albint
            lnd
  );

  if($model_name eq "LAPS"){
  @display_choice = (
            'All',
            'Terrain elevation',
            'Land use',
            'Land fraction',
            'Max snow albedo',
            'Max snow albedo interpolated',
  );

  @render_choice = qw(
            meta
            avg
            use
            ldf
            alb
            albint
  );
  } 

  # Choose domain id number.
  $id_label=$dg_output->Label(-text => "  Choose domain ID",
                              -font => "Helvetica -10")
               ->grid(-row => 0, -column => 1, 
                      -sticky => 'sw');

  # Choose graphic.
  $dg_output->Label(-text => "  Choose graphic type",
                              -font => "Helvetica -10")
               ->grid(-row => 0, -column => 3, 
                      -sticky => 'sw');

  # Choose what you want to do.
  $dg_output->Label(-text => "  Choose what you want to do",
                              -font => "Helvetica -10")
               ->grid(-row => 0, -column => 5, 
                      -sticky => 'sw');

  $select_id_mb=$dg_output->Menubutton(-indicator => 1,
                           -tearoff => 0,
                           -relief => 'raised',
                           -bd => 2,
                           -width => 6,
                           -anchor => 'c',
                           -activebackground => $update_color,
                          )
                   ->grid( -row => 1, -column => 1);

  configure_domain_id_mb();

  # List of graphics.
  $d_type="";
  $view_graphics_mb=$dg_output->Menubutton(-indicator => 1,
                           -tearoff => 0,
                           -relief => 'raised',
                           -bd => 2,
                           -width => 25,
                           -anchor => 'c',
                           -activebackground => $update_color,
                           -textvariable => \$d_type)
             ->grid( -row => 1, -column => 3);

  my $i=0;
  until ($i > scalar(@display_choice)-1 ) {
      $view_graphics_mb->command(-label => $display_choice[$i],
                                 -command => [ \&set_graphics_type, $i]);
      $i++;
  }

  $b_create_graphics=$dg_output->Button(
                           -width => 16,
                           -command => [\&create_graphics],
                           -text => 'Create graphics')
                   ->grid( -row => 1, -column => 5);

  $b_view_graphics=$dg_output->Button(
                           -width => 16,
                           -command => [\&view_graphics],
                           -text => 'View graphics')
                   ->grid( -row => 3, -column => 5);

  $b_delete_graphics=$dg_output->Button(
                           -width => 16,
                           -command => [\&delete_graphics],
                           -text => 'Delete graphics')
                   ->grid( -row => 5, -column => 5);

  $b_save_graphics_cmd=$dg_output->Button(
                           -width => 16,
                           -command => [\&write_graphics_cmd],
                           -text => 'Save Command',
                           -state => 'disabled')
                   ->grid( -row => 7, -column => 5);

  # Add balloon messages to inform user.
  #$static_select
  $balloon->attach($b_create_graphics, -msg => 
    "Create static geographical data images from
static/static.$window_domain_arg
netCDF file using NCARGraphics ncl (4.3 or higher)."); 
  $balloon->attach($b_view_graphics, -msg => 
    "View static geographical data images\nusing NCARGraphics idt."); 
  $balloon->attach($b_delete_graphics, -msg => 
    "Delete static geographical data images."); 
  $balloon->attach($b_save_graphics_cmd, -msg => 
    "Save generate_images.pl with args to a file."); 
  

  #--------- View Domain Output Label -------------
  $dg_output->Label(-text => "Geographical Data")
           ->place( -x => 5, -y => -7);

}

# ----------------------------------
# configure_domain_id_mb
#
# $select_id_mb=$dg_output->Menubutton(-indicator => 1,
# Configure the domain id menubar.
# ----------------------------------
sub configure_domain_id_mb {
  if (defined $ROOT_NCARG) {
    # Clear list. 
    if (defined $select_id_mb && $select_id_mb->cget(-menu)) {
       $select_id_mb->cget(-menu)->delete(0, 'end'); }

    # Build list. 
    my $i;
    for ($i=1; $i <= $num_nests; $i++) {
       $select_id_mb->command(-label => $i,
                              -command => [ \&set_domain_id_num, $i]);
    }
  }
}

# ----------------------------------
# reset_domain_graphics_vars
#
# Domain d01.
# Generate all graphics types.
# ----------------------------------
sub reset_domain_graphics_vars {

    if (!defined $view_graphics_mb) { return; }

    # Reset all widgets.
    set_graphics_but_state(0); 
    set_button_state(0,$view_graphics_mb);
    $d_type="Choose";
    $g_type="";
    $select_id_mb->configure(-text => "Choose");

    
    set_domain_id_num(1); 
    set_graphics_type(0); 
    #if($model_name eq "LAPS"){ 
    if($num_nests < 2){ 
       $id_label->gridForget();
       $select_id_mb->gridForget();
    } else {
       $id_label->grid();
       $select_id_mb->grid();
    }

    # Check for localization files. 
    if(!-e $static_select) {
       # Don't continue.
       $select_id_mb->configure(-bg => $bg_color);
       set_button_state(0,$select_id_mb);
       return;

    } else {
       # Continue and let user select values.
       $select_id_mb->configure(-bg => $colorY2);
       set_button_state(1,$select_id_mb);
    }
    
    # Force a choice for the user the first time. 2/14/06
    if(-e $static_select) {
       set_domain_id_num(1); # 'd01'
       set_graphics_type(0); # 'All'
    }
}

# ----------------------------------
# set_domain_id_num 
#
# $d_num is used throughout and 
# $d_index is used in create_graphics(). 
# ----------------------------------
sub set_domain_id_num {
   ($d_num)="d0@_";
   ($d_index)=@_;
    $select_id_mb->configure(-text => $d_num, -bg => $bg_color);

    set_button_state(1,$view_graphics_mb);
    if ($d_type eq "Choose") { 
	$view_graphics_mb->configure(-bg => $colorY2); 
	$hint_msg="Select a graphics type.";
        set_graphics_but_state(0); 
    } else {
        $hint_msg="Create, View, or Delete $d_num graphic images from localization output.";
        if($model_name eq "LAPS"){
           $hint_msg="Create, View, or Delete graphic images from localization output.";}
    
        # Continue if static_select exists.
        set_graphics_action_state();
    }
}

sub set_graphics_type {
    my ($iii)=@_[0];
    $d_type=$display_choice[$iii];
    $g_type=$render_choice[$iii];
    $view_graphics_mb->configure(-bg => $bg_color);
    set_graphics_action_state();
    $hint_msg="Create, View, or Delete $d_num graphic images from localization output.";
    if($model_name eq "LAPS"){
       $hint_msg="Create, View, or Delete graphic images from localization output.";}
}

# ----------------------------------
# set_graphics_action_state 
#
# If NCAR ncgm file exist, change the state of buttons.
# ----------------------------------
sub set_graphics_action_state {

    # Check for presence of domain static file.
    #if ($g_type ne "" && $d_num ne "") { 
    if(-e $static_select) {
       set_button_state(1,$b_create_graphics); 
       set_button_state(1,$b_save_graphics_cmd);
    } else {
       # Failed to find file.
       $hint_msg="Cannot create graphics images without localization NetCDF output file, 
$static_select.";
       return(1);
    }

    # -- Create image filename.
    $ncg_file="$dataroot_select/$domain_select/static/$g_type.$d_num.ncgm";
    if($model_name eq "LAPS"){
     $ncg_file="$dataroot_select/$domain_select/static/$g_type.ncgm";
    }

    # Check for presence of NCAR ncgm file(s).
    # They must exist in order to view or delete them.
    if(-e $ncg_file) {
       set_button_state(1,$b_view_graphics);
       set_button_state(1,$b_delete_graphics);
       $b_view_graphics->configure(-bg => $colorY2); 
       $b_create_graphics->configure(-bg => $bg_color); 
    } else {
       set_button_state(0,$b_view_graphics);
       set_button_state(0,$b_delete_graphics);
       $b_create_graphics->configure(-bg => $colorY2); 
       $b_view_graphics->configure(-bg => $bg_color); 
    }

    # Success
    return(0);
}

# ----------------------------------
# set_graphics_but_state 
#
# Only If NCARG && NCL_COMMAND are set, then 
# set graphics butttons state.
# ----------------------------------
sub set_graphics_but_state {
   my ($g_state)=@_;

   if (defined $ROOT_NCARG) { 
      set_button_state($g_state,$b_create_graphics);
      set_button_state($g_state,$b_view_graphics);
      set_button_state($g_state,$b_delete_graphics);
      $b_create_graphics->configure(-bg => $bg_color); 
      $b_view_graphics->configure(-bg => $bg_color); 
      $b_delete_graphics->configure(-bg => $bg_color); 
   }

}

# ----------------------------------
# check_staticFile 
#
# ----------------------------------
sub check_staticFile {

    # Check for presence of static file.
    if(!-e $static_select) {
       info_dbox("File Not Found", "No NetCDF file found, the domain
has not been successfully localized.");
       # Fail
       return(1);
    }

    # Success
    return(0);
}

# ----------------------------------
# create_graphics 
#
# Only If NCARG && NCL_COMMAND are set, then 
# create a domain ncgm output file.
# ----------------------------------
sub create_graphics {

    # -- Check for presence of static file.
    if(check_staticFile()) { return(1); }

    # -- Check for presence of image file.
    if(-e "$ncg_file") {
       my $ans = yesno_dbox("Graphics Exists", 
"Geographical data images already exist for
$ncg_file.
Re-create $d_num graphics anyway?","No"); 
       if($ans eq "No"){ return(1); } 
    }

    if ($d_index eq "") {
       # Fail
       info_dbox("Domain Not Selected", "    Choose domain.                    ");
       return(1);
    }

    # -- Create ncl graphics.
    watch_cursor(1);
    my $ncl_dir="$ROOT_INSTALL/graphics/ncl";
    my $img_cmd=
    "$ncl_dir/generate_images.pl -domain=$dataroot_select/$domain_select -grid $d_index";
    if($model_name eq "LAPS"){
       $ncl_dir="$ROOT_INSTALL/etc";
       $img_cmd= "$ncl_dir/generate_images.pl -domain=$dataroot_select/$domain_select";
    }

    # -- Create single image or mega image.
    # -------------------------------------
    if($g_type eq "meta") {
       $img_cmd="$img_cmd -mode=$g_type";
    } else {
       $img_cmd="$img_cmd -type=$g_type";
    }


    $hint_msg="Creating $d_num images from localization output. This will take a FEW minutes...";
    # -- Update idletasks.
    $mw->update();
    $mw->idletasks();

    chdir $ncl_dir     or warn "Cannot change dir to $ncl_dir: $!";
    print "Create ncl graphics command: $img_cmd\n";
    my $my_catch=Tk::catch { `$img_cmd`; };
    watch_cursor(0);

    $hint_msg="";
    if ($my_catch eq "") { 
       # -- There is a major problem with the executable file.
       fail_dbox("Script Problem", 
        "Problem with data script.\n\nLook at log file $logFile for more information."); 
       run_sys::run_sys("$img_cmd",1);
   
       return(1);

    } elsif ($my_catch =~ m/fail|error|not work/i ) { 
       fail_dbox("Runtime Error", 
        "Error running command: $img_cmd: $my_catch");
       print "Error running command: $img_cmd: $my_catch"; # for log_file, too.
       return(1);
    }


    # -- Display graphics.
    if(-e $ncg_file) {
       # Success
       view_graphics();
       sleep(1);
       
       # Enable view and delete buttons.
       set_graphics_but_state(1); 
       return(0);
    } else {
       # Fail
       return(1);
    }
}

# ----------------------------------
# delete_graphics 
#
# ----------------------------------
sub delete_graphics {

  # Check for presence of image file.
  if(-e $ncg_file) {
     my $ans=okcancel_dbox("Delete Graphics", 
"Geographical data images will be deleted for $d_num:
$ncg_file."); 
     if($ans eq "Cancel"){ return(1); } 
  } else {
     fail_dbox("File Not Found", 
     "File '$ncg_file'\nnot found. Cannot view image graphics."); 
     return(1);
  }

  # Delete file(s).
  if ($g_type eq "meta") {
    my $i=0;
    until ($i > scalar(@display_choice)-1 ) {
       $ncg_file=
       "$dataroot_select/$domain_select/static/$render_choice[$i].$d_num.ncgm";
       if ($model_name eq 'LAPS') {
          $ncg_file=
          "$dataroot_select/$domain_select/static/$render_choice[$i].ncgm";
       }
       if(-e $ncg_file) { unlink $ncg_file;
print "removing $ncg_file\n"; }
       $i++;
    }

  } else {
    unlink $ncg_file;
  }

  # Enable/disable graphics action buttons.
  set_graphics_action_state();

  return(0);
}

# ----------------------------------
# view_graphics
#
# idt - NCARG image display tool.
# ----------------------------------
sub view_graphics {

  my $ncg_cmd="$ROOT_NCARG/bin/idt -bg gray85";
  print "Create ncl graphics command: $ncg_cmd\n";
  if (!-e $ncg_cmd ) {
     info_dbox("Problem with path", "'$ROOT_NCARG/bin/idt' does not exist.\n
Variable 'NCARG_ROOT' is set to '$ROOT_NCARG' 
in file '$ROOT_INSTALL/config_paths'
this may need to be reset.\n");
  }

  # Display graphics.
  if(-e $ncg_file) {
    #my $my_result=system ("$ncg_cmd $ncg_file &");
    #if($my_result != 0) {
    #   info_dbox("Problem with NCARG.", "Cannot run '$ncg_cmd $ncg_file'.");
    #}

    my $my_catch=Tk::catch { `$ncg_cmd $ncg_file &`; };

    if ($my_catch eq "") { 
       # -- There is a major problem with the executable file.
       run_sys::run_sys("$ncg_cmd",1);
       fail_dbox("Script Problem", 
        "Problem with image display command \n\nDoes $ncg_cmd exist? Look at log file $logFile for more information."); 
   
       return(1);

    } elsif ($my_catch =~ m/fail|error|not work/i ) { 
       fail_dbox("Runtime Error", 
        "Error running command: $ncg_cmd: $my_catch");
       print "Error running command: $ncg_cmd: $my_catch"; # for log_file, too.
       return(1);
    }

  } else {
     fail_dbox("File Not Found", 
     "File '$ncg_file'\nnot found. Cannot view image graphics."); 
     return(1);
  }

}

# ----------------------------------
# write_graphics_cmd 
#
# Write the command line to a file in $ROOT_INSTALL/user_scripts.
# ----------------------------------

sub write_graphics_cmd {

    # -- Create ncl graphics.
    my $ncl_dir="$ROOT_INSTALL/graphics/ncl";
    my $img_cmd=
    "$ncl_dir/generate_images.pl -domain=$dataroot_select/$domain_select -grid $d_index";
    if($model_name eq "LAPS"){
       $ncl_dir="$ROOT_INSTALL/etc";
       $img_cmd= "$ncl_dir/generate_images.pl -domain=$dataroot_select/$domain_select";
    }

    # -- Create single image or mega image.
    # -------------------------------------
    if($g_type eq "meta") {
       $img_cmd="$img_cmd -mode=$g_type";
    } else {
       $img_cmd="$img_cmd -type=$g_type";
    }


  # Create directory, if necessary.
  my $user_scripts="$ROOT_INSTALL/user_scripts";
  my $user_scripts_file="$user_scripts/generate_images_$domain_select.sh";

  if (!-d $user_scripts) {
      mkdir $user_scripts, 0777
      or (&info_dbox("Make Directory Error", "Cannot make $user_scripts.
Try to run chmod in another terminal window, then press 'Write File' again.")); 
  }

  # Write command.
  open(ICMD,">$user_scripts_file") or &fail_dbox("Failure","Can't open file: $user_scripts_file."), return(1);
  print ICMD "# The following command runs generate_images.pl to process input data.\n";
  print ICMD "$img_cmd\n";
  close(ICMD);
  chmod 0775, $user_scripts_file;

  # Let user know file was written.
  &info_dbox("Wrote File", "Wrote file 'generate_images_$domain_select.sh' located in directory\n$user_scripts containing the command:\n\n$img_cmd"); 

  return(0);
}


# Not currently used. There is only ncgm image output.

# ----------------------------------
# create_iviewer 
#
# Create image viewer window to view output images 
# from the Domain Localization.
# ----------------------------------
sub create_iviewer { 
    
   my $WIDTH =850;
   my $HEIGHT=850;
   $iviewer_mw = $mw->Toplevel(-title => "Image Viewer",
                       -bg => '#d9d9d9',
                       -visual => 'truecolor');
   $iviewer_mw->geometry("-100+100");
 
   $i_can_scroll = $iviewer_mw->Scrolled('Canvas', 
        -width          => 0.99 * $WIDTH, 
        -height         => 0.99 * $HEIGHT,
        -relief         => 'sunken',
        -borderwidth    => 2,
        -scrollbars     => 'ose',
        -background     => 'black',
        -scrollregion   => [0, 0, $WIDTH, $HEIGHT],
    )->pack(-fill => 'both');
    $i_can_scroll->Subwidget('xscrollbar')->configure(-width => 10);
    $i_can_scroll->Subwidget('yscrollbar')->configure(-width => 10);


   $i_can = $i_can_scroll->Subwidget("canvas");
   $i_photo = $i_can->Photo;
   $iviewer_mw->Button(-text => "Close", -command => [\&forget_iviewer])->pack;

}

# ----------------------------------
# forget_iviewer 
#
# Destroy the iviewer window.
# ----------------------------------
sub forget_iviewer { 
    $iviewer_mw->destroy if Tk::Exists($iviewer_mw);
    #$view_graphics_mb->configure(-text => 'Choose variable');
}

# ----------------------------------
# view_loc_graphics
#
# Load gif or png format images of localization by
# first calling ncl (NCAR Graphics Command Line) to 
# generate the image.
# ----------------------------------
sub view_loc_graphics {
    return;
    my($i_choice)=@_;

    watch_cursor(1);

    # Change the label on the selector menubar.
    $view_graphics_mb->configure(-text => $display_choice[$i_choice]);
    #print "Calling ncl to display >$render_choice[$i_choice].ncgm<.\n";

    # Create image destroy viewer, else clear out old image.
    if (!Tk::Exists($iviewer_mw)) {
        create_iviewer();
    } else {
        # Clean up
        if ($i_can->find('withtag', 'geog_image')) {
            $i_can->delete("geog_image");
            $i_photo->blank;
        }
    }

    # Load image in viewer.
    my $img="$dataroot_select/$domain_select/static/$render_choice[$i_choice].png";
    $i_photo->read($img, -shrink);
    $i_can->createImage( 0, 0, 
           -anchor => 'nw',
           -image => $i_photo, 
           -tags  => "geog_image");

    # Change viewer title.
    $iviewer_mw->configure(-title => "Image Viewer - $img", -fg => $colorW);

    $view_graphics_mb->focus;
    watch_cursor(0);
}

### Return 1 to the calling use statement ###
1;
