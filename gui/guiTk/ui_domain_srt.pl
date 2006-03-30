
# ---------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# ui_domain_srt.pl
# 	Domain Selection set-up and relocation tool.
#       graphical user interface routine where most frames, 
#	notebook panels, widgets, and buttons are created.
#
# Author: Paula McCaslin   30 Jan 2003  Original Version
# ---------------------------------------------------------------------------

#use warnings;
#use strict;
use strict 'subs';
use strict 'refs';

# ----------------------------------------
# create_domain_mw 
#
# Create Domain Selection Tool.
# ----------------------------------------
sub create_domain_mw {
  no strict 'refs';

  #-- Domain Set Up & Relocatable Frame (dom_frame) -- 
  $domain_mw=$sys_frm->Frame();
  $domain_mw->Frame(-height => 6)->pack();
 
  $nb_frame =$domain_mw->Frame()
                  ->pack(-expand => 1, -fill => 'both');

  $dom_frame =$domain_mw->Frame(-relief => 'flat', -bd => 2)
                  ->pack(-expand => 0, -fill => 'x');

  @panel_tag= ('',
               'Domain',
               'Horizontal Grid',
               'Vertical Grid',
               'Namelist Parms',
               'Localize Domain',
               'Domain Graphics',
               );

  $ipanel_min=1;
  my $tabpad=20;

  if (defined $ROOT_NCARG) { 
     $ipanel_max=scalar(@panel_tag)-1;
  } else {
     # Count should not include panel7.
     $ipanel_max=scalar(@panel_tag)-2;
  }


  #---- Create domain noteBook & panels ----
  $nb = $nb_frame->NoteBook(-disabledforeground => 'gray53', 
                           #-inactivebackground => $bg_color,
                           ) #-tabpadx => $tabpad)
         ->pack(-expand => 1, -fill => 'both', 
                -padx => 0, -pady => 0);


     $ival=1;
     $panel1 = $nb->add($ival,
                       -label => $panel_tag[$ival],
                       -state => 'disabled',
                       -raisecmd => [\&set_panel_number, $ival]);
     $ival++;
     $panel2 = $nb->add($ival,
                       -label => $panel_tag[$ival],
                       -state => 'disabled',
                       -raisecmd => [\&set_panel_number, $ival]);
     $ival++;
     $panel3 = $nb->add($ival,
                       -label => $panel_tag[$ival],
                       -state => 'disabled',
                       -raisecmd => [\&set_panel_number, $ival]);
     $ival++;
     $panel4 = $nb->add($ival,
                       -label => $panel_tag[$ival],
                       -state => 'disabled',
                       -raisecmd => [\&set_panel_number, $ival]);
     $ival++;
     $panel5 = $nb->add($ival,
                       -label => $panel_tag[$ival],
                       -state => 'disabled',
                       -raisecmd => [\&set_panel_number, $ival]);

     if (defined $ROOT_NCARG) { 
        $ival++;
        $panel6 = $nb->add($ival,
                       -label => $panel_tag[$ival],
                       -state => 'disabled',
                       -raisecmd => [\&set_panel_number, $ival]);
     }


  if (0) {
  # Add balloon messages to the tabs to help inform the users.
  $balloon->attach($nb,
              -balloonmsg => \$msg,
              -motioncommand => sub {
                  my($nb,$x,$y)=@ARG;
                  $x-=$nb->rootx; # adjust screen to widget coords
                  $y-=$nb->rooty;
                  my $tab = $nb->identify($x,$y);
                  if (defined $tab) {
                      $num="tpanel$tab";
                     #$num=join("tpanel", $tab);
                      $msg="Set variables related to the ${$num}";
                      0; # show balloon
                  } else {     1; # cancel balloon
                  }
              });
  }

  #--- Instantiate all the domain panels and their respective widgets.
  create_next_previous();
  create_domainSelectors_panel();
  create_horizGrid_panel();
  create_vertGrid_panel();
  create_locParms_panel();
  create_localizeDomain_panel();
  if (defined $ROOT_NCARG) { create_domainGraphics_panel(); }
}

# ----------------------------------
# create_next_previous;
#
# Create Back and Next navigation buttons.
# ----------------------------------
sub create_next_previous { 

  #--------- Action Selectors ----------
  $w_cntrls=$dom_frame->Frame() ->pack(); 

  $cntl_back=$w_cntrls->Button(-text => "<Back",
                      -command => [\&raise_panel, -1],
                      -width => 7,
                      -state => 'disabled',
                      -relief => 'groove') 
                         ->pack(-side => 'left', -anchor => 'n', 
                                -padx => 1, -pady => 3);
  $cntl_next=$w_cntrls->Button(-text => "Next>",
                      -command => [\&raise_panel, 1],
                      -width => 7,
                      -state => 'disabled',
                      -relief => 'groove') 
                         ->pack(-side => 'left', -anchor => 'n',
                                -padx => 1, -pady => 3);

  # Add balloon message to inform user
  $balloon->attach($cntl_back, -msg => "Go back to the previous panel."); 
  $balloon->attach($cntl_next, -msg => "Advance to the next panel\nonce all values are selected."); 
}

# ----------------------------------
# create_domainSelectors_panel 
#
# 'Domain' panel.   (d_frame)
# Choose what you want to do:
#      "Create new domain",
#      "Load existing domain",
#      "Copy existing domain", or
#      "Delete existing domain",
# ----------------------------------
sub create_domainSelectors_panel {

  $d_frame=$panel1->Frame()
                   ->pack(-expand => 1, -fill => 'both',
                          -padx => 10, -pady => 6, -anchor => 'nw');

  $d_enter=$d_frame->Frame(-relief => 'groove', -bd => 2)
                    ->pack(-expand => 0, -fill => 'y',
                           -padx => 0, -pady => 3, -anchor => 'nw');

  $d_exists=$d_frame->Frame();

#--------- Frame: Choose domain mode.

  #--------- Grid Information -------
  $d_enter->gridRowconfigure(0, -minsize => 13);
  $d_enter->gridRowconfigure(3, -minsize =>  7);
  $d_enter->gridRowconfigure(6, -minsize =>  7);
  $d_enter->gridRowconfigure(9, -minsize => 13);
  $d_enter->gridColumnconfigure(0, -minsize => 30);
  $d_enter->gridColumnconfigure(2, -minsize => 30);
  $d_enter->gridColumnconfigure(4, -minsize =>  5);
  $d_enter->gridColumnconfigure(8, -minsize => 30);


  # --- Choose mode.
  $row_count=1;
  $d_enter->Label(-text => "  Choose what you want to do",
                              -font => "Helvetica -10")
               ->grid(-row => $row_count, -column => 1, 
                      -columnspan => 1, -sticky => 'sw');

  $row_count++;
  @domain_choices = ("Choose Mode",
            "Create new domain",
            "Load existing domain",
            "Copy existing domain",
            "Delete existing domain",
  );

  $domain_mb=$d_enter->Menubutton(-indicator => 1,
                           -tearoff => 0,
                           -relief => "raised",
                           -bd => 2,
                           -width => 23,
                           -anchor => 'c',
                           -fg => $colorN,
                           -activebackground => $update_color,
                           -text => $domain_choices[0])
             ->grid( -row => $row_count,
                     -column => 1,
                     -sticky => 'w',
                     -padx => 5);
  $i=1;
  until ($i>4) {
      $domain_mb->command(-label => $domain_choices[$i],
                          -command => [\&switch_domains, $i]);
      $i++;
  }
  $domain_mb->focus;

  # --- Name of domain.
  $row_count=1;
  $domain_lab=$d_enter->Label(-text => "Enter domain name (e.g. \"Alaska\")",
                              -font => "Helvetica -10")
               ->grid(-row => $row_count, -column => 3, 
                      -columnspan => 1, -sticky => 'sw');

  $row_count++;
  $domain_ent=$d_enter->Entry(-width => 20,
                            -justify => 'left',
                            -textvariable => \$domain_select)
               ->grid( -row => $row_count, -column => 3, -sticky => 'w');


  # --- Dataroot
  $row_count++;
  $dataroot_lab=$d_enter->Label(-text => 
                    "Accept default or Enter path to your data (i.e. DATAROOT)",
                    -font => "Helvetica -10")
               ->grid(-row => $row_count, -column => 3, 
                      -columnspan => 2, -sticky => 'sw');
  $row_count++;
  $dataroot_ent=$d_enter->Entry(-width => 55,
                            -justify => 'left',
                            -textvariable => \$dataroot_select)
               ->grid( -row => $row_count, -column => 3, 
                       -columnspan => 2, -sticky => 'w');

  $dataroot_ent->bind('<Return>', [\&add_checkmark]);

#  # Highlight the domain name, indicating it should be changed.
#  my $numindex=($domain_ent->index('end') - length $domain_select);
#  $dataroot_ent->selectionRange($numindex, 'end');
    
  $dataroot_browse_button=$d_enter->Button(-text => "Browse...", 
                            -width => 7,
                            #-disabledforeground => $bg_color,
                            -command => [\&browse_dataroot])
               ->grid( -row => $row_count, -column => 5, -sticky => 'sw');


  # --- Simulation description label
  $row_count++;
  $domaindesc_lab=$d_enter->Label(-text => "Simulation description",
                                  -font => "Helvetica -10")
               ->grid(-row => $row_count, -column => 3, 
                      -columnspan => 1, -sticky => 'sw');

  # --- User description label
  $userdesc_lab=$d_enter->Label(-text => "User",
                                -font => "Helvetica -10");

  # --- Simulation description entry
  $domaindesc_ent=$d_enter->Entry(-width => 40,
                            -justify => 'left',
                            -textvariable => \$nl_var{SIMULATION_NAME}[0]);

  # --- User description entry
  $userdesc_ent=$d_enter->Entry(-width => 25, 
                            -justify => 'left',
                            -textvariable => \$nl_var{USER_DESC}[0]);

 if ($model_name eq "WRF") {

  $userdesc_lab->grid(-row => $row_count, -column => 4, 
                      -columnspan => 2, -sticky => 'sw');

  $row_count++;
  $domaindesc_ent->grid(-row => $row_count, -column => 3, 
                      -columnspan => 1, -sticky => 'w');

  $userdesc_ent->grid(-row => $row_count, -column => 4, 
                      -columnspan => 2, -sticky => 'w', -padx => 3);

 } elsif ($model_name eq "LAPS") {

  $row_count++;
  $domaindesc_ent->grid(-row => $row_count, -column => 3, 
                      -columnspan => 3, -sticky => 'w');

  $domaindesc_ent->configure(-width => 65,
                             -textvariable => \$nl_var{C80_DESCRIPTION}[0]);

 }

#--------- Frame: Existing Domains --------

#  # --- Template dir 
#  my $hname=hostname(); 
#  if ($hname ne "") {
#     $hname="Domains listed in $hname:$ROOT_TEMPLATES";
#  } else { 
#     $hname="Domains listed in $ROOT_TEMPLATES";
#  }
#  $d_exists->Label(-text => "$hname",
#                   -fg => $colorN)
#                  #-font => "Helvetica -10")
#           ->pack( -expand => 0, -fill => 'y',
#                   -anchor => 'w');

  # --- Simulation description entry

  # Domains were found in data directories 
  # then stored in array @domain_found.

  $domain_lb=$d_exists->Scrolled('HList', 
                         # HList
                         (-indicator => 1,
                          -separator => ",",
                          #-command => \&load_domainSelection,    # double click
                          #-browsecmd => \&load_domainSelection, # single click
                         ),
                         # TList
                         #(-height => 100),
                         #-height => 25,
                         -width => 25,
                         -scrollbars => "ose") 
              ->pack(-expand => 0, -fill => 'y',
                     -side => 'left', -anchor => 'w');
  $domain_lb->Subwidget('xscrollbar')->configure(-width => 10);
  $domain_lb->Subwidget('yscrollbar')->configure(-width => 10);

  my $hname=hostname(); 
  $domain_balloon=$mw->Balloon(-background => "#ffffbb",
                        -initwait => 500,
                        -balloonposition => 'mouse');
  $domain_balloon->attach($domain_lb, -msg => 
"List of domains found in:
$hname:$ROOT_TEMPLATES;

Check mark indicates a localized domain.");

  #--------- Bind domain selections --------
  $domain_lb->bind('<ButtonRelease-1>', [\&load_domainSelection]);

  $preview_can=$d_exists->Scrolled('Canvas', 
                                  (-bg => $bg_canvas, 
                                   -width => 450,
                                   -relief => 'groove', -bd => 2),
                        -scrollbars => 'ose')
              ->pack(-expand => 1, -fill => 'both',  
                     -side => 'left', -anchor => 'w');
  $preview_can->Subwidget('xscrollbar')->configure(-width => 10);
  $preview_can->Subwidget('yscrollbar')->configure(-width => 10);

  $p_can=$preview_can->Subwidget('canvas');

  # Preview of existing domain image.
  $p_can->createImage(0, 0, -anchor => 'nw',
                            -tags => ["preview_img", "ptags"]);

  $p_label=$preview_can->Label(-textvariable => \$etext,
                   -fg => $colorW, 
                   -bg => $bg_canvas, 
                   -font => $error_font)
           ->place(-x => 10, -y => 10);

  $p_photo=$p_can->Photo();
}

# ----------------------------------
# browse_dataroot
#
# Browse button was pressed. Call browseDialog with initial
# arguments for setting ENV variable directory path values.
# ----------------------------------
sub browse_dataroot {
    no strict 'refs';


    # Check for a good directory path.
    if (!-d $dataroot_select) {
       my $ans = info_dbox("Dataroot Problem", 
                 "Directory '$dataroot_select' does not exist!"); 
       if($ans eq "Ok"){ return;}
    }
    
    my $my_dir=$dataroot_select;
    # Browse for directory to set a new path.
    $browseDialog=$mw->FileDialog(-Title =>"Select data output directory",
 		                  -Create => 0,
                                  -PathEntryLabel => "Data Root:",
                                  -FileEntryLabel => "Namelist:",
                                  -Path => $my_dir,
                                  -File => $namelist_arg,
                                  ##-FPat => '*.*',
 	 	                  -ShowAll => 'NO');
    $browseDialog->Show(-SelDir => 1, -Horiz => 1);
                       #-DisableShowAll => 1 );

    # Check again for a good directory path.
    $my_dir = $browseDialog->{Configure}{-Path};
    if (!-d $my_dir) {
       my $ans = yesno_dbox("Dataroot Problem", 
                 "Directory '$my_dir' does not exist!\nAccept anyway?","No"); 
       if($ans eq "Yes"){ return; } 
    } 

    #Strip off backslash at end of var
    $my_dir =~ s/\/$//  if (substr($my_dir, -1, 1) eq '/');
    $dataroot_select=$my_dir;

    # Destroy when finished.
    $browseDialog->destroy if Tk::Exists($browseDialog);
}

# ----------------------------------
# set_domain_widget_state
#
# Set state of domain entry widgets based on selected domain mode.
# ----------------------------------
sub set_domain_widget_state {
   my ($arg)=@ARG;


   if ($arg == 1) {
        
      # Enable widgets.
      $domain_lab->configure  (-fg => $normal_color,
                               -text => "Enter domain name (e.g. \"Alaska\")");
      $dataroot_lab->configure(-fg => $normal_color,
         -text => "Accept default or Enter path to your data (i.e. DATAROOT)");
      $domaindesc_lab->configure(-fg => $normal_color);
      $userdesc_lab->configure(-fg => $normal_color);
      set_button_state(1,$dataroot_browse_button);
      set_button_list_state(1,$domain_ent,$dataroot_ent,$domaindesc_ent,$userdesc_ent);

   } elsif ($arg == 0) {

      # Disable widgets.
      $domain_lab->configure  (-fg => $disabled_color, -text => "Domain name");
      $dataroot_lab->configure(-fg => $disabled_color, -text => "Path to data");
      $domaindesc_lab->configure(-fg => $disabled_color);
      $userdesc_lab->configure(-fg => $disabled_color);
      set_button_state(0,$dataroot_browse_button);
      set_button_list_state(0,$domain_ent,$dataroot_ent,$domaindesc_ent,$userdesc_ent);

   } elsif ($arg == 2) {

      # Disable widgets.
      $domain_lab->configure  (-fg => $disabled_color, -text => "Domain name"),
      set_button_list_state(0,$domain_ent);

      # Enable widgets.
      $dataroot_lab->configure(-fg => $normal_color,
         -text => "Accept default or Enter path to your data (i.e. DATAROOT)");
      $domaindesc_lab->configure(-fg => $normal_color);
      $userdesc_lab->configure(-fg => $normal_color);
      set_button_state(1,$dataroot_browse_button);
      set_button_list_state(1,$dataroot_ent,$domaindesc_ent,$userdesc_ent);
 
   }
}


# --- Subroutines to switch between new and existing domains. ---


# -------------------------------------------
# switch_domains
#
# New_domain 1=New,  2=Existing
# Change the interface by showing or disabling the 
# list of existing domains (e.g. templates).
# Then, load a namelist into the application.
# -------------------------------------------
sub switch_domains {
    my ($idx)=@ARG;
     
    # Set $domain_mode, 0 or 1.
    $domain_mode=$idx;
     
    if ($domain_mode == 1 && $panel_index == 1) { 
         # Reset domain name only if 'Domain' panel is raised.
         $domain_select=$domain_current=$hints_label=""; 
         $dataroot_select=$ROOT_DATA;
    }

    # Change the label on the selector menubar: Create, Edit, Copy, or Delete.
    $domain_mb->configure(-text => $domain_choices[$idx],
                          -fg => $normal_color);

    # Update access to tabs.
    access_to_tabs(0);
    
    # Set state of domain entry widgets.
    $hint_msg="";
    reset_vars();
    $domain_lb->configure(-selectbackground => 'cornsilk');

    # Clear selected domain, if any.
    $domain_lb->selectionClear();
    $domain_lb->anchorClear();

    # For consistency reset fine scale editing to bbox(=1).
    $grid_val_restrict=1; $rb_edit_bbox->invoke(); #restrict_grid_var_calc();

    # Reset vertical widgets to default.
    reset_vert_wdgets();
    
    if ($domain_mode == 1){

      # New Domain =1;

      # Enable 'Next' button.
      set_button_state(1,$cntl_next);

      $d_exists->packForget();
      
      # Read install root namelist.
      load_namelist();

      # Reset 'Horizontal Grid'.
      present_grid_editor(1);

      # Important to reset these vars.
      $nl_var{XDIM}[0]=$nx_dim=$NXmax;

      # Set state of domain entry widgets.
      set_domain_widget_state(1);
      $domain_ent->focus;

    } elsif ($domain_mode >= 2){

      # Set state of domain entry widgets.
      set_domain_widget_state(0);

      # Existing Domain =2;
      $d_exists->pack(-expand => 1, -fill => 'both',
                      -padx => 0, -pady => 3, -anchor => 'nw');
   
      # Move focus.
      switch_focus_main();

      # If user presses clear, raise 'Domain' panel.
      # Otherwise, domain label will be out of sync with 
      # user selected domain.
      $panel_index=1;
      raise_panel(0);
       
      # Clear entry box values.
      $domain_select="";
      $dataroot_select="";
      $nl_var{SIMULATION_NAME}[0]=""; #WRF
      $nl_var{USER_DESC}[0]="";       #WRF
      $nl_var{C80_DESCRIPTION}[0]=""; #LAPS
   
      # Clear image.
      $etext="";
      load_image_preview(0);

      # Disable 'Next' button.
      set_button_state(0,$cntl_next);

      if ($domain_mode == 3){
         # Delete domain should highlight tan.
         $domain_lb->configure(-selectbackground => 'tan'); 

      } elsif ($domain_mode == 4){
         # Delete domain should highlight red.
         $domain_lb->configure(-selectbackground => 'red'); 
      }

    }

    # Update idletasks.
    $mw->update();
    $mw->idletasks();

    ## When the value for domain_select changes: 
    
    # Update "Interpolate Data" Tool.
    interp_data_but_state(0);
}

# ----------------------------------
# delete_domain 
#
# Use needs to confirm delete. 
# ----------------------------------
sub delete_domain {

   if($domain_select eq 'default') {
        $hint_msg= "You cannot delete the 'default' domain.";
        return(1);
   }

   my $delete_it="$dataroot_select/$domain_select"; 
   if ( -e $delete_it) {
      my $ok=okcancel_dbox("Deleting files", 
"Are you sure that you want to recursively delete $root_env_label[2]:\n$delete_it."); 
      if($ok eq "Ok"){ 
         if (-e $delete_it) { 
            my $stat=system("rm -rf $delete_it"); 
            if ($stat == 0) { 
               # Remove $domain_select's checkmark.
               if ($domain_lb->indicatorExists($domain_select) ) { 
                   $domain_lb->indicatorDelete($domain_select); }

            } else {
               fail_dbox("Problem removing directory", 
                          "Problem removing directory"); 
               $hint_msg= "Problem removing directory.";
               return(1); 
            }

         }
      } else {
        return(1);
      }
   }

   $delete_it="$ROOT_TEMPLATES/$domain_select"; 
   if ( -e $delete_it) {
      my $ok=okcancel_dbox("Deleting files", 
      "Are you sure that you want to recursively delete:\n$delete_it."); 

      if($ok eq "Ok"){ 
         if ( -e $delete_it) { 
            my $stat=system("rm -rf $delete_it"); 
            if (!$stat) { remove_domain_entry(); }
         }
      } else {
        return(1);
      };

   }

   # Success.
   return(0);
}

# ----------------------------------
# remove_domain_entry
#
# Remove $domain_select from display list - $domain_lb,
#                   and from array list - @domain_found.
# ----------------------------------
sub remove_domain_entry {

    # Remove $domain_select from display list - $domain_lb.
    $domain_lb->delete('entry',$domain_select);

    # Remove $domain_select from array list - @domain_found.
    my $offset=0;
    foreach (@domain_found) { 
       if ($domain_select eq $_) { 
          splice(@domain_found, $offset, 1);
          @dfound_previous=@domain_found;
          last; 
       }
       $offset++;
    }
}

# ----------------------------------
# load_image_preview 
#
# Called from subroutines:
#     load_domainSelection 
#     switch_domains and 
# ----------------------------------
sub load_image_preview {
   my ($arg)=@ARG;

   
   # Clear image - reduce possibility of memory leak. 
   if ($p_can->find('withtag', 'preview_img')) {
      $p_can->delete("preview_img");
      $p_photo->blank;
   }
   if ($arg eq 0) { return; } # Done.


   # Load domain image label & file.
   $p_can->createImage( 0, 0, 
           -anchor => 'nw',
           -image => $p_photo, 
           -tags  => "preview_img");

   #$p_can->itemconfigure("preview_img", -image => $p_photo);
   $p_can->configure(-scrollregion => [ $p_can->bbox("preview_img") ]);
   $p_can->yviewMoveto(0.15);
   $p_label->raise;

   #$p_can->move("preview_img", 15, 0);
#  $cen_y1=$p_can->height()/2;
#  $cen_y2=$preview_img_ymax/2;
#  $yoff=$cen_y1-$cen_y2;
#  print "  $yoff\n";
#  $p_can->yviewMoveto($yoff);
#  # Center bbox in canvas.
#  #$p_can->move("preview_img", 0, $yoff);

}

# ----------------------------------
# load_domainSelection
#
# Call sub load_domain, if success press 'Next'
# otherwise disable 'Next'.
# ----------------------------------
sub load_domainSelection {

   # If empty selection, return.
   # T.Tom if($domain_lb->info('selection') eq "") { return(1); }
   @domain_select_array=$domain_lb->info('selection');
   $domain_select_element=@domain_select_array[0];
   if($domain_select_element eq "") { return(1); }

   # Update access to tabs.
   access_to_tabs(0);
    
   # Disable 'Next' button.
   set_button_state(0,$cntl_next);

   # Update "Interpolate Data" button.
   interp_data_but_state(0);

   # Clear image.
   $etext="";
   load_image_preview(0);

   if (!load_domain()) {
      # Set state of Horizontal Grid widgets.
      if ($domain_mode == 2){ set_domain_widget_state(2); }

      # Enable 'Next' button.
      set_button_state(1,$cntl_next);

      # Update access to tabs.
      access_to_tabs(1);
      if ($domain_mode == 3) { access_to_tabs(0); }

      # determine if domain_select is localized.
      if (-e $static_select) {
         # Update "Interpolate Data" button.
         interp_data_but_state(1);

      } else {

         # Disable graphics buttons.
         set_graphics_but_state(0);
      }

   } else {
      $hints_label="";
   }

   # Fill 'Localization Parms' text window editor with namelist
   nl_text_fill();
}

# ----------------------------------
# load_domain
#
# Load existing domain template and snapshot image
# via user mouse B-1 selection.
# Then, load a namelist into the application.
# ----------------------------------
sub load_domain {
   # Domain exists.

   # Clear data root label.
   reset_vars();
   @_=$domain_lb->info('selection');
   #T.Tom if(@_[0] eq "") { return(1); } #If empty, timing prob.
   #T.Tom $domain_select=@_[0];
   @domain_select_array=$domain_lb->info('selection');
   $domain_select_element=@domain_select_array[0];
   if($domain_select_element eq "") { return(1); }
   $domain_select=$domain_select_element;

   $domain_lb->anchorClear();
   $dataroot_select="";
   $hints_label="- $domain_choices[$domain_mode]: $domain_select";

   # If problem with selection, return.
   if(defined $nl_var{$domain_select}) { # ne "") {
      info_dbox("Dataroot Error", $nl_var{$domain_select}); 
      return(1); 
   }

   # Reset 'Horizontal Grid' to a normal size.
   present_grid_editor(1);


   # Check for dataroot file (dataroot.txt).
   # -------------------------------
   my $data_path_file="$ROOT_TEMPLATES/$domain_select/dataroot.txt";

   if (-e $data_path_file) { 
      # Fill dataroot from contents of dataroot.txt file.
      $dataroot_select=`cat $data_path_file`;
      chomp($dataroot_select);
  
      # Cannot resolve contents of file.
      if (! -e $dataroot_select) { return(1); }
  
   } elsif ($domain_mode != 4) { #$domain_mode ==2,or 3.

      # Cannot find file.
      fail_dbox("Path to data Not Found", 
        "Cannot find file containing a one line entry of the directory path to data: 
         $data_path_file\nCannot continue!"); 
      if (! -e "$ROOT_TEMPLATES/$domain_select") { remove_domain_entry(); }
  
      $domain_select="";
      $dataroot_select="";
      return(1); 
   }

   # ----

   if ($domain_mode == 4) {

      # Delete domain.
      delete_domain();

      # Clear selected domain, if any.
      $domain_lb->selectionClear();
      $domain_lb->anchorClear();

      # Clear all.
      $dataroot_select="";
      $domain_select="";
      $etext="";

      # Disable 'Next' button.
      set_button_state(0,$cntl_next);

      # Success =0 but need to exit on error status.
      return(1);
   }

 
   # Check for namelist file (wrfsi.nl -or- nest7grid.parms).
   # -------------------------------
   my $nl="$ROOT_TEMPLATES/$domain_select/$namelist_arg";
   if (-e $nl) { 
      # Read template namelist and 'merge' with install namelist.
      load_namelist();
      if ($nmm) { assign_nl_vars_nmm(); }
   } else {
      info_dbox("Domain Error", "Domain namelist $nl cannot be found."); 
      return(1); 
   }

   # Check for image domain file (domain.gif).
   # Load domain image label & file.
   # -------------------------------
   $my_photo="$ROOT_TEMPLATES/$domain_select/domain.gif";
   if (-e $my_photo) {
      $etext="Image Preview";
      $p_photo->read($my_photo, -shrink);
      $preview_img_xmax=$p_photo->width();
      $preview_img_ymax=$p_photo->height();

      # Load image.
      load_image_preview(1);

   } else {
      # Clear image.
      $etext="No Image Found";
      load_image_preview(0);
   }


   # Make a bounding box from namelist (without dynamic user input).
   make_bbox_from_nlist();
 
   # For consistency reset fine scale editing to bbox(=1).
   $grid_val_restrict=1; $rb_edit_bbox->invoke(); #restrict_grid_var_calc();

   # If a file containing the original bbox exists, load bbox.
   ##load_orig_bbox();

   set_button_state(1,$b_startover);

   # Build case name from user selected entry.
   # -------------------------------
   if ($domain_mode == 3) {
      set_button_state(0,$b_startover);
  
      # Build case name from user selected entry.
      $domain_select=make_save_case($domain_select);
      # Set state of domain entry widgets.
      set_domain_widget_state(1);
      $hints_label="- $domain_choices[$domain_mode]: $domain_select";
   $hint_msg= "This domain name is '$domain_select', enter a different domain name if you would like to.";
   }

   # Automatically update localization_command.
   # -------------------------------
   if($domain_select ne $domain_previous) { localization_command(); }

   # Set static_select (DATAROOT/$domain_select/static/static.wrfsi.d01).
   # -------------------------------
   $static_select=
      "$dataroot_select/$domain_select/static/static.$window_domain_arg";

   # Set graphics vars and disable graphics buttons.
   # -------------------------------
   reset_domain_graphics_vars();


   return(0); 
}


# ----------------------------------
# make_save_case
#
# Given an existing domain, add a case number from 1 to 999.
# ----------------------------------
sub make_save_case {
  my ($domain_name) = @_;
  my ($case_num, $new_case_num, $new_domain_name, $n);

  if ($domain_name =~ m/_\d{3}$/) { 
     # Domain name ends in exactly 3 digits preceded by an
     # underscore: _nnn, strip this off.
     $n=length($domain_name);
     $domain_name=substr($domain_name,0,$n-4);
  }

  foreach (@domain_found){
    next if( m/^\./);
    if ( m/$domain_name+_/ ) { #$domain_name plus underscore

       $n=length($_);
       $case_num=substr($_,$n-3,$n);
       if($case_num < 999){
          $new_case_num=$case_num+1;
       }else{
          my $msg="Currently, system is limited to no more than 1000 cases.\nThe current count for $domain_name is $case_num.";
          my $ans = info_dbox("Too many cases.", $msg);
          if($ans eq "Ok"){ return;}
       }
    }
  }

  if(!defined $case_num || !defined $new_case_num){ $new_case_num = "001"; }
  if(length($new_case_num)==1){$new_case_num = "00".$new_case_num;}
  if(length($new_case_num)==2){$new_case_num =  "0".$new_case_num;}
  #$new_case_num=sprintf ("%.2f", $new_case_num);

  $new_domain_name = join("_", $domain_name, $new_case_num);
  return $new_domain_name;
}

#______________________________________________
#
# Controls for the Domain Selection Tool
#______________________________________________


# -------------------------------------------
# access_to_tabs
#
# The tabs need to be managed. Sometimes allowing
# access other times not.
# -------------------------------------------
sub access_to_tabs {
   my ($access)=@_;

   my $ibegin=2;
   if (!$access) { 
     if ($panel_index == 2){$ibegin=3;}
     for $i ($ibegin .. $ipanel_max) {$nb->pageconfigure($i, -state => 'disabled');}
   } else {
     for $i ($ibegin .. $ipanel_max) {$nb->pageconfigure($i, -state => 'normal');}
   }


}

# -------------------------------------------
# set_panel_number 
#
# Command that is invokes when a panel is presented 
# (either manually or via the 'Next', 'Back' buttons).
# -------------------------------------------
sub set_panel_number {
   no strict 'refs';

   # Set Application panel index every time a a new panel is selected.
   if($domain_mb->cget(-text) ne "Choose Mode") {
      $panel_index=$ARG[0];
      set_button_state(1,$cntl_back);
      set_button_state(1,$cntl_next);
      $cntl_next->focus;

      change_in_panel(1,$panel_index);
   }

   #raise_panel(0);
}

# -------------------------------------------
# raise_panel
#
# Raise the panel.
# -------------------------------------------

sub raise_panel {
   my ($change)=@_;

   # Increment/decrement panel_index.
   if($change == 1) { 
      # When moving from one panel to the next checks
      # are made and values are updated.
      $panel_index++;
      if(change_in_panel($panel_index,$panel_index) ) {$panel_index--; return;}

   } elsif($change == -1) { 
      $panel_index--;
      $hint_msg="";
   }

   #---- 

   # Enable current tab.
   $nb->pageconfigure($panel_index, -state => 'normal');

   # Raise command => set_panel_number($panel_index)
   if($panel_index <= $ipanel_max) {$nb->raise($panel_index);

   } else {
     # When done stepping forward, advance to 'Initial Data' Tool.
     $bc_tool->invoke();
   }

   #---- 

   # Disable "Back/Next" button as needed.
   if($panel_index >= $ipanel_max) { 
      # Disable "Next>" button on the last panel.
      $panel_index=$ipanel_max; 
###   set_button_state(0,$cntl_next);
###   switch_focus_main();

   } elsif($panel_index <= 1) { 
      # Disable "<Back" button on the first panel.
      $panel_index=$ipanel_min; 
      set_button_state(0,$cntl_back);
   }


}


# -------------------------------------------
# change_in_panel
#
# When moving from one panel to the next checks
# are made and values are updated.
# -------------------------------------------
sub change_in_panel {
 my ($pi_begin,$pi_end)=@ARG;
   
 if(check_domain_name()) { return(1); }

 #print "\n--> sub change_in_panel; $pi_begin thru $pi_end;\t";

 for $pi ($pi_begin .. $pi_end) {

   ##print "$pi ";
   if($pi <= 2) {

      # 'Domain' & 'Horizontal Grid'
      # ---------------------------
      $hints_label="- $domain_choices[$domain_mode]: $domain_select";

      if($globalMap == 1) {
         $hint_msg="To create and edit a new domain, use mouse button 1. For more information, press the 'i' above the map.";
      } else { 
         $hint_msg=""; 
      }

   } elsif($pi == 3) {
      #'Vertical Grid'
      # ------------------

      ##print "'Vertical Grid' ";

      ## See if 'Update Map' button is disabled.
      if($b_update->cget(-state) ne 'normal') {

         # Cannot continue unless a map exists.
         info_dbox("No Domain Box",
          "A domain box has NOT been drawn or\nloaded from an existing domain.");
         return(1);
      }

      ## See if parameters have changed.
      if($b_update_highlight) {
         if($globalMap == 1) {
            ## The map needs to be updated.
            my $ans=okcancel_dbox("Press 'Update Map'", 
"The domain map parameters have changed. Press 'Update Map' to continue.
By pressing 'Ok' the 'Update Map' button will be invoked automatically.");

            # Update the map, after calling sync_nl_vars.
            if($ans eq 'Ok'){ wrap_render_new_map(); }

        } else { 
            info_dbox("Image Not Saved",
"The domain's parameters have changed but have not been updated.
Either 'Update Map' or 'Reset Values' and step forward via 'Next'.")
         }
         return(1);
      }

      # Save info, automatically.
      if($domain_select ne $domain_previous) { localization_command(); }
      write_domain_files_and_image();

      # Update access to tabs.
      if ($domain_mode > 1) { access_to_tabs(1); }

      # Set hint_msg with 'Vertical Grid' msg.
      vert_hint_msg();

   } elsif($pi == 4) {
      #'Namelist Parms'
      # ------------------

      nl_text_readit();
      write_namelist();
      $hint_msg="Variables are used as input to gridgen_model, check paths to geographical files.
See documentation on parameters SILAVWT_PARM_WRF and TOPTWVL_PARM_WRF.";
      if ($model_name eq "LAPS" ) { 
      $hint_msg="Variables are used as input to gridgen_model, check paths to geographical files."; }

   } elsif($pi == 5) {
      #'Localize Domain'
      # ------------------

      #if (sync_geog_path()) {return(1);}
      update_force_loc_mb();
      write_namelist();
      $hint_msg="Press 'Run Localization' if you are ready to run gridgen_model and process geographical data.";

      # Don't let user move beyond 'Domain Graphics', if restrict_ui is set.
      if(!$NCARG_ROOT or !$NCL_COMMAND) {
         set_button_state(0,$cntl_next);
         switch_focus_main();
      }

   } elsif($pi == 6) {
      #'Domain Graphics'
      # ------------------

      # Set graphics vars and disable graphics buttons.
      #reset_domain_graphics_vars();
      if(!-e $static_select) {
         $hint_msg="Cannot create graphics images without localization output files ($static_select).";
      } else {
         $hint_msg="Create or view graphics of the domain's geographical data, if the domain has been localized.";
      }

      # Don't let user move beyond 'Domain Graphics', if restrict_ui is set.
      if($restrict_ui) {
         set_button_state(0,$cntl_next);
         switch_focus_main();
      }

   }

 }
 return(0);

}

# --------------------------------------
# check_domain_name 
#
# --------------------------------------
sub check_domain_name {
    
   if ($domain_select eq "") {
      # Blank name is NOT accepted.

      my $response=info_dbox("Domain Name Problem", 
      "Select domain name other than \"$domain_select\"."); 
      return(1);
   };

   if($domain_select eq $domain_current) { return(0); }

   unless ($domain_mode == 2) {  
      # Make sure new domain name is NOT found in the list of existing 
      # domain names. If found, then ask user to confirm choice.
      # -----------------------

      foreach (@domain_found) {
         if ($domain_select eq $_){
             my $resp=yesno_dbox("Domain Name Problem", 
                "Domain name '$domain_select' already exists.\nWould you like to delete '$domain_select'.", "No"); 
             if ($resp eq "Yes"){
                # Delete domain.
                if(delete_domain()){$domain_select=""; return(1);}
             } else {
                $domain_select="";
                return(1);
             }
         };
      }

   };

   # Success.
   return(0);
}


### Return 1 to the calling use statement ###
1;


