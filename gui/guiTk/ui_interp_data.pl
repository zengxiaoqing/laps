# ---------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# ui_interp_data.pl
# 	Graphical user interface routine for running wrfprep.pl
#       to interpolate data.
#
# Author: Paula McCaslin   30 May 2003  Original Version
# ---------------------------------------------------------------------------

#use warnings;
#use strict;
use strict 'subs';
use strict 'refs';

# ----------------------------------
# create_interp_data_mw
#
# Setup Model Background Files.
# ----------------------------------
sub create_interp_data_mw {

  ## Create Interp Data frame.
  $interp_mw = $sys_frm->Frame();
  $interp_mw->Frame(-height => 4)->pack();
#  $interp_mw->Label(-text => "$model_name SI Interpolate Data",
#                    -font => $italics_font,
#                    -fg => $colorN)
#                  ->pack(-expand => 0, -fill => 'x', -pady => 0);


  #---- Create NoteBook & Panels ----
  $d_nb = $interp_mw->NoteBook(-disabledforeground => 'gray53', 
                           -inactivebackground => $bg_color,
                           -tabpadx => 20,
                           )
         ->pack(-expand => 1, -fill => 'both', 
                -padx => 1, -pady => 0);

  my @panel_tag= ('',
               'Controls',
               'Script',
               );

  $d_ival=1;
  $d_panel1=$d_nb->add($d_ival, -label => $panel_tag[$d_ival],
                                -raisecmd => sub {
$hint_msg="Interpolate Data - controls the execution of wrfprep which interpolates gridded data to a specific domain.\n";
                                                 } );


  $d_ival=2;
  $d_panel2=$d_nb->add($d_ival, -label => $panel_tag[$d_ival],
                                -raisecmd => [\&sync_all_si_vars]);

  #---
  &create_initControls_panel();
  &create_interp_data_editor();
}

# ----------------------------------
# create_interp_data_editor
#
# ----------------------------------
sub create_interp_data_editor {
 
  my $interp_frame=$d_panel2->Frame()
               ->pack(-expand => 1, -fill => 'both',
                      -padx => 15, -pady => 15);

#--------- Frame: GRIB Prep --------

  $interp_data=$interp_frame->Frame(-relief => 'groove', -bd => 2)
               ->pack(-expand => 0, -fill => 'y', -anchor => 'nw', 
                      -ipadx => 33, -ipady => 35);

  #--------- Grid Information -------
  $interp_data->gridRowconfigure(6, -minsize => 20);
  $interp_data->gridRowconfigure(9, -minsize => 30);
 #$interp_data->gridRowconfigure(13, -minsize => 10);
  $interp_data->gridColumnconfigure(1, -minsize => 440);
  $interp_data->gridColumnconfigure(1, -minsize => 440);
  $interp_data->gridColumnconfigure(4, -minsize => 0);
  $interp_data->gridColumnconfigure(6, -minsize => 5);
  my $my_width=8; 
  my $len=36;
  my $hour=3;

  #---
  my $roww=1;

  $roww+=2;

  $interp_data->Label(-text => "Data files:",
                            -fg => $colorN) 
                   ->grid( -row => $roww++, -column => 1,
                           -columnspan => 5, -sticky => 'w');

  $wrfp_idfiles=$interp_data->Scrolled('Text', 
                                       (-height => 4,
                                        -width => 95,
                                        -padx => 5,
                                        -wrap => 'word',
                                        -highlightthickness => 0, 
                                        -spacing1=> 2,
                                        -selectbackground => 'gray80',
                                        #-exportselection => 1,
                                        -relief => 'groove'),
                             -takefocus => '0',
                             -scrollbars => 'oe')
              ->grid(-row => $roww++, -column => 1, -columnspan => 5, -sticky => 'w');
  $wrfp_idfiles->Subwidget('xscrollbar')->configure(-width => 10);
  $wrfp_idfiles->Subwidget('yscrollbar')->configure(-width => 10);
  $wrfp_idfiles->tagConfigure('blue_flag', -font => $bold_font,
                                           -foreground => $colorN);
  $wrfp_idfiles->tagConfigure('red_flag',  -font => $bold_font,
                                           -foreground => $colorR);
  $wrfp_idfiles->tagConfigure('orange_flag',  -font => $bold_font,
                                           -foreground => $colorO);

  $wrfprep_but_list=$interp_data->Button(-text => 'List', 
                                 -width => 8, -justify => 'right',
                                 -command => [\&find_initial_data_files])
              ->grid(-row => $roww, -column => 5, 
                     -sticky => 'e', -pady => 5, -ipady => 2);
  #---
  $roww+=2;
  $interp_data->Label(-text => "Command line Options ( wrfprep.pl --help):", 
                            -fg => $colorN) 
                   ->grid( -row => $roww++, -column => 1,
                           -columnspan => 5, -sticky => 'w');

  $wrfp_help=$interp_data->Scrolled('Text', 
                                       (-height => 7,
                                        -width => 95,
                                        -padx => 5,
                                        -wrap => 'none',
                                        -highlightthickness => 0, 
                                        -spacing1=> 2,
                                        -selectbackground => 'gray80',
                                        #-exportselection => 1,
                                        -relief => 'groove'),
                             -takefocus => '0',
                             -scrollbars => 'ose')
              ->grid(-row => $roww++, -column => 1, -columnspan => 5, -sticky => 'w');
  $wrfp_help->Subwidget('xscrollbar')->configure(-width => 10);
  $wrfp_help->Subwidget('yscrollbar')->configure(-width => 10);


  my $help_msg=`$ROOT_INSTALL/etc/wrfprep.pl -h`;
  $wrfp_help->insert('end', $help_msg);
  $wrfp_help->configure(-state => 'disabled');
  $wrfp_help->see("16.0"); # View row 12 of help cmd.

  #---
  $roww+=2;
  $interp_data->Label(-text => "Command line for wrfprep.pl:", 
                            -fg => $colorN) 
                   ->grid( -row => $roww, -column => 1, 
                           -columnspan => 5, -sticky => 'w');
  $roww++;
  $wrfp_text=$interp_data->Scrolled('Text', 
                                       (-height => 3,
                                        -width => 95,
                                        -padx => 5,
                                        -wrap => 'none',
                                        -highlightthickness => 0, 
                                        -spacing1=> 2,
                                        -exportselection => 1,
                                        ),
                             -takefocus => '0',
                             -scrollbars => 'ose')
              ->grid(-row => $roww++, -column => 1, -columnspan => 5, 
                     -sticky => 'w');
  $wrfp_text->Subwidget('xscrollbar')->configure(-width => 10);
  $wrfp_text->Subwidget('yscrollbar')->configure(-width => 10);

  $roww++;
  $wrfprep_but_run=$interp_data->Button(-text => 'Run', 
                                 -width => 8, -justify => 'right',
                                 -command => [\&run_wrfprep_pl])
              ->grid(-row => $roww, -column => 3, 
                     -sticky => 'e', -pady => 5, -ipady => 2);

  my $wrfprep_but_save=$interp_data->Button(-text => 'Save', 
                                 -width => 8, -justify => 'right',
                                 -command => [\&write_wrfprep_cmd])
              ->grid(-row => $roww, -column => 4, 
                     -sticky => 'e', -pady => 0, -ipady => 2);

  my $wrfprep_but_restore=$interp_data->Button(-text => 'Restore', 
                                 -width => 8, -justify => 'right',
                                 -command => [\&update_wrfprep_cmd])
              ->grid(-row => $roww, -column => 5, 
                     -sticky => 'e', -pady => 5, -ipady => 2);

  # Add balloon message to inform user
  $balloon->attach($wrfprep_but_list, -msg => 
"List the Initial data files created by running grib_prep.pl.
Processed files are in the (ANALPATH, LBCPATH, LSMPATH, and 
CONSTANS_PATH) directories."); 
  $balloon->attach($wrfprep_but_run, -msg => 
    "Run wrfprep.pl.\nThis may takes several minutes."); 
  $balloon->attach($wrfprep_but_save, -msg => 
    "Write the wrfprep.pl command to a script located in\n$ROOT_INSTALL/user_scripts\n to be run at a later time or added to a cronfile."); 
  $balloon->attach($wrfprep_but_restore, -msg => 
    "Restore basic wrfprep.pl command with arguments."); 

  #---
  $interp_data->Label(-text => "Configure script - wrfprep.pl")
       ->place(-x => 5, -y => -7);
}

# ----------------------------------
# update_wrfprep_cmd
#
# As the arguments are changed, update the command line.
# ----------------------------------
sub update_wrfprep_cmd {

  $wrfp_text->delete('1.0', 'end');
  $wrfp_text->tagConfigure('blue_flag', -font => $bold_font,
                                        -foreground => $colorN);

  # Create cmd using current selected_doamin.
  $d_cmd=
  "$ROOT_INSTALL/etc/wrfprep.pl -d $dataroot_select/$domain_select -f 24 -t 3 -s YYYYMMDDHH";
  $wrfp_text->insert('end', $d_cmd);
  
  # Show far right side of scrolled box.
  $wrfp_text->xviewMoveto(1.0);

  # Highlight wrfprep.pl.
  my $d_idx2=length("$ROOT_INSTALL/etc/wrfprep.pl");
  my $d_idx1=$d_idx2-length("wrfprep.pl");
  $wrfp_text->tagAdd('blue_flag', "1.$d_idx1", "1.$d_idx2");

  # Highlight YYYYMMDDHH.
  $d_idx2=length("$d_cmd");
  $d_idx1=$d_idx2-length("-s YYYMMDDHH");
  $wrfp_text->tagAdd('blue_flag', "1.$d_idx1", "1.$d_idx2");

}

# ----------------------------------
# get_wrfprep_cmd
#
# Get the command and parse the line.
# ----------------------------------
sub get_wrfprep_cmd {

  # Get command.
  my $d_cmd=$wrfp_text->get("1.0", "end");
  $d_cmd =~ s/-s YYYYMMDDHH//;

  # Test for null command line.
  #if ($d_cmd =~ m/^\s$/) {
  if ($d_cmd eq "") {
      &fail_dbox("Command Line Null", 
        "The command line has been left blank.\n Please select a data source."); 
      $wrfp_text->delete('1.0', 'end');
      return(1);
  }

  # Show far right side of scrolled box.
  $wrfp_text->xviewMoveto(1.0);

  # Highlight wrfp_prep.pl.
  my $g_idx2=length("$ROOT_INSTALL/etc/wrfp_prep.pl");
  my $g_idx1=$g_idx2-length("wrfp_prep.pl");
  $wrfp_text->tagAdd('blue_flag', "1.$g_idx1", "1.$g_idx2");

  return($d_cmd);

}

# ----------------------------------
# write_wrfprep_cmd
#
# Write the command line to a file in $ROOT_INSTALL/user_scripts.
# ----------------------------------

sub write_wrfprep_cmd {

  # Get command and model name.
  my ($d_cmd)=&get_wrfprep_cmd();
  if ($d_cmd == 1) { return; } 

  # Create directory, if necessary.
  my $user_scripts="$ROOT_INSTALL/user_scripts";
  my $user_scripts_file="$user_scripts/wrfprep_$domain_select.sh";

  if (!-d $user_scripts) {
      mkdir $user_scripts, 0777
      or (&info_dbox("Make Directory Error", "Cannot make $user_scripts.
Try to run chmod in another terminal window, then press 'Write File' again.")); 
  }

  # Write command.
  open(WCMD,">$user_scripts_file") or &fail_dbox("Failure","Can't open file: $user_scripts_file."), return(1);
  print WCMD "# The following command runs wrfprep.pl to process input data.\n";
  print WCMD "$d_cmd\n";
  close(WCMD);
  chmod 0775, $user_scripts_file;

  # Let user know file was written.
  &info_dbox("Wrote File", "Wrote file 'wrfprep_$domain_select.sh' located in directory\n$user_scripts containing the command:\n\n$d_cmd"); 

  return(0);
}

# ----------------------------------
# run_wrfprep_pl
#
# Run the wrfprep.pl.
# ----------------------------------

sub run_wrfprep_pl {

  # Get command and model name.
  my ($d_cmd)=&get_wrfprep_cmd();
  if ($d_cmd == 1) { return; } 
 
  # Application cannot currently handle switch -r.
  if ($d_cmd=~ m/-r/) { 
       my $my_ans0=&yesno_dbox("Command Line Option Error", 
       "Without an updated wrf.nl, wrfprep.pl should not run with -r switch.
See documentation for further information.\n
\t\tWould you like to continue?", "Yes");
       if($my_ans ne "No"){ return(1); } 
  }

  # Command will be executed.
  &watch_cursor(1);
  $hint_msg= "Running 'wrfprep.pl'.  This may take SEVERAL minutes. 
On success, interpolated data files are written to: $ROOT_DATA/$domain_select/siprd";
#On success, interpolated data files are written to: $dataroot_select/$domain_select/siprd";
  $mw->update();
  $mw->idletasks();

  my $my_catch=Tk::catch { `$d_cmd`; };
  &watch_cursor(0);

  if ($my_catch eq "") { 
     # There is a major problem with the executable file.
     &fail_dbox("Script Problem", "Problem with data script.\n
Look at log file $logFile for more information."); 
     &run_sys::run_sys("$d_cmd",1);
     $hint_msg=$my_catch;
     return(1);
  }

  # ---- Was script successful? ----
  my $dir="$dataroot_select/$domain_select/log";
  my @most_recent=`ls -t $dir/\*wrfprep`;
  my $fn=$most_recent[0];
  chomp ($fn);

  if (-z $fn) {
     $hint_msg="Failure: $fn has zero size.";
  } elsif (open(FN,$fn)) {
     my @fn=<FN>;
     close(FN);
     foreach(@fn){
        next if !/fail/i && !/normal/i && !/died/i && !/no match/i;
        if (m/fail|died/i){$_="$_ Look at log file $fn.";}
        if (m/match/i){chomp($_); $_="$_  Press 'Commit' if INIT_ROOT, LSM_ROOT, etc values have recently changed.\nLook at log file $fn.";}
        $hint_msg=$_;
     }
  } else {
     $hint_msg="Could not determine success of wrfprep.pl.
Look at log file $fn.";
  }


}

# ----------------------------------
# interp_data_but_state
#
# Set the button state of the Interpolate Data
# when the value for domain_select changes. 
# ----------------------------------

sub interp_data_but_state {
   my ($change)=@_;

   if ($change) {
      # Activate "Interpolate Data" tool
      # and update the wrfprep.pl command.
      &update_wrfprep_cmd();
   }

   &set_button_state($change,$id_tool);
}

# --------------------------------------
# sync_all_si_vars 
#
# Sync the directory path to sfcfiles (i.e. EXT_DATAROOT) 
# if the entered value differs from orig value.
#
# First check si_path.
# Fill 'Data files' box.
# Update namelist and DATAROOT namelist.
# --------------------------------------

sub sync_all_si_vars {
#print "\nSUB sync_all_si_vars\n";

      # First check si_path.
      # -------------------
      if ($ROOT_EXT ne $si_path) { 
         my $ans=&yesno_dbox("Update SI Path",
"The EXT_DATAROOT: $ROOT_EXT,
differs from the path to external data root: $si_path.\n 
By pressing 'Yes' the external path will become 
equal to EXT_DATAROOT=$ROOT_EXT?","No");
         if($ans eq 'Yes'){&set_path_to_si();}
      };

      # Fill 'Data files' box.
      # -------------------
      &find_initial_data_files();

      # Update DATAROOT namelist.
      # -------------------
      write_namelist();
      write_dataroot_namelist();
} 

# --------------------------------------
# localization_command_check
#
# --------------------------------------
sub localization_command_check {

      # Clear listing of data files.
      $wrfp_idfiles->configure(-state => 'normal');
      $wrfp_idfiles->delete('1.0', 'end');
      $wrfp_idfiles->configure(-state => 'disabled');

      # Confirm paths.
      if ($si_path ne $ROOT_EXT) {
         my $resp=&yesno_dbox("Path Discrepancy", 
"ANALPATH = $nl_var{ANALPATH}[0]
LBCPATH = $nl_var{LBCPATH}[0]
LSMPATH = $nl_var{LSMPATH}[0]
CONSTANTS_PATH = $nl_var{CONSTANTS_PATH}[0]\n
Would you like to set these paths to equal
EXT_DATAROOT=$ROOT_EXT?", "Yes");
         if($resp eq "Yes"){ 
            &set_path_to_si();
         } else {
            if (!-d $si_path) { 
                # Clear listing of data files.
                $wrfp_idfiles->configure(-state => 'normal');
                #$wrfp_idfiles->delete('1.0', 'end');
                $wrfp_idfiles->insert('end', " *** NO FILES ***", 'red_flag');
                $wrfp_idfiles->configure(-state => 'disabled');
                
                # Disable 'Run'. 
                &set_button_state(0,$wrfprep_but_run);

                return(1); 
            }

         }
      }
 
      # Enable 'Run'. 
      &set_button_state(1,$wrfprep_but_run);

      return(0); 
}

# ----------------------------------
# find_initial_data_files
#
# Find initial data files for display on Interpolate Data panel.
# ----------------------------------
sub find_initial_data_files {

      # Enable 'Run'. 
      &set_button_state(1,$wrfprep_but_run);

      # Set vars.
      my $a=$nl_var{ANALPATH}[0];
      my $b=$nl_var{LBCPATH}[0];
      my $c=$nl_var{LSMPATH}[0];
      my $d=$nl_var{CONSTANTS_PATH}[0];

      # Reduce the number of listings, if possible.
      my @plist;
      if ($a eq $b && $b eq $c && $c eq $d) {
         # All paths are identical.
         @plist=($a);

      } elsif (($a ne $b && $b eq $c && $c eq $d) ||
               ($a ne $b && $a eq $c && $c eq $d) ||
               ($a ne $b && $a eq $c && $b eq $d)) {
         @plist=($a, $b);

      } elsif (($a eq $b && $b ne $c && $c eq $d) ||
               ($a eq $b && $b ne $c && $a eq $d)) {
         @plist=($a, $c);

      } elsif (($a eq $b && $b eq $c && $c ne $d) ||
               ($a eq $b && $b eq $c && $c ne $d)) {
         @plist=($a, $d);

      } elsif ($a ne $b && $b ne $c && $c eq $d) {
         @plist=($a, $b, $c);

      } elsif ($a ne $b && $b eq $c && $c ne $d) {
         @plist=($a, $b, $d);

      } elsif ($a eq $b && $b ne $c && $c ne $d) {
         @plist=($a, $c, $d);

      } else {
         # All paths are unique.
         @plist=($a, $b, $c, $d);
      } 

      # Assume data files can be found in EXT_DATAROOT/extprd.
      if ($plist[0] eq "") {
         @plist=("$EXT_DATAROOT/extprd");
         $myhead="ANALPATH, etc";
         print "****** ASSUMPTION ******\n"; 
      }

      # List data files in text widget.
      #-------------------------------
      $wrfp_idfiles->configure(-state => 'normal');
      $wrfp_idfiles->delete('1.0', 'end');

      foreach $path (@plist) {
         $hpath = $path;
         $path =~ s/'//g;  # Strip off quotes.
         #@dfiles=`ls -t $path`;

         # Find label.
         if ($hpath eq $a) {$myhead="ANALPATH";
           if (scalar(@plist) == 1) {$myhead="All paths";}
         }elsif ($hpath eq $b) {$myhead="LBCPATH";
         }elsif ($hpath eq $c) {$myhead="LSMPATH";
         }elsif ($hpath eq $d) {$myhead="CONSTANTS_PATH";}
         $wrfp_idfiles->insert('end', "$myhead ", 'blue_flag');
         $wrfp_idfiles->insert('end', "$path\n", 'blue_flag');

         # Find files.
         if (opendir(THATDIR, "$path")) { 
            @dafiles=grep 
{$_ ne "CVS" and $_ ne "README" and $_ ne '.' and $_ ne '..'} readdir(THATDIR);
            closedir(THATDIR);
            chomp (@dafiles);

            # No files found.
            if (scalar(@dafiles) < 1) {
              $wrfp_idfiles->insert('end', " *** NO FILES ***", 'red_flag');
              # Disable 'Run'. 
              &set_button_state(0,$wrfprep_but_run);
            }

            # Files found.
            foreach (@dafiles) { $wrfp_idfiles->insert('end', "$_\t"); }
         } else { 
            # Directory problem.
            $wrfp_idfiles->insert('end', " *** INVALID DIRECTORY ***", 
            'orange_flag');
         }
         $wrfp_idfiles->insert('end', "\n\n");
      } 
      $wrfp_idfiles->delete('end - 2 chars', 'end');

      $wrfp_idfiles->configure(-state => 'disabled');
      #$wrfp_idfiles->see("end"); # View last row.
}

### Return 1 to the calling use statement ###
1;
