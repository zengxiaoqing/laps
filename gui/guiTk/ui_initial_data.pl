# ---------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# ui_initial_data.pl
# 	Graphical user interface routine to create grib_prep.pl 
#       which configures initial boundary conditions.
#
# Author: Paula McCaslin   20 May 2003  Original Version
# ---------------------------------------------------------------------------

#use warnings;
#use strict;
use strict 'subs';
use strict 'refs';

# ----------------------------------
# create_initial_data_mw
#
# Setup model background files.
# ----------------------------------
sub create_initial_data_mw {

  ## Create Input Data frame.
  $bkgnd_dat_mw = $sys_frm->Frame();
  $bkgnd_dat_mw->Frame(-height => 6)->pack();
#  $bkgnd_dat_mw->Label(-text => "$model_name SI Input Data",
#                    -font => $italics_font,
#                    -fg => $colorN)
#                  ->pack(-expand => 0, -fill => 'x', -pady => 0);

 
  #---- Create NoteBook & Panels ----
  $g_nb = $bkgnd_dat_mw->NoteBook(-disabledforeground => 'gray53', 
                           -inactivebackground => $bg_color,
                           -tabpadx => 20,
                           )
         ->pack(-expand => 1, -fill => 'both', 
                -padx => 1, -pady => 0);

  my @panel_tag= ('',
               'Sources',
               'Script',
               );

  my $g_ival=1;
  $g_panel1=$g_nb->add($g_ival, -label => $panel_tag[$g_ival],
                                -raisecmd => sub {
$hint_msg="Initial Data - controls the execution of grib_prep which decodes GRIB files.\n";
                                                 } );
  $g_ival=2;
  $g_panel2=$g_nb->add($g_ival);
  $g_nb->raise($g_ival);
  $g_nb->pageconfigure($g_ival, -label => $panel_tag[$g_ival],
                                -raisecmd => [\&write_grib_namelist]);
  $g_nb->raise(1);

  #---
  &create_mdl_source_editor();
  &create_grib_prep_editor(); 
}

# ----------------------------------
# create_mdl_source_editor
#
# ----------------------------------
sub create_mdl_source_editor {
 
  my $g1_frame=$g_panel1->Frame()
              ->pack(-expand => 1, -fill => 'both',
                     -padx => 15, -pady => 15);

  my $row_count=1;
  my $m_grib=$g1_frame->Frame(-relief => 'groove', -bd => 2)
              ->pack(-expand => 0, -fill => 'y', -anchor => 'nw', 
                     -ipadx => 35, -ipady => 20);

  $m_grib->Label(-text => "Configure namelist - grib_prep.nl  ")
# to control execution of grib_prep (which decodes GRIB files).")
       ->place(-x => 5, -y => -7);

  #--------- Grid Information -------
  $m_grib->gridRowconfigure(100, -minsize => 40);
  $m_grib->gridColumnconfigure(1, -minsize => 420);


  # --- Matrix frame.
  $grib_matrix=$m_grib->Frame()
              ->grid(-row => 0, -column => 1, -columnspan => 4, 
                     -sticky => 'nw',
                     -pady => 20);

  #--------- Grid Information -------
#  $grib_matrix->gridColumnconfigure(3, -minsize => 20);

  #--

  $row_count=0;
  my $col=1;

  $m_grib->Label( -textvariable => \$grib_prep_nl,
                  -justify => 'left',
                  -fg => $colorN) 
       ->place(-x => 200, -y => -7);

  # --- SI Data Source labels.
  my $g_src_lab=$grib_matrix->Label(-text => "GRIB source\nname\n", 
                           -justify => 'left',
                           -relief => 'flat', -width => 13,
                           -fg => $colorN) 
                   ->grid( -row => $row_count, 
                           -column => $col++, 
                           -sticky => 'w');
  my $g_vtab_lab=$grib_matrix->Label(-text => 
                                   "GRIB Vtable\nused to extract\nvariables", 
                           -justify => 'left',
                           -relief => 'flat', -width => 13,
                           -fg => $colorN) 
                   ->grid( -row => $row_count, 
                           -column => $col++, 
                           -sticky => 'w');
  my $g_path_lab=$grib_matrix->Label(-text => "Path to GRIB source\n\n", 
                           -justify => 'left',
                           -relief => 'flat', -width => 40,
                           -fg => $colorN) 
                   ->grid( -row => $row_count, 
                           -column => $col++, 
                           -sticky => 'w');
  my $g_cycle_lab=$grib_matrix->Label(-text => "Cycle - hours\nbetween runs\n", 
                           -justify => 'left',
                           -relief => 'flat', -width => 13,
                           -fg => $colorN) 
                   ->grid( -row => $row_count, 
                           -column => $col++, 
                           -sticky => 'w');
  my $g_delay_lab=$grib_matrix->Label(-text => 
                                    "Delay - hours\nafter inital\nvalid time", 
                           -justify => 'left', 
                           -relief => 'flat', -width => 13,
                           -fg => $colorN) 
                   ->grid( -row => $row_count, 
                           -column => $col++, 
                           -sticky => 'w');
  # Add balloon message to inform user
  $balloon->attach($g_src_lab, -msg => 
    "GRIB source data sets you wish to process.\nAlso, prefix for the output files."); 
  $balloon->attach($g_vtab_lab, -msg => 
    "Suffix used to obtain Vtable from the \n$ROOT_EXT/static\nto extract variables."); 
  $balloon->attach($g_path_lab, -msg => 
    "Directory to GRIB data."); 
  $balloon->attach($g_cycle_lab, -msg => 
    "The number of hours between new model runs."); 
  $balloon->attach($g_delay_lab, -msg => 
    "This defines the number of hours after\nthe initial valid time that you expect\n the entire run to be available."); 


  # --- SI Initial Data - Sources entries table.

  my $add_grib_nl_but=$m_grib->Button(-width => 8, 
                   -text => "Add",
                   -justify => 'left',
                   -command => sub {
                       &add_mdl_source_entries($matrix_rows,$matrix_rows+1);
                   })
          ->grid(  -row => 3, -column => 2, 
                   -sticky => 'e', -pady => 5, -ipady => 2);

  my $save_grib_nl_but=$m_grib->Button(-width => 8, 
                   -text => "Save",
                   -justify => 'left',
                   -command => [\&write_grib_namelist]) 
          ->grid(  -row => 3, -column => 3, 
                   -sticky => 'e', -ipady => 2);

  my $reload_grib_nl_but=$m_grib->Button(-width => 8, 
                   -text => "Reload",
                   -justify => 'left',
                   -command => [\&read_grib_namelist]) 
          ->grid(  -row => 3, -column => 4, 
                   -sticky => 'e', -ipady => 2);

  # Add balloon message to inform user
  $balloon->attach($add_grib_nl_but, -msg => 
    "Add a new line entry to add\nadditional source dataset(s)."); 
  $balloon->attach($save_grib_nl_but, -msg => 
    "Write grib_prep.nl to\n$ROOT_EXT/static.
(This is invoked automatically when 'Script' is pressed.)"); 
  $balloon->attach($reload_grib_nl_but, -msg => 
    "Reload grib_prep.nl from \n$ROOT_EXT/static."); 

}

# ----------------------------------
# add_mdl_source_entries
#
# Create model grib background file entries.
# ----------------------------------
sub add_mdl_source_entries {
  my ($my_min,$my_max)=@_;

  my $row_count=($my_min+1)*3;
  my $col=1;

  for ($ii=$my_min; $ii < $my_max; $ii++) {
     $g_entry[0][$ii]=$grib_matrix->Entry(-width => 13, 
                        -justify => 'left',
                        -textvariable => \$nl_var{SRCNAME}[$ii]) 
                ->grid( -row => $row_count,
                        -column => $col++, 
                        -sticky => 'w');
     $g_entry[1][$ii]=$grib_matrix->Entry(-width => 13, 
                        -justify => 'left',
                        -textvariable => \$nl_var{SRCVTAB}[$ii]) 
                ->grid( -row => $row_count,
                        -column => $col++, 
                        -sticky => 'w');
     $g_entry[2][$ii]=$grib_matrix->Entry(-width => 40, 
                        -justify => 'left',
                        -textvariable => \$nl_var{SRCPATH}[$ii]) 
                ->grid( -row => $row_count,
                        -column => $col++, 
                        -sticky => 'w');
     $g_entry[3][$ii]=$grib_matrix->Entry(-width => 13, 
                        -justify => 'right',
                        -textvariable => \$nl_var{SRCCYCLE}[$ii]) 
                ->grid( -row => $row_count,
                        -column => $col++, 
                        -sticky => 'w');
     $g_entry[4][$ii]=$grib_matrix->Entry(-width => 13, 
                        -justify => 'right',
                        -textvariable => \$nl_var{SRCDELAY}[$ii]) 
                ->grid( -row => $row_count,
                        -column => $col++, 
                        -sticky => 'w');
     $row_count=($ii+2)*3;
     $col=1;
  }
  if ($my_max > $matrix_rows) {$matrix_rows=$my_max;}
  if ($my_min != 0) {$g_entry[0][$my_min]->focus;} # Help user with focus.
}

# ----------------------------------
# create_grib_prep_editor 
#
# Create the girb executable command.
# ----------------------------------
sub create_grib_prep_editor {

  my $g2_frame=$g_panel2->Frame()
              ->pack(-expand => 1, -fill => 'both',
                     -padx => 15, -pady => 15);

#--------- Frame: GRIB Prep --------

  my $s_grib=$g2_frame->Frame(-relief => 'groove', -bd => 2)
                   ->pack(-expand => 0, -fill => 'y', -anchor => 'nw', 
                          -ipadx => 33, -ipady => 35);

  #--------- Grid Information -------
  $s_grib->gridRowconfigure(2, -minsize => 50);
  $s_grib->gridRowconfigure(5, -minsize => 30);
  $s_grib->gridRowconfigure(10, -minsize => 35);
  $s_grib->gridColumnconfigure(1, -minsize => 520);

if(0){
  # Added feature to create a cronfile times.
  $s_grib->gridColumnconfigure(8, -minsize => 5);
  $s_grib->gridColumnconfigure(10, -minsize => 5);
  $s_grib->gridColumnconfigure(12, -minsize => 350);
  $s_grib->gridColumnconfigure(20, -minsize =>  40);

  $mm=0; $hh='4,10,16,22'; $dd='*'; $MM='*'; $DD='*';

  my $roww=1;
  my $my_font="Helvetica -10";
  my $mm_lab=$s_grib->Label(-text => "Min(0-59)", 
                            -font => $my_font,
                            -fg => $normal_color) 
                   ->grid( -row => $roww, 
                           -column => 3, 
                           -sticky => 'w');
  my $hh_lab=$s_grib->Label(-text => "Hour(0-23)", 
                            -font => $my_font,
                            -fg => $normal_color) 
                   ->grid( -row => $roww, 
                           -column => 5, 
                           -sticky => 'w');
  my $dd_lab=$s_grib->Label(-text => "Date(1-31)", 
                            -font => $my_font,
                            -fg => $normal_color) 
                   ->grid( -row => $roww, 
                           -column => 7, 
                           -sticky => 'w');
  my $MM_lab=$s_grib->Label(-text => "Month(1-12)", 
                            -font => $my_font,
                            -fg => $normal_color) 
                   ->grid( -row => $roww, 
                           -column => 9, 
                           -sticky => 'w');
  my $DD_lab=$s_grib->Label(-text => "WeekDay(0-6)", 
                            -font => $my_font,
                            -fg => $normal_color) 
                   ->grid( -row => $roww, 
                           -column => 11, 
                           -sticky => 'w');
  $roww++; 
  my $cron_sched_lab=$s_grib->Label(-text => "Crontab Schedule:", 
                            -fg => $colorN) 
                   ->grid( -row => $roww, 
                           -column => 1, 
                           -sticky => 'w');

  my $my_width=6; 
  my $mm_ent=$s_grib->Entry(-width => $my_width, 
                           -justify => 'right',
                           -textvariable => \$mm) 
                   ->grid( -row => $roww,
                           -column => 3, 
                           -sticky => 'w');
  my $hh_ent=$s_grib->Entry(-width => $my_width, 
                           -justify => 'right',
                           -textvariable => \$hh) 
                   ->grid( -row => $roww,
                           -column => 5, 
                           -sticky => 'w');
  my $dd_ent=$s_grib->Entry(-width => $my_width, 
                           -justify => 'right',
                           -textvariable => \$dd) 
                   ->grid( -row => $roww,
                           -column => 7, 
                           -sticky => 'w');
  my $MM_ent=$s_grib->Entry(-width => $my_width, 
                           -justify => 'right',
                           -textvariable => \$MM) 
                   ->grid( -row => $roww,
                           -column => 9, 
                           -sticky => 'w');
  my $DD_ent=$s_grib->Entry(-width => $my_width, 
                           -justify => 'right',
                           -textvariable => \$DD) 
                   ->grid( -row => $roww++,
                           -column => 11, 
                           -sticky => 'w');
}
  $my_width=8; 
  $len=36;
  $hour=3;

  $roww=1;
  $s_grib->Label(-text => "Choose one Data Source for each command line entry:",
                           -fg => $colorN) 
              ->grid( -row => $roww, -column => 1, -sticky => 'w'); 

  $grib_src_mb=$s_grib->Menubutton(-indicator => 1,
                           -tearoff => 0,
                           -relief => "raised",
                           -bd => 2,
                           -width => 10,
                           -anchor => 'c',
                           -fg => $normal_color,
                           -activebackground => $update_color)
              ->grid( -row => $roww++, -column => 2, 
                      -columnspan => 2, -sticky => 'e'); 

if(0){
  $s_grib->Label(-text => "Forecast Length (hr):", 
                            -fg => $colorN) 
                   ->grid( -row => $roww, -column => 1, -sticky => 'w');

  $s_grib->Entry(-width => $my_width, 
                           -justify => 'right',
                           -textvariable => \$len) 
                   ->grid( -row => $roww++, -column => 3, -sticky => 'e');

  $s_grib->Label(-text => "Forecast Interval (hr):", 
                            -fg => $colorN) 
                   ->grid( -row => $roww, -column => 1, -sticky => 'w');

  $s_grib->Entry(-width => $my_width, 
                           -justify => 'right',
                           -textvariable => \$hour) 
                   ->grid( -row => $roww++, -column => 3, -sticky => 'e');
   
  my $grib_prep_but_update=$s_grib->Button(-text => 'Create Command', 
                                 -width => 15, -justify => 'right',
                                 -command => [\&create_grib_prep_cmd])
              ->grid(-row => $roww++, -column => 3, -columnspan => 1, 
                     -sticky => 'e',
                     -ipady => 2);
}

  #---
  $roww++;
  $s_grib->Label(-text => "Command line Options ( grib_prep.pl --help):", 
                           -fg => $colorN) 
                   ->grid( -row => $roww++, -column => 1, -sticky => 'w');

  $grib_help=$s_grib->Scrolled('Text', 
                                       (-height => 10,
                                        -width => 95,
                                        -padx => 5,
                                        -wrap => 'none',
                                        -highlightthickness => 0, 
                                        -spacing1=> 2,
                                        -exportselection => 1,
                                        -selectbackground => 'gray80',
                                        -relief => 'groove'),
                             -takefocus => '0',
                             -scrollbars => 'ose')
              ->grid(-row => $roww++, -column => 1, -columnspan => 3,
                     -sticky => 'w', -pady => 3);
  $grib_help->Subwidget('xscrollbar')->configure(-width => 10);
  $grib_help->Subwidget('yscrollbar')->configure(-width => 10);

  my $help_msg=`$ROOT_INSTALL/etc/grib_prep.pl -h`;
  $grib_help->insert('end', $help_msg);
  $grib_help->configure(-state => 'disabled');
  $grib_help->see("28.0"); # View row 12 of help cmd.

  #---
  $roww++;
  $s_grib->Label(-text => "Command line to run grib_prep.pl:", 
                            -fg => $colorN) 
                   ->grid( -row => $roww++, -column => 1, -sticky => 'w');

  $grib_text=$s_grib->Scrolled('Text', 
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
              ->grid(-row => $roww++, -column => 1, -columnspan => 3,
                     -sticky => 'w', -pady => 3);
  $grib_text->Subwidget('xscrollbar')->configure(-width => 10);
  $grib_text->Subwidget('yscrollbar')->configure(-width => 10);

  $roww++;
  my $grib_prep_but_run=$s_grib->Button(-text => 'Run', 
                                 -width => 8, -justify => 'right',
                                 -command => [\&run_grib_prep_pl])
              ->grid(-row => $roww, -column => 2, 
                     -sticky => 'e', -padx => 2, -pady => 5, -ipady => 2);

  my $grib_prep_but_save=$s_grib->Button(-text => 'Save', 
                                 -width => 8, -justify => 'right',
                                 -command => [\&write_grib_prep_cmd])
              ->grid(-row => $roww, -column => 3, 
                     -sticky => 'e', -pady => 0, -ipady => 2);

  # Add balloon message to inform user
  $balloon->attach($grib_prep_but_run, -msg => 
    "Runs grib_prep.pl command with arguments."); 
  $balloon->attach($grib_prep_but_save, -msg => 
    "Writes grib_prep.pl command with arguments\nto a script located in\n$ROOT_INSTALL/user_scripts."); 

  #---
  $s_grib->Label(-text => "Configure script - grib_prep.pl")
       ->place(-x => 5, -y => -7);

}

# ----------------------------------
# update_grib_sources
#
# Fill grib sources menubutton with user's (new) input.
# ----------------------------------
sub update_grib_sources {

  &clear_grib_prep_cmd();
  $grib_src_mb->cget(-menu)->delete(0, 'end');
  for ($ii=0; $ii < $nl_var_max{SRCNAME}; $ii++) {
      $grib_src_mb->command(-label => $nl_var{SRCNAME}[$ii],
                            -command => [\&update_grib_prep_mb, $ii]);
      $_=$nl_var{SRCPATH}[$ii];
      s/\'//g;
      if (-d $_) {
         $g_entry[2][$ii]->configure(-fg => $normal_color);
      } else {
         $g_entry[2][$ii]->configure(-fg => $colorO);
      }
  }

  # Change the label on the selector menubar.
  $grib_src_mb->configure(-text => "Choose", -bg => $colorY2);
  $grib_src_mb->focus;
}

# ----------------------------------
# update_grib_prep_mb
#
# Change the label on the selector menubar.
# ----------------------------------
sub update_grib_prep_mb {
  my($idx)=@_;

  $g_src = $nl_var{SRCNAME}[$idx];
  $g_src =~ s/\'//g;      # Strip quotes from variable.
  $grib_src_mb->configure(-text => $g_src, -bg => $bg_color);

  $_=$nl_var{SRCPATH}[$idx];
  s/\'//g;
  if (! -d $_) {
     # Dir does not exist, skip creating command.
     &info_dbox("Directory Not Found", "Path to GRIB source,
$nl_var{SRCPATH}[$idx], 
does not exist for $g_src."); 
     return;

  } else {
     # Create command.
     &create_grib_prep_cmd();
  } 

}

# ----------------------------------
sub clear_grib_prep_cmd { $grib_text->delete('1.0', 'end'); }

# ----------------------------------
# create_grib_prep_cmd
#
# As the arguments are changed, update the command line.
# ----------------------------------
sub create_grib_prep_cmd {

  &clear_grib_prep_cmd();
  $grib_text->tagConfigure('blue_flag', -font => $bold_font,
                                        -foreground => $colorN);
  # Get the label on the selector menubar.
  my $g_src=$grib_src_mb->cget(-text);

  #$g_cmd="$mm $hh $dd $MM $DD $ROOT_INSTALL/etc/grib_prep.pl -c comp -d $ROOT_EXT -l $len -t $hour -q \"01:00:00\" -u wrf-laps $g_src > /dev/null 2>&1";

  # Create command.
  my $g_opts="-l $len -t $hour";
  my $g_cmd="$ROOT_INSTALL/etc/grib_prep.pl -d $ROOT_EXT $g_opts -s YYYYMMDDHH $g_src";
  $grib_text->insert('end', $g_cmd);
  
  # Show far right side of scrolled box.
  $grib_text->xviewMoveto(1.0);

  # Highlight grib_prep.pl.
  my $g_idx2=length("$ROOT_INSTALL/etc/grib_prep.pl");
  my $g_idx1=$g_idx2-length("grib_prep.pl");
  $grib_text->tagAdd('blue_flag', "1.$g_idx1", "1.$g_idx2");

  # Highlight YYYYMMDDHH.
  $g_idx2=length($g_cmd)-length($g_src);
  $g_idx1=$g_idx2-length("-s YYYMMDDHH");
  $grib_text->tagAdd('blue_flag', "1.$g_idx1", "1.$g_idx2");
}

# ----------------------------------
# get_grib_prep_cmd
#
# Get the command and parse the line.
# ----------------------------------
sub get_grib_prep_cmd {
 
  # Get command and model name.
  my $g_cmd=$grib_text->get("1.0", "end");
  $g_cmd =~ s/-s YYYYMMDDHH//;
  my $g_src=$grib_src_mb->cget(-text);

  # Test for null command line.
  #if ($g_cmd eq "") {
  if ($g_src eq "Choose" || $g_cmd =~ m/^\s$/) {
      &fail_dbox("Command Line Null", 
        "The command line has been left blank.\n Please select a data source."); 
      &clear_grib_prep_cmd();
      return(1);
  }

  # Model name is needs to be the last element in the string.
  if ($g_cmd =~ m/$g_src$/) {
  } else {
     $g_cmd =~ s/$g_src//;    # Strip out model.
     chomp ($g_cmd);         
     $g_cmd="$g_cmd $g_src";  # Add model to end of line.
     &clear_grib_prep_cmd();
     $grib_text->insert('end', $g_cmd);

     # Show far right side of scrolled box.
     $grib_text->xviewMoveto(1.0);

     # Highlight grib_prep.pl.
     my $g_idx2=length("$ROOT_INSTALL/etc/grib_prep.pl");
     my $g_idx1=$g_idx2-length("grib_prep.pl");
     $grib_text->tagAdd('blue_flag', "1.$g_idx1", "1.$g_idx2");
  }

  return($g_cmd,$g_src);

 }

# ----------------------------------
# write_grib_prep_cmd
#
# Write the command line to a file in $ROOT_INSTALL/user_scripts.
##  if (!-d "$ROOT_EXT/log")    { mkdir "$ROOT_EXT/log", 0777
##       or warn("Could not mkdir $ROOT_EXT/log\n"), return(1); }
# ----------------------------------

sub write_grib_prep_cmd {

  # Get command and model name.
  my ($g_cmd,$g_src)=&get_grib_prep_cmd();
  if ($g_cmd == 1) { return; } 

  # Create directory, if necessary.
  my $user_scripts="$ROOT_INSTALL/user_scripts";
  my $user_scripts_file="$user_scripts/grib_prep_$g_src.sh";

  if (!-d $user_scripts) {
      mkdir $user_scripts, 0777
      or (&info_dbox("Make Directory Error", "Cannot make $user_scripts.
Change write permissions on INSTALLROOT (chmod), then press 'Write File' again.")), 
         return(1);
  }

  # Write command.
  open(GCMD,">$user_scripts_file") or 
    &fail_dbox("Failure","Can't open file: $user_scripts_file."), return(1);
  print GCMD "# The following command runs grib_prep.pl to process initial input data.\n";
  print GCMD "$g_cmd\n";
  close(GCMD);
  chmod 0775, $user_scripts_file;

  # Let user know file was written.
  &info_dbox("Wrote File", "Wrote file 'grib_prep_$g_src.sh' located in directory\n$user_scripts containing the command:\n\n$g_cmd"); 

  return(0);
}
# ----------------------------------
# run_grib_prep_pl
#
# Run the grib_prep.pl.
# ----------------------------------

sub run_grib_prep_pl {

  # Get command and model name.
  my ($g_cmd,$g_src)=&get_grib_prep_cmd();
  if ($g_cmd == 1) { return; } 
      
  # Command will be executed.
  &watch_cursor(1);
  $hint_msg= "Running 'grib_prep.pl' ...processing files in: $ROOT_EXT/work/$g_src. This may take SEVERAL minutes. 
Output data and log files are written to $ROOT_EXT/: extprd/ and log/.";
  $mw->update();
  $mw->idletasks();
 
  my $my_catch=Tk::catch { `$g_cmd`; };
  &watch_cursor(0);

  if ($my_catch eq "") { 
     # There is a major problem with the executable file.
     &fail_dbox("Script Problem", "Problem with data script.\n
Look at log file $logFile for more information."); 
     &run_sys::run_sys("$g_cmd",1);
     $hint_msg=$my_catch;
     return(1);
  }


  # ---- Was external script successful? ----
  my $dir="$ROOT_EXT/log";
  my @most_recent=`ls -t $dir/\*log`;
  my $fn=$most_recent[0];
  chomp ($fn);

  if (-z $fn) {
     $hint_msg="Failure: $fn has zero size.";
  } elsif (open(FN,$fn)) {
     my @fn=<FN>;
     close(FN);
     foreach(@fn){
        next if !/not found/i && !/normal/i && !/decod/i && !/prob/i;
        if (m/not found/ || m/prob/ || m/decod/) {$_="$_ Look at log file $fn.";}
        $hint_msg=$_;
     }
  } else {
     $hint_msg="Could not determine success of grib_prep.pl. 
Look at log file $fn.";
  }

if(0){
  # ---- Was script successful? ----
  my $testDir="$ROOT_EXT/extprd";
  my @dfiles;

  if (opendir(DIR, "$testDir")) {
     @dfiles = grep { /^[$g_src]/ } readdir(DIR);
     closedir(DIR);

     foreach (@dfiles) {
        s/'//g;    # Strip off quotes.
        $wtime=(stat($_))[9];
        print "file: $_, write: $wtime\n";
     }
    
     my $d_writetime=(stat($testDir))[9];
     my $diff=$d_writetime-$wtime;
     print "dir: $testDir, write: $d_writetime\n";
     print "file: $_, write: $wtime\n";
     print "Time difference: $diff\n";


     if ($diff < 1000) { 
       #$hint_msg="Wrote initial data files to: $testDir";
       print "Wrote initial data files to: $testDir\n";
     } else {
       #$hint_msg="Possible problem with initial data files found in: $testDir";
       print "Possible problem with initial data files found in: $testDir\n";
     }
   
     &watch_cursor(0);

  } else {
     return(1);
  }
}

  #--- 
#  # Fill grib sources (grib_src_mb).
#  &update_grib_sources;

}

# ----------------------------------
# wrap_set_ext_dataroot
#
# Browse for a new EXT_DATAROOT. 
# ----------------------------------
sub wrap_set_ext_dataroot {

    $ROOT_EXT=&browse4dir($ROOT_EXT);
    #&set_ext_dataroot();
}

# ----------------------------------
# set_ext_dataroot
#
# Set new EXT_DATAROOT. If necessary files exist, 
# then great. If not, then copy them from the ROOT_INSTALL.
# ----------------------------------
sub set_ext_dataroot {

    my $err_stat=0;
    if (!-d "$ROOT_EXT") {
       &info_dbox("Directory Error", "The EXT_DATAROOT entered does not exist.
It will be replaced with $ROOT_INSTALL/extdata.");
       $ROOT_EXT="$ROOT_INSTALL/extdata"; 
       $err_stat=1;
    }

    if ((-d "$ROOT_EXT/extprd") && (-d "$ROOT_EXT/work") &&
        (-d "$ROOT_EXT/static") && (-e "$ROOT_EXT/static/grib_prep.nl") ){
       # Necessary dirs and files exist.

    } else {

       # Necessary dirs and/or file(s) do not exist - copy them.
       my $rsp=&yesno_dbox("Missing Directory and Files", 
"The EXT_DATAROOT entered does not contain the expected 
subdirectories and/or file grib_prep.nl.\n
Would you like to have the necessary directories created and the 
necessary files coped to:
   $ROOT_EXT?\n
This involves making directories:
   $ROOT_EXT/log, 
   $ROOT_EXT/extprd,
   $ROOT_EXT/work
and, 
   cp -fR $ROOT_INSTALL/extdata/static to:
   $ROOT_EXT", "Yes"); 
       if ($rsp eq "Yes"){ 
          # Command will be executed.
          if (!-d "$ROOT_EXT/log")    { mkdir "$ROOT_EXT/log", 0777
            or warn("Could not mkdir $ROOT_EXT/log\n"), return; }
          if (!-d "$ROOT_EXT/extprd") { mkdir "$ROOT_EXT/extprd", 0777
            or warn("Could not mkdir $ROOT_EXT/extprd\n"), return; }
          if (!-d "$ROOT_EXT/work")   { mkdir "$ROOT_EXT/work", 0777
            or warn("Could not mkdir $ROOT_EXT/work\n"), return; }
          &watch_cursor(1);
          my $my_result=system ("cp -fR $ROOT_INSTALL/extdata/static $ROOT_EXT\n");
          &watch_cursor(0);
          $err_stat=1;

       } else {
          # Force a value for ROOT_EXT.
          $ROOT_EXT="$ROOT_INSTALL/extdata"; 
          $err_stat=1;
       }
    }
    
    # Re-assign filename value.
    $grib_prep_nl_NEW="$ROOT_EXT/static/grib_prep.nl";

    return($err_stat);
}

#__________________________________________________
#
# Read and write grib_namelist file.
#__________________________________________________


# --------------------------------------
# read_grib_namelist 
#
# Read grib_prep.nl namelist.
# --------------------------------------
sub read_grib_namelist {

  if ($model_name eq "LAPS") { return(1); }
  open(NL,"$grib_prep_nl") or 
        warn("Can't open grib_prep.nl: $grib_prep_nl.\n"), return(1);
  my @lines = <NL>; 
  close(NL);

  if ($model_name eq "WRF") {
     # Turn on store vars flags.
     $store_gnl_keys=1;
     $install_nl=1;

     &get_namelist_array(@lines);

#     # Because grib_prep.nl and wrfsi.nl both have identical sections
#     # (filetimespec) and variables (listed below), we store each 
#     # in separate array, ie. $nl_var, $nl_var2.
#     foreach $key (@sched_vars) {
#        $nl_var2{$key}[0] = $nl_var{$key}[0]; 
#        print "$key $nl_var{$key}[0]\n"; 
#     }
#

     # Turn off store vars flags.
     $store_gnl_keys=0;
     $install_nl=0;
  }

  # Create model grib background file entries.
  &add_mdl_source_entries(0,$nl_var_max{SRCNAME});

  # Clean up matrix.
  for ($ii=$matrix_rows; $ii >= $nl_var_max{SRCNAME}; $ii--) {

     # Erase remaining entries.
     foreach $key (@gsrc_vars) {
        $nl_var_orig{$key}[$ii]=$nl_var{$key}[$ii]=""; }

     # Delete blank widgets in columns 0-5, inclusively.
     for ($ll=0; $ll < 5; $ll++) {
        $xx=$g_entry[$ll][$ii];
        if ( Exists ($xx) ){ $xx->destroy; }
     }
  }

  # Fill grib sources (grib_src_mb).
  &update_grib_sources();

}

# --------------------------------------
# sort_grib_namelist
#
# Clean up grib namelist.
# --------------------------------------
sub sort_grib_namelist {

  # Sort all entries, eliminating blank entries.
  my $new_count=0;
  my $old_count=0;
  for ($old_count=0; $old_count <= $matrix_rows; $old_count++) {
     foreach $key (@gsrc_vars) {
        # Bubble sort.
        $nl_var{$key}[$new_count]     =$nl_var{$key}[$old_count];
        $nl_var_orig{$key}[$new_count]=$nl_var_orig{$key}[$old_count];
        $nl_var_max{$key}=$new_count;
     }
#    if ($nl_var{SRCNAME}[$old_count] eq "" ) {$new_count--;}
     if ($nl_var{SRCNAME}[$old_count] =~ m/^\s*$|\'\'/) {
     print "old: $nl_var{SRCNAME}[$old_count]\n";
     $new_count--;} 
     else {print "new: $nl_var{SRCNAME}[$old_count]\n";}
     $new_count++; 
  }

  # Clean up matrix.
  for ($ii=$matrix_rows; $ii >= $new_count; $ii--) {

     # Erase remaining entries.
     foreach $key (@gsrc_vars) {
        $nl_var_orig{$key}[$ii]=$nl_var{$key}[$ii]="";
     }

     # Delete blank widgets in columns 0-5, inclusively.
     for ($ll=0; $ll < 5; $ll++) {
        $xx=$g_entry[$ll][$ii];
        if ( Exists ($xx) ){ $xx->destroy; }
     }
  }
  $matrix_rows=$nl_var_max{SRCNAME}=$new_count;
print "\nLIMIT1: $new_count\n";

}

# --------------------------------------
# write_grib_namelist 
#
# Write grib_prep.nl namelist.
# --------------------------------------
sub write_grib_namelist {

  # Clean up grib namelist.
  &sort_grib_namelist();

  # Compare all variables to the original values.
  # If differences occur, then write all vars.
  # ---------------------------------
  my $count=0;
# my $jj=0;
# for ($ii=0; $ii < $nl_var_max{SRCNAME}; $ii++) {
#    $jj++;
#    foreach $key (@gsrc_vars) {

# print " no. sections: $num_gnl_sections\n";
  for $ii (1 .. $num_gnl_sections) {

#    print " no. entries: $num_gnl_entries[$ii]\n";
     for $jj (1 .. $num_gnl_entries[$ii]) {
        my $key = $gnl_var_array[$ii][$jj];

#       print " key: $key\t";
        my $limit=$nl_var_max{$key};
        if($nl_var_max_orig{$key} > $nl_var_max{$key}) {$limit=$nl_var_max_orig{$key};}

print "LIMIT2: $limit\n";

        for ($kk=0; $kk < $limit; $kk++) {

           if ($nl_var_orig{$key}[0] =~ m/^'/ ) {    # If orig entry has quote,
              $nl_var{$key}[$kk] =~ s/'//g;          # then strip off existing 
              $nl_var{$key}[$kk] ="\'$nl_var{$key}[$kk]\'"; # and add new quotes.
           }
   
           # Compare.
           if ($nl_var{$key}[$kk] ne $nl_var_orig{$key}[$kk]) {
              print "$kk: $nl_var{$key}[$kk] ne $nl_var_orig{$key}[$kk]\n";
              $nl_var_orig{$key}[$kk]=$nl_var{$key}[$kk];
              $count++;
              last; # Count the $key only once regardless of how
                    # many elements or its $nl_var_max{$key} change.
           } else {
#             print "$kk: $nl_var{$key}[$kk] eq $nl_var_orig{$key}[$kk]\n";

           }
        } #$kk
     }; #$jj
  }; #$ii


  # If changes is non-zero write all vars to grib prep namelist file (GNL).
  # ---------------------------------
# if ($count > 0) {

     if ($model_name eq "LAPS") { return(1); }
     open(GNL,">$grib_prep_nl") or 
         warn("Can't open grib_prep.nl: $grib_prep_nl.\n"), return;

     # Consider entry widgets for these entries.
     print GNL "\&filetimespec\n"; #section header
     print "\&filetimespec\n"; #section header
     foreach $key (@sched_vars) {
        # Deal with wrfsi.nl AND grib_prep.nl namelist files 
        # containing the @sched_vars; add indicator 
        # (eg. START_YEAR becomes XSTART_YEAR). See sub get_namelist_array().
        my $tmp_key="X$key";

        print GNL " $key = $nl_var{$tmp_key}[0]\n"; 
        print " $key = $nl_var{$tmp_key}[0]\n"; 
     }
     print GNL "/\n"; #section close

     #--
     print GNL "\&gpinput_defs\n"; #section header
     print "\&gpinput_defs\n"; #section header

     foreach $key (@gsrc_vars) {
        #print "key: $key num: $nl_var_max{$key}\n";
        for ($k=0; $k < $nl_var_max{$key}; $k++) {
           $val=$nl_var{$key}[$k]; 

           if ($val eq '' ) { last; }
           if ($k<1) { 
              print GNL " $key = $val"; 
              print " $key = $val"; 
           } elsif ($val =~ m/\// ) {  # If entry has a slash, ie. a directory
              print GNL ", \n\t\t$val"; 
              print ", \n\t\t$val"; 
           } else {
              print GNL ", $val"; 
              print ", $val"; 
           }
        }
        print GNL "\n"; 
        print "\n"; 
     }

     print GNL "/\n"; #section close
     close(GNL);


     # Fill grib sources (grib_src_mb).
     &update_grib_sources();

     $hint_msg="Successfully wrote '$grib_prep_nl'.\n And, updated Data Source pull-down menu on the 'Scripts' panel.";
# }

  #&set_button_state(0,$reload_grib_nl_but);
}

### Return 1 to the calling use statement ###
1;

