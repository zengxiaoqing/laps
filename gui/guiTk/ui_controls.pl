
# ---------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# ui_controls.pl 
# 	Collection of subroutines with error checking, 
#	dialog boxes, logic and navigation routines.
#
# Author: Paula McCaslin   30 Jan 2003  Original Version
# ---------------------------------------------------------------------------


#use warnings;
#use strict;
use strict 'subs';
use strict 'refs';
use vars qw($geog_path_button);
 
# -------------------------------------------
# increment_value
#
# Increment the value of the variable passed (operand) by amount
# e.g. print "$operand ($$operand) $amount\n";
#      increment stdlon (44.0) by -0.1.
# -------------------------------------------

sub increment_value {
    no strict 'refs';
    my ($operand,$amount,$win)=@ARG;

    # Bring focus to the entrybox associated with the up/down arrow.
    $win->focus;

    # Do the math.
    $$operand = $$operand + $amount;
}

# -------------------------------------------
# Message dialog boxes pop up to assist the
# user to confirm their decisions and/or 
# accuracy of entered information.
# 
# okcancel_dbox - warning to user
# info_dbox - information for user
# yesno_dbox - confirm user request 
# fail_dbox - request is not possible
# not_yet_imp - request has not yet been implemented
# -------------------------------------------

sub info_dbox {
    my ($mytitle,$mymsg)=@ARG;
    $mw->messageBox( -type => 'ok',
                     -image => $info_im,
                     -title => $mytitle,
                     -message => $mymsg,
                     -font=> $balloon_font,
                     -wraplength=>'7i',
                     -default => 'Ok');
}
sub okcancel_dbox {
    my ($mytitle,$mymsg)=@ARG;
    $mw->messageBox( -type => 'okcancel',
                     -image => $warning_im,
                     -title => $mytitle,
                     -message => $mymsg,
                     -font=> $balloon_font,
                     -wraplength=>'7i',
                     -default => 'Ok');
}
sub yesno_dbox {
    my ($mytitle,$mymsg,$default)=@ARG;
    $mw->messageBox( -type => 'yesno',
                     -image => $questhead_im,
                     -title => $mytitle,
                     -message => $mymsg,
                     -font=> $balloon_font,
                     -wraplength=>'7i',
                     -default => $default);
}
sub yesnocancel_dbox {
    my ($mytitle,$mymsg)=@ARG;
    $mw->messageBox( -type => 'yesnocancel',
                     -image => $questhead_im,
                     -title => $mytitle,
                     -message => $mymsg,
                     -font=> $balloon_font,
                     -wraplength=>'7i',
                     -default => 'Yes');
}
sub fail_dbox {
    my ($mytitle,$mymsg)=@ARG;
    $mw->messageBox( -type => 'ok',
                     -image => $warning_im,
                     -title => $mytitle,
                     -message => $mymsg,
                     -font=> $balloon_font,
                     -wraplength=>'7i',
                     -default => 'Ok');
}

#---------------------------------------------
# nyi 
# This feature has not yet been implemented.
#---------------------------------------------
sub nyi {
    $mw->messageBox( -type => 'ok',
                     -image => $warning_im,
                     -title => "Not Yet Implemented",
                     -message => "This feature has not yet been implemented",
                     -font=> $balloon_font,
                     -wraplength=>'7i',
                     -default => 'Ok');
}

# ----------------------------------------
# confirm_exit
#
# Confirm exit of application
# ----------------------------------------

sub confirm_exit {
    my $disp = &yesno_dbox("Confirm Exit", 
               "Are you sure you want to exit?", "Yes");
    if($disp eq "Yes"){ 
       $exit_now=1;

       if (!$globalMap && $domain_select ne "") {

          my $disp2 = &yesno_dbox("Write File", 
             "Do you want to save current namelist?", "No");
          # If the namelist exists, just update the
          # namelist else write_domain_files_and_image.
          if($disp2 eq "Yes"){
             if (&write_namelist()){ &write_domain_files_and_image();}
          }


       }
       exit;  
    } 
}

# -------------------------------------------
# about_gui
#
# Lets user know release version of gui.
# -------------------------------------------

sub about_gui {

    $about_gui=`cat "$GUI_TK/version"`; 
    &info_dbox("About this Tool", $about_gui);
}

# -----------------------------------
# show_balloon_tips
#
# Change cursor to a watch when processing
# data, return cursor when finished.
# -----------------------------------
sub show_balloon_tips {

      if ($cb_balloon) {
        $balloon->configure(-state => 'balloon');
        $domain_balloon->configure(-state => 'balloon');
      } else {
        $balloon->configure(-state => 'none');
        $domain_balloon->configure(-state => 'none');
      }

}


#----------------------------------
# launch_webpage
#
# Invoke a mozilla process to display on-line WRF SI help.
#----------------------------------

sub launch_webpage {
    my($my_url)=@_;

    # Set the cursor to look like a watch since the 
    # mozilla process may take awhile...
    &watch_cursor(1);
 
    my $result=system("mozilla -remote 'openURL($my_url,new-window)'");
    if($result != 0) {
      $result=system("mozilla -geometry =700x800+550+150 $my_url&");
      $hint_msg="Trying to start a Mozilla process";
    }
    if($result != 0) {
      &info_dbox("Problem with Mozilla.", "Cannot start a Mozilla process.");
      $hint_msg="Cannot start a Mozilla process";
    }

    &watch_cursor(0);
}


# -----------------------------------
# set_button_highlight
#
# Set da button state to 'active' or 'gray' 
# to indicate to the user an action needs to 
# be taken.
# -----------------------------------
sub set_button_highlight {
   my ($state,$button)=@ARG;

      if ($state) {
         $button->configure(-bg => $update_color);
      } else {
         $button->configure(-bg => $bg_color);
      }
}

# -----------------------------------
# set_button_state
#
# Set da button state to 'normal' or 'disabled' 
# to prevent the user from advancing until
# all criteria is met.
# -----------------------------------
sub set_button_state {
   my ($state,$button)=@ARG;

      if ($state) {
         $button->configure(-state => 'normal');
      } else {
         $button->configure(-state => 'disabled');
      }
}

# -----------------------------------
# set_button_list_state
#
# Set da button list state to 'normal' or 'disabled' 
# to prevent the user from advancing until
# all criteria is met.
# -----------------------------------
sub set_button_list_state {
   my ($state,@button_list)=@ARG;

   foreach my $button (@button_list) {
      if ($state) {
         $button->configure(-state => 'normal', -fg => $normal_color);
      } else {
         $button->configure(-state => 'disabled', -fg => $disabled_color);
         &switch_focus_main();
      }
   }

}

# -----------------------------------
# watch_cursor
#
# Change cursor to a watch when processing
# data, return cursor when finished.
# -----------------------------------
sub watch_cursor {
    my ($delay)=@ARG;

    if($delay == 1) {
      # Change cursor icon to a wristwatch.
      $mw_cursor=$mw->cget(-cursor);
      $can_cursor=$can->cget(-cursor);

      $mw->configure(-cursor => 'watch');
      $can->configure(-cursor => 'watch');

    }elsif($delay == 0) {

      # Return cursor icon to an arrow.
      if($mw_cursor eq 'watch') {
         # Don't reset active cursor to a watch this could happen if
         # watch_cursor is called from several subroutines at once.
         $mw->configure(-cursor => "left_ptr");
         $can->configure(-cursor => "crosshair");

      } else {

         $mw->configure(-cursor => $mw_cursor);
         $can->configure(-cursor => $can_cursor);
      } 

    }elsif($delay == 2) {

      #$mw->configure(-cursor => $cursor_Diag2);
      #$can->configure(-cursor => $cursor_Diag2);
      $mw->configure(-cursor => 'watch');
      $can->configure(-cursor => 'watch');

    }elsif($delay == 3) {

      $mw->configure(-cursor => 'watch');
      $can->configure(-cursor => 'watch');

    }

    $mw->update;
    $mw->idletasks;
     
}

# ----------------------------------------
# switch_focus_main
#
# Sometimes the focus needs to be taken away 
# from a current widget. This does the trick.
# ----------------------------------------
sub switch_focus_main { $mw->focus; }

# ----------------------------------------
# delete_hints
#
# Clear the contents of the Hints and Info area.
# ----------------------------------------
sub delete_hints { $hint_msg=""; }

# ----------------------------------------
# splash_val
#
# Sometimes processes take much time.
# This does the trick to inform user.
# ----------------------------------------
sub splash_val { 
   my ($val)=@_;

   if ($val == 1) {
      $pctComplete=0; 
      $splash->deiconify;
      $splash->raise;
   } elsif ($val == 0) {
      $splash->withdraw;
   } else {
      $pctComplete+=$val; 
   }
   $splash->update();

}

# ----------------------------------
# browse4dir
#
# Browse button was pressed to browse directories.
# ----------------------------------

sub browse4dir {
    no strict 'refs';
    my ($my_path)=@_;
    print "browse4dir: $my_path \n";

    # Check for a good directory path.
    if (!-d $my_path) {
       if ($my_path eq "" || $my_path=~ m/ /) {
          $my_path='/'; 
       } else {
          my $ans = &info_dbox("Directory Problem", 
                    "Directory '$my_path' does not exist!"); 
          if($ans eq "Ok"){ return;}
       }
    }
    
    my $my_dir=$my_path;
    # Browse for directory to set a new path.
    $browseDialog2=$mw->FileDialog(-Title =>"Select data output directory",
 		                  -Create => 0,
                                  -PathEntryLabel => "Path:",
                                  -Path => $my_dir,
 	 	                  -ShowAll => 'NO');
    $browseDialog2->Show(-SelDir => 1, -Horiz => 1);

    # Check again for a good directory path.
    $my_dir = $browseDialog2->{Configure}{-Path};
    if (!-d $my_dir) {
       my $ans = &yesno_dbox("Directory Problem", 
                 "Directory '$my_dir' does not exist!\nAccept anyway?","No"); 
       if($ans eq "Yes"){ return; } 
    } 

    #Strip off backslash at end of var
    $my_dir =~ s/\/$//  if (substr($my_dir, -1, 1) eq '/');

    # Destroy when finished.
    $browseDialog2->destroy if Tk::Exists($browseDialog2);

    return($my_dir);
}


# ----------------------------------------
# selectColor
#
# ----------------------------------------
sub selectColor {
    if ($globalMap == 1) {
       &info_dbox("No Lines to Color",
         "Select new line colors after pressing 'Update Map'.");
       return;
    }
    my $tag = '$_[0]';
    #my $tag = $_[0];
    #my $tag2="$tag";
    my $cColor = $can->chooseColor(
        #-initialcolor => $can->itemcget($tag, -fill),
        -popover => $can, 
        -popanchor => 'w',
        -overanchor => 'e', 
        -title => "Select Color");
    $can->itemconfigure($tag, -fill => $cColor) if defined $cColor;
}

# ----------------------------------------
# how_long
# ----------------------------------------
sub how_long {
    my ($string)=@_;
    my ($user, $system, $cuser, $csystem) = times();
    $total_time=$user+ $system+ $cuser+ $csystem;
    $current_time=$total_time-$prev_time;
    $prev_time=$total_time;
    printf "$string has consumed %.3f seconds\n", $current_time;
}
 
# ----------------------------------
# create_tsplash
# 
# Display the message window with the message 
# that this is a Time Consuming Process.
# ----------------------------------
sub create_tsplash {

  ## Create time splash window.
  $tsplash = $mw->Toplevel(-title => "Time Consuming Process",
                              -bg => '#d9d9d9',
                              -visual => 'truecolor');
  $tsplash->geometry("+348+100");
  $tsplash->focus;
  $tsplash->protocol('WM_DELETE_WINDOW', [\&confirm_exit]);


  $tframe=$tsplash->Frame(-relief => 'flat', -bd => 2)
              ->pack(-anchor => 'c');

#  $tframe->Label(
#                  # Image route
#                  -image => $info_im,)
#                  ->pack(-side => 'left', -padx => 25, -pady => 25);

  $tframe->Label(
                  # Ascii route
                  -font => $bold_font, 
                  -text => "This could take several minutes...\nas your domain is $nx_dim by $ny_dim.")
                  ->pack(-padx => 15, -pady => 25);

  $tsplash->bind ('<Visibility>',  sub {$tsplash->raise(); $mw->idletasks} );
}

# ----------------------------------
# forget_tsplash 
#
# Destroy the tsplash window
# ----------------------------------
sub forget_tsplash { $tsplash->destroy if Tk::Exists($tsplash); }

### Return 1 to the calling use statement ###
1;
