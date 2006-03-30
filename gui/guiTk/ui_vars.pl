# ---------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# srt_vars.pl 
# 	Initializes the global variables such as min, max, 
#	file names. Contains file read and write subroutines.
#
# Author: Paula McCaslin   30 Jan 2003  Original Version
# ---------------------------------------------------------------------------


#use warnings;
#use strict;
use strict 'subs';
use strict 'refs';

use vars qw($geog_path_orig $geog_path); 

#__________________________________________________
#
# Define all variables.
#__________________________________________________

# --------------------------------------
# define_env_vars
#
# Application variables.
# --------------------------------------
sub define_mdl_vars {

   if ($model_name eq "WRF") {
     
     $window_domain_arg ="wrfsi.d01";
     $window_domain_arg2="wrfsi";
     if ($nmm) {$window_domain_arg ="wrfsi.rotlat";}
     if ($nmm) {$window_domain_arg2="wrfsi.rotlat";}
     $namelist_arg="wrfsi.nl";
     @root_env_label = qw(SOURCE_ROOT INSTALLROOT MOAD_DATAROOT);

     # Set grib_prep_nl.
     $grib_prep_nl="$ROOT_EXT/static/grib_prep.nl";
     #$wrf_nl="$SOURCE_ROOT/data/static/wrf.nl";
     $wrf_nl="$ROOT_TEMPLATES/default/wrf.nl"; 
     $readmeURL="http://wrfsi.noaa.gov/wrfsiREADME.html";


   } elsif ($model_name eq "LAPS") {

     $window_domain_arg ="nest7grid";
     $window_domain_arg2="laps";
     $namelist_arg="nest7grid.parms"; 
     @root_env_label = qw(LAPSSRCROOT INSTALLROOT LAPSDATAROOT);
     $readmeURL="http://laps.noaa.gov/software/README_0-27-19.html";
   }

   $MDL_NAMELIST="$ROOT_TEMPLATES/default/$namelist_arg"; 
   $dataroot_select=$ROOT_DATA;
   
   $helpURL="http://wrfsi.noaa.gov/gui/users_guide/V3/users_guide.html";
}

sub set_app_variables {

  $Debug=0;
  $d_num="d01";
  $store_nl_keys=1;
  $firstTime=1; 
  $startover_vars_init=0;
  $previous_args="";
  $NXmax=101;
  $perform_calc=1;
  $calc_centerLatLon=1;
  $deg2rad = 0.017453293;
  $polar_cap=0;

  # Nesting vars
  $grid_buffer=4;
  $max_nests=9;
  $nestLab[0]="d01";

  # Check for executables.
  $gridgen_mdl="$ROOT_INSTALL/bin/gridgen_model.exe";
  #$gridgen_mdl="$ROOT_INSTALL/bin/gridgen_nmm.exe" if $projection_type eq 'RL';
  $gen_map_exe="$GUI_EXE/gen_map_bkgnd.exe";
  $pwrap_ll_xy_convert_exe="$GUI_EXE/pwrap_ll_xy_convert.exe"; 
  if (!-e $gridgen_mdl) { 
      die "Executable $gridgen_mdl doesn't exist and is mandatory.\n"; }
  if (!-e $gen_map_exe) { 
      die "Executable $gen_map_exe doesn't exist.\n"; }
  if (!-e $pwrap_ll_xy_convert_exe) { 
      die "Executable $pwrap_ll_xy_convert_exe doesn't exist.\n"; }

  # Temporary text output files.
  $mapVectorFile="$tempdir/vector_instructions.tk";
  $setsupData="$tempdir/setsup.dat";
  $logFile="$tempdir/wrf_tools.log";
  if ($model_name eq "LAPS") { $logFile="$tempdir/ui_tools.log"; }

  # Find for map background files (bcd).
  find_bcdFiles();

  #GUI label, LAPS variable, WRFSI var, NCARGraphics supmap num, map_setup var
  @proj_choices = (
            ["Polar stereographic", "plrstr", "polar",    "1", "PS"],
            ["Lambert Conformal",   "lambrt", "lambert",  "3", "LC"],
            ["Mercator",            "merctr", "mercator", "9", "ME"],
            ["Rot Lat-Lon (NMM)",   "lambrt", "rotlat",   "3", "RL"],
  );
  # Counting from zero, three elements are =2.
  $max_proj_choices=2;

  if ($nmm) {
  #GUI label, LAPS variable, WRFSI var, NCARGraphics supmap num, map_setup var
     @proj_choices = (
            ["Rot Lat-Lon (NMM)",   "lambrt", "rotlat",   "3", "RL"],
     );
     $max_proj_choices=0;
  }

  # Two independant namelist files contain the following variables,
  # wrfsi.nl and grib_prep.nl.
  @sched_vars=('START_YEAR', 'START_MONTH', 'START_DAY', 
               'START_HOUR', 'START_MINUTE', 'START_SECOND', 
               'END_YEAR', 'END_MONTH', 'END_DAY', 
               'END_HOUR', 'END_MINUTE', 'END_SECOND', 
               'INTERVAL');
 
  @gsrc_vars= ('SRCNAME', 'SRCVTAB', 'SRCPATH', 'SRCCYCLE', 'SRCDELAY');

  @nesting_vars= qw(PARENT_ID RATIO_TO_PARENT
		DOMAIN_ORIGIN_LLI DOMAIN_ORIGIN_LLJ 
		DOMAIN_ORIGIN_URI DOMAIN_ORIGIN_URJ
		);
  
}

# --------------------------------------
# gui_vars
#
# Initialize gui specific application variables
# such as image file names, colors, fonts, etc
# --------------------------------------
sub gui_vars {

# --- set color & misc options ---

  $colorN="navy";
  $colorW="white";
  $colorG="gray50";
  $colorR="red";
  $colorY="#ffffbb";
  $colorY="#ffffaa"; # little brighter yellow
  $colorY2="#ffff33";
  $colorL="#d9d9ef";
  $colorL2="#d9d9ff";
  $colorg="#c0d0d0";
  $colorO="orangered";
  $normal_color="black";
  $disabled_color="gray50";
  $update_color="cornsilk";
  $update_color=$colorY;
  $bg_color="#d9d9d9";        #$bg_color=$mw->cget(-bg) "gray65";
  $bg_canvas="#27ef645aaa7e"; #$bg_canvas="darkslateblue";
  $mw_cursor=$mw->cget(-cursor);
# ---

  #$normal_font = "Helvetica -12"; 
  #$thin_font = 'Helvetica 11';
  $thin_font = "Helvetica -12"; 
  $bold_font = "Helvetica -12 bold"; 
  $balloon_font = "-*-helvetica-bold-n-*-*-12-*-*-*-*-*-*-*";
  $error_font   = "-*-helvetica-bold-*-*-*-14-*-*-*-*-*-*-*";
  $small_font = "Helvetica 9 bold";
  $italics_font = "Helvetica 9 italic bold";
  $preview_font = "helvetica 16 italic bold";
  $large_font = "Times 24";
  $legend_font  = 'helvetica 11 bold';
  $grayfg=$mw->ItemStyle('text', -fg => 'gray50', 
                                 -bg => 'gray85',
                                 -font => $bold_font,
                                 -selectforeground => 'gray40'
                        );
  if($GUI_TK) {
    $small_up   = "\@$GUI_TK/icons/small-up.xbm";
    $small_down = "\@$GUI_TK/icons/small-down.xbm";
    $my_watch   = "\@$GUI_TK/icons/watch.xbm";
    $my_watch_msk= "$GUI_TK/icons/watch.mask.xbm";

    $info_im=$mw->Photo(-file => "$GUI_TK/icons/info.ppm");
    $warning_im=$mw->Photo(-file => "$GUI_TK/icons/warning.ppm");
    $questhead_im=$mw->Photo(-file => "$GUI_TK/icons/questhead.ppm");

    my $CenterPoint="$GUI_TK/icons/CenterPoint.xpm";
    my $CenterPoint_dot="$GUI_TK/icons/CenterPoint_dot.xpm";
    my $CenterPoint_dotY="$GUI_TK/icons/CenterPoint_dotY.xpm";
    my $ModelLabel="$GUI_TK/icons/wrf_label.gif";
    if ($model_name eq "LAPS") { $ModelLabel="$GUI_TK/icons/laps_banner.gif"; }
    my $DefaultMap="$GUI_TK/icons/cylindrical_map.gif";
    my $null_domain="$GUI_TK/icons/Null_Domain.gif";
    my $ZoomOut="$GUI_TK/icons/ZoomOut.gif";
    my $ZoomIn="$GUI_TK/icons/ZoomIn.gif";
    my $NewNest="$GUI_TK/icons/NewNest.gif";

    $centerPoint=$mw->Photo(-file => $CenterPoint);
    $centerPoint_dot=$mw->Photo(-file => $CenterPoint_dot);
    $centerPoint_dotY=$mw->Photo(-file => $CenterPoint_dotY);
    $banner_image=$mw->Photo(-file => $ModelLabel);
    $dmap_image=$mw->Photo(-file => $DefaultMap);
    $zoom_out=$mw->Photo(-file => $ZoomOut);
    $zoom_in=$mw->Photo(-file => $ZoomIn);
    $newNest=$mw->Photo(-file => $NewNest);

    $img_xmax=$dmap_image->width();
    $img_ymax=$dmap_image->height();
  } 

  set_resource_options();

}

# --- set resource options ---
sub set_resource_options {

  $mw->optionAdd(qw/Idt*canvas.width:   400/);
  $mw->optionAdd(qw/Idt*canvas.height:  200/);
  $mw->optionAdd(qw/*foreground                        black/);
  $mw->optionAdd(qw/*selectForeground                  black/);
  $mw->optionAdd(qw/*selectColor                       yellow/);
  $mw->optionAdd(qw/*background                        gray85/);
  $mw->optionAdd(qw/*activeBackground                  gray83/);
  $mw->optionAdd(qw/*blinkingHighlightColor            CornSilk/);
  $mw->optionAdd(qw/*disabledForeground                gray48/);
  $mw->optionAdd(qw/*displacedTimeBackground           black/);
  $mw->optionAdd(qw/*displacedTimeForeground           yellow/);
  $mw->optionAdd(qw/*highlightBackground               gray85/);
  $mw->optionAdd(qw/*lightBackground                   white/);
  $mw->optionAdd(qw/*productAvailableForeground        springgreen1/);
 #$mw->optionAdd(qw/*shadowBackground                  \#818181/);
  $mw->optionAdd(qw/*titleColor                        blue/);
  $mw->optionAdd(qw/*troughColor                       gray75/);
  $mw->optionAdd(qw/*urgentActiveBackground            red/);
  $mw->optionAdd(qw/*urgentBackground                  red2/);
  $mw->optionAdd(qw/*Checkbutton.selectColor           yellow/);
  $mw->optionAdd(qw/*Entry.background                  WhiteSmoke/);
  $mw->optionAdd(qw/*Entry.background                  gray85/);
  $mw->optionAdd(qw/*Entry.foreground                  Black/);
  $mw->optionAdd(qw/*Entry.selectBackground            CornSilk/);
  $mw->optionAdd(qw/*HList.background                gray85/);
  $mw->optionAdd(qw/*HList.foreground                Black/);
  $mw->optionAdd(qw/*HList.selectBackground          CornSilk/);
  $mw->optionAdd(qw/*HList.selectForeground          Black/);
  $mw->optionAdd(qw/*HList.exportSelection           no/);
  $mw->optionAdd(qw/*Listbox.background                gray85/);
  $mw->optionAdd(qw/*Listbox.foreground                Black/);
  $mw->optionAdd(qw/*Listbox.selectBackground          CornSilk/);
  $mw->optionAdd(qw/*Listbox.selectForeground          Black/);
  $mw->optionAdd(qw/*Listbox.exportSelection           no/);
  $mw->optionAdd(qw/*Menu.selectColor                  yellow/);
  $mw->optionAdd(qw/*Radiobutton.selectColor           yellow/);
  $mw->optionAdd(qw/*Text.background                   gray85/);
  $mw->optionAdd(qw/*Text.foreground                   Black/);
  $mw->optionAdd(qw/*Text.selectBackground             CornSilk/);
  $mw->optionAdd(qw/*Text.font  -*-helvetica-medium-r-*-*-*-120-*-*-*-*-*-*/);
  $mw->optionAdd(qw/*HList.font  -*-helvetica-medium-r-*-*-*-120-*-*-*-*-*-*/);
  $mw->optionAdd(qw/*Listbox.font  -*-helvetica-medium-r-*-*-*-120-*-*-*-*-*-*/);
  $mw->optionAdd(qw/*fixedFont  -*-lucidatypewriter-medium-*-*-*-*-120-*-*-*-*-*-*/);
  $mw->optionAdd(qw/*font       -*-helvetica-*-r-*-*-12-120-*-*-*-*-*-*/);
  $mw->optionAdd(qw/*italicFont -*-helvetica-*-o-*-*-12-*-*-*-*-*-*-*/);
  $mw->optionAdd(qw/*statusFont -*-helvetica-medium-r-*-*-12-120-*-*-*-*-*-*/);
  $mw->optionAdd(qw/*urgentFont -*-helvetica-bold-r-*-*-12-140-*-*-*-*-*-*/);
  $mw->optionAdd(qw/*TkFDialog.icons.canvas.background WhiteSmoke/);
  $mw->optionAdd(qw/*menubarPadY                       2/);
  $mw->optionAdd(qw/*Button.borderWidth                2/);
  $mw->optionAdd(qw/*Button.highlightThickness         1/);
  $mw->optionAdd(qw/*Button.padX                       7/);
  $mw->optionAdd(qw/*Button.padY                       2/);
  $mw->optionAdd(qw/*Canvas.highlightThickness         1/);
  $mw->optionAdd(qw/*Checkbutton.borderWidth           2/);
  $mw->optionAdd(qw/*Checkbutton.highlightThickness    1/);
  $mw->optionAdd(qw/*Entry.borderWidth                 2/);
  $mw->optionAdd(qw/*Entry.highlightThickness          1/);
  $mw->optionAdd(qw/*HList.borderWidth               2/);
  $mw->optionAdd(qw/*HList.highlightThickness        1/);
  $mw->optionAdd(qw/*HList.selectBorderWidth         1/);
  $mw->optionAdd(qw/*Listbox.borderWidth               2/);
  $mw->optionAdd(qw/*Listbox.highlightThickness        1/);
  $mw->optionAdd(qw/*Listbox.selectBorderWidth         1/);
  $mw->optionAdd(qw/*NoteBook.borderWidth              2/);
  $mw->optionAdd(qw/*NoteBook.tabPadX                  10/);
  $mw->optionAdd(qw/*NoteBook.backPageColor            gray85/);
  $mw->optionAdd(qw/*Menu.borderWidth                  1/);
  #$mw->optionAdd(qw/*Menubutton.activeBackground       cornsilk/);
  $mw->optionAdd(qw/*Menubutton.highlightThickness     1/);
  $mw->optionAdd(qw/*Radiobutton.borderWidth           2/);
  $mw->optionAdd(qw/*Radiobutton.highlightThickness    1/);
  $mw->optionAdd(qw/*Scale.borderWidth                 2/);
  $mw->optionAdd(qw/*Scale.highlightThickness          1/);
  $mw->optionAdd(qw/*Scrollbar.borderWidth             2/);
  $mw->optionAdd(qw/*Scrollbar.highlightThickness      1/);
  $mw->optionAdd(qw/*Text.borderWidth                  2/);
  $mw->optionAdd(qw/*Text.highlightThickness           1/);
  $mw->optionAdd(qw/*TkFDialog.borderWidth             10/);
  $mw->optionAdd(qw/*TkFDialog.relief                  flat/);

  if (defined $ROOT_NCARG) {
    $mw->optionAdd(qw/idt*geometry:               300x300+0+0/);
    $mw->optionAdd(qw/idt*TopLevelShell*geometry: +0+536/);
  }

}    

#__________________________________________________
#
# I/O.
#__________________________________________________

# --------------------------------------
# load_namelist
#
# Select a namelist to be read. In the case of a
# domain template start again with the Install Root
# namelist, because a domain template only show the
# variables that have been modified and is not a 
# complete namelist.
# --------------------------------------
sub load_namelist {

  # Read Install Root namelist  
  $install_nl=1;
  $mdl_namelist=$MDL_NAMELIST;
  read_namelist();

  if ($domain_mode != 1) {

      # Additionally, read Dataroot namelist -- merge these lists.
      $install_nl=0;
      $mdl_namelist="$ROOT_TEMPLATES/$domain_select/$namelist_arg"; 
      read_namelist(); 
  }

  if($model_name eq "WRF"){ 
     $num_nests=$nl_var{NUM_DOMAINS}[0];
     if ($num_nests > 1) { 
        fill_nl_disp_ary(); 
        reduce_nl_maxs_nest(); 
     }
  }

  # Load pressure levels for LAPS from pressures.nl.
  if ($model_name eq "LAPS") {
       if ($domain_mode == 1) { 
          $install_nl=1;
          $mdl_namelist="$ROOT_TEMPLATES/default/pressures.nl"; 
          read_namelist(); 
          $install_nl=0;
       } else {
          $levels_index=0;
          $install_nl=1;
          $mdl_namelist="$ROOT_TEMPLATES/$domain_select/pressures.nl"; 
          read_namelist(); 
          $install_nl=0;
       }
     $num_nests=1;
  }

  # Fill application variables with model namelist variables
  assign_nl_vars(0);
  set_nx(); 
  set_dx(); 

  # For new domain always reset nx and projection label.
  if ($domain_mode == 1) {
     $nl_var{XDIM}[0]=$nl_var{NX_L_CMN}[0]=$nx_ddim=$NXmax;
     $proj_label="None";
  }

  # Load and display vertical cooridinate levels.
  if (Exists($ve_mb)) { set_vert_scheme_current(); }

}

# ---------------------------------------
# read_namelist 
#
# Read namelist, then parse contents.
# ---------------------------------------
sub read_namelist {

  if (-z "$mdl_namelist") {
     fail_dbox("Namelist Size Error", "The file $mdl_namelist\nhas zero size.");
     die "The file $mdl_namelist\nhas zero size.";
  }

  open(NL,"$mdl_namelist") or warn "Can't open mdl_namelist: $mdl_namelist\n", 
  $hint_msg="Can't open mdl_namelist: $mdl_namelist", return(1);
  my @lines = <NL>; 
  close(NL);

  get_namelist_array(@lines);
}

#__________________________________________________
#
# Parse namelist.
#__________________________________________________


# --------------------------------------
# get_namelist_array 
#
# Parse info from WRF (wrfsi.nl) 
# lines will contain one of:
#       1) SECTION NAME
#       2) SOLITARY BACKSLASH
#       3) KEY/VALUE PAIR
#       4) MULTI-LINE ENTRY
#
# Store the namelist entries into:
#       1) store key/value pairs in a hash (%nl_var) 
#       2) store keys into an array (@nl_var) to organize printing,
#          as hashes loose order in which entries are stored.
# --------------------------------------

#use vars qw(%nl_var_orig %nl_var_max
#            @nl_var @nl_var_array $nl_var_max 
#            @num_nl_entries $num_nl_sections @nl_section
#            $nl_sig_levels $geog_path $geog_path_orig);
#
sub get_namelist_array {
  no strict 'refs';
  #use strict;

  my (@orig_lines) = @ARG;
  my $levels_flag=0; 
  my @values;
  my $valstring;
  #my $key;

  my $i=0;
  my $j=0;
  my $k=0;

  foreach (@orig_lines) {

    # ---- parse and skip solitary backslash at end of section ----

    s/^\s//g;                       # clear out beginning whitespace
    if (/^\/$/) { next; }           # skip backslash at end of last value
    if (/^#/)   { next; }           # ditch comments
    if (/^c\s/) { next; }           # ditch comments
    if (/^C\s/) { next; }           # ditch comments

    # ---- clean up line -----------------------
                                    # s/\s+//g;  #clear out all whitespace
    s/^\s//g;                       # clear out beginning whitespace
    s/\!.*//;                       # strip off comments
    s/\s*$//g;                      # clear out ending whitespace
    s/\,*$//g;                      # strip off ending commas

    # ---- parse and store section name -----------------------


    if (/^&/) {                     # line begins w/ '&' 
      s/&//g;

      $i++;
      $j=0;
      if ($store_nl_keys) {
          # because this is the full namelist!
          $nl_section[$i] = $ARG;

      } elsif ($store_gnl_keys) {
          # because this is the full grib prep namelist!
          $gnl_section[$i] = $ARG;
      }

      next; 
    }              

    # ---- parse and store key/value pairs ----------------------

    if (/=/) {    

      s/\/$//;                         # strip backslash at end of last value
     #s/\'//g;                         #strip quotes
                                       ##s/\'//g;   #strip quotes

      $j++;
      ($key,$valstring) = split /\s*=\s*/, $ARG; #equal sign with any whitespace
      $valstring =~ s/^\s//g;          # clear out beginning whitespace
      $key = uc $key;                  # upper case

      if($key eq "SIMULATION_NAME" || $key eq "USER_DESC" || 
         $key eq "C80_DESCRIPTION"){
         @values = $valstring;
      } else {
         $valstring =~ s/\s+//g;       # clear out all whitespace
         @values = split /,/, $valstring;
      }
      #print "\t$j: $key = $valstring\n";


      # Process vars for wrfsi.nl/nest7grid.parms AND grib_prep.nl differently.
      # ----
      if ($store_nl_keys) {
         # Full namelist, so $key is new (eg. XDIM).
         # ----
         $nl_var_array[$i][$j] = $key;       

      } elsif ($store_gnl_keys) {
         # Make filespec vars unique (eg. START_YEAR becomes XSTART_YEAR).
         # ----
         if ($key =~ m/^START|^END|^INTERVAL/) {$key="X$key";}

         # Grib Prep namelist, so $key is new (eg. SRCNAME).
         # ----
         $gnl_var_array[$i][$j] = $key;       
         $nl_var_max_orig{$key} = scalar(@values); 
      }
 
      # Store the variables.
      # ----
      $nl_var_max{$key} = scalar(@values); 
      for ($k=0; $k < $nl_var_max{$key}; $k++) {
        $nl_var{$key}[$k] = $values[$k]; 
        $nl_var_prev{$key}[$k] = $values[$k]; 
        if ($install_nl) { $nl_var_orig{$key}[$k] = $values[$k]; };
      }

      # For var ACTIVE_SUBNESTS, treat element=0 as a list.
      # ----
      if ($key eq 'ACTIVE_SUBNESTS') { 
	$nl_var{$key}[0]=$valstring; 
	$nl_var_max{$key}=1;
      }

    } else {     

    # ---- parse multi-line key entries, i.e. LEVELS & SRCPATH -------------
      #print "\t$j: $key = $valstring\n";
      if ($key =~ m/LEVELS|PRESSURES/) { 
      # LEVELS variable is multi-line in wrfsi.nl, clean up variable.

         if (!$levels_flag) {
           # When code gets here, this manipulates the line read previously.
           $levels_flag=1; 

           @nl_sig_levels=();
           for ($k=0; $k < (scalar(@values)); $k++) {
              $nl_var{$key}[$k] = $values[$k]; 
              if ($install_nl) { $nl_var_orig{$key}[$k] = $values[$k]; };
              push(@nl_sig_levels, $values[$k]);
              $levels_index=$k+1;
           }

         }

         # Manipulate the current line.
         s/\s+//g;                  #clear out all whitespace
         @values = split /,/, $ARG;

         $kk=$k=0;
         for ($kk=$levels_index; $kk < ($levels_index + scalar(@values)); $kk++) {
            $nl_var{$key}[$kk] = $values[$k]; 
            if ($install_nl) { $nl_var_orig{$key}[$kk] = $values[$k]; };
            #print  "$levels_index $values[$k]\n";
            push(@nl_sig_levels, $values[$k]);
            $k++;
         }
         $nl_var_max{$key}=$levels_index=$kk;
         $nl_var_max_orig{$key}=$kk;

      }

      if ($key =~ m/SRCPATH/) { 
      # SRCPATH variable is multi-line in grib_prep.nl, clean up variable.
         
         s/\s+//g;                  #clear out all whitespace
         @values = split /,/, $ARG;
         $nl_var{$key}[$k] = $values[0]; 
	 if ($install_nl) { $nl_var_orig{$key}[$k] = $values[0]; };
         $k++;

         $nl_var_max{$key}=$k;
         #print "ARG $k is $ARG\n"; # ($values[0])\n";
      }


    } #end key/value pair

    if ($store_nl_keys)       { $num_nl_entries[$i] = $j; 
    } elsif ($store_gnl_keys) { $num_gnl_entries[$i] = $j; 
}

  } #end foreach

  if ($store_nl_keys)       { 
     # Model namelist.
     $num_nl_sections=$i;

  } elsif ($store_gnl_keys) { 
     # Grib_prep namelist.
     $num_gnl_sections=$i; 

  }
 

  #----

  # For geog and si paths, only. Create a common PATH variable 
  # by stripping off its suffix leaving only its prefix.
  if (1) {
     $ARG=$nl_var{TOPO_30S}[0];
     if ($model_name eq "LAPS") {$ARG=$nl_var{PATH_TO_TOPT30S}[0];}
     s/\/[^\/]*$//;                # Strip off directory path.
     s/'//g;                       # Strip off quote(s).
     $geog_path=$ARG;
  #--
     $ARG=$nl_var{LBCPATH}[0];
     s/\/[^\/]*$//;                # Strip off directory path.
     s/'//g;                       # Strip off quote(s).
     $si_path=$ARG;
     if ($si_path ne $ROOT_EXT) {
        $si_msg="Replace the string: 
$si_path
with: 
$ROOT_EXT?
Use MouseButton1 to make this change,
use MouseButton2 to undo this change.";
     } else {
        $si_msg="Nothing to update.";
     }
  }

  # Once initial namelist variables and hash keys are stored 
  # never change them, store_nl_keys=0.
  $store_nl_keys=0;

  # WRF SI versions 1.3.2 (and earlier) had multiple dims
  # current versions necessitate only one.
  if($model_name eq "WRF"){
     $nl_var_max{XDIM}[0]=$nl_var_max{YDIM}[0]=1;
  }

}

#__________________________________________________
#
# Fill display variables with model variables and visa versa.
#__________________________________________________

# ---------------------------------------
# assign_nl_vars 
#
# Set app variables from model namelist variables
# ---------------------------------------
sub assign_nl_vars {
  my ($return_flag)=@_;
  no strict 'refs';

  if($model_name eq "WRF"){
     
     # Set all the app variables from model namelist.
     $num_nests=$nl_var{NUM_DOMAINS}[0];

     $grid_cen_lat_cmn=$nl_var{MOAD_KNOWN_LAT}[0];
     $grid_cen_lon_cmn=$nl_var{MOAD_KNOWN_LON}[0];
     $stdlon=$nl_var{MOAD_STAND_LONS}[0];
     $truelat1=$nl_var{MOAD_STAND_LATS}[0];
     $truelat2=$nl_var{MOAD_STAND_LATS}[1];

     $grid_spacing_m_cmn=$nl_var{MOAD_DELTA_X}[0];
     $dx_orig=$grid_spacing_dkm=($grid_spacing_m_cmn/1000.);
     $nx_orig=$nx_ddim=$nl_var{XDIM}[0];
     $ny_orig=$ny_dim=$nl_var{YDIM}[0];
     $ptop_mb=$nl_var{PTOP_PA}[0]/100;

     # Don't set projection for existing domains
     # e.g. those called from 'sub reinstate_gen_map_vars'.
     if ($return_flag) {return};

     # Set proj related variables.
     set_proj_widget($nl_var{MAP_PROJ_NAME}[0]);

  } elsif($model_name eq "LAPS"){

     # Set all the app variables from model namelist.
     $grid_cen_lat_cmn=$nl_var{GRID_CEN_LAT_CMN}[0];
     $grid_cen_lon_cmn=$nl_var{GRID_CEN_LON_CMN}[0];
     $stdlon=$nl_var{STANDARD_LONGITUDE}[0];
     $truelat1=$nl_var{STANDARD_LATITUDE}[0];
     $truelat2=$nl_var{STANDARD_LATITUDE2}[0];

     $grid_spacing_m_cmn=$nl_var{GRID_SPACING_M_CMN}[0];
     #$grid_spacing_km=$grid_spacing_dkm=($grid_spacing_m_cmn/1000.);
     $dx_orig=$grid_spacing_dkm=($grid_spacing_m_cmn/1000.);
     $nx_orig=$nx_ddim=$nl_var{NX_L_CMN}[0];
     $ny_orig=$ny_dim=$nl_var{NY_L_CMN}[0];
     $pr_levels=$nl_var{NK_LAPS}[0]-1;
     $pr_levels=$nl_var_max{PRESSURES}-1;

     $v::pbot_mb=$nl_var{PRESSURES}[0];
     $ptop_mb=$nl_var{PRESSURES}[$pr_levels];

     # Don't set projection for existing domains
     # e.g. those called from 'sub reinstate_gen_map_vars'.
     if ($return_flag) {return};

     # Set proj related variables.
     set_proj_widget($nl_var{C6_MAPROJ}[0]);
  }

}

# --------------------------------------
# set_proj_widget
#
# Set the projection widget's proj_label to the 
# correct value when a namelist is read in and
# when user chooses new projection type.
# --------------------------------------
sub set_proj_widget {
    my ($pro_arg) = @ARG;


    $pro_arg =~ s/\'//g;      #strip quotes to test projection var

    my $ix=0;
    until ($ix>$max_proj_choices) {
      if ($pro_arg eq $proj_choices[$ix][1] ||    # Test for LAPS
          $pro_arg eq $proj_choices[$ix][2]) {    # Test for WRF

          if (Exists($maproj_mb)) { set_proj_vars($ix); }

          last;                                  #if found then exit loop
      }
      $ix++;
    }


    
}
# ---------------------------------------
# sync_nl_vars 
#
# Set model namelist vars from the app variables.
# ---------------------------------------
sub sync_nl_vars {

  if($model_name eq "WRF"){

     $nl_var{NUM_DOMAINS}[0]=$num_nests;
     $nl_var{MOAD_KNOWN_LAT}[0]=$grid_cen_lat_cmn;
     $nl_var{MOAD_KNOWN_LON}[0]=$grid_cen_lon_cmn;
     $nl_var{MOAD_STAND_LONS}[0]=$stdlon;
     $nl_var{MOAD_STAND_LATS}[0]=$truelat1;
     $nl_var{MOAD_STAND_LATS}[1]=$truelat2;
     $nl_var{MOAD_DELTA_X}[0]=($grid_spacing_dkm*1000);
     $nl_var{MOAD_DELTA_Y}[0]=($grid_spacing_dkm*1000);
     $nl_var{XDIM}[0]=$nx_ddim;
     $nl_var{YDIM}[0]=$ny_dim;
     $nl_var{PTOP_PA}[0]=$ptop_mb*100;

     if ($num_nests==1){
        $nl_var{DOMAIN_ORIGIN_URI}[0]=$nx_ddim;
        $nl_var{DOMAIN_ORIGIN_URJ}[0]=$ny_dim;
     }

     # Set proj widget sets MAP_PROJ_NAME
     $nl_var{MAP_PROJ_NAME}[0]="$proj_choices[$p_index][2]";

# RAR - The following lines will have to be modified as more dynamic cores
#       are added to the system.
#
     $nl_var{OUTPUT_COORD}[0]="\'ETAP\'";
     $nl_var{OUTPUT_COORD}[0]="\'NMMH\'" if $projection_type eq 'RL';
    #$nl_var{OUTPUT_COORD}[0]="\'NMMH\'" if $nl_var{MAP_PROJ_NAME}[0] =~ /rotlat/i;
#
# RAR END
  
     # Deal with nests.
     if ($num_nests > 1) { 
         
        #sync_nl_vars_nest();  # Applies to only one of the indices, not all.
        reduce_nl_maxs_nest();
        set_nest_index($num_nests);
        sync_nl_maxs_nest();
     }
     configure_domain_id_mb();

  } elsif($model_name eq "LAPS"){

     # Set all the app variables from model namelist.
     $nl_var{GRID_CEN_LAT_CMN}[0]=$grid_cen_lat_cmn;
     $nl_var{GRID_CEN_LON_CMN}[0]=$grid_cen_lon_cmn;
     $nl_var{STANDARD_LONGITUDE}[0]=$stdlon;
     $nl_var{STANDARD_LATITUDE}[0]=$truelat1;
     $nl_var{STANDARD_LATITUDE2}[0]=$truelat2;
     $nl_var{GRID_SPACING_M_CMN}[0]=$grid_spacing_m_cmn;
     $nl_var{GRID_SPACING_M_CMN}[0]=($grid_spacing_dkm*1000);
     $nl_var{NX_L_CMN}[0]=$nx_dim;
     $nl_var{NY_L_CMN}[0]=$ny_dim;
     $nl_var{NK_LAPS}[0]=$nl_var_max{PRESSURES};

     # Set proj widget sets MAP_PROJ_NAME
     $nl_var{C6_MAPROJ}[0]="\'$proj_choices[$p_index][1]\'";
  }
}

# ---------------------------------------
# assign_nl_vars_nmm 
#
# Set app variables from NMM model namelist variables,
# particularly change grid_spacing from degrees to km. 
# ---------------------------------------
sub assign_nl_vars_nmm {

     my $x=$nl_var{MOAD_DELTA_X}[0] / 0.006500;
     $x=sprintf ("%.1f",$x);
     $dx_orig=$grid_spacing_dkm=$x;
     $nl_var{MOAD_DELTA_Y}[0]=$nl_var{MOAD_DELTA_X}[0]=$grid_spacing_m_cmn=$x*1000.;
     set_dx();

     display_degrees();

     # Turn off the button that indicates that a value has changed,
     # we forced a change don't want to be prevented from advancing.
     $b_update_highlight=0;
}

# ---------------------------------------
# sync_nl_vars_nmm 
#
# Set NMM model namelist vars from the app variables,
# particularly change grid_spacing from km to degrees. 
# ---------------------------------------
sub sync_nl_vars_nmm {

     $nl_var{MOAD_DELTA_X}[0]=$nx_degrees;
     $nl_var{MOAD_DELTA_Y}[0]=$ny_degrees;

}


# --------------------------------------
# assign_nl_vars_nest {
#
# Set app variables from model namelist vars.
# --------------------------------------
sub assign_nl_vars_nest {

  if($model_name eq "WRF"){
     
     $parent_of_nest=$nl_var{PARENT_ID}[$nest_index];
     set_nest_ratio($nl_var{RATIO_TO_PARENT}[$nest_index]);
     $nest_lli=$nl_var{DOMAIN_ORIGIN_LLI}[$nest_index];
     $nest_llj=$nl_var{DOMAIN_ORIGIN_LLJ}[$nest_index];
     $nest_uri=$nl_var{DOMAIN_ORIGIN_URI}[$nest_index];
     $nest_urj=$nl_var{DOMAIN_ORIGIN_URJ}[$nest_index];
  
     # Activate the nest widgets, too.
     activate_nest_widgets();
     
  } elsif($model_name eq "LAPS"){
  } 

}

# --------------------------------------
# sync_nl_vars_nest {
#
# Set model namelist vars from the app variables.
# --------------------------------------
sub sync_nl_vars_nest {

  if($model_name eq "WRF"){
     
    $nl_var{PARENT_ID}[$nest_index]=$parent_of_nest;
    $nl_var{RATIO_TO_PARENT}[$nest_index]=$nest_ratio;
    $nl_var{DOMAIN_ORIGIN_LLI}[$nest_index]=$nest_lli;
    $nl_var{DOMAIN_ORIGIN_LLJ}[$nest_index]=$nest_llj;
    $nl_var{DOMAIN_ORIGIN_URI}[$nest_index]=$nest_uri;
    $nl_var{DOMAIN_ORIGIN_URJ}[$nest_index]=$nest_urj;

    fill_nest_table_entry();
    fill_nl_disp_ary();

  } elsif($model_name eq "LAPS"){
  }

}

# --------------------------------------
# fill_nl_disp_ary 
#
# For display purposes make a single variable, 
# nl_var{$key}[0], contain all array information,
# kinda bogus way to make multiple values appear 
# in one entrybox on Localization Parms panel 
# but it works.
# --------------------------------------
sub fill_nl_disp_ary {

    if ($model_name eq "WRF") {

      # Fill Nest Declaration variables.
      #-------------------------
      for $key (@nesting_vars) {

        # Fill first element of array.
	#------
        if ($key eq 'DOMAIN_ORIGIN_URI') { 
           # Set element [0] = xdim.
           @nl_disp_ary{$key}=$nl_var{$key}[0]=$nl_var{XDIM}[0];
        } elsif ($key eq 'DOMAIN_ORIGIN_URJ'){ 
           # Set element [0] = ydim.
           @nl_disp_ary{$key}=$nl_var{$key}[0]=$nl_var{YDIM}[0];
        } else {
           # All remaining element [0]'s are=1.
           @nl_disp_ary{$key}=1; 
        }

        # Fill remaining array element(s).
	#------
        for $count (1 .. $num_nests) {
          if($nl_var{$key}[$count] ne "") {
            @nl_disp_ary{$key}="$nl_disp_ary{$key}, $nl_var{$key}[$count]";
          }
        }
        $nl_var{$key}[0]=$nl_disp_ary{$key}; 
      }

      # Filter Subnests input.
      #-------------------------
      filter_active_subnest_input();
    } 

} 

# --------------------------------------
# hide_nl_disp_ary
#
# Hide nl_disp_ary.
# --------------------------------------
sub hide_nl_disp_ary {
    my ($filter_flag)=@_;

      # Fill Nest Declaration variables.
      #-------------------------
      for $key (@nesting_vars) {
         @nl_disp_ary{$key}=$nl_var{$key}[0]; 
         if ($key eq 'DOMAIN_ORIGIN_URI') { 
            $nl_var{$key}[0]=$nl_var{XDIM}[0];
         } elsif ($key eq 'DOMAIN_ORIGIN_URJ'){ 
            $nl_var{$key}[0]=$nl_var{YDIM}[0];
         } else {
            $nl_var{$key}[0]=1; 
         }
      }

      # Filter Subnests input.
      #-------------------------
      if ($filter_flag) {
         filter_active_subnest_input();
      }
}

# --------------------------------------
# show_nl_disp_ary
#
# Show nl_disp_ary.
# --------------------------------------
sub show_nl_disp_ary {

      # Fill Nest Declaration variables.
      #-------------------------
      for $key (@nesting_vars) { $nl_var{$key}[0]=$nl_disp_ary{$key}; } 
}

# --------------------------------------
# filter_active_subnest_input
#
# Sort and filter the nl_disp_ary array.
# Fill Active Subnests variable.
# --------------------------------------
sub filter_active_subnest_input {

    my $key='ACTIVE_SUBNESTS';
    my @input;

    if ($nl_var{ACTIVE_SUBNESTS}[0] =~ m/,/) {
       $nl_var{ACTIVE_SUBNESTS}[0] =~ s/\s+//g; #clear out all whitespace
       $nl_var{ACTIVE_SUBNESTS}[0] =~ s/,/ /g;
       @input = split / /, $nl_var{ACTIVE_SUBNESTS}[0];
    } else {
       @input = $nl_var{ACTIVE_SUBNESTS}[0];
    }

    # Sort list into descending order.
    # ---
    @active_subnests=sort @input;

    # Clean up list.
    # ---
    @input=@active_subnests;

    my $i=0;
    my $j=0;
    my $my_max=scalar(@active_subnests);
    
    # Check the size of array.
    if ($my_max > 1) {
       for my $ith (0 .. $my_max-1) {
          # Note: $num_nests=$nl_var{NUM_DOMAINS}[0];
          if ($active_subnests[$ith] > $num_nests || 
              $active_subnests[$ith] < 2 ) {

             # Delete input less than 0 or greater than 1.
             splice(@input, $i-$j, 1);
             $j++;

          } elsif ($active_subnests[$ith] == $active_subnests[$ith-1] ) {

             # Delete duplicate input.
             splice(@input, $i-$j, 1);
             $j++;
          }
          $i++;
       } 

    } elsif ($my_max == 1) {
       # Single value.
       @input=@active_subnests;

    } 

    if ($input[0] == "" || $input[$ith] > $num_nests || $input[$ith] < 2 ) {
       # Default value=2, if otherwise null.
       @input=2; 
    }


    # Assign filtered values to var.
    # ---
    my $storage_bin;
    foreach (@input) {
	$storage_bin = "$storage_bin, $_";
    }
    $storage_bin =~ s/^\s+//;
    $storage_bin =~ s/^,\s+//;
    $nl_var{$key}[0] = $storage_bin;

} 

# --------------------------------------
# reduce_nl_maxs_nest {
#
# Set the model namelist MAX from the nest index variables.
# --------------------------------------
sub reduce_nl_maxs_nest {

    # If array contents are deleted, then reduce num_nest by one.
    my $my_max=$num_nests-1;
    for ($i=$my_max; $i > 1; $i--) {
       if ($nl_var{PARENT_ID}[$i] eq "" ||
           $nl_var{DOMAIN_ORIGIN_LLI}[$i] eq "" ||
           $nl_var{DOMAIN_ORIGIN_LLJ}[$i] eq "" || 
           $nl_var{DOMAIN_ORIGIN_URI}[$i] eq "" || 
           $nl_var{DOMAIN_ORIGIN_URJ}[$i] eq "" ||

           $nl_var{DOMAIN_ORIGIN_LLI}[$i] <= $grid_buffer ||
           $nl_var{DOMAIN_ORIGIN_LLJ}[$i] <= $grid_buffer || 
           $nl_var{DOMAIN_ORIGIN_URI}[$i] <= $grid_buffer || 
           $nl_var{DOMAIN_ORIGIN_URJ}[$i] <= $grid_buffer) { 
          $num_nests=$i;
       }
    }

  #  sync_nl_maxs_nest();
}

# --------------------------------------
# sync_nl_maxs_nest {
#
# Set the model namelist MAX from the nest index variables.
# --------------------------------------
sub sync_nl_maxs_nest {

    if($num_nests <= $nest_id_num) {
       $num_nests=$nest_id_num;
       for $key (@nesting_vars) {
          $nl_var_max{$key}=$nest_id_num;
       }
       print " called sync_nl_MAXs_nest to set new NUM_DOMAINS $nest_id_num\n\n"
	if $Debug;
    }
    update_nest_id_mb();
}

#__________________________________________________
#
# Write files and create directories.
#__________________________________________________


# --------------------------------------
# write_domain_files_and_image 
#
# Write the two domain files to the templates dir:
#          wrfsi.nl -or- nest7grid.parms 
#          dataroot.txt
# --------------------------------------
sub write_domain_files_and_image {

   # Write wrfsi.nl -or- nest7grid.parms 

   if($domain_select eq 'default') {
        $hint_msg= "You CANNOT write to the 'default' domain.";
        return;
   }

   # Create template if first time.
   create_template_dir();      

   $hint_msg="...writing namelist.";
   if (write_namelist()) { return(1); }

   # Create template if first time.
   write_domain_image();

   # Add domain_select to list of domains.
   add_domain_to_list();

   # Set $domain_current=$domain_select;
   $domain_current=$domain_select;

   #fill_nest_table_entry();
   if ($model_name ne "LAPS") { 
      write_nest_table(); 
   } else {
      write_pressures_namelist(1); 
   }

}

# --------------------------------------
#  write_dataroot_txt
#
# Write dataroot file to store the domain's DATAROOT.
# --------------------------------------
sub write_dataroot_txt  {

   my $data_path_file="$ROOT_TEMPLATES/$domain_select/dataroot.txt";
   # Don't ask if (-e $data_path_file)
   open(DPF,">$data_path_file") or 
      warn("Can't open file: $data_path_file.\n"), return(1);
   print DPF "$dataroot_select\n";
   close(DPF);

}

# --------------------------------------
#  add_domain_to_list
#
# Add domain_select to list of domains, if domain is brand new 
# (and insert alphabetically) otherwise, return.
# --------------------------------------

sub add_domain_to_list  {

   my $x1=0;
   $after_entry=""; # Global.
   foreach (@domain_found) {
      if ($domain_select eq $_){ 
         return; 
      } elsif ($domain_select gt $_){ 
         # Alphabetical sorting.
         $after_entry=$_;
      } else {
         last;
      }
      $x1++;
   }

   # Similar to find_existingDomains 
   push(@dfound_previous,$domain_select);
   flag_existingLocalizations($domain_select);
   @domain_found=sort(@dfound_previous);
   @dfound_previous=@domain_found;
  
   # Move scrollbar to show domain_select.
   my $x2=scalar(@domain_found);
   $domain_lb->yview(moveto => ($x1/$x2));
}

# --------------------------------------
# create_template_dir
#
# Create template directory in which to write namelist for domain.
# --------------------------------------
sub create_template_dir {
  no strict 'refs';

  my $droot="$ROOT_TEMPLATES/$domain_select";
  if (-d "$droot") {return;}

  # Create new directory.
  if (-d "$ROOT_TEMPLATES") {
     if (!mkdir $droot, 0777) {
        info_dbox("Make Directory Error", 
        "Cannot make $droot.\nPermission denied?!\nRun chmod in another terminal window, then press 'Next>' again."); 
        return(1);
     }

  } else {
     info_dbox("Domain Error", "Domain path $ROOT_TEMPLATES does not exist."); 
     return(1);
  }
}

# --------------------------------------
# write_namelist 
#
# Write namelist for domain in TEMPLATE.
# See also below, write_dataroot_namelist.
# --------------------------------------
sub write_namelist {
  no strict 'refs';

  # Check for valid domain_select.
  if($domain_select eq "") {
     $hint_msg= "You CANNOT write $namelist_arg to the a 'null' domain.";
     return(1);
  } elsif($domain_select eq 'default') {
     $hint_msg= "You CANNOT write to $namelist_arg the 'default' domain.";
     return(1);
  }

  # Sync the variables.
  sync_nl_vars();
  if ($nmm) {sync_nl_vars_nmm();}

  # Compare all variables to the original values.
  # Store new values in array "changes" for printing.
  # ---------------------------------
  my @changes;
  my $count=0;
  for $i (1 .. $num_nl_sections) {

     for $j (1 .. $num_nl_entries[$i]) {
        my $key = $nl_var_array[$i][$j];
 
        #print " key: $key\t";
        #print  " $nl_var_max{$key}\n";

        for ($k=0; $k < $nl_var_max{$key}; $k++) {

           if ($nl_var_orig{$key}[$k] =~ m/^'/ ) { # If orig entry has quote,
             $nl_var{$key}[$k] =~ s/'//g;          # then strip off existing
             $nl_var{$key}[$k] ="\'$nl_var{$key}[$k]\'"; # and add new quotes.
           }

           # Compare.
           if ($nl_var{$key}[$k] ne $nl_var_orig{$key}[$k]) {
        ###      print "$k: $nl_var{$key}[$k] ne $nl_var_orig{$key}[$k]\n";
              ($changes[$count][0], $changes[$count][1])=($i, $j);
              $count++;
              last; # Count the $key only once regardless of how
                    # many elements or its $nl_var_max{$key} change.
        ###   } else {
        ###      print "$k: $nl_var{$key}[$k] eq $nl_var_orig{$key}[$k]\n";
           }
        } #$k
     } #$j
  }; #$i

  # Hide all nest values by calling sub nl_disp_ary.
  if (!$nmm and $model_name eq "WRF") { hide_nl_disp_ary(1); }

  # Write "changes" to domain template namelist file (TNL).
  # ---------------------------------
  my $val; 
  if ($count <= 0) { 
      $reason_to_update_loc=0;
      return(1);
  } else {
      $reason_to_update_loc=1;

      $template_nl="$ROOT_TEMPLATES/$domain_select/$namelist_arg"; 
      open(TNL,">$template_nl") or 
          warn("First create dir for $template_nl.\n"), return(1);

      print TNL "&$nl_section[ $changes[0][0] ]\n"; #section header
      for ($c=0; $c < $count; $c++) {
        ($i, $j)=($changes[$c][0], $changes[$c][1]); 
    
        if ($changes[$c][0] ne $changes[$c-1][0]) {
           if ($c>0) { 
              print TNL "/\n"; 
              print TNL "&$nl_section[$i]\n"; #section header
           }
        }

        $key = $nl_var_array[$i][$j]; 
        for ($k=0; $k < $nl_var_max{$key}; $k++) {
           $val=$nl_var{$key}[$k]; 
           if ($k<1) { 
              print TNL " $key = $val"; 
           } elsif ($k%5 == 0 && $key eq 'LEVELS') {
              print TNL ", \n\t$val"; 
           } else {
              print TNL ", $val"; 
           }
     
        } 
        print TNL "\n"; 
      } 
      print TNL "/\n"; #section close
      close(TNL);

  } #$count


  # Show all nest values by calling sub nl_disp_ary.
  if ($nmm) {
     assign_nl_vars_nmm();
  } elsif ($model_name eq "WRF") { 
     show_nl_disp_ary(); 
  }

  $hint_msg= "Successfully wrote namelist to $template_nl.";

  # Write dataroot.txt.
  write_dataroot_txt();

  # Write pressures.nl for LAPS.
  if ($model_name eq "LAPS") { write_pressures_namelist(0); }


 # Success.
 return(0);
}

# --------------------------------------
# write_full_namelist 
#
# Write namelist for domain in TEMPLATE.
# See also below, write_dataroot_namelist.
# --------------------------------------
sub write_FULL_namelist {
  no strict 'refs';

  # Check for valid domain_select.
  if($domain_select eq "") {
     $hint_msg= "You CANNOT write $namelist_arg to the a 'null' domain.";
     return(1);
  } elsif($domain_select eq 'default') {
     $hint_msg= "You CANNOT write to $namelist_arg the 'default' domain.";
     return(1);
  }

  # Write all vars to domain template namelist file (TNL).
  # ---------------------------------
  my $val;
      $reason_to_update_loc=1;

      $domain_nl_dir="$ROOT_DATA/$domain_select/static"; 
      if (!-d $domain_nl_dir) { mkdir -p "$domain_nl_dir", 0777 
         or warn("Could not mkdir $domain_nl_dir\n"), 
         $hint_msg="Could not mkdir $domain_nl_dir\n", return(1); }


      $domain_nl="$ROOT_DATA/$domain_select/static/$namelist_arg"; 
      open(TNL,">$domain_nl") or 
          warn("First create dir for $domain_nl.\n"), return(1);

      for $i (1 .. $num_nl_sections) {
       print TNL "&$nl_section[$i]\n"; #section header
       for $j (1 .. $num_nl_entries[$i]) {

        my $key = $nl_var_array[$i][$j];
        for ($k=0; $k < $nl_var_max{$key}; $k++) {


         $key = $nl_var_array[$i][$j]; 
         for ($k=0; $k < $nl_var_max{$key}; $k++) {
           $val=$nl_var{$key}[$k]; 
           if ($k<1) { 
              print TNL " $key = $val"; 
           } elsif ($k%5 == 0 && ($key eq 'LEVELS')) {
              print TNL ", \n\t$val"; 
           } else {
              print TNL ", $val"; 
           }
     
         } 
         print TNL "\n"; 
        } 
       } 
       print TNL "/\n"; #section close
      } 
      close(TNL);


  $hint_msg= "Successfully wrote FULL namelist to $domain_nl.";

 # Success.
 return(0);
}

# --------------------------------------
# write_textWindow_namelist 
#
# Write namelist for domain in DATAROOT.
# --------------------------------------
sub write_textWindow_namelist {

      # Write template namelist file (TNL).
      my $usr_input= $nl_text->get("1.0", "end");
      $template_nl="$ROOT_TEMPLATES/$domain_select/$namelist_arg"; 
      open(TNL,">$template_nl") or 
          warn("First create dir for $template_nl.\n"), return(1);
      print TNL "$usr_input";
      close(TNL);
      $hint_msg= "Successfully wrote namelist to $template_nl.";

      # Write dataroot namelist file (DNL).
      $domain_nl="$ROOT_DATA/$domain_select/static/$namelist_arg"; 
      if (-e "$domain_nl") {
        open(DNL,">$domain_nl") or 
          warn("First create dir for $domain_nl.\n"), return(1);
        print DNL "$usr_input";
        close(DNL);
        $hint_msg= "$hint_msg \nAlso wrote namelist to $domain_nl.";
      }
}

# --------------------------------------
# write_dataroot_namelist 
#
# Write namelist for domain in DATAROOT.
# --------------------------------------
sub write_dataroot_namelist {
  no strict 'refs';

  # Check for valid domain_select.
  if($domain_select eq "") {
     $hint_msg= "You CANNOT write $namelist_arg to the a 'null' domain.";
     return(1);
  } elsif($domain_select eq 'default') {
     $hint_msg= "You CANNOT write to $namelist_arg the 'default' domain.";
     return(1);
  }

  # Sync the variables.
  sync_nl_vars();
  if ($nmm) {sync_nl_vars_nmm();}

  # Store new values in array "changes" for printing.
  # ---------------------------------
  my $count=0;
  my @changes;
  for $i (1 .. $num_nl_sections) {

     for $j (1 .. $num_nl_entries[$i]) {
        my $key = $nl_var_array[$i][$j];
        for ($k=0; $k < $nl_var_max{$key}; $k++) {

           if ($nl_var_orig{$key}[$k] =~ m/^'/ ) { # If orig entry has quote,
             $nl_var{$key}[$k] =~ s/'//g;          # then strip off existing
             $nl_var{$key}[$k] ="\'$nl_var{$key}[$k]\'"; # and add new quotes.
           }
           ($changes[$count][0], $changes[$count][1])=($i, $j);
           $count++;
           last; # Count the $key only once regardless of how
                 # many elements or its $nl_var_max{$key} change.

        } #$k
     } #$j
  }; #$i

  # Hide nl_disp_ary.
  if (!$nmm and $model_name eq "WRF") { hide_nl_disp_ary(0); }

  # Write "changes" to domain template namelist file (DNL).
  # ---------------------------------
  my $val;
  if ($count <= 0) { 
      $reason_to_update_loc=0;
      return(1);
  } else {
      $reason_to_update_loc=1;

      $domain_nl="$dataroot_select/$domain_select/static/$namelist_arg"; 
      open(DNL,">$domain_nl") or 
          warn("First create dir for $domain_nl.\n"), return(1);

      print DNL "&$nl_section[ $changes[0][0] ]\n"; #section header
      for ($c=0; $c < $count; $c++) {
        ($i, $j)=($changes[$c][0], $changes[$c][1]); 
    
        if ($changes[$c][0] ne $changes[$c-1][0]) {
           if ($c>0) { 
              print DNL "/\n"; 
              print DNL "&$nl_section[$i]\n"; #section header
           }
        }

        $key = $nl_var_array[$i][$j]; 
        for ($k=0; $k < $nl_var_max{$key}; $k++) {
           $val=$nl_var{$key}[$k]; 
           if ($k<1) { 
              print DNL " $key = $val"; 
           } elsif ($k%5 == 0 && $key eq 'LEVELS') {
              print DNL ", \n\t$val"; 
           } else {
              print DNL ", $val"; 
           }
     
        } 
        print DNL "\n"; 
      } 
      print DNL "/\n"; #section close
      close(DNL);
  } #$count

  # Show all nest values by calling sub nl_disp_ary.
  if ($nmm) {
     assign_nl_vars_nmm();
  } elsif ($model_name eq "WRF") { 
     show_nl_disp_ary(); 
  }

  $hint_msg="$hint_msg\nand to DATAROOT $domain_nl.";

 # Success.
 return(0);
}

# -----------------------------------
# write_domain_image
#
# -----------------------------------
sub write_domain_image {


      if($domain_select eq 'default') {
           $hint_msg= "You CANNOT write to the 'default' domain.";
           return;
      }

      # Is the canvas viewable to the screen.
      if (!$can->viewable && $exit_now) {
          present_tool(1);
          #print "TRYING TO PRINT IMAGE\n";
          $panel_index=2;
          raise_panel(0);

          # Update idletasks.
          $mw->update();
          $mw->idletasks();
      }

      # Is the canvas mapped to the screen.
      if (!$can->ismapped) {return(1);}

      my $ext="gif";
      my $img_output="$ROOT_TEMPLATES/$domain_select/domain.$ext"; 

      # If horizontal variables change or we are copying a domain, 
      # then save a new image.
      my $must_write_img=0;
   
      if($domain_mode ==3 && $domain_select ne $domain_previous) { 
         $domain_previous=$domain_select;
         $must_write_img=1;
   
      } else {
         # If the variables to generate a map change, then
         # must 'Update Map' in order to write_domain_image.
         $must_write_img=check_gen_map_vars(1);
      }

      if (!$must_write_img) { return; };

      # Write image.
      watch_cursor(1);
      $hint_msg="...writing namelist & '$img_output'."; 
      
      # Make bbox white.
      $grid_val_restrict=1;
      restrict_grid_var_calc();

      hide_tags(); 
      hide_tags_nest(); 
      remove_parentbox();

      $mw->update;
      my $win_img=$can->Photo(#-format=>'Window', -width=>$xmax, -height=>$ymax,
                              -data => oct($can->id));
      my $bd=3;
      my $xmax=$win_img->width  - $bd;
      my $ymax=$win_img->height - $bd;
      $bd=2;
      $win_img->write($img_output, 
                               -format => $ext, 
                               -from => $bd,$bd,$xmax, $ymax);
      watch_cursor(0);
      update_tags(); 
      $hint_msg="Successfully wrote '$img_output'.";

      #---

      # A change was made to namelist. 
      # If localization is present, it could be out of sync - so replace checkmark.
      if ($domain_mode == 2 && $domain_lb->indicatorExists($domain_select) ) {
         $domain_lb->indicatorDelete($domain_select);
	 add_checkmark();
      }
  

if(0){
      # Write Cylindrical Equidistant projection bbox, only for a new
      # domain or an existing domain with existing orig_bbox.txt file.
      if ($domain_mode == 1 || $orig_bbox_exists) { 
         my $orig_bbox="$ROOT_TEMPLATES/$domain_select/.orig_bbox.txt"; 
         open(CYL,">$orig_bbox") or warn("Can't open file: $orig_bbox\n"), return;
         print CYL "$orig_bbox_L\n";
         print CYL "$orig_bbox_T\n";
         print CYL "$orig_bbox_R\n";
         print CYL "$orig_bbox_B\n";
         close(CYL);
      }
}

}

# --------------------------------------
# write_pressures_namelist 
#
# Write pressures namelist for LAPS domain in TEMPLATE.
# --------------------------------------
sub write_pressures_namelist {
  my ($count) = @ARG;
  no strict 'refs';

  # Compare all variables to the original values.
  # Store new values in array "changes" for printing.
  # ---------------------------------
  my $key="PRESSURES";

  for ($k=0; $k < $nl_var_max{$key}; $k++) {

     # Compare level totals.
     if ($nl_var_max_orig{$key} ne $nl_var_max{$key}) {
        $nl_var_max_orig{$key}=$nl_var_max{$key};
        $count++;
     } elsif ($nl_var{$key}[$k] ne $nl_var_orig{$key}[$k]) {
        # Compare level values.
        $count++;
        last; # Count the $key only once regardless of how
              # many elements or its $nl_var_max{$key} change.
     }

  }

  # Write "changes" to domain template namelist file (PNL).
  # ---------------------------------
  my $val;
  if ($count <= 0) { 
      $reason_to_update_loc=0;
      return(1);
  } else {
      $reason_to_update_loc=1;

      $press_nl="$ROOT_TEMPLATES/$domain_select/pressures.nl"; 
      open(PNL,">$press_nl") or 
          warn("First create dir for $press_nl.\n"), return(1);

      print PNL "&pressures_nl\n"; #section header
      print PNL "pressures=\n"; 
    
      for ($k=0; $k < $nl_var_max{$key}; $k++) {
         $val=$nl_var{$key}[$k]; 
         print PNL "$val, \n"; 
      } 
      print PNL "/\n"; #section close
      close(PNL);

  } #$count


  $hint_msg= "Successfully wrote pressures_namelist to $press_nl.";


 # Success.
 return(0);
}


# -----------------------------------
# load_orig_bbox
#
# A original Cylindrical Equidistant projection bbox
# is saved to a file, to allow a user to re-edit the 
# center point lat/lon values.
# -----------------------------------
sub load_orig_bbox {
return;

      # Deactivate startover_button.
      set_button_state(0,$b_startover);
 
      $orig_bbox_exists=0; 

      # Load Cylindrical Equidistant projection bbox.
      my $orig_bbox="$ROOT_TEMPLATES/$domain_select/.orig_bbox.txt"; 

      if (-e $orig_bbox) {
         # orig_bbox file exist.
         open(CYL,"$orig_bbox") or warn("Can't open file: $orig_bbox\n"), return;
         my @line = <CYL>; 
         close(CYL);

         $orig_bbox_exists=1; 
         chomp (@line);
         ($bbox_L, $bbox_T, $bbox_R, $bbox_B)=@line;
         if ($bbox_L == $bbox_R || $bbox_T == $bbox_B) {return;}

         # Activate startover_button.
         set_button_state(1,$b_startover);
      }
}

# -----------------------------------
# check_gen_map_vars
#
# If the variables to generate a map change, then
# must 'Update Map' in order to write_domain_image.
# -----------------------------------
sub check_gen_map_vars{
      #my($write_nl)=@_;
      my($write_nl)=0; # This should always be 0, not 1 like I thought.
      my $change=0;

  # Sync the variables.
  sync_nl_vars();

   #--
      
      if ($model_name eq "WRF") {
         @img_parms=qw(MOAD_KNOWN_LAT 
                    MOAD_KNOWN_LON
                    MOAD_STAND_LONS 
                    MOAD_STAND_LATS
                    MOAD_DELTA_X 
                    MOAD_DELTA_Y
                    XDIM YDIM 
                    MAP_PROJ_NAME

                    NUM_DOMAINS
                    DOMAIN_ORIGIN_LLI DOMAIN_ORIGIN_LLJ 
                    DOMAIN_ORIGIN_URI DOMAIN_ORIGIN_URJ
                    PARENT_ID RATIO_TO_PARENT
                    );


      } elsif ($model_name eq "LAPS") {
         @img_parms=qw(GRID_CEN_LAT_CMN
             GRID_CEN_LON_CMN
             STANDARD_LONGITUDE
             STANDARD_LATITUDE
             STANDARD_LATITUDE2
             GRID_SPACING_M_CMN
             NX_L_CMN NY_L_CMN
             C6_MAPROJ);
      }

   #--

      # Loop thru all parms.
      foreach $kkey (@img_parms) {
       for ($kk=0; $kk < $nl_var_max{$kkey}; $kk++) { # Deal with array vars.
         if ($nl_var_prev{$kkey}[$kk] ne $nl_var{$kkey}[$kk]) {
             #print "$kkey, $nl_var_prev{$kkey}[$kk] NE $nl_var{$kkey}[$kk]\n";
             if($write_nl) {$nl_var_prev{$kkey}[$kk]=$nl_var{$kkey}[$kk];}
             $change=1;
         }
       }
      }

       
   #--

      # Test for change in bcd file.
      if ($bcd_prev ne $bcd_choice) { 
          if($write_nl) {$bcd_prev=$bcd_choice;}
          $change=1; 
      }
      
      return $change;
}

#__________________________________________________
#
# Trace variables in application.
#__________________________________________________


# ---------------------------------------
# trace_vars 
#
# Set a trace on textvariables.
# 'r' -fetch
# 'w' -store
# 'u' -destroy (But char_checker has: "return if $op eq 'u'").
#
# Note: everytime these values change the "Update Map" button
# is highlighted.
# ---------------------------------------
sub trace_vars {

$op='w';

$mw->traceVariable(\$nx_ddim, $op=> [\&char_checker, 'nx_ddim']);
$mw->traceVariable(\$ny_dim, $op=> [\&char_checker, 'ny_dim']);
$mw->traceVariable(\$stdlon, $op=> [\&char_checker, 'stdlon']);
$mw->traceVariable(\$truelat1, $op=> [\&char_checker, 'truelat1']);
$mw->traceVariable(\$truelat2, $op=> [\&char_checker, 'truelat2']);
$mw->traceVariable(\$geog_path, $op=> [\&char_checker, 'geog_path']);
$mw->traceVariable(\$domain_select, $op=> [\&char_checker, 'domain_select']);
$mw->traceVariable(\$dataroot_select, $op=> [\&char_checker, 'dataroot_select']);
$mw->traceVariable(\$grid_spacing_dkm, $op=> [\&char_checker, 'grid_spacing_dkm']);
$mw->traceVariable(\$grid_cen_lat_cmn, $op=> [\&char_checker, 'grid_cen_lat_cmn']);
$mw->traceVariable(\$grid_cen_lon_cmn, $op=> [\&char_checker, 'grid_cen_lon_cmn']);
$mw->traceVariable(\$nl_var{TOPTWVL_PARM_WRF}[0], $op=> [\&char_checker, 'TOPTWVL_PARM_WRF']);

#$main::mw->traceVariable(\$temp_k, $op=> [\&char_checker], 'temp_k');
#$mw->traceVariable(\$ptop_mb, $op=> [\&char_checker, 'ptop_mb']);
#$main::mw->traceVariable(\$v::pbot_mb, $op=> [\&char_checker], 'v::pbot_mb');

}

# -----------------------------------
# char_checker
#
# Traced variables are checked every time a user 
# presses a key. Only value that are valid are
# accepted.
# -----------------------------------
sub char_checker {
    use strict;
    no strict 'refs';

    my($index, $value, $op, $textvar) = @_;
    return if $op eq 'u';
    my($keypress) = ($Tk::event->A); # get event object key press.

    # Check for valid alpha numeric string.
    # ------------------------------------------

    if ($textvar eq "domain_select" || 
        $textvar eq "dataroot_select" ||
        $textvar eq "geog_path" ||
        $textvar eq "si_path") {

      if ($keypress =~ /\/|-|_|\.|\w/) { # slash, hyphen, dash, period, alpha numeric - okay 
      } elsif ($keypress eq "") {
      } elsif ($keypress eq "^") {
          $value =~ s/\^//; 
      } elsif ($keypress eq "\\") {
          $value =~ s/\\//; 
      } else {
          $value =~ s/[$keypress]//;   # Strip remaining non-numerics.
      }
 

    } else {

    # Check for valid string with one or more digits.
    # ------------------------------------------

      if ($value !~ /^([-+]?|[-+]?\d+\.?\d*|[-+]?\.?\d*)$/o) {
          # Matches a string staring with an optional +/- sign,
          # followed by one or more digits,
          # followed by an optional dot,
          # followed by zero or more digits.


         if ($value =~ /\D$/) {

            # strip off last char, it is NOT a numeric value. 
            $value =~ s/.$//; 

         } else {

            # find & strip off last char typed, it is NOT a numeric value. 
            if ($keypress eq "^") {
                $value =~ s/\^//; 
            } elsif ($keypress eq "\\") {
                $value =~ s/\\//; 
            } elsif($keypress ne "") {
                $value =~ s/[$keypress]//;   # Strip remaining non-numerics.
            }

         }
            

      };

      # Check bounds of certain variables.
      if ($textvar eq "grid_cen_lat_cmn") {
          if ($value >   90.0) { $value=  90.0; }
          if ($value <  -90.0) { $value= -90.0; }
          ${$textvar}=$value;
      } elsif ($textvar eq "grid_cen_lon_cmn" || $textvar eq "stdlon") {
          if ($value >  180.0) { $value= 180.0; }
          if ($value < -180.0) { $value=-180.0; }
          ${$textvar}=$value;
      } elsif ($textvar eq "nx_ddim" || $textvar eq "ny_dim") {
          if ($value >  99999) { $value= 99999; }
          if ($value <  0    ) { $value= 10; }
          ${$textvar}=$value;
      } elsif ($textvar eq "grid_spacing_dkm") {
          if ($value >  99999) { $value= 99999; }
          if ($value <  0    ) { $value= 0.001; }
      } elsif ($textvar eq "TOPTWVL_PARM_WRF") {
          #$value=int($value); 
          #if ($value <=  1   ) { $value= 1; 
          #} elsif ($value ==3) { $value= 2; 
          #} elsif ($value > 4) { $value= 4; 
	  #} 
      }

      # Update needs to be pressed, a parameter has changed.
      highlight_update_button(1);
      
    };

    return $value;
}

#__________________________________________________
#
# Find files in system directories.
#__________________________________________________


# --------------------------------------
# find_existingDomains 
#
# Find the existing domains by searching
# install_dir/templates for namelist files
# i.e. wrfsi.nl -or- nest7grid.parms.  
# --------------------------------------
sub find_existingDomains {
  
    my @temp;

    # Get model info from namelist.
    if(-e "$ROOT_TEMPLATES/$namelist_arg"){unlink "$ROOT_TEMPLATES/$namelist_arg";}
    #@domain_found=`find $ROOT_TEMPLATES -follow -name $namelist_arg`;
    @domain_found=`find $ROOT_TEMPLATES -name $namelist_arg`;
  
    if ( scalar(@domain_found) > 0) {
          # Strip off the full pathname and the sort list of domains.
          foreach (@domain_found) {
            s/\s+//g;
            s{$ROOT_TEMPLATES/}//g;
            s/\/$namelist_arg//g;
            push(@temp,$ARG);
          }
          @domain_found=sort(@temp);

    } else {
  
       info_dbox("Existing Domain Error", 
           "No existing domain files ($namelist_arg) were found in dir 
            $ROOT_TEMPLATES."); 
    }


    $after_entry="";
    if (@domain_found ne @dfound_previous) {
       flag_existingLocalizations(@domain_found);
       @dfound_previous=@domain_found;
    }
}

# --------------------------------------
# flag_existingLocalizations
#
# Find the existing domain localizations by 
# searching dataroot/domain for static files
# i.e. static.wrfsi -or- static.laps.
# --------------------------------------
sub flag_existingLocalizations {

   use vars qw/$after_entry/;
   my (@my_list)=@_;

   my $iii=0;
   foreach (@my_list) {
     if ($after_entry eq "") {
        # New list of templates to flag.
        $domain_lb->add($_, -text => $_, -at => $iii);
        $iii++;
     } else {
        # Individual template to flag.
        $domain_lb->add($_, -text => $_, -after => $after_entry);
     }

     # Does dataroot.txt file exist (containing path to dataroot).
     my $idata_path_file="$ROOT_TEMPLATES/$_/dataroot.txt";
     if (-e $idata_path_file) { 

         # Does dataroot path exists.
         chop ($idataroot_select=`cat $idata_path_file`);
         if (-e $idataroot_select) { 

             # Does localization exists (is static file presence).
             my $istatic_file_d=
                "$idataroot_select/$_/static/static.$window_domain_arg";

             if (-e $istatic_file_d) { 
                # File static.wrfsi.d01 exists.
               
                    # Valid and localized domain.
                    $domain_lb->indicator('create', $_,
                              -itemtype => 'image', 
                              -image => $mw->Bitmap(-data => $check_blk) );

             } else {
                # Valid but domain not localized.
             }

         } else {
             # Contents of dataroot.txt file are an invalid dataroot.
             $domain_lb->itemConfigure($_, 0, -style=>$grayfg);
             $nl_var{$_}=
"Dataroot file:
>$ROOT_TEMPLATES/$_/dataroot.txt<
contains invalid value: 
>$idataroot_select<.";

         }

     } else {
       # No dataroot.txt file.
       $domain_lb->itemConfigure($_, 0, -style=>$grayfg);
       $nl_var{$_}=
"The dataroot file does not exist
>$ROOT_TEMPLATES/$_/dataroot.txt<";

     } #if     
   } #foreach

}

# --------------------------------------
# find_bcdFiles
#
# Find the map background files by searching
# install_dir/gui/data/maps for .bcd (binary 
# cartographic data) files. e.g. Latlon10.bcd.
# --------------------------------------
sub find_bcdFiles {

  my @temp;

  # Get model info from namelist.
  @bcd_found=`find $GUI_MAP -name \*.bcd`;

  if ( scalar(@bcd_found) > 0) {
     # Strip off the full pathname and the sort list of domains.
     foreach (@bcd_found) {
       s/\s+//g;
       s{$GUI_MAP/}//g;
       s/\.bcd//g;
       push(@temp,$ARG);
     }
     @bcd_found=sort(@temp);

  } else {

     fail_dbox("Map Background Files Error", 
           "No bcd (binary cartographic data) files were found in dir 
            $GUI_MAP."); 
     return(1);
  }

  $bcdUSCounties="US_Counties";
  $bcdFileDefault="US_States";
  $bcdFileName= join ('', $GUI_MAP, '/', $bcdFileDefault, '.bcd');
  if (!-e $bcdFileName) { die "Background file $bcdFileName doesn't exist.\n"; }

  return(0);
}

# --------------------------------------
# add_checkmark
#
#
# --------------------------------------
sub add_checkmark {

    # (Re)set static_select (DATAROOT/$domain_select/static/static.xxxx).
    # -------------------------------
    $static_select=
      "$dataroot_select/$domain_select/static/static.$window_domain_arg";

    # Check for successful domain localization.
    # -------------------------------
    if(-e $static_select) {
       # Add checkmark.
       $domain_lb->indicator('create', $domain_select,
                             -itemtype => 'image', 
                             -image => $mw->Bitmap(-data => $check_blk) );
    } elsif ($domain_mode == 2 && $domain_lb->indicatorExists($domain_select)) {
       # Remove checkmark.
       $domain_lb->indicatorDelete($domain_select);
    }
}

### Return 1 to the calling use statement ###
1;

#__DATA__

$check_blk=<<'END';
/* XBM */
#define check_mrk_width 8
#define check_mrk_height 11
static unsigned char check_bits[] = {
   0x80, 0x80, 0x40, 0x40, 0x22, 0x23, 0x16, 0x16, 0x1c, 0x0c, 0x00};
END

$check_box=<<'END';
/* XBM */
#define check_box_width 9
#define check_box_height 12
static unsigned char check_box_bits[] = {
   0xff, 0x00, 0x81, 0x00, 0xc1, 0x00, 0xc1, 0x00, 0xa1, 0x00, 0xa3, 0x00,
   0x97, 0x00, 0x95, 0x00, 0x9d, 0x00, 0x89, 0x00, 0xff, 0x00, 0x00, 0x00};
END

$check_blu=<<'END';
/* XPM */
static char * check_blk_xpm[] = {
"8 11 3 1",
"       c None",
"+      c gray85",
".      c #1111ff",
"+++++++.",
"+++++++.",
"++++++.+",
"++++++.+",
"+.+++.++",
"..+++.++",
"+..+.+++",
"+..+.+++",
"++...+++",
"++..++++",
"++++++++"};
END

$check_blk2=<<'END';
/* XPM */
static char * check_blk2_xpm[] = {
"9 11 3 1",
" 	c None",
".	c gray85",
"+	c black",
".......+.",
"......++.",
"......+..",
".....++..",
".....+...",
".....+...",
".+..++...",
".++.+....",
"..+++....",
"...++....",
"....+...."};
END

$check_red=<<'END';
/* XPM */
static char * star_blk_xpm[] = {
"8 11 3 1",
" 	c None",
"+	c gray85",
".	c #FF0000",
"++++++++",
"++++++++",
".+++++.+",
"..+++.++",
"+..+.+++",
"++..++++",
"++...+++",
"+.++..++",
".++++..+",
"++++++++",
"++++++++"};
END

$check_red3=<<'END';
/* XPM */
static char * check_red_xpm[] = {
"9 11 3 1",
" 	c None",
".	c gray85",
"+	c #ff0000",
".......+.",
"......++.",
"......+..",
".....++..",
".....+...",
".....+...",
".+..++...",
".++.+....",
"..+++....",
"...++....",
"....+...."};
END

$my_watch=<<'END';
#define watch_width 16
#define watch_height 16
#define watch_x_hot 7
#define watch_y_hot 7
static char watch_bits[] = {
   0xf8, 0x07, 0xf8, 0x07, 0xf8, 0x07, 0xfc, 0x0f, 0x86, 0x18, 0x83, 0x30,
   0x81, 0xe0, 0xc1, 0xe1, 0xc1, 0xe1, 0x21, 0xe0, 0x13, 0x30, 0x06, 0x18,
   0xfc, 0x0f, 0xf8, 0x07, 0xf8, 0x07, 0xf8, 0x07};
END

$my_watch_msk=<<'END';
#define watch_width 16
#define watch_height 16
static char watch_bits[] = {
   0xf8, 0x07, 0xf8, 0x07, 0xf8, 0x07, 0xfc, 0x0f, 0xfe, 0x1f, 0xff, 0x3f,
   0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x3f, 0xfe, 0x1f,
   0xfc, 0x0f, 0xf8, 0x07, 0xf8, 0x07, 0xf8, 0x07};
END

