#!/bin/csh

# routine to localize satellite usage

# define Constants.inc

if ($1 == 8) then

cat g08coef.dat GOES8_IMAGER_ASCII_COEFF.DAT >! coef.dat

/usr/nfs/lapb/parallel/laps/bin/binary_coeff.x

cat g08_dry_control_file.dat GOES8_IMAGER_DRY_CONTROL_FILE.DAT >! Dry_Control_File.dat
cat g08_wet_control_file.dat GOES8_IMAGER_WET_CONTROL_FILE.DAT >! Wet_Control_File.dat
cat g08_ozo_control_file.dat GOES8_IMAGER_OZO_CONTROL_FILE.DAT >! Ozo_Control_File.dat



else if ($1 == 10) then

cat g10coef.dat GOES10_IMAGER_ASCII_COEFF.DAT >! coef.dat

/usr/nfs/lapb/parallel/laps/bin/binary_coeff.x

cat g10_dry_control_file.dat GOES10_IMAGER_DRY_CONTROL_FILE.DAT >! Dry_Control_File.dat
cat g10_wet_control_file.dat GOES10_IMAGER_WET_CONTROL_FILE.DAT >! Wet_Control_File.dat
cat g10_ozo_control_file.dat GOES10_IMAGER_OZO_CONTROL_FILE.DAT >! Ozo_Control_File.dat

endif
