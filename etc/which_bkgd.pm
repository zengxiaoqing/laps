package which_bkgd;
sub which_bkgd
{
    @lines=@_;

# Determine background model.

    $modeltype = "missing";

#Newer way with using 'get_modelfg_3d'
    $searchstring = "Successfully obtained: ";
    foreach (@lines)
    {
       if (/$searchstring(.*)/)
       {
           ($modeltime,$modeltype) = split(" ",$1);
           last;
       }
    }

    if ($modeltype eq "missing") {print "Can't determine model background.\n"}

    if ($modeltype eq 'fua' || $modeltype eq 'fsf')
    {
        $searchstring = "Reading 3d";
        $searchstring = "Reading 2d" unless ($modeltype eq 'fua');

        foreach (@lines)
        {
           if(/$searchstring(.*)/)
           {
               @string_elems=split(" ",$1);
               $bkgdpath = @string_elems[0];
               @string_elems = split("/",$bkgdpath);
               $modelid = @string_elems[$#string_elems];
               last;
           }
        }
    }

    print "BACKGROUND FIELDS:  ";
    if ($modeltype eq 'fua' || $modeltype eq 'fsf')
    {
# $modelid = "local";
    }
    else
    {
        $bgmodelfile = $filename;
        $bgmodelfile =~ s/$logname/lga/;
        $modelid = "???";
        if(open(BGMODEL,$bgmodelfile))
        {
           foreach(<BGMODEL>)
           {

#     if(/cmodel (\w*)/) {$modelid = $1;}

              if (/cmodel: (.*)/)           {(@modelid) = split(" ", $1); $modelid = @modelid[0];}
              if (/Reading - (.*)/)         {$pathname = $1; last; close BGMODEL; }
              if (/reading cdfname\: (.*)/) {$pathname = $1; last; close BGMODEL; }
           }
           $basename = $1 if ($pathname =~ /([^\/]*)$/);
        }
        if($basename =~ /(\d\d\d\d\d)(\d\d\d\d)(\d\d)(\d\d)/)
        {
           $runtime = $2/100;
           $fcsthr = $4;
#   print "Using $fcsthr hr fcst from $modelid started at $runtime UTC \n";
        }
        elsif($basename =~ /(\d\d\d\d)(\d\d)(\d\d)\_(\d\d)(\d\d)/)
        {
              $runtime = $4;
#   print "Using grids from $modelid started at $runtime UTC \n";
        }

    }

    if($modeltime =~ /(\d\d)(\d\d\d)(\d\d\d\d)(\d\d\d\d)/ )
    {
       $runtime = $3;
       $fcsthr = $4;
       $lengthmdlid=length $modelid;
#      print "Using $fcsthr fcst from $modelid model started at $runtime UTC \n";
#      print "$modeltype $modelid $1$2$3$4 \n";
    }

    return ($modelid,$modeltype,$runtime,$fcsthr,$1,$2);

}
1;
