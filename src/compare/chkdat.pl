#!/usr/local/bin/perl

$ah = shift || ($ah = 1);

$oldlapsdata = "/data/mdle/laps/nest7grid/lapsprd";
$newlapsdata = "/data/mdla/ROC_laps/data/lapsprd";

opendir(ND,$oldlapsdata);
@olddata = readdir(ND);
closedir(ND);

foreach $dir (@olddata){
    next if($dir =~ /^\.\.*$/);
    next unless($dir =~ /lga/);
#              || $dir =~ /vrc/
#              || $dir =~ /lvd/ );

    chdir($oldlapsdata);    
    stat($dir);
    next unless (-d _);
#    print "checking dir $dir\n";

#    next unless($dir =~ /lga/);

    chdir("$oldlapsdata/$dir");
    opendir(DD,".");
    @files = readdir(DD);
    close(DD);
    foreach $file (@files){
        next if($file =~ /^\.\.*$/);
        
        stat($file);
        $age = -M _;
#        print "Here $age $file\n";
	next unless ($age < $ah/24.);
        print "Comparing in $dir\n";
        
        $val = `cmp $file $newlapsdata/$dir/$file`;
        if($val != 0){
	    system("/data/mdla/newlaps/src/compare/compare.exe << EOF
                   1\n$file\n$newlapsdata/$dir/\n$oldlapsdata/$dir/\nEOF");
	    print "$val Data for $dir/$file differ\n" if($val !=0);
	}
    }
}







