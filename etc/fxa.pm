package fxa;
sub Get_env'fxa{ #'
    local($fxa_env_file) = "/awips/fxa/.environs";
    open(FXA,"$fxa_env_file") || return 0;
    @fxaenv = <FXA>;
    close(FXA);
    $fxa_env_file.=".".`hostname`;
    chomp($fxa_env_name);
    open(FXA,"$fxa_env_file");
    push(@fxaenv,<FXA>);
    close(FXA);
    foreach(@fxaenv){
	next if /^#/;
	next unless(/\s*([^\s]+)\s+([^\s]+)\s*$/);
        $evar = $1;
        $eval = $2;
        $eval =~ s#\$\{(.*)\}#$ENV{$1}#g;
	$eval =~ s#\$([^\s\/]+)#$ENV{$1}#g;    
        $ENV{$evar} = $eval;
    }
    close(FXA);
    return 1;
}

1;
sub Set_logdir'fxa{
    
    if( -d $ENV{LOG_DIR}){
#	local($yymmdd) = `date -u +%y%m%d`;
#	chomp($yymmdd);
        local(@gmtime) = gmtime;
        $gmtime[4]++;
        for($i=3;$i<=5;$i++){
          $gmtime[$i]="0".$gmtime[$i] if(length($gmtime[$i]<2));
        }
        local($yymmdd) = $gmtime[5].$gmtime[4].$gmtime[3];
#        print "hera $yymmdd\n"; chomp($yymmdd); print $yymmdd;


	$LAPS_LOG_PATH= "$ENV{LOG_DIR}/$yymmdd";

	if(! -d $LAPS_LOG_PATH){
	    mkdir $LAPS_LOG_PATH, 0777 || 
                die "Could not create Log directory $LAPS_LOG_PATH";
	}
	$LAPS_LOG_PATH .= "/laps";
	if(! -d $LAPS_LOG_PATH){
	    mkdir $LAPS_LOG_PATH, 0777 || 
                die "Could not create Log directory $LAPS_LOG_PATH";
	}
    }else{
	die "Could not find log dir $ENV{LOG_DIR}";
    }
    return $LAPS_LOG_PATH;	
}

1;
