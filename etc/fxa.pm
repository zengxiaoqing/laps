package fxa;

sub Get_env'fxa{
    local($fxa_env_file) = "/awips/fxa/.environs";
    open(FXA,"$fxa_env_file") || return 0;

    while(<FXA>){
	next if /^#/;
	/\s*([^\s]+)\s+([^\s]+)\s*$/;
        $evar = $1;
        $eval = $2;
        $eval =~ s#\$\{(.*)\}#$ENV{$1}#;
	$eval =~ s#\$([^\s\/]+)#$ENV{$1}#;    
        $ENV{$evar} = $eval;
    }
    close(FXA);
    return 1;
}

1;

sub Set_logdir'fxa{
    
    if( -d $ENV{LOG_DIR}){
	local($yymmdd) = `date -u +%y%m%d`;
	chomp($yymmdd);
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
	
}
1;
