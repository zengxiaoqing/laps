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
