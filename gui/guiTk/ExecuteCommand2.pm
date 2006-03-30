
# ----------------------------------------------------------
# Modified code originally:
#      Proc::Killfam
#      Tk::ExecuteCommand
# customized as srt_execute_cmd.pl
# ----------------------------------------------------------
$Proc::Killfam::VERSION = '1.0';

package Proc::Killfam;

use Exporter;
use base qw/Exporter/;
use subs qw/get_pids/;
use vars qw/@EXPORT @EXPORT_OK $ppt_OK/;
use strict;

@EXPORT = qw/killfam/;
@EXPORT_OK = qw/killfam/;

# We need Proc::ProcessTable to work properly.  If it's not available,
# then we act like Perl's builtin kill() command.

BEGIN {
    $ppt_OK = 1;
    eval "require Proc::ProcessTable";
    if ($@) {
	$ppt_OK = 0;
#	warn "Proc::ProcessTable missing, can't kill sub-children.";
    }
}

sub killfam {

    my($signal, @pids) = @_;

    if ($ppt_OK) {
	my $pt = Proc::ProcessTable->new;
	my(@procs) =  @{$pt->table};
	my(@kids) = get_pids \@procs, @pids;
	@pids = (@pids, @kids);
    }

    kill $signal, @pids;

} # end killfam

sub get_pids {

    my($procs, @kids) = @_;

    my @pids;
    foreach my $kid (@kids) {
	foreach my $proc (@$procs) {
	    if ($proc->ppid == $kid) {
		my $pid = $proc->pid;
		push @pids, $pid, get_pids $procs, $pid;
	    } 
	}
    }
    @pids;

} # end get_pids

$Tk::ExecuteCommand::VERSION = '1.1';

package Tk::ExecuteCommand;

use IO::Handle;
#require Proc::Killfam;
use Tk::widgets qw/LabEntry ROText/;
use base qw/Tk::Frame/;
use strict;

Construct Tk::Widget 'ExecuteCommand';

sub Populate {

    my($self, $args) = @_;


    $self->SUPER::Populate($args);

#    my $frame0 = $self->Frame()->pack(-anchor => 'w');

#    my $frame2 = $frame0->Frame->pack(-side => 'left', -padx => 5, -anchor => 'n');
#    $frame2->Label(
#        -text => 'Executing Command:', -fg => 'navy')
#        ->pack(-anchor => 'w');


    #my $frame1 = $frame0->Frame()->pack(-side => 'top');
    #my $frame2 = $frame0->Frame()->pack(-side => 'bottom');

#    my $scroll = $frame1->Scrollbar(-orient => "horizontal" );
#    my $jobprocess = $frame1->Entry(
#        -textvariable => \$self->{-command},
#        -exportselection => 1,
#        -xscrollcommand => [ 'set' => $scroll ]);
##        -width => 70,
##        ->pack(-expand => 1, -fill => 'x');
##    $scroll->pack(-expand => 1, -fill => 'x');
#    $scroll->configure(-command => [ $jobprocess => 'xview' ]);
#
#    my $sys_cmd=$jobprocess->cget(-text); 

    my $sys_cmd=\$self->{-command};
    my $text = $self->Scrolled('ROText', -wrap => 'none', -height => 12);
    $text->pack(qw/-side top -expand 1 -fill both/); 

    $self->Advertise('text' => $text);
    $self->OnDestroy([$self => 'kill_command']);

    $self->{-finish} = 0;

    $self->ConfigSpecs(
        -command    => [qw/METHOD command Command/, 'sleep 5; pwd'],
        -display_command    => [qw/METHOD display_command Command/],
        -scrollbars => [$text, qw/scrollbars Scrollbars se/],
    );
  
    if (0) {
      my $stop = $self->Button(-text => 'Stop', -width => 8, -bg => 'cornsilk',
      )->pack(qw/-side left -pady 5 -padx 30 -anchor n/);
      $stop->configure(
        -text    => 'Stop', -width => 8,
        -relief  => 'raised',
        -state   => 'normal',
        -command => [\&kill_command, $self],
      );
      $self->Advertise('stop' => $stop);
    }

    my $doit = $self->Button(-text => 'Localize', -width => 8, 
                             -bg => 'cornsilk')
       ->pack(qw/-pady 5 -padx 5 -anchor e/);

#    $self->Frame->pack(qw/pady 10/);
#    $self->Label(-text => ' Standard output and error (if available):',
#                    -fg => 'navy')->pack(-anchor => 'w');

    #$self->Advertise('jobprocess' => $jobprocess);
    $self->Advertise('doit' => $doit);
    $self->_reset_doit_button;

} # end Populate

sub command {

    my($self, $command) = @_;
    $self->{-command} = $command;

} # end command

sub display_command {

    my($self, $command) = @_;
    my $tx = $self->Subwidget('text');
    $tx->delete('1.0','end'); 
    $tx->insert('end', 'The system command to execute is: '); 
    $tx->insert('end', $command); 
    $self->Subwidget('doit')->configure(-background => '#ffffaa'); # yellow

} # end display_command

sub _flash_doit {

    # Flash "Do It" by alternating its background color.

    my($self, $option, $val1, $val2, $interval) = @_;

    if ($self->{-finish} == 0) {
	$self->Subwidget('doit')->configure($option => $val1);
	$self->idletasks;
	$self->after($interval, [\&_flash_doit, $self, $option, $val2,
            $val1, $interval]);
    }

} # end _flash_doit

sub _read_stdout {

    # Called when input is available for the output window.  Also checks
    # to see if the user has clicked Cancel.

    my($self) = @_;

    my $t;

    if ($self->{-finish}) {
	$self->kill_command;
    } else {
	my $h = $self->{-handle};
	if ( sysread $h, $_, 4096 ) {
	    $t = $self->Subwidget('text');
	    $t->insert('end', $_);
	    $t->yview('end');
	    $self->idletasks;
       	    #$self->after(500, [ sub {$t->insert('end', ".\n") }] );
	} else {
	    $t = $self->Subwidget('text');
            $t->insert('end', "\n\n\t\t\t--- Process has finished --- \n");
	    $self->{-finish} = 1;
	    $t->yview('end');

            # Call srt_localize_domain subroutine.
            &system_tool::localization_done;
	}
    }
	
} # end _read_stdout

sub _reset_doit_button {

    # Establish normal "Do It" button parameters.

    my($self) = @_;

    my $doit = $self->Subwidget('doit');
    my $doit_bg = ($doit->configure(-background))[3];
    $doit->configure(
        -text       => 'Run Localization', -width => 14,
        -relief     => 'raised',
        -background => $doit_bg,
        -state      => 'normal',
        -command    => [sub {
            &system_tool::localization_began;
	    my($self) = @_;
            $self->{-finish} = 0;
            $self->Subwidget('doit')->configure(
                -text   => 'Working ...',
                -relief => 'sunken',
                -state  => 'disabled'
            );
            $self->execute_command;
        }, $self],
    );

} # end _reset_doit_button

# Public methods.

sub execute_command {

    # Execute the command and capture stdout/stderr.

    my($self) = @_;
    
    my $h = IO::Handle->new;
    die "IO::Handle->new failed." unless defined $h;
    $self->{-handle} = $h;

    $self->{-pid} = open $h, $self->{-command} . ' 2>&1 |';
    if (not defined $self->{-pid}) {
	$self->Subwidget('text')->insert('end',
            "'" . $self->{-command} . "' : $!\n");
	$self->kill_command;
	return;
    }
    $h->autoflush(1);
    $self->fileevent($h, 'readable' => [\&_read_stdout, $self]);

    my $doit = $self->Subwidget('doit');
    $doit->configure(
        -text    => 'Running...', -width => 14,
        -relief  => 'flat',
        -state   => 'normal',
        -command => sub { print ""; },
        # Make inactive.
        #-text    => 'Cancel', -width => 8,
        #-command => [\&kill_command, $self],
    );

    my $doit_bg = ($doit->configure(-background))[3];
    $self->_flash_doit(-background => $doit_bg, qw/cornsilk 750/);

    $self->waitVariable(\$self->{-finish});
    
} # end execute_command

sub get_status {

    # Return a 2 element array of $? and $! from last command execution.

    my($self) = @_;

    my $stat = $self->{-status};
    return (defined $stat ? @$stat : undef);

} # end get_status

sub kill_command {
    
    # A click on the blinking Cancel button resumes normal operations.

    my($self) = @_;

    $self->{-finish} = 1;
    my $h = $self->{-handle};
    return unless defined $h;
    $self->fileevent($h, 'readable' => ''); # clear handler
    Proc::Killfam::killfam 'TERM', $self->{-pid} if defined $self->{-pid};
    close $h;
    $self->{-status} = [$?, $!];
    $self->_reset_doit_button;

} # end kill_command


1;

__END__

=head1 NAME

Proc::Killfam - kill a list of pids, and all their sub-children

=head1 SYNOPSIS

 use Proc::Kilfam;
 killfam $signal, @pids;

=head1 DESCRIPTION

B<killfam> accepts the same arguments as the Perl builtin B<kill> command,
but, additionally, recursively searches the process table for children and
kills them as well.

=head1 EXAMPLE

B<killfam 'TERM', ($pid1, $pid2, @more_pids)>;

=head1 KEYWORDS

kill, signal

=cut

1;

__END__

=head1 NAME

Tk::ExecuteCommand - execute a command asynchronously (non-blocking).

=for pm Tk/ExecuteCommand.pm

=for category Widgets

=head1 SYNOPSIS

S<    >I<$exec> = I<$parent>-E<gt>B<ExecuteCommand>;

=head1 DESCRIPTION

Tk::ExecuteCommand runs a command yet still allows Tk events to flow.  All
command output and errors are displayed in a window.

This ExecuteCommand mega widget is composed of an LabEntry widget for
command entry, a "Do It" Button that initiates command execution, and
a ROText widget that collects command execution output.

While the command is executing, the "Do It" Button changes to a "Cancel"
Button that can prematurely kill the executing command. The B<kill_command>
method does the same thing programmatically.

=over 4

=item B<-command>

The command to execute asynchronously.

=back

=head1 METHODS

=over 4

=item C<$exec-E<gt>B<execute_command>;>

Initiates command execution.

=item C<$exec-E<gt>B<get_status>;>

Return a 2 element array of $? and $! from last command execution.

=item C<$exec-E<gt>B<kill_command>;>

Terminates the command.  This subroutine is called automatically via an
OnDestroy handler when the ExecuteCommand widget goes away.

=back

=head1 EXAMPLE

I<$exec> = I<$mw>-E<gt>B<ExecuteCommand>;

=head1 KEYWORDS

exec, command, fork, asynchronous, non-blocking, widget

=head1 COPYRIGHT

Copyright (C) 1999 - 2001 Stephen O. Lidie. All rights reserved.

This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut

