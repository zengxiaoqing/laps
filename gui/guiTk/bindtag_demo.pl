 use Tk;

my $MW = tkinit;

my $Display = $MW->Text(-wrap => "word")->pack();
my @btags = $Display->bindtags;
$Display->bindtags( [ @btags[1,0,2,3] ] );
$Display->bind("<KeyPress>", \&processKey);

MainLoop;

sub processKey {
 shift->markSet('insert','end');
} 
