#!/usr/bin/perl -T
# this simply cats files.txt and adds an immediate expiration date

#for security
$ENV{'PATH'}="";

use HTTP::Date;
my $stringGMT = time2str(time());   # Format as GMT ascii time
print "Content-type: text/plain\n".
      "Expires: $stringGMT\n\n";

open FILES, "/bin/cat files.txt|";
while (<FILES>) {
    print;
}
