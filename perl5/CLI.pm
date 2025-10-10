package CLI;

use strict;
use warnings;
use Data::Dumper;
use Params::Check;

#.......................................
sub new {
  my($class, %args) = @_;
  return bless \%args, $class;
}
#.......................................
sub usage {
  my($self) = @_;
  print "USAGE HERE\n";
}

#.......................................
sub dump {
  my($self) = @_;
  print Dumper($self);
}

1;

