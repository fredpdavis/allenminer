=head1 NAME

allenminer.pl - a tool for identification of genes expressed in patterns
                of interest using the Allen Brain Atlas

=head1 AUTHOR

Fred P. Davis (fredpdavis@gmail.com)

=head1 LICENCE AND COPYRIGHT

Copyright 2008,2012 Fred P. Davis (fredpdavis@gmail.com).
See the file COPYING for copying permission.

This file is part of ALLENMINER.

ALLENMINER is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ALLENMINER is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ALLENMINER.  If not, see <http://www.gnu.org/licenses/>.

=cut

use lib 'perl_api' ;
use alnmnr;

use strict;
use warnings;

main() ;

sub main {

   alnmnr::run({ARGV => \@ARGV}) ;

}
