=head1 NAME

alnmnr::getinp.pm - routines to get input for allenminer runs

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

package alnmnr::getinp ;
use strict;
use warnings;
use Cwd;
use File::Path qw/mkpath/ ;
use File::Temp qw/tempfile tempdir/ ;
use POSIX qw/floor ceil/ ;
use alnmnr::io ;
use alnmnr::runs ;


=head2 _getinp_roiop()

   Title:       _getinp_roip()
   Function:    getinput() routine for run mode roiop
   Args:        $_->{options} = any command line options parsed by getinput()
                $_->{options}->{roidef_fn} = ROI definition file.
   Returns:     $_->{roidef} = parsed ROI definition
                  (readin_roi_file() structure)

=cut

sub _getinp_roiop {

   my $in = shift ;
   my $options = $in->{options} ;

   my $type2fields = {
      'peel'                    => { roidef_fn  => 1,
                                     layers     => 1 },

      'flip'                    => { roidef_fn  => 1 },

      'subtract'                => { roidef_fn  => 1,
                                     roi1       => 1,
                                     roi2       => 1,
                                     new_name   => 1 },

      'rename'                  => { roidef_fn  => 1,
                                     new_name   => 1,
                                     old_name   => 1 },

      'separate_adjacent'       => { roidef_fn  => 1,
                                     layers     => 1 },

      'partition'               => { roidef_fn  => 1,
                                     axes       => 1,
                                     numbins    => 1 },
   } ;

   if (!defined $options->{type} ||
       !exists $type2fields->{$options->{type}}) {
      die "must specify -type, one of: ".
          join(", ", sort keys %{$type2fields})."\n"; }

   foreach my $field (keys %{$type2fields->{$options->{type}}}) {
      if (!exists $options->{$field}) {
         die "$field must be specified\n";}}

# Defaults to adult if no age specified.
   if (!defined $options->{age} || $options->{age} eq 'adult') {
      $options->{age} = 'P56' ; }

# load roi definition file
   my $out ;
   $out->{roidef} = alnmnr::io::readin_roi_file({
      fn => $options->{roidef_fn},
      age => $options->{age}
   }) ;

   return $out ;

}



=head2 _getinp_search()

   Title:       _getinp_search()
   Function:    getinput() routine for run mode roi_expr
   Args:        $_->{options} = any command line options parsed by getinput()
                $_->{options}->{roidef_fn} = ROI definition file.
                $_->{options}->{xpr_list_fn} = optional list of XPR files
                  to process; if not specified, parses full database

   Returns:     $_->{roidef} = parsed ROI definition
                  (readin_roi_file() structure)
                $_->{xpr_list} = [file1, 2, ...n] - if xpr_list_fn is specified

=cut

sub _getinp_search {

   my $in = shift ;
   my $options = $in->{options} ;
   my $fields = {
      numjobs => 0,
      cluster_fl => 0,
      roidef_fn => 1,
      xpz_list_fn => 0,
      age => 0,
      lrcheck => 0,
      xyz_details => 0,
      xyz_details_compress => 0,
      xyz_details_out_suffix => 0,
      xyz_details_out_dir => 0,
      xyz_details_out_subdir => 0,
      compress_fl => 0,
      out_fn => 0,
   } ;

#      xpr_list_fn => 0,
#      use_sva => 0,

   my $out ;
   foreach my $field (keys %{$fields}) {
      if ($fields->{$field} == 1 && !exists $options->{$field}) {
         die "$field must be specified\n" ; } }

# Defaults to adult if no age specified.
   if (!defined $options->{age} || $options->{age} eq 'adult') {
      $options->{age} = 'P56' ; }

# load roi definition file
   $out->{roidef} = alnmnr::io::readin_roi_file({
      fn => $options->{roidef_fn},
      age => $options->{age}
   }) ;

## If SVA, convert ROI to 200 micron resolution
#   if (exists $options->{sva_list_fn} || $options->{use_sva}) {
#      my $roidef_200 = convert_roidef_100_to_200({ roidef => $out->{roidef} });
#      $out->{roidef} = {} ;
#      $out->{roidef} = $roidef_200 ;
#      $out->{use_sva} = 1;
#   }

#   if (exists $options->{xpr_list_fn}) {
#      $out->{xpr_list} = readin_1col_arr({fn => $options->{xpr_list_fn}}) ; }

#   if (exists $options->{sva_list_fn}) {
#      $out->{xpr_list} = readin_1col_arr({fn => $options->{sva_list_fn}}) ; }
#

#fpd111122_0838  - note all expression files are basically XPZ format now...
   if (exists $options->{xpz_list_fn}) {
      $out->{xpz_list} = alnmnr::io::readin_1col_arr({
                           fn => $options->{xpz_list_fn}}) ; }

   return $out ;

}


=head2 _getinp_simsearch()

   Title:       _getinp_simsearch()
   Function:    getinput() routine for run mode roi_expr
   Args:        $_->{options} = any command line options parsed by getinput()
                $_->{options}->{roidef_fn} = ROI definition file (optional).
                $_->{options}->{sva_query_fn} = SVA file to use as query
                $_->{options}->{sva_query_list_fn} = list of query SVA files
                  either sva_query_fn or sva_query_list_fn must be specified
                $_->{options}->{sva_target_list_fn} = optional list of SVA files
                  to process; if not specified, parses full database
                $_->{options}->{numjobs} = number of jobs to run on cluster

   Returns:     $_->{roidef} = parsed ROI definition
                  (readin_roi_file() structure)
                $_->{xpr_list} = [file1, 2, ...n] - if xpr_list_fn is specified

=cut

sub _getinp_simsearch {

   my $in = shift ;
   my $options = $in->{options} ;
   my $fields = {
      numjobs => 0,
      cluster_fl => 0,
      roidef_fn => 0,
      xpz_target_list_fn => 0,
      xpz_query_list_fn => 0,
      xpz_query_fn => 0,
      plane => 0,
      calc_pval => 0,
   } ;

   my $out ;
   foreach my $field (keys %{$fields}) {
      if ($fields->{$field} == 1 && !exists $options->{$field}) {
         die "$field must be specified" ; } }

# Defaults to adult if no age specified.
   if (!defined $options->{age} || $options->{age} eq 'adult') {
      $options->{age} = 'P56' ; }

# load roi definition file
   if (exists $options->{roidef_fn}) {
      $out->{roidef} = alnmnr::io::readin_roi_file({
                           fn => $options->{roidef_fn},
                           age => $options->{age}}) ;
   }

# Make sure any XPZ-labeled things are also available as SVA options
#   if (exists $options->{xpz_query_list_fn}) {
#      $options->{sva_query_list_fn} = $options->{xpz_query_list_fn} ; }
#
#   if (exists $options->{xpz_query_fn}) {
#      $options->{sva_query_fn} = $options->{xpz_query_fn} ; }
#
#   if (exists $options->{xpz_target_list_fn}) {
#      $options->{sva_target_list_fn} = $options->{xpz_target_list_fn} ; }

# load xpz query information
   if (exists $options->{xpz_query_list_fn}) {
      $out->{xpz_query_list} = alnmnr::io::readin_1col_arr({
         fn => $options->{xpz_query_list_fn}});
   } elsif ( exists $options->{xpz_query_fn}) {
      $out->{xpz_query_list} = [$options->{xpz_query_fn}] ;
   } else {
      die "ERROR: must specify either xpz_query_fn or xpz_query_list_fn\n";
   }

# load xpz target list if specified
   if (exists $options->{xpz_target_list_fn}) {
      $out->{xpz_target_list} = alnmnr::io::readin_1col_arr({
         fn => $options->{xpz_target_list_fn}}); }

   return $out ;

}


=head2 _getinp_convert()

   Title:       _getinp_convert)
   Function:    getinput() routine for run mode convert
   Args:        $_->{options} = any command line options parsed by getinput()
                $_->{options}->{roidef_fn} = ROI definition file.
   Returns:     $_->{roidef} = parsed ROI definition
                  (readin_roi_file() structure)

=cut

sub _getinp_convert {

   my $in = shift ;
   my $options = $in->{options} ;

   my $type2fields = {
      'roi2pdb'         => { roidef_fn  => 1,
                             pdb_fn     => 1 },

      'roi3d_to_2d'     => { roidef_fn  => 1,
                             xpz_fn     => 1 },


      'xpr2pdb'         => { xpr_fn     => 1,
                             pdb_fn     => 1 },

      'xpz2pdb'         => { xpz_fn     => 1,
                             pdb_fn     => 1 },

      'sva2pdb'         => { sva_fn     => 1,
                             pdb_fn     => 1 },

      'xprarchive2pdb'  => { out_dir    => 1 },

      'xpr2txt'         => { xpr_fn     => 1,
                             txt_fn     => 1 },

      'xpz2txt'         => { xpz_fn     => 1,
                             txt_fn     => 1 },

      'sva2txt'         => { sva_fn     => 1,
                             txt_fn     => 1 },
   } ;

   if (!defined $options->{type} ||
       !exists $type2fields->{$options->{type}}) {
      die "must specify -type, one of: ".
          join(", ", sort keys %{$type2fields})."\n"; }

   foreach my $field (keys %{$type2fields->{$options->{type}}}) {
      if (!exists $options->{$field}) {
         die "$field must be specified\n";}}

# Defaults to adult if no age specified.
   if (!defined $options->{age} || $options->{age} eq 'adult') {
      $options->{age} = 'P56' ; }

# load roi definition file
   my $out ;
   if (exists $options->{roidef_fn}) {
      $out->{roidef} = alnmnr::io::readin_roi_file({
         fn => $options->{roidef_fn},
         age => $options->{age}
      }) ;
   }
   return $out ;

}


sub _getinp_fastsearch {

   my $in = shift ;
   my $options = $in->{options} ;
   my $fields = {
      query_specs_fn => 1,
      aba_results_fn => 0,
      outfile_prefix => 0,
      out_fn => 0,
      age => 0,
   } ;

   my $queries ;

# Defaults to adult if no age specified.
   if (!defined $options->{age} || $options->{age} eq 'adult') {
      $options->{age} = 'P56' ; }

# Call create-data if specified
   if (exists $options->{prep_data}) {
      alnmnr::runs::create_fastsearch_data() ;
      exit;
   }

   my $out ;
   foreach my $field (keys %{$fields}) {
      if ($fields->{$field} == 1 && !exists $options->{$field}) {
         die "$field must be specified\n" ; } }

   return $queries ;

}



=head2 _getinp_make_report()

   Title:       _getinp_make_report()
   Function:    routine to format -search and -simsearch results in tex/html
   Args:        $_->{options} = any command line options parsed by getinput()
                $_->{options}->{roidef_fn} = ROI definition file.
                $_->{options}->{xpr_list_fn} = optional list of XPR files
                  to process; if not specified, parses full database

   Returns:     $_->{roidef} = parsed ROI definition
                  (readin_roi_file() structure)
                $_->{xpr_list} = [file1, 2, ...n] - if xpr_list_fn is specified

=cut

sub _make_report {

   my $in = shift ;
   my $options = $in->{options} ;
   my $fields = {
      roidef_fn => 1,   # name of ROI definition file
      search_results_fn => 1,
      age => 0,         # defaults to P56 (adult)
      out_prefix => 1,  # prefix of output file names
      query_plane => 0, # only show hits to query imageseries in this PLANE (if simsearch) 
      hit_plane => 0,   # only show hits in this PLANE (if simserach)
      plane => 0,       # only show hits in this PLANE (if not simserach)
      num_images => 0,  # images to retrieve; defaults to 10
      sortby_field => 0,     # defaults to enrichment score, similarity score
      sortby_record => 0,    # for results with multiple record types, ie patterned
                             # search, defaults to patterning records
   } ;

   my $out ;
   foreach my $field (keys %{$fields}) {
      if ($fields->{$field} == 1 && !exists $options->{$field}) {
         die "$field must be specified" ; } }

# Defaults to adult if no age specified.
   if (!defined $options->{age} || $options->{age} eq 'adult') {
      $options->{age} = 'P56' ; }

# load roi definition file
   $out->{roidef} = alnmnr::io::readin_roi_file({
      fn => $options->{roidef_fn},
      age => $options->{age}
   }) ;

#fpd111122_0838  - note all expression files are basically XPZ format now...
   return $out ;

}


1;
