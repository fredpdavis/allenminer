=head1 NAME

alnmnr.pm - routines to mine the Allen Brain Atlas.

=head1 DESCRIPTION

AllenMiner.pm contains routines to parse and analyze expression data from
the adult and developing Allen Brain Atlases and spinal cord atlas.

See docs/allenminer_users_manual.pdf for detailed usage instructions.

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

=head1 SUBROUTINES

=cut

package alnmnr ;
use strict;
use warnings;
use SGE ;
use Exporter ;
our @ISA = qw/Exporter/ ;

use File::Copy ;
use File::Path ;
use File::Temp qw/tempfile tempdir/ ;
use Sys::Hostname qw/hostname/ ;
use POSIX qw/floor ceil/ ;

use alnmnr::spinalcord ;
use alnmnr::download;
use alnmnr::getinp ;
use alnmnr::runs ;

sub getspecs {

   my $specs = {"version" => "2.1"} ;

#----- EDIT THIS SECTION TO SPECIFY YOUR LOCAL CONFIGURATION ------------
#
   $specs->{allen_data_dir} =
      '/gpfs/gsfs4/users/davisfp/github.repos/allenminer/data' ;

   $specs->{cluster}->{cluster_mode} = 0 ; # Use cluster by default? 1=yes, 0=no
   $specs->{cluster}->{head_node} = 'login-eddy' ; # Name of head node
   $specs->{cluster}->{numjobs} = 300 ;              # Number of nodes to use
   $specs->{cluster}->{qstat_sleep} = 60 ;         # Polling frequency (sec)
#
#------------------------------------------------------------------------

   $specs->{allen_spinalcord_fn} = #spinal expression data.
      $specs->{allen_data_dir}.'/allen_spinalcord_data.txt' ;
   $specs->{allen_xpz_dir} = #Developing atlas 3D expression files
      $specs->{allen_data_dir}.'/xpz_files' ;
   $specs->{allen_devel_atlas_dir} = #Developing atlas files
      $specs->{allen_data_dir}.'/devel_atlas' ;
   $specs->{allen_adult_atlas_dir} = #Developing atlas files
      $specs->{allen_data_dir}.'/adult_atlas' ;
   $specs->{allen_fastsearch_dir} = #Developing atlas 3D expression files
      $specs->{allen_data_dir}.'/fastsearch' ;

# {atlas_downsample} is like a pre-emptive 'conversion_factor' (used later
#  to switch betwen tissue-volume and view-volume coordinates) done while
#  parsing the atlas Annotation file; done only for adult, because file is
#  so big (~150MB binary)

   $specs->{atlas_downsample} = {
      "P56" => {from => 25, to => 200} #file is just too big to load to memory.
   } ;

#tissue-volume      "old_adult_xpr"   => {x => 133, y =>  81, z => 115, spacing => 100 },
   $specs->{atlas_dims} = { #tissue_volume
      "P56"   => {x => 528, y => 320, z => 456, spacing =>  25 },
      "E11.5" => {x => 111, y => 119, z =>  64, spacing =>  50 },
      "E13.5" => {x => 177, y => 216, z => 137, spacing =>  50 },
      "E15.5" => {x => 226, y => 315, z => 155, spacing =>  50 },
      "E18.5" => {x => 187, y => 119, z => 112, spacing =>  50 },
      "P4"    => {x => 243, y => 136, z => 160, spacing =>  50 },
      "P14"   => {x => 267, y => 158, z => 196, spacing =>  50 },
      "P28"   => {x => 290, y => 161, z => 210, spacing =>  50 },
   } ;

   $specs->{dims} = { #tissue_volume (after any {atlas_downsample} application)
      "P56"   => {x =>  67, y =>  41, z =>  58, spacing => 200 }, #after downsample
      "E11.5" => {x => 111, y => 119, z =>  64, spacing =>  50 },
      "E13.5" => {x => 177, y => 216, z => 137, spacing =>  50 },
      "E15.5" => {x => 226, y => 315, z => 155, spacing =>  50 },
      "E18.5" => {x => 187, y => 119, z => 112, spacing =>  50 },
      "P4"    => {x => 243, y => 136, z => 160, spacing =>  50 },
      "P14"   => {x => 267, y => 158, z => 196, spacing =>  50 },
      "P28"   => {x => 290, y => 161, z => 210, spacing =>  50 },
   } ;

#view-volume      "old_adult_xpr"   => {x => 133, y =>  81, z => 115, spacing => 100 },
   $specs->{view_volume} = {
      "P56"   => {x =>  67, y =>  41, z =>  58, spacing => 200 },
      "E11.5" => {x =>  70, y =>  75, z =>  40, spacing =>  80 },
      "E13.5" => {x =>  89, y => 109, z =>  69, spacing => 100 },
      "E15.5" => {x =>  94, y => 132, z =>  65, spacing => 120 },
      "E18.5" => {x =>  67, y =>  43, z =>  40, spacing => 140 },
      "P4"    => {x =>  77, y =>  43, z =>  50, spacing => 160 },
      "P14"   => {x =>  68, y =>  40, z =>  50, spacing => 200 },
      "P28"   => {x =>  73, y =>  41, z =>  53, spacing => 200 },
   } ;

#   $specs->{images}->{zoomout2scale} = {
#      1 => 80,
#      2 => 40,
#      3 => 20,
#      4 => 10,
#      5 => 5,
#      6 => 2.5,
#   } ;

   $specs->{report} = {
      num_images          => 10,
      hits_per_page       => 50,
      image_zoom          => 3,
      image_margins_perc  => 0.1,
      sortby       => {
         search =>                { record      => 'E',
                                    field       => 'specificity_energy', },

         search_patterned =>      { record      => 'P',
                                    field       => 'infocontent_nl', },

         simsearch =>             { record      => '',
                                    field       => 'similarity_score', }
      },
   } ; 

# Default split/sort by fields
   $specs->{report}->{rec2specs} = {
      FE  => {search  => 'fastsearch_enrichment',
              sortby  => 'enrichment_energy',
              splitby => ['ROI1', 'ROI2']},

      F   => {search  => 'fastsearch',
              sortby  => 'enrichment_energy',
              splitby => ['ROI']},

      P   => {search  => 'search_patterned',
              sortby  => 'energy_infocontent_nl',
              splitby => ['ROI_basename']},

      S   => {search  => 'search_patterned',
              sortby  => 'energy_bin_enrichment',
              splitby => ['ROI']},

      SIM => {search  => 'simsearch',
              sortby  => 'similarity_score',
              splitby => ['ROI']},

      E   => {search  => 'search',
              sortby  => 'specificity_energy',
              splitby => ['ROI']},
   };
   $specs->{report}->{search2rec} = {
      search                    => 'energy',
      fastsearch                => 'F',
      fastsearch_enrichment     => 'FE',
      simsearch                 => 'SIM',
      search_patterned          => 'P',
   } ;


   $specs->{download} = {
      adult_xpz_search_chunk => 200,
      devel_xpz_search_chunk => 100,
      URL => {
         interactive_imageseries => 'http://mouse.brain-map.org/experiment/show/',
# http://mouse.brain-map.org/experiment/show/75492683

         image_info => 'http://mouse.brain-map.org/aba/api/image/info?path=',
         image_part1 => 
            'http://mouse.brain-map.org/aba/api/image?zoom=',
         image_part2 => '&path=',

# http://mouse.brain-map.org/experiment/siv?imageId=71661837&imageType=ish,expression&showSubImage=y&coordSystem=pixel&x=4584&y=4126&z=50
         interactive_image_part1 =>
            'http://mouse.brain-map.org/experiment/siv?imageId=',
         interactive_image_part2 =>
            '&imageType=ish,expression&showSubImage=y&coordSystem=pixel&',
#            x=4584&y=4126&z=50

         adult_search_part1 =>
            'http://mouse.brain-map.org/grid_data/v1/search/gene?term=',
         adult_search_part2 => '&startRow=',
         adult_search_part3 => '&numRows=',

         devel_search_part1 => 'http://developingmouse.brain-map.org/api/v2/data/Gene/query.xml?query=',
         devel_search_part2 => '&start_row=',
         devel_search_part3 => '&num_rows=',
         devel_search_part4 => '&criteria=rma::include,data_sets[delegate$eq%27true%27][failed$eq%27false%27](products[id$eq3],specimen(donor(age[id$in1,2,4,5,7,11,14]))),rma::options[order$eq%27genes.acronym%27]',

         adult_xpz => "http://mouse.brain-map.org/grid_data/v1/visualize/",
         devel_xpz =>
            "http://developingmouse.brain-map.org/grid_data/v1/visualize/",
         devel_atlas => "http://www.brain-map.org/BrainExplorer2/".
                     "Atlases/Developing_Mouse_6.zip",
         adult_atlas => "http://www.brain-map.org/BrainExplorer2/".
                     "Atlases/Mouse_Brain_7.zip",
         spinalcord_part1 =>
            'http://mousespinal.brain-map.org/imageseries/list/',
         spinalcord_part2 => 
            '.html?clear=true&contains=false&'.
            'excludeFail=on&gene_term=&searchSym=t',
         fastsearch_prefix => 
            "https://zenodo.org/record/29427/files/v".
            $specs->{version}."_fastsearch_",
      }
   } ;
   $specs->{download}->{fn}->{adult_xpz_list} =
      $specs->{allen_xpz_dir}."/adult_xpz_list.txt" ;
   $specs->{download}->{fn}->{devel_xpz_list} =
      $specs->{allen_xpz_dir}."/devel_xpz_list.txt" ;

   foreach my $age (keys %{$specs->{atlas_dims}}) {
      $specs->{fastsearch_fn}->{$age} =
         $specs->{allen_fastsearch_dir}."/fastsearch_$age.txt.gz"; }

   $specs->{getinp_routines} = {
      search            =>      \&alnmnr::getinp::_getinp_search,
      fastsearch        =>      \&alnmnr::getinp::_getinp_fastsearch,
      simsearch         =>      \&alnmnr::getinp::_getinp_simsearch,
      roiop             =>      \&alnmnr::getinp::_getinp_roiop,
      convert           =>      \&alnmnr::getinp::_getinp_convert,
      make_report       =>      \&alnmnr::getinp::_make_report,
   } ;

   $specs->{run_routines} = {
      search            =>      \&alnmnr::runs::_run_search,
      fastsearch        =>      \&alnmnr::runs::_run_fastsearch,
      simsearch         =>      \&alnmnr::runs::_run_simsearch,
      roiop             =>      \&alnmnr::runs::_run_roiop,
      convert           =>      \&alnmnr::runs::_run_convert,
      spinesearch       =>      \&alnmnr::spinalcord::_run_spinesearch,
      download          =>      \&alnmnr::download::_run_download_data,
      make_report       =>      \&alnmnr::runs::_make_report,
   } ;

   return $specs ;

}


=head2 run()

   Title:       run()
   Function:    Reads input and runs the proper AllenMiner routine.
   Args:        $_ = hash of raw flags (eg ARGV command line)
   Returns:     DBI database handle to pibaes

=cut

sub run {

   my $raw_in = shift ;
   my $specs = getspecs() ;

# get input
   my $in = getinp({
      specs => $specs,
      ARGV => $raw_in->{ARGV}
   }) ;

   foreach my $raw_key (keys %{$raw_in}) {
      $in->{$raw_key} = $raw_in->{$raw_key} ; }

   $in->{specs} = $specs ;

# run the routine for this mode
   $specs->{run_routines}->{$in->{mode}}->($in) ;

}


=head2 getinp()

   Title:       getinp()
   Function:    General get input routine that calls a specific getinp routine
                  appropriate for the current AllenMiner run mode
   Args:        $_->{ARGV} = [ARGV[0], [1],...]
                  o arrayref of command line options
                $_->{specs} = allenminer specs [optional]-
                  o if not specified will call getspecs()]

   Returns:     hash of input parameters that includes both the ARGV options
                  and anything read in by the mode-specific getin routine
                  ($specs->{getin_routines}->{run_mode}):

                $_->{key} = value for input parameter key

=cut

sub getinp {

   my $raw_in = shift ;
   my $specs ;
   if (!exists $raw_in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $raw_in->{specs} ;
   }

   my $j = 0 ;
   my $options ;
   while ($j <= $#{$raw_in->{ARGV}}) {
      my $type = $raw_in->{ARGV}->[$j] ; $type =~ s/^\-// ;
      my $value = $raw_in->{ARGV}->[($j + 1)] ;
      $options->{$type} = $value ;
      $j+=2 ;
   }

   if (!defined $options->{mode}) {
      die "ERROR: Run mode not specified (-mode XXX)\n   specify one of: ".
         join(', ', sort keys %{$specs->{run_routines}})."\n";}

   my $in ;
   if (exists $specs->{getinp_routines}->{$options->{mode}}) {
      $in = $specs->{getinp_routines}->{$options->{mode}}->({
         options => $options}) ;
   }

   foreach my $key (keys %{$options}) {
      $in->{$key} = $options->{$key} ; }

   my @keys = keys %{$in} ;

   return $in ;

}


1;
