=head1 NAME

alnmnr::core.pm - core ALLENMINER routines (non-IO/-ROI)

=head1 AUTHOR

Fred P. Davis (fredpdavis@gmail.com)

=head1 LICENCE AND COPYRIGHT

Copyright 2008,2011 Fred P. Davis (fredpdavis@gmail.com).
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

package alnmnr::core;
use strict;
use warnings;
use Cwd;
use File::Path qw/mkpath/ ;
use File::Temp qw/tempfile tempdir/ ;
use POSIX qw/floor ceil/ ;



=head2 calc_expr()

   Title:       calc_expr()
   Function:    Calculates the expression level in a region of interest
                  in an ABA 3D expression (XPZ) file - roi defined by either
                  (1) cuboid x/y/zmin/max bounds on all axes
                  (2) ABA brain structure id
                  (3) voxel x,y,z list

   Args:        $_->{xpz}        - xpz file [or on command line]
                $_->{roidef}    - region of interest definitions
                $_->{age}        - developing stages (E11.5,...) or 'adult'
                $_->{display_fl} - flag to print results to STDOUT
                $_->{out_fn}     - file for output
                $_->{out_fh}     - file handle for output

                $_->{xyz_details_out_fn} - file name for XYZ details
                $_->{xyz_details_compress_fl} - flag to compress XYZ details

   Files in:    XPZ file - per ABA format ($_->{xpr} or command line -xpz)
                ROI definition (from readin_roid_file())

   Returns:     $_->{roi_name}->{num_points} number of defined points in ROI
                $_->{roi_name}->{total_expression_level}
                $_->{roi_name}->{mean_expression_level}
                $_->{roi_name}->{total_num_expressors}
                $_->{roi_name}->{mean_num_expressors}

=cut

sub calc_expr {

# changed roidef from old full info list to the simpler:
#  point2roi->{xyz}->{roiname}->{roi_subname}++ ;
#  roi_meta->{roiname}->{roi_subnmae}->{roi_type} = boxbounds|brainstrx|point;
#  roi_meta->{roiname}->{roi_subnmae}->{roi_spec} = line from roidef file;

   my $in = shift ;

   my $out_fh ;
   if (exists $in->{out_fh}) {
      $out_fh = $in->{out_fh} ;
   } elsif (exists $in->{out_fn}) {
      open($out_fh, ">".$in->{out_fn}) ;
   } else {
      open($out_fh, ">-") ;
   }

   my $raw_measures = ['intensity', 'density', 'energy'] ;

   my $details_out_fh ;
   my @details_headers = qw/roi x y z/;
   push @details_headers, @{$raw_measures} ;
   if (exists $in->{xyz_details_out_fn}) {
      my $dirbase = dirname($in->{xyz_details_out_fn}) ;
      if (!-s $dirbase) { mkpath($dirbase); }
      open($details_out_fh, ">".$in->{xyz_details_out_fn}) ;
      print {$details_out_fh} "#".join("\t", @details_headers)."\n" ;
   }

# 3 possible types of output records; last 2 if -patterning set
#  E. expression levels, whole brain enrichment/specificity
#  S. enrichment/specificity for partition vs a whole structure.
#  P. patterning and gradient score within each partitioned ROI.

   my $headers = {};
   $headers->{E} = [] ;

#120208_1510 changed from roi_name -> ROI. for make_report() consistency
   @{$headers->{E}} = ( qw/E gene_name plane xpz_name ROI/,
                        qw/num_roidef_points numpoints/) ;
   map {push @{$headers->{E}}, $_, "nl_".$_; } @{$raw_measures} ;
   map {push @{$headers->{E}}, "specificity_".$_, "specificity_nl_".$_; }
      @{$raw_measures} ;
   push @{$headers->{E}}, qw/specificity_numpoints/ ;

   if (exists $in->{patterning} && $in->{patterning} == 1) {
#120208_1510 changed from roi_name -> ROI
      @{$headers->{S}} = (qw/S gene_name plane xpz_name ROI/) ;

      @{$headers->{P}} = (qw/P gene_name plane xpz_name ROI_basename/,
                          qw/total_roi_energy/);
      foreach my $measure (qw/numpoints density intensity energy/) {
         push @{$headers->{P}}, (
            $measure."_infocontent",
            $measure."_infocontent_nl",
            $measure."_gradient_num",
            $measure."_gradient_den",
            $measure."_gradient_score") ;

         push @{$headers->{S}}, (
            $measure."_bin_specificity",
            $measure."_bin_enrichment",
            $measure."_bin_value",
            $measure."_total_value") ;
      }
   }

#just send back headers - for output file procesing
   if (exists $in->{return_headers_only}) {
      return $headers ; }

   my $display_fl = 0 ;
   if (exists $in->{display_fl} && $in->{display_fl} == 1) {
      $display_fl = 1 ; }

# do point enumeration before hand - only pass in points roi;
# preferably all in one hash to speed things up.

# 1. read in ROI definition: cuboid definition of {x|y|z}{min|max}
   my $roidef = $in->{roidef};

   my $xpz_data = $in->{xpz} ;

   my $roi_totals ;
   foreach my $roi_name (keys %{$roidef->{roi_meta}}) {
      map { $roi_totals->{$roi_name}->{$_} = 0; } @{$raw_measures} ;

      $roi_totals->{$roi_name}->{numpoints} = 0 ;
      $roi_totals->{$roi_name}->{num_roidef_points} =
         $roidef->{num_roidef_points}->{$roi_name} ;
   }
   

#Improvements: smarter data structure for holding these points?; 
#  to enable faster point in roi check....
# or otherway around: roi -> point list...
   map {$roi_totals->{internal_background}->{$_} = 0;}
      (qw/numpoints num_expressors/) ;

   my $roi_points ;
   foreach my $point (@{$xpz_data->{xyz}}) {
      my $xyz = $point->{x}.','.$point->{y}.','.$point->{z} ;
      $roi_totals->{internal_background}->{numpoints}++ ;

      map { $roi_totals->{internal_background}->{$_} += $point->{$_}; }
            @{$raw_measures} ;

      if (!exists $roidef->{point2roi}->{$xyz}) {next;}
      foreach my $roi_name_subname (keys %{$roidef->{point2roi}->{$xyz}}) {
         my ($roi_name, $roi_subname) = split(/\t/, $roi_name_subname) ;
         push @{$roi_points->{$roi_name}}, $point ;
         $roi_totals->{$roi_name}->{numpoints}++ ;

         map { $roi_totals->{$roi_name}->{$_} += $point->{$_}; }
               @{$raw_measures} ;

         if (defined $details_out_fh) {
            my @details_outvals = () ;
            foreach my $details_f (@details_headers) {
               if ($details_f eq 'roi') {
                  push @details_outvals, $roi_name ;
               } else {
                  push @details_outvals, $point->{$details_f} ;
               }
            }
            print {$details_out_fh} join("\t", @details_outvals)."\n" ;
         }
      }
   }

   if (exists $in->{display_headers} && $in->{display_headers} == 1) {
      foreach my $type (keys %{$headers}) {
         print $out_fh '#'.join("\t", @{$headers->{$type}})."\n" ; } }


   foreach my $roi_name ( 'internal_background',
                          (sort keys %{$roidef->{roi_meta}})) {

      if ($roi_name eq 'internal_background') {
         $roi_totals->{$roi_name}->{num_roidef_points} =
            $roi_totals->{$roi_name}->{numpoints} ; }

      my $outvals = {
         E                     => 'E',
         gene_name             => $xpz_data->{meta}->{gene_name},
         xpz_name              => $xpz_data->{meta}->{xpz_name},
         plane                 => $xpz_data->{meta}->{plane},
         ROI                   => $roi_name,
         num_roidef_points     => $roi_totals->{$roi_name}->{num_roidef_points},
         numpoints             => $roi_totals->{$roi_name}->{numpoints},
      } ;

      foreach my $type (@{$raw_measures}) {

         if ($roi_totals->{$roi_name}->{numpoints} == 0) {

            $roi_totals->{$roi_name}->{"nl_".$type} = 'NA' ;
            $roi_totals->{$roi_name}->{$type}       = 'NA' ;

         } else {

            $roi_totals->{$roi_name}->{"nl_".$type} = sprintf("%.3f",
               $roi_totals->{$roi_name}->{$type} /
               $roi_totals->{$roi_name}->{numpoints}) ;

            $roi_totals->{$roi_name}->{$type} = sprintf("%.3f",
               $roi_totals->{$roi_name}->{$type}) ;

         }

         $outvals->{$type}       = $roi_totals->{$roi_name}->{$type} ;
         $outvals->{'nl_'.$type} = $roi_totals->{$roi_name}->{"nl_".$type} ;

      }

      my @tf ;
      map {push @tf, $_;
           push @tf, "nl_".$_;} @{$raw_measures} ;
      push @tf, 'numpoints' ;

      foreach my $type (@tf) {
         if ($roi_totals->{$roi_name}->{$type} ne 'NA' &&
             $roi_totals->{'internal_background'}->{$type} ne 'NA' &&
             $roi_totals->{'internal_background'}->{$type} > 0) {

            $roi_totals->{$roi_name}->{"specificity_$type"} = 
               sprintf("%.3f",
               $roi_totals->{$roi_name}->{$type} / 
               $roi_totals->{'internal_background'}->{$type}) ;

         } else {
            $roi_totals->{$roi_name}->{"specificity_$type"} = 'NA' ;
         }

         $outvals->{"specificity_$type"} = 
            $roi_totals->{$roi_name}->{"specificity_$type"} ;
      }

      print {$out_fh} join("\t", @{$outvals}{@{$headers->{E}}})."\n"
         if $display_fl ;
   }

   if (exists $in->{patterning} && $in->{patterning} == 1 && $display_fl) {
      my $pattern = calc_patterning({ roi_totals     => $roi_totals,
                                      raw_measures   => $raw_measures });

      foreach my $roibase (keys %{$pattern}) {

# CHANGE so all 'data_types' shown on same line

         my $outvals = {} ;

# Print P records (patterning/gradient scores per base ROI)

         $outvals->{P} = {
            P                  => 'P',
            gene_name          => $xpz_data->{meta}->{gene_name},
            xpz_name           => $xpz_data->{meta}->{xpz_name},
            plane              => $xpz_data->{meta}->{plane},
            ROI_basename       => $roibase,
            total_roi_energy   => $pattern->{$roibase}->{total}->{energy},
         } ;

         foreach my $data_type (sort keys %{$pattern->{$roibase}->{bindata}->[0]}){
            foreach my $measure (qw/infocontent infocontent_nl gradient_num/ ,
                                 qw/gradient_den gradient_score/) {
               $outvals->{P}->{$data_type."_".$measure} =
                  $pattern->{$roibase}->{summary}->{$data_type}->{$measure}; }}

         print {$out_fh} join("\t", @{$outvals->{P}}{@{$headers->{P}}})."\n" ;

# Print S records (specificity per bin)
# bindata, summary, total.
         foreach my $i (0 .. $#{$pattern->{$roibase}->{bindata}}) {

            my $roiname = $roibase."_part$i" ;

            $outvals->{S} = {
               S                  => 'S',
               gene_name          => $xpz_data->{meta}->{gene_name},
               xpz_name           => $xpz_data->{meta}->{xpz_name},
               plane              => $xpz_data->{meta}->{plane},
               ROI                => $roiname,
            } ;

            foreach my $data_type (sort keys %{$pattern->{$roibase}->{
                                               "bindata"}->[0]}){
               $outvals->{S}->{$data_type."_bin_specificity"} =
                  $pattern->{$roibase}->{bindata}->[$i]->{
                        $data_type}->{specificity};

               $outvals->{S}->{$data_type."_bin_enrichment"} =
                  $pattern->{$roibase}->{bindata}->[$i]->{
                        $data_type}->{enrichment};

               $outvals->{S}->{$data_type."_bin_value"} =
                  $pattern->{$roibase}->{bindata}->[$i]->{
                        $data_type}->{val};

               $outvals->{S}->{$data_type."_total_value"} =
                  $pattern->{$roibase}->{total}->{$data_type} ;
            }

            print {$out_fh} join("\t", @{$outvals->{S}}{@{$headers->{S}}})."\n";
         }
      }
   }
      
#fpd111221_1209 - no one seems to use $roi_points...
#   if (exists $in->{return_roi_points_fl})      { return $roi_points;}
   if (exists $in->{xyz_details_out_fn})        { close($details_out_fh);}
   if (exists $in->{xyz_details_compress_fl}) {
      system("gzip ".$in->{xyz_details_out_fn}) ; }

# no one uses this hash...
#   return $roi_totals ;
   return 1 ;

}



sub calc_expr_sim {
# purpose: given two XPZ files and optionally ROI definitions, calc correlation
#
# 0. parse XPZ files for voxels registered in ROI
# 1. foreach ROI, calculate pearson's correlation a la NeuroBLAST

   my $in = shift ;
   my $xpz_fn = [$in->{xpz1_fn}, $in->{xpz2_fn}] ;
   my $roidef = $in->{roidef} ;
   my $age = $in->{age} ;

   my $pval_numshuffle = 0 ;
   if (exists $in->{pval_numshuffle}) {
      $pval_numshuffle = $in->{pval_numshuffle} ; }

   my $roidef_mode = '100'; #by default expect roi definitions in 100um mode
   if (exists $in->{roidef_mode}) {
      $roidef_mode = $in->{roidef_mode} ; }

   my $roidef_200 = $in->{roidef};

# fpd111122_0853  - don't think new XPZ format requires ROI down-sampling
#   if (defined $roidef) {
#      if ($roidef_mode == 100 && $age eq 'adult') {
#         $roidef_200 = convert_roidef_100_to_200({
#            roidef => $roidef }) ;
#      } else {
#         $roidef_200 = $in->{roidef}
#      }
#   }

# Parse XPZ files if needed:
   my ($xpz1_data, $xpz2_data) ;
   if (exists $in->{xpz1_data}) {
      $xpz1_data = $in->{xpz1_data} ;
   } else {
      $xpz1_data = alnmnr::io::aba_parse_xpz({
         fn             => $in->{xpz1_fn},
         roidef         => $roidef_200,
      }) ;
   }

   if (exists $in->{xpz2_data}) {
      $xpz2_data = $in->{xpz2_data} ;
   } else {
      $xpz2_data = alnmnr::io::aba_parse_xpz({
         fn             => $in->{xpz2_fn},
         roidef         => $roidef_200,
      }) ;
   }

# if ROIdef is defined, calc expression similarity over each ROI,
#  otherwise, over all shared voxels.

   my $roi_similarity ;
   if (defined $roidef_200) {

      my $roi_energies = {};
      my $num_shared_points = {} ;
      foreach my $point (sort keys %{$roidef_200->{point2roi}}) {

         foreach my $roi_name (keys %{$roidef_200->{point2roi}->{$point}}) {
            if (!exists $num_shared_points->{$roi_name}) {
               $num_shared_points->{$roi_name} = 0 ; } }

# make sure point exists in both SVA files
         if (!exists $xpz1_data->{$point} ||
             !exists $xpz2_data->{$point}) {next;}

         map {push @{$roi_energies->{$_}->[0]}, $xpz1_data->{$point} ;
              push @{$roi_energies->{$_}->[1]}, $xpz2_data->{$point} ;
              $num_shared_points->{$_}++ ; }
          (keys %{$roidef_200->{point2roi}->{$point}}) ;
      }


# 2. Calculate voxel-level pearsons correlation (if any shared points)
      foreach my $roi (keys %{$roi_energies}) {
         if ($num_shared_points->{$roi} < 2) {
            $roi_similarity->{$roi} = {error_fl => 'less than 2 shared points'};
         } else {
            $roi_similarity->{$roi} = calc_correlation_coeff({
               mode     => 'pearson',
               a        => $roi_energies->{$roi}->[0],
               b        => $roi_energies->{$roi}->[1],
               pval_numshuffle => $pval_numshuffle,
            }) ;
         }
      }

   } else { #otherwise over all voxels shared between both SVA files

      my $energies = [] ;
      my $num_shared_points = 0 ;
      foreach my $point (sort keys %{$xpz1_data}) {
         if (!exists $xpz2_data->{$point}) {next;}
         $num_shared_points++ ;
         push @{$energies->[0]}, $xpz1_data->{$point} ;
         push @{$energies->[1]}, $xpz2_data->{$point} ;
      }

      if ($num_shared_points < 2) {
         $roi_similarity->{ALLVOXELS} = {error_fl=>'less than 2 shared points'};
      } else {
         $roi_similarity->{ALLVOXELS} = calc_correlation_coeff({
            mode            => 'pearson',
            a               => $energies->[0],
            b               => $energies->[1],
            pval_numshuffle => $pval_numshuffle,
         }) ;
      }
   }

   return $roi_similarity;
}


sub DEPRECATED_postproc_results_xpr_roi_entropy {
#fpd111222_1250 rolled into -mode search -patterning 1

   my $in = shift ;

   my $out_fh ;
   if (exists $in->{out_fn}) {
      open($out_fh, '>'.$in->{out_fn}) ;
   } else {
      open($out_fh, '>-') ;
   }

   my $f2i ;
   foreach my $j ( 0 .. $#{$in->{headers}->{headers}}) {
      $f2i->{$in->{headers}->{headers}->[$j]} = $j ; }

   print {$out_fh} '#'.join("\t", qw/xpr_name gene roi entry_type_summary total_roi_expression_level data_type infocontent infocontent_nl gradient_num gradient_den gradient_score/)."\n";
   print {$out_fh} '#'.join("\t", qw/xpr_name gene roi entry_type_bindata bin_number data_type bin_specificity bin_enrichment bin_value total_value/)."\n" ;

   my $xpr2roi_data ;
   my $xpr_metainfo ;
   my $last_xpr_roi ; my $last_roibase ; my $last_xprname ;
   my $cur_xprroi_info ;
   open(RESULTSF, $in->{results_fn}) || die "results_fn not specified";
   while (my $line = <RESULTSF>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      my $entry = {};
      my @t = split(/\t/, $line) ;
      foreach my $j (0 .. $#t) {
         $entry->{$in->{headers}->{headers}->[$j]} = $t[$j] ; }

      my $roi = $entry->{roi_name} ;
      if ($roi !~ /part/) {next;}
      my ($roibase, $partno) = ($roi =~ /(.+)_part([0-9]+)/) ;

      if ($entry->{plane} ne 'sagittal' &&
          $entry->{plane} ne 'coronal') { next; }
      my $slice = $entry->{plane} ;
      my $xpr_name = $entry->{xpr_name} ;

      my $lastline_fl = 0 ;
      if (eof(RESULTSF)) {
         foreach my $f (sort keys %{$f2i}) {
            if ($entry->{$f} eq 'NA') {$entry->{$f} = 0;}
            $cur_xprroi_info->{$roibase}->{bins}->[$partno]->{$f} =
               $entry->{$f} ;
            if ($f eq 'expression_level' || $f eq 'numpoints' ||
                $f eq 'num_roidef_points') {
               $cur_xprroi_info->{$roibase}->{total}->{$f} += $entry->{$f} ; }
         }
         $lastline_fl = 1 ;
      }

      if (defined $last_xpr_roi && $last_xpr_roi ne $xpr_name."\t".$roibase ||
          $lastline_fl == 1) { #LAST XPR ; do IC, IC_NL calcs

         foreach my $data_type (qw/numpoints expression_level/) {
            my $infocontent = { raw => 0,
                                nl => 0 } ;

            my $gradient_score = {num => 0, den=>0, score=> 0};

            my ($bin2p, $bin2p_ref, $bin2enrichment) = ([],[],[]);
            if ($cur_xprroi_info->{$last_roibase}->{total}->{expression_level}
               > 0) {
            foreach my $part_no (0 ..
               $#{$cur_xprroi_info->{$last_roibase}->{bins}}) {

               if ($cur_xprroi_info->{$last_roibase}->{bins}->[$part_no]->{
                     $data_type} == 0) {
                  $bin2p->[$part_no]            = 0 ;
                  $bin2p_ref->[$part_no]        = 0 ;
                  $bin2enrichment->[$part_no]   = 0 ;
                  next;
               }

   my $p = $cur_xprroi_info->{$last_roibase}->{bins}->[$part_no]->{$data_type} /
      $cur_xprroi_info->{$last_roibase}->{total}->{$data_type} ;
               $bin2p->[$part_no] = $p ;

   my $p_ref =
    $cur_xprroi_info->{$last_roibase}->{bins}->[$part_no]->{num_roidef_points} /
      $cur_xprroi_info->{$last_roibase}->{total}->{num_roidef_points} ;
               $bin2p_ref->[$part_no] = $p_ref ;

               $bin2enrichment->[$part_no] = $bin2p->[$part_no] /
                                             $bin2p_ref->[$part_no] ;

               $infocontent->{raw} += $p * (log($p) / log(2)) ;
               $infocontent->{nl} += $p * (log($p / $p_ref) / log(2)) ;
            }

            foreach my $part_no (0 .. $#{$bin2p}) {
               print {$out_fh} join("\t", $last_xprname,
                  $cur_xprroi_info->{$last_roibase}->{bins}->[0]->{gene_name},
                  $roibase,
                  "BINDATA",
                  $part_no,
                  $data_type,
                  $bin2p->[$part_no],
                  $bin2enrichment->[$part_no],
       $cur_xprroi_info->{$last_roibase}->{bins}->[$part_no]->{$data_type},
                  $cur_xprroi_info->{$last_roibase}->{total}->{$data_type},
               )."\n" ;
            }

            foreach my $part_no (1 .. $#{$bin2p_ref}) {
               $gradient_score->{num} += ($bin2enrichment->[$part_no] -
                                          $bin2enrichment->[($part_no - 1)]) ;
               $gradient_score->{den} += abs($bin2enrichment->[$part_no] -
                                             $bin2enrichment->[($part_no - 1)]);
            }
            if ($gradient_score->{den} > 0) {
            $gradient_score->{score} = $gradient_score->{num} /
                                       $gradient_score->{den} ;
            } else {
            $gradient_score->{score} = 'UNDEF' ;
            }
            }
            $infocontent->{raw} = -1 * $infocontent->{raw} ;


            print {$out_fh} join("\t", $last_xprname,
               $cur_xprroi_info->{$last_roibase}->{bins}->[0]->{gene_name},
               $roibase,
               "SUMMARY",
               $cur_xprroi_info->{$last_roibase}->{total}->{expression_level},
               $data_type,
               $infocontent->{raw},
               $infocontent->{nl},
               $gradient_score->{num},
               $gradient_score->{den},
               $gradient_score->{score},
            )."\n" ;
         }

         $cur_xprroi_info = {} ; #reset all info
      }

      foreach my $f (sort keys %{$f2i}) {
         if ($entry->{$f} eq 'NA') {$entry->{$f} = 0;}
         $cur_xprroi_info->{$roibase}->{bins}->[$partno]->{$f} =
            $entry->{$f} ;
         if ($f eq 'expression_level' || $f eq 'numpoints' ||
             $f eq 'num_roidef_points') {
            $cur_xprroi_info->{$roibase}->{total}->{$f} += $entry->{$f} ; }
      }
      $last_xprname = $xpr_name ;
      $last_xpr_roi = $xpr_name."\t".$roibase ;
      $last_roibase = $roibase ;
   }

}


# input: calc_expr() values for a single XPZ file
# output: gradient scores for any ROI with defined partitions.
sub calc_patterning {

   my $in               = shift ;
   my $roi_totals       = $in->{roi_totals} ;
   my $raw_measures     = $in->{raw_measures} ;

   my $roi2numbins = {};
   foreach my $roi (keys %{$roi_totals}) {
      if ($roi !~ /_part[0-9]/) {next;}
      my ($roibase, $partno) = ($roi =~ /(.+)_part([0-9]+)/) ;
      if (!exists $roi2numbins->{$roibase} ||
          ($partno + 1) > $roi2numbins->{$roibase}) {
         $roi2numbins->{$roibase} = $partno + 1 ;
      }
   }

   my $out = {};
   foreach my $roibase (keys %{$roi2numbins}) {
#      print STDERR "ROIBASE $roibase, number of bins: ".$roi2numbins->{$roibase}."\n";

      my $cur_roi_info = { bins         => [],
                           total        => {} };

      foreach my $data_type ('num_roidef_points', 'numpoints',
                             @{$raw_measures}) {
         foreach my $i (0 .. $roi2numbins->{$roibase} - 1) {
            my $roi = $roibase."_part$i" ;
            my $curval = 0 ;
            if ($roi_totals->{$roi}->{$data_type} ne 'NA') {
               $curval = $roi_totals->{$roi}->{$data_type}; }

            $cur_roi_info->{bins}->[$i]->{$data_type}  = $curval ;
            $cur_roi_info->{total}->{$data_type}      += $curval ;
         }
         $out->{$roibase}->{total}->{$data_type} =
            sprintf("%.3f",$cur_roi_info->{total}->{$data_type});
      }

      foreach my $data_type ('numpoints', @{$raw_measures}) {
         my $infocontent        = { raw         => 0,
                                    nl          => 0 };
         my $gradient_score     = { num         => 0,
                                    den         => 0,
                                    score       => 0 };

         my ($bin2p, $bin2p_ref, $bin2enrichment) = ([],[],[]);

         if ($cur_roi_info->{total}->{energy} > 0) {
         foreach my $i (0 .. $roi2numbins->{$roibase} - 1) {
            if ($cur_roi_info->{bins}->[$i]->{$data_type} == 0) {
               $bin2p->[$i]               = 0 ;
               $bin2p_ref->[$i]           = 0 ;
               $bin2enrichment->[$i]      = 0 ;
            } else {
               my $p =     $cur_roi_info->{bins}->[$i]->{$data_type} /
                           $cur_roi_info->{total}->{$data_type} ;

               my $p_ref = $cur_roi_info->{bins}->[$i]->{num_roidef_points} /
                           $cur_roi_info->{total}->{num_roidef_points} ;

               $bin2p->[$i]          = $p ;
               $bin2p_ref->[$i]      = $p_ref ;
               $bin2enrichment->[$i] = $p / $p_ref ;

               $infocontent->{raw} += $p * (log($p) / log(2)) ;
               $infocontent->{nl}  += $p * (log($p / $p_ref) / log(2)) ;
            }

            $out->{$roibase}->{bindata}->[$i]->{$data_type} = {
               specificity   => sprintf("%.3f",$bin2p->[$i]),
               enrichment    => sprintf("%.3f",$bin2enrichment->[$i]),
               val           => sprintf("%.3f",
                                 $cur_roi_info->{bins}->[$i]->{$data_type}),
            };
         }

         foreach my $i (1 .. $#{$bin2p_ref}) {
            $gradient_score->{num} += ($bin2enrichment->[$i] -
                                       $bin2enrichment->[($i - 1)]) ;

            $gradient_score->{den} += abs($bin2enrichment->[$i] -
                                          $bin2enrichment->[($i - 1)]);
         }

         if ($gradient_score->{den} > 0) {
            $gradient_score->{score} = sprintf("%.3f",$gradient_score->{num} /
                                                      $gradient_score->{den}  );
         } else {
            $gradient_score->{score} = 'UNDEF' ;
         }
         }
         $infocontent->{raw} = -1 * $infocontent->{raw} ;

         $out->{$roibase}->{summary}->{$data_type} = {
            infocontent         => sprintf("%.3f",$infocontent->{raw}),
            infocontent_nl      => sprintf("%.3f",$infocontent->{nl}),
            gradient_num        => sprintf("%.3f",$gradient_score->{num}),
            gradient_den        => sprintf("%.3f",$gradient_score->{den}),
            gradient_score      => $gradient_score->{score},
         } ;
      }
   }

   return $out ;
}


=head2 val_is_in_arr()

   Title:       val_is_in_arr()
   Function:    Checks an array for presence of a value
   Args:        $_->{array} arrayref
                $_->{value} value to check for

   Returns:     1 if present
                0 if absent

=cut

sub val_is_in_arr {

   my $in = shift ;
   foreach my $j ( 0 .. $#{$in->{array}}) {
      if ($in->{array}->[$j] eq $in->{value}) {return 1;} }

   return 0 ;
}


=head2 calc_correlation_coeff()

   Title:       calc_correlation_coeff()
   Function:    Routine to compare probe score between two datasets
   Args:        ->{mode} = 'spearman' or 'pearson'

   Returns:     ->{r} = correlation coefficient (spearman or pearson)
                ->{full_stats} = full correlation coefficient info returned from
                  _correlation_coefficient_spearman() or
                  _correlation_coefficient_pearson()

   brought in from fpdgenlib::chip::chip

=cut

sub calc_correlation_coeff {

   my $in = shift;
   my $stats ;
   if (!exists $in->{mode}) { $in->{mode} = 'spearman'; }
   $stats->{mode} = $in->{mode} ;

   if ($in->{mode} eq 'pearson') {
      $stats->{full_stats} = _calc_correlation_pearson({ #pass in score arrays
         a => $in->{a},
         b => $in->{b},
      }) ;
      $stats->{r} = $stats->{full_stats}->{pearson_r} ;

      if (exists $in->{pval_numshuffle} && $in->{pval_numshuffle} > 0) {
         my @randa = @{$in->{a}} ;
         my $n_shuffle_better = 0;
         foreach my $i (1 .. $in->{pval_numshuffle}) {
            fy_shuffle(\@randa) ;
            my $shuffle_stat = _calc_correlation_pearson({ a => \@randa,
                                                           b => $in->{b} }) ;
            if ($shuffle_stat->{pearson_r} >= $stats->{r}) {
               $n_shuffle_better++ ; }
         }
         $stats->{sig}->{n_shuffle_better} = $n_shuffle_better ;
         $stats->{sig}->{n_shuffle} = $in->{pval_numshuffle} ;
         $stats->{sig}->{pval} = sprintf("%.8g", $n_shuffle_better /
                                                 $in->{pval_numshuffle}) ;
      } else {
         $stats->{sig}->{pval} = 'NA' ;
      }

   } elsif ($in->{mode} eq 'spearman') {
      $stats->{full_stats} = _calc_correlation_spearman({
         a => $in->{a},
         b => $in->{b},
      }) ;
      $stats->{r} = $stats->{full_stats}->{spearman_r} ;

      if (exists $in->{pval_numshuffle} && $in->{pval_numshuffle} > 0) {
         my @randa = @{$in->{a}} ;
         my $n_shuffle_better = 0;
         foreach my $i (1 .. $in->{pval_numshuffle}) {
            fy_shuffle(\@randa) ;
            my $shuffle_stat = _calc_correlation_spearman({ a => \@randa,
                                                            b => $in->{b} }) ;
            if ($shuffle_stat->{spearman_r} >= $stats->{r}) {
               $n_shuffle_better++ ; }
         }
         $stats->{sig}->{n_shuffle_better} = $n_shuffle_better ;
         $stats->{sig}->{n_shuffle} = $in->{pval_numshuffle} ;
         $stats->{sig}->{pval} = sprintf("%.8g", $stats->{n_shuffle_better} /
                                                 $in->{pval_numshuffle}) ;
      } else {
         $stats->{sig}->{pval} = 'NA' ;
      }

   }

   return $stats ;

}


=head2 _calc_correlation_pearson()

   Title:       _calc_correlation_pearson()
   Function:    Routine to calculate pearson correlation between 2 arrays
   Args:        ->{a} = arrayref of data set 1
                ->{b} = arrayref of data set 2

   Returns:     ->{a|b}->{total} = total of scores in each dataset
                ->{a|b}->{n} = number of points in each dataset
                ->{a|b}->{mean} = mean score in each dataset

                ->{a|b}->{ss}, ->{ab}->{ss} = pieces of the Pearson score eqn

                ->{pearson_r} = Pearson correlation coefficient

=cut

sub _calc_correlation_pearson {

   my $in = shift ;

   my $stats ;
   foreach my $side (qw/a b/) {
      foreach my $j (0 .. $#{$in->{$side}}) {
         $stats->{$side}->{total} += $in->{$side}->[$j] ;
      }
      $stats->{$side}->{n} = $#{$in->{$side}} + 1 ;
      $stats->{$side}->{mean}= $stats->{$side}->{total} / $stats->{$side}->{n};
   }

   foreach my $j (0 .. $#{$in->{a}}) {
      my $cur_diffs ;
      foreach my $side (qw/a b/) {
         $cur_diffs->{$side} = $in->{$side}->[$j] - $stats->{$side}->{mean};
         $stats->{$side}->{ss} += $cur_diffs->{$side} * $cur_diffs->{$side} ;
      }
      $stats->{ab}->{ss} += $cur_diffs->{a} * $cur_diffs->{b} ;
   }

# fpd111206_0950 - new check, if either one of vectors has 0 variance, return 0.
# TODO - double check, is this ok, or return an UNDEF?
   if ($stats->{a}->{ss} == 0 ||
       $stats->{b}->{ss} == 0) {
      $stats->{pearson_r} = 0 ;
   } else {
      $stats->{pearson_r} = $stats->{ab}->{ss} /
         sqrt($stats->{a}->{ss} * $stats->{b}->{ss})  ;
   }

   return $stats ;

}


=head2 _calc_correlation_spearman()

   Title:       _calc_correlation_spearman()
   Function:    Routine to calculate spearman correlation between 2 arrays
                  Returns Pearson R if there are tied ranks.

   Args:        ->{a} = arrayref of data set 1
                ->{b} = arrayref of data set 2

   Returns:     ->{spearman_r} = Spearman correlation coefficient
                  (if there were tied ranks, then a Pearson's R is returned)

=cut

sub _calc_correlation_spearman {

   my $in = shift ;

   my $score_ordered_probes ;
   my $ranks ; my $tied_ranks_fl = 0 ;
   foreach my $side (qw/a b/) {
      @{$score_ordered_probes->{$side}} = sort {
         $in->{$side}->[$b] <=> $in->{$side}->[$a]} (0 .. $#{$in->{$side}}) ;
# break ties
      my $ties ;
      $ranks->{$side}->[$score_ordered_probes->{$side}->[0]] = 0 ;
      foreach my $j ( 1 .. $#{$score_ordered_probes->{$side}}) {
         my $cur_ind = $score_ordered_probes->{$side}->[$j] ;
         my $last_ind = $score_ordered_probes->{$side}->[($j - 1)] ;
         $ranks->{$side}->[$cur_ind] = $j ;
         my $curscore = $in->{$side}->[$cur_ind] ;
         if ($curscore == $in->{$side}->[$last_ind]) {
            $ties->{$curscore}->{$cur_ind} = $j ;
            $ties->{$curscore}->{$last_ind} = $j - 1 ;
         }
      }

      foreach my $tied_score (keys %{$ties}) {
         $tied_ranks_fl = 1 ;
         my $rank_sum = 0; my $num_points = 0;
         foreach my $id (keys %{$ties->{$tied_score}}) {
            $num_points++ ;
            $rank_sum += $ties->{$tied_score}->{$id} ;
         }
         my $rank_mean = $rank_sum / $num_points ;
         foreach my $id (keys %{$ties->{$tied_score}}) {
            $ranks->{$side}->[$id] = $rank_mean ; }
      }
   }


   my $stats ;
   if ($tied_ranks_fl) {
      my $pearson_stats = _calc_correlation_pearson({
         a => $ranks->{a},
         b => $ranks->{b},
      }) ;
      $stats->{spearman_r} = $pearson_stats->{pearson_r} ;
   } else { #otherwise can use shortcut
      my $ss_rank_diff = 0;
      my $num_points = $#{$ranks->{a}} + 1 ;
      foreach my $j (0 .. $#{$ranks->{a}}) {
         my $cur_rank_diff = $ranks->{a}->[$j] - $ranks->{b}->[$j];
         $ss_rank_diff += $cur_rank_diff * $cur_rank_diff ;
      }
      $stats->{spearman_r} = 1 - ((6 * $ss_rank_diff) /
                                  ($num_points * $num_points - 1)) ;
   } 

   return $stats ;

}


=head2 safe_copy()

   Function:    Safely copy a file to a directory (using File::Copy::copy),
                  retries 14 times, and prints an error if it didnt work
   Return:      nothing
   Args:        $_[0] = source filename
                $_[1] = target directory

=cut

sub safe_copy {

   my $file = shift ;
   my $dest = shift ;
   my $tries = 15 ;

   if (!-s $file) {
      print STDERR "ERROR: couldnt find $file\n"; return}

   my $res = 0 ;
   while (($tries > 0) && ($res == 0 )) {
      $res = File::Copy::copy($file, $dest) ;
      $tries-- ;
   }

   if (!-s $dest) {
      print STDERR "ERROR: couldnt copy $file to $dest\t$!\n"; }

   return ;

}


sub convert_3d_to_2d_coords {
# in: parsed XPZ, 3D x,y,z coordinates
# out: 2D image id, x,y coordinates
#
# purpose: convert 3D x,y,z ref atlas coordinates to 2D image x,y,z
# 1. apply the image-series' TRV (reference to specimen volume)
# 2. choose closest z-slice image.
# 3. apply the image's TVS (specimen volume to sample image) 

   my $in = shift ;
   my $specs ;
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   my $xpz = $in->{xpz} ;

# hold in 4D homogeneous coordinate to simplify transforms
   my $x_3d = [[ @{$in->{x_3d}}, 1]] ;
#   print STDERR "original 3D coordinates:\n";
#   print STDERR print_matrix_octave($x_3d)."\n";

   my $trv = [] ;
   {
      my @t = split(',', $xpz->{meta}->{trv}) ;
      push @{$trv}, [$t[0], $t[3], $t[6], 0] ;
      push @{$trv}, [$t[1], $t[4], $t[7], 0] ;
      push @{$trv}, [$t[2], $t[5], $t[8], 0] ;
      push @{$trv}, [$t[9], $t[10], $t[11], 1] ;
   }
#   print STDERR "TRV:\n".print_matrix_octave($trv)."\n\n";

# Transform from reference to specimen volume space
   my $x_specvol = mmult($x_3d, $trv) ;
#   print STDERR "coordinates in specimen volume:\n";
#   print STDERR print_matrix($x_specvol)."\n";

   my $est_zval = $x_specvol->[0]->[2];

# Find image with closest z_value
   my @zval = sort {$a <=> $b} keys %{$xpz->{meta}->{images}} ;
# need something like a binary search
   my ($closest_zval, $closest_zdist) ;
   foreach my $zval ( @zval ) {
      if (defined $closest_zdist &&
          abs($zval - $est_zval) > $closest_zdist) {next;}
      $closest_zval = $zval ;
      $closest_zdist = abs($zval - $est_zval) ;
   }

# Transform from reference to specimen volume space
   my $tvs = [] ;
   {
      my @t = split(',', $xpz->{meta}->{images}->{$closest_zval}->{tvs}) ;
      $tvs = [ [ $t[0], $t[2], 0,  0],
               [ $t[1], $t[3], 0,  0],
               [     0,     0, 1,  0],
               [ $t[4], $t[5], 0,  1] ] ;
   }


   my $x_image  = mmult($x_specvol, $tvs) ;
#   print STDERR "coordinates in image (z = $closest_zval):\n";
#   print STDERR print_matrix_octave($x_image)."\n\n";

   return {
      z => $closest_zval,
      x => $x_image->[0],
   } ;

}

sub mmult {
# code from: http://rosettacode.org/wiki/Matrix_multiplication#Perl
   our @a; local *a = shift;
   our @b; local *b = shift;
   my @p = [];
   my $rows = @a;
   my $cols = @{ $b[0] };
   my $n = @b - 1;
   for (my $r = 0 ; $r < $rows ; ++$r) {
      for (my $c = 0 ; $c < $cols ; ++$c) {
         $p[$r][$c] += $a[$r][$_] * $b[$_][$c] foreach 0 .. $n; } }
   return [@p];
}


sub print_matrix {

   my $mat = shift ;
   my $out ;
   foreach my $j ( 0 .. $#{$mat}) {
      foreach my $k ( 0 .. $#{$mat->[$j]}) {
         $out .= $mat->[$j]->[$k]."\t";
      }
      $out .= "\n" ;
   }
   return $out ;

}


sub print_matrix_octave {

   my $mat = shift ;
   my $out = '[';
   foreach my $j ( 0 .. $#{$mat}) {
      foreach my $k ( 0 .. $#{$mat->[$j]}) {
         $out .= $mat->[$j]->[$k].',';
      }
      $out = substr($out,0,length($out) - 1) ;
      $out .= ";" ;
   }
   $out = substr($out,0,length($out) - 1) ;
   $out .= ']' ;
   return $out ;

}


=head2 fy_shuffle()

   Title:       fy_shuffle()
   Function:    Shuffles an array in place
   Args:        $_ = arrayref
   Returns:     Nothing

=cut

sub fy_shuffle {
   my $deck = shift;  # $deck is a reference to an array
   my $i = @$deck;
   while ($i--) {
      my $j = int rand ($i+1);
      @$deck[$i,$j] = @$deck[$j,$i];
   }
}


1;
