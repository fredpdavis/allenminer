=head1 NAME

alnmnr::runs.pm - interface routines for different ALLENMINER run modes

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

package alnmnr::runs ;
use strict;
use warnings;
use Cwd qw/getcwd abs_path/;
use File::Path qw/mkpath/ ;
use File::Temp qw/tempfile tempdir/ ;
use File::Basename ;
use POSIX qw/floor ceil/ ;

use alnmnr::core ;
use alnmnr::io ;
use alnmnr::roiops ;


=head2 _run_search()

   Title:       _run_search()
   Function:    Calculates expression values within a user defined region
                 of interest

   Args:        $_->{options} = any command line options parsed by getinp()
                $_->{cluster_fl} = 1 if to be dispatched to cluster
                $_->{specs}->{roidef_fn} - name of file specifying ROI(s)

                $_->{format_fl} [optional] - 1 to make plots of expression
                 level distributions in each specified ROI; requires R

                $_->{xpr_list} = [xpr file path 1, 2, ...n] - [optional]
                $_->{specs} = allenminer specs [optional]


   Returns:     to STDOUT or file output (if cluster)
                1. xpz file basename
                2. gene name
                3. Slice (coronal or sagittal)
                4. Number of bits set = number of ABA regions with expressing
                     voxel

=cut

sub _run_search {

   my $in = shift ;
   my $age = $in->{age} ;
   my $roidef = $in->{roidef} ;
   my $specs = $in->{specs} ;
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   if (!exists $in->{display_headers_fl}) {
      $in->{display_headers_fl} = 1 ; }

   if (!exists $in->{patterning}) {
      $in->{patterning} = 0 ; }

   if (!exists $in->{lrcheck}) {
      $in->{lrcheck} = 0 ; }
   if ($in->{lrcheck}) {
      $in->{patterning} = 1 ; }

# If the directory is not a full path, change it to one.
   if (exists $in->{xyz_details_out_dir}) {
      $in->{xyz_details_out_dir} = abs_path($in->{xyz_details_out_dir}) ; }

# Default to current directory; set this so cluster nodes get a common
#  output directory even if not specified.
   if (!exists $in->{xyz_details_out_dir}) {
      $in->{xyz_details_out_dir} = getcwd; }

#NOTE 111206_1532  Add this behavior to all run modes
#  - ie compress_fl option, and if out_fn specified with a gz$ file, compress it

   my ($results_fn, $results_fh) ;
   if (!exists $in->{out_fn}) {
      ($results_fh, $results_fn)=
         tempfile("results_search.XXXXX", SUFFIX => ".out") ;
   } else {
      $results_fn = $in->{out_fn} ;
      if ($results_fn =~ /gz$/) {
         $results_fn =~ s/\.gz$// ;
      }
      open($results_fh, ">$results_fn") ;
   }

#1. get list of expression files to process
# (1) run full-blown on everything
# (2) run full-blown on a specified list of xpr files

   my $xpz_prelist ;
   if (!exists $in->{xpz_list}) {
      @{$xpz_prelist} = glob($specs->{allen_xpz_dir}.'/'.$age.'/*/*/*') ;
   } else {
      $xpz_prelist = $in->{xpz_list} ;
   }
#   print STDERR "TARGETS: ".join(", ", @{$xpz_prelist})."\n";


   if (!exists $in->{cluster_fl}) {
      $in->{cluster_fl} = $specs->{cluster}->{cluster_mode} ; }

# If not specified, set to default cluster mode
   if (!exists $in->{cluster_fl}) {
      $in->{cluster_fl} = $specs->{cluster}->{cluster_mode} ; }

# If cluster flag set, set specs
   if ($in->{cluster_fl}) {
      foreach my $t_key (keys %{$specs->{cluster}}) {
         if ($t_key eq 'cluster_mode') {next;}
         if (!exists $in->{cluster}->{$t_key}) {
            $in->{cluster}->{$t_key} = $specs->{cluster}->{$t_key} ; }
      }
   }

   my $xpz_list = $xpz_prelist;

#2. dispatch equal numbers to cluster nodes; merge aftewards
   my $search_headers = alnmnr::core::calc_expr({
                           return_headers_only  => 1,
                           patterning           => $in->{patterning}
                        }) ;

   if ($in->{display_headers_fl}) {
      foreach my $type (sort keys %{$search_headers}) {
         print {$results_fh} '#'.join("\t",
               @{$search_headers->{$type}})."\n";
      }
   }

   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1 ) { # master node:
# 1. split input, 2. cluster run self with xpz_list_fn specified, cluster_fl = 0
# 3. merge STDOUT and STDERR output

      print "* _run_search() ".localtime() if(!exists $in->{quiet_fl});
      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{run_search_in},
       $temp_fn->{run_search_in}) =
       tempfile("splits_run_search_input.XXXXX") ;

      ($temp_fh->{run_search_err},
       $temp_fn->{run_search_err}) =
       tempfile("splits_run_search_SGEerr.XXXXX") ;

      foreach my $xpz_fn (@{$xpz_list}) {
         print {$temp_fh->{run_search_in}} $xpz_fn."\n" ; }
      close($temp_fh->{run_search_in}) ;

      my $split_dir = tempdir("splits_run_search.XXXXX") ;
      my $splits = SGE::_clust_split_ins({
         fn => $temp_fn->{run_search_in},
         dir => $split_dir,
         numjobs => $in->{cluster}->{numjobs},
      }) ;

      my ($perlscript_fh, $perlscript_fn) =
         tempfile("am.run_search.XXXXX", SUFFIX =>'.ami.pl') ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use alnmnr ;

main() ;

sub main {

   alnmnr::run({
      cluster_fl        => 0,
      patterning        => ".$in->{patterning}.",
      lrcheck           => ".$in->{lrcheck}.",
      roidef_fn         => \"".$in->{roidef_fn}."\",
      ARGV              => \\\@ARGV,
      out_fn            => \"-\",
      display_fl        => 1,
   }) ;

}
" ;

      my $roidef_fnbase = basename($in->{roidef_fn}) ;

      my ($sgescript_fh, $sgescript_fn)  =
         tempfile("am.run_search.XXXXX", SUFFIX => ".SGE.sh") ;
      my $sge_outdir = tempdir("SGEOUT.run_search.XXXXX") ;

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $in->{cluster}->{nodespecs}) {
         print {$sgescript_fh} $in->{cluster}->{nodespecs}."\n" ;}

      my $am_command = "perl $perlscript_fn -mode search ".
                       "-age $age -cluster_fl 0 ".
                       "-xpz_list_fn \$input1 ".
                       "-roidef_fn $roidef_fnbase ".
                       "-display_headers_fl 0";

      foreach my $option_type (qw/xyz_details xyz_details_compress xyz_details_out_suffix xyz_details_out_dir xyz_details_out_subdir/) {
         if (exists $in->{$option_type}) {
            $am_command .= ' -'.$option_type."  ".$in->{$option_type} ; } }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1 = \$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir
cp ".$in->{roidef_fn}." \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
$am_command

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn \$scratchdir/$roidef_fnbase
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = SGE::_clust_qsub({
         host => $in->{cluster}->{head_node},
         sgescript_fn => $sgescript_fn,
      }) ;
      print STDERR "      job $qsub_job_id\n" ;

      while (1) {
         sleep $in->{cluster}->{qstat_sleep} ;
         my $job_status = SGE::_clust_qstat({
            host => $in->{cluster}->{head_node},
            job_id => $qsub_job_id
         });
         if ($job_status) {last;}
      }

#3. get the results and piece them back together
      close($results_fh) ;

      my $header_line_strings = [];
      foreach my $type (keys %{$search_headers}) {
         push @{$header_line_strings},
               '#'.join("\t", @{$search_headers->{$type}}); }
      my $header_line = join("\n", @{$header_line_strings}) ;

      SGE::_clust_merge_outs({
         header         => $header_line,
         script_fn      => $sgescript_fn,
         out_fn         => $results_fn,
         err_fn         => $temp_fn->{run_search_err},
         job_id         => $qsub_job_id,
         outdir         => $sge_outdir,
         numjobs        => $splits->{numjobs}
      }) ;

   } else {

# If -lrcheck specified, create -LR_part0 and -_part1 ROI for all ROI
#     and set -patterning option
      if ($in->{lrcheck}) { #do it once and print the ROI to a temp file for calcs?
         print STDERR "Prepare for -lrcheck:\n" ;
         alnmnr::roiops::create_lrcheck_roi({
            roidef      => $roidef,
            age         => $age,
         }) ;
         $in->{patterning} = 1 ;

         if ($in->{lrcheck} == 2) {
            my ($t_fh, $t_fn) = tempfile("lrcheck.XXXXX", SUFFIX=>".roi");
            close($t_fh) ;
            alnmnr::io::display_roidef({
               age    => $age,
               out_fn => $t_fn,
               roidef => $roidef,
            }) ;
         }
      }

#3. actually run the roidefs against each xpz file now
      foreach my $j (0 .. $#{$xpz_list}) {
         my $xpz_fn = $xpz_list->[$j] ;
#         print STDERR "now on $xpz_fn\n" ;
         if (!-s $xpz_fn) {
            print STDERR "SKIPPING MISSING OR EMPTY FILE: $xpz_fn\n" ;
            next;
         }
         my $xpz_data = alnmnr::io::aba_parse_xpz({fn => $xpz_fn}) ;

         if (exists $xpz_data->{error_fl}) {
            print STDERR "ERROR: skipping $xpz_fn - not read properly: ".
                         $xpz_data->{error_fl}."\n" ;
            next;
         }

         my $options_search = {
            display_headers     => 0,
            out_fh              => $results_fh,
            xpz                 => $xpz_data,
            roidef              => $roidef,
            display_fl          => 1,
            age                 => $age,
            patterning          => $in->{patterning}
         } ;

# if the xyz_details flag is set, make sure it's output file is specified
         if (exists $in->{xyz_details} && $in->{xyz_details} == 1) {
            my $xpz_base = basename($xpz_fn) ;
            my $details_out_fn = $xpz_base ;
            $details_out_fn =~ s/\.xpz$// ;
            if (exists $in->{xyz_details_out_suffix}) {
               $details_out_fn .= '.'.$in->{xyz_details_out_suffix} ;
            } else {
               $details_out_fn .= ".$$.roi_xyz_details" ;
            }

            if (exists $in->{xyz_details_out_subdir}) {
               $details_out_fn = substr($details_out_fn,0,1).'/'.
                                 $details_out_fn ;
            }

            if (exists $in->{xyz_details_out_dir}) {
               $details_out_fn = $in->{xyz_details_out_dir}."/".
                                 $details_out_fn ;
            }

            $options_search->{xyz_details_out_fn} =
               $details_out_fn;

            if (exists $in->{xyz_details_compress} &&
                $in->{xyz_details_compress} == 0) {
               $options_search->{xyz_details_compress_fl} = 1;}
         }
         alnmnr::core::calc_expr($options_search);
      }
      close($results_fh) ;
   }


   if ((exists $in->{out_fn} && $in->{out_fn} =~ /gz$/) ||
       (exists $in->{compress_fl} && $in->{compress_fl} == 1)) {
      my $gzip_tcom = "gzip $results_fn" ;
      system($gzip_tcom) ;
   }

}


=head2 _run_simsearch()

   Title:       _run_simsearch()
   Function:    Given query XPZ file(s), calculates expression similarity
                 to other XPZ files; optionally within defined ROI.

   Args:        $_->{options} = any command line options parsed by getinp()
                $_->{cluster_fl} = 1 if to be dispatched to cluster
                $_->{specs}->{roidef_fn} - name of file specifying ROI(s)

                $_->{format_fl} [optional] - 1 to make plots of expression
                 level distributions in each specified ROI; requires R

                $_->{xpz_query_list} = [xpz file path 1, 2, ...n] - [optional]

                $_->{xpz_target_list} = [xpz file path 1, 2, ...n] - [optional]
                $_->{specs} = allenminer specs [optional]

   Returns:     to STDOUT or file output (if cluster)
                1. query XPZ file basename
                2. query gene name
                3. query slice (coronal or sagittal)
                4. target XPZ file basename
                5. target gene name
                6. target slice (coronal or sagittal)
                7. roi
                8. pearson's correlation
                9. number of voxels compared

=cut

sub _run_simsearch {

#1. get list of xpz files

   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $age = $in->{age} ;
   my $specs = $in->{specs};
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   if (!exists $in->{display_headers_fl}) {
      $in->{display_headers_fl} = 1 ; }

   if (!exists $in->{lrcheck}) {
      $in->{lrcheck} = 0 ; }

   my $pval_numshuffle = 0 ;
   if (exists $in->{calc_pval}) {
      $pval_numshuffle = $in->{calc_pval} ; }

# If not specified, set to default cluster mode
   if (!exists $in->{cluster_fl}) {
      $in->{cluster_fl} = $specs->{cluster}->{cluster_mode} ; }

# If cluster flag set, set specs
   if ($in->{cluster_fl}) {
      foreach my $t_key (keys %{$specs->{cluster}}) {
         if ($t_key eq 'cluster_mode') {next;}
         if (!exists $in->{cluster}->{$t_key}) {
            $in->{cluster}->{$t_key} = $specs->{cluster}->{$t_key} ; }
      }
   }

# Set up output file
   my ($results_fn, $results_fh) ;
   if (!exists $in->{out_fn}) {
      ($results_fh, $results_fn)=
         tempfile("results_simsearch_query.XXXXX", SUFFIX => ".out");
   } else {
      $results_fn = $in->{out_fn} ;
      open($results_fh, ">$results_fn") ;
   }

# set up query and target XPZ file lists
   my $xpz_query_list = $in->{xpz_query_list};
   my $xpz_target_list = [];
   if (!exists $in->{xpz_target_list}) {

      if (exists $in->{plane}) {
         @{$xpz_target_list} = glob($specs->{allen_xpz_dir}.'/'.$age.
                                    '/*/*/*'."_".$in->{plane}.'_*') ;
      } else {
         @{$xpz_target_list} = glob($specs->{allen_xpz_dir}.'/'.$age.
                                    '/*/*/*');
      }

   } else {

      if (exists $in->{plane}) {
         foreach my $t_fn (@{$in->{xpz_target_list}}) {
            if ($t_fn !~ /_$in->{plane}_/) {next;}
            push @{$xpz_target_list}, $t_fn ;
         }
      } else {
         $xpz_target_list = $in->{xpz_target_list} ;
      }
   }


   my $headers = [(qw/SIM query_file query_gene query_plane hit_file hit_gene/,
                   qw/hit_plane ROI similarity_score pval num_shared_cubes/,
                   qw/query_mean_cube_level hit_mean_cube_level/)] ;
   my $header_line = '#'.join("\t", @{$headers}) ;


   if ($in->{display_headers_fl}) {
      print {$results_fh} '#'.join("\t", @{$headers})."\n" ; }

#2. dispatch equal numbers to cluster nodes; merge aftewards
   if ( exists $in->{cluster_fl} && $in->{cluster_fl} == 1 ) { # master node:
# 1. split input, 2. cluster run self with xpz_target_list_fn specified,
#    cluster_fl = 0
# 3. merge STDOUT and STDERR output

      print "* _run_simsearch() ".localtime() if(!exists $in->{quiet_fl});
      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{target_list_in}, $temp_fn->{target_list_in}) =
       tempfile("splits_run_simsearch.target_list.XXXXX") ;

      ($temp_fh->{query_list_in}, $temp_fn->{query_list_in}) =
       tempfile("splits_run_simsearch.query_list.XXXXX") ;

      ($temp_fh->{SGEerr}, $temp_fn->{SGEerr}) =
       tempfile("splits_run_simsearch.XXXXX") ;

      foreach my $xpz_target_fn (@{$xpz_target_list}) {
         print {$temp_fh->{target_list_in}} $xpz_target_fn."\n";}
      close($temp_fh->{target_list_in}) ;

      foreach my $xpz_query_fn (@{$xpz_query_list}) {
         print {$temp_fh->{query_list_in}} $xpz_query_fn."\n";}
      close($temp_fh->{query_list_in}) ;

      my $split_dir = tempdir("splits_run_simsearch.XXXXX") ;
      my $splits = SGE::_clust_split_ins({
         fn => $temp_fn->{target_list_in},
         dir => $split_dir,
         numjobs => $in->{cluster}->{numjobs},
      }) ;

      my ($perlscript_fh, $perlscript_fn) =
         tempfile("am.run_simsearch_query.XXXXX", SUFFIX =>'.ami.pl') ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use alnmnr;

main() ;

sub main {

   alnmnr::run({
      cluster_fl        => 0,
      xpz_query_list_fn => \"".$temp_fn->{query_list_in}."\",
      calc_pval         => ".$pval_numshuffle.",
      lrcheck           => ".$in->{lrcheck}.",
      ARGV              => \\\@ARGV,
      out_fn            => \"-\",
      display_fl        => 1,
   }) ;

}
" ;

      my ($sgescript_fh, $sgescript_fn)  =
         tempfile("am.run_simsearch_query.XXXXX", SUFFIX => ".SGE.sh") ;
      my $sge_outdir = tempdir("SGEOUT.run_simsearch_query.XXXXX") ;

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $in->{cluster}->{nodespecs}) {
         print {$sgescript_fh} $in->{cluster}->{nodespecs}."\n" ;}

      my $am_command = "perl $perlscript_fn -mode simsearch ".
                       "-display_headers_fl 0 ".
                       "-xpz_target_list_fn \$input1 -age $age -cluster_fl 0 ".
                       "-xpz_query_list_fn ".$temp_fn->{query_list_in} ;
      my $roidef_fnbase = basename($in->{roidef_fn}) ;
      if (exists $in->{roidef_fn}) {
         $am_command .= ' -roidef_fn '.$roidef_fnbase ; }

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1 = \$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir
cp ".$temp_fn->{query_list_in}." \$scratchdir\n" ;

      if (exists $in->{roidef_fn}) {
         print {$sgescript_fh} "cp ".$in->{roidef_fn}." \$scratchdir\n\n"; }

      print {$sgescript_fh} "cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
$am_command

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn \$scratchdir/".$temp_fn->{query_list_in}."\n";

      if (exists $in->{roidef_fn}) {
         print {$sgescript_fh} "rm -f \$scratchdir/$roidef_fnbase\n"; }

      print {$sgescript_fh} "cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = SGE::_clust_qsub({
         host => $in->{cluster}->{head_node},
         sgescript_fn => $sgescript_fn,
      }) ;
      print STDERR "      job $qsub_job_id\n" ;

      while (1) {
         sleep $in->{cluster}->{qstat_sleep} ;
         my $job_status= SGE::_clust_qstat({
            host => $in->{cluster}->{head_node},
            job_id => $qsub_job_id
         });
         if ($job_status) {last;}
      }

#3. get the results and piece them back together
      close($results_fh) ;
      SGE::_clust_merge_outs({
         header         => $header_line,
         script_fn      => $sgescript_fn,
         out_fn         => $results_fn,
         err_fn         => $temp_fn->{SGEerr},
         job_id         => $qsub_job_id,
         outdir         => $sge_outdir,
         numjobs        => $splits->{numjobs}
      }) ;

   } else {

#3. actually run the roidefs against each xpr file now

#NOTE 111122_0845 - not sure but probably don't have to downsample adult rois, since XPZ is full res
#      my $roidef_200 ;
#      if ($age eq 'adult') {
#         $roidef_200 = convert_roidef_100_to_200({ roidef => $roidef }) ;
#      } else { #otherwise, no conversion; XPZ already synced with ROI coordinates
#         $roidef_200 = $roidef ;
#      }


# If -lrcheck specified, create -LR_part0 and -_part1 ROI for all ROI
      if ($in->{lrcheck}) { #do it once and print the ROI to a temp file for calcs?
         print STDERR "Prepare for -lrcheck:\n" ;
         alnmnr::roiops::create_lrcheck_roi({
            roidef      => $roidef,
            age         => $age,
         }) ;
         $in->{patterning} = 1 ;
      }

      my $roidef_200 = $roidef;

      print STDERR "Starting search\n";
      my ($query_xpz, $query_meta) ;
      foreach my $j ( 0 .. $#{$xpz_query_list}) {
         my $xpz_query_fn = $xpz_query_list->[$j];
         print STDERR "ON QUERY: $xpz_query_fn\n";

         my $cur_query_xpz = alnmnr::io::aba_parse_xpz({
                                 fn             => $xpz_query_fn,
                                 energy_only    => 1,
                                 roidef         => $roidef_200 })  ;

         if (exists $cur_query_xpz->{error_fl}) {
            print STDERR "ERROR: skipping QUERY $xpz_query_fn ".
                         "- not read properly: ".
                         $cur_query_xpz->{error_fl}."\n" ;
            next;
         }

         my $num_query_points = keys %{$cur_query_xpz->{energy}} ;
         print STDERR "NUM QUERY POINTS: $num_query_points\n";
         if ($num_query_points == 0) {
            print STDERR "NOTE: skipping QUERY $xpz_query_fn ".
                         "- no expressing voxels (in ROI if specified)\n";
            next;
         }

         $query_xpz->[$j] = $cur_query_xpz->{energy} ;

         $query_meta->[$j]->{basename} = basename($xpz_query_fn) ;
         $query_meta->[$j]->{plane} = $cur_query_xpz->{meta}->{plane} ;
         $query_meta->[$j]->{gene_name} = $cur_query_xpz->{meta}->{gene_name} ;
      }

      foreach my $j (0 .. $#{$xpz_target_list}) {
         my $xpz_target_fn = $xpz_target_list->[$j] ;
         my $target_xpz = alnmnr::io::aba_parse_xpz({
                              fn                => $xpz_target_fn,
                              energy_only       => 1,
                              roidef            => $roidef_200}) ;

         if (exists $target_xpz->{error_fl}) {
            print STDERR "ERROR: skipping TARGET $xpz_target_fn ".
                         "- not read properly: ".
                         $target_xpz->{error_fl}."\n" ;
            next;
         }

         my $num_target_points = keys %{$target_xpz->{energy}} ;
         if ($num_target_points == 0) {
            print STDERR "NOTE: skipping TARGET $xpz_target_fn ".
                         "- no expressing voxels (in ROI if specified)\n";
            next;
         }

         my $target_meta = {} ;
         $target_meta->{basename} = basename($xpz_target_fn) ;
         $target_meta->{plane} = $target_xpz->{meta}->{plane} ;
         $target_meta->{gene_name} = $target_xpz->{meta}->{gene_name} ;

         foreach my $j (0 .. $#{$query_xpz}) {
            my $options_simsearch = {
               out_fh           => $results_fh,
               xpz1_data        => $query_xpz->[$j],
               xpz2_data        => $target_xpz->{energy},
               roidef_mode      => 200,
               roidef           => $roidef_200,
               age              => $age,
               pval_numshuffle  => $pval_numshuffle,
            } ;

# also display the average per shared voxel in the ROI and
# average per all voxels in the ROI.
# helps to rank the genes to get not only similar correlation tracking
#  but also similar levels.

            my $expr_sim = alnmnr::core::calc_expr_sim(
               $options_simsearch);

            foreach my $roi_name (keys %{$expr_sim}) {
               if (exists $expr_sim->{$roi_name}->{error_fl}) {
                  $expr_sim->{$roi_name}->{r} =
                     $expr_sim->{$roi_name}->{error_fl} ;
                  $expr_sim->{$roi_name}->{sig}->{pval} = 'NA' ;
                  $expr_sim->{$roi_name}->{full_stats}->{a}->{n} = 0 ;
                  $expr_sim->{$roi_name}->{full_stats}->{a}->{mean} = 'NA' ;
                  $expr_sim->{$roi_name}->{full_stats}->{b}->{mean} = 'NA' ;
               } else {
                  $expr_sim->{$roi_name}->{r} =
                     sprintf("%.5f", $expr_sim->{$roi_name}->{r}) ;

                  $expr_sim->{$roi_name}->{full_stats}->{a}->{mean} =
                     sprintf("%.5f",
                        $expr_sim->{$roi_name}->{full_stats}->{a}->{mean}) ;

                  $expr_sim->{$roi_name}->{full_stats}->{b}->{mean} =
                     sprintf("%.5f",
                        $expr_sim->{$roi_name}->{full_stats}->{b}->{mean}) ;

               }

               print {$results_fh} join("\t",
                  "SIM",
                  $query_meta->[$j]->{basename},
                  $query_meta->[$j]->{gene_name},
                  $query_meta->[$j]->{plane},
                  $target_meta->{basename},
                  $target_meta->{gene_name},
                  $target_meta->{plane},
                  $roi_name,
                  $expr_sim->{$roi_name}->{r},
                  $expr_sim->{$roi_name}->{sig}->{pval},
                  $expr_sim->{$roi_name}->{full_stats}->{a}->{n},
                  sprintf("%.5f",
                     $expr_sim->{$roi_name}->{full_stats}->{a}->{mean}),
                  sprintf("%.5f",
                     $expr_sim->{$roi_name}->{full_stats}->{b}->{mean}),
               )."\n";
            }
         }
      }
      close($results_fh) ;
   }

}


=head2 _run_roi_flip()

   Title:       _run_roi_flip()
   Function:    Flips ROI(s) across the mid-sagittal plane
   Args:        $_->{roidef} - region of interested definitions

=cut

sub _run_roi_flip {

   my $in = shift ; #should come in with roidefs
   my $roidef = $in->{roidef} ;

#2. iterate over ROI and peel layers
   my $flipped_roi = alnmnr::roiops::flip_roi({
      age    => $in->{age},
      roidef => $in->{roidef},
   }) ;

#3. print out all ROI to file
   if (exists $in->{out_fn}) {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         out_fn => $in->{out_fn},
         roidef => $flipped_roi
      }) ;
   } else {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         roidef => $flipped_roi }) ;
   }

}


=head2 _run_roiop()

   Title:       _run_roiop()
   Function:    Parent routine for running ROI operations
   Args:        $_ - necessary options for the actual operation routine

=cut

sub _run_roiop {

   my $in = shift ; #should come in with roidefs
   if (   $in->{type} eq 'peel')          { _run_roi_peel($in); }
   elsif ($in->{type} eq 'flip')          { _run_roi_flip($in); }
   elsif ($in->{type} eq 'subtract')      { _run_roi_subtract($in); }
   elsif ($in->{type} eq 'rename')        { _run_roi_rename($in); }
   elsif ($in->{type} eq 'separate_adjacent') {_run_roi_separate_adjacent($in);}
   elsif ($in->{type} eq 'partition')     { _run_roi_partition($in); }

}


=head2 _run_convert()

   Title:       _run_convert()
   Function:    Parent routine for converting file formats
   Args:        $_ - necessary options for the actual operation routine

=cut

sub _run_convert {

   my $in = shift ; #should come in with roidefs
   if (   $in->{type} eq 'roi2pdb')             { _run_roi2pdb($in); }
   elsif ($in->{type} eq 'xpr2pdb')             { _run_xpr2pdb($in); }
   elsif ($in->{type} eq 'xpz2pdb')             { _run_xpz2pdb($in); }
   elsif ($in->{type} eq 'sva2pdb')             { _run_sva2pdb($in); }
   elsif ($in->{type} eq 'xprarchive2pdb')      { _run_xprarchive2pdb($in); }
   elsif ($in->{type} eq 'xpr2txt')             { _run_xpr2txt($in); }
   elsif ($in->{type} eq 'xpz2txt')             { _run_xpz2txt($in); }
   elsif ($in->{type} eq 'sva2txt')             { _run_sva2txt($in); }
   elsif ($in->{type} eq 'roi3d_to_2d')         { _run_roi3d_to_2d($in); }

}


=head2 _run_roi_peel()

   Title:       _run_roi_peel()
   Function:    Peels off edge layers from ROI(s)
   Args:        $_->{roidef} - region of interestd definitions
                $_->{layers} - number of edge layers to remove

=cut

sub _run_roi_peel {

   my $in = shift ; #should come in with roidefs
   my $roidef = $in->{roidef} ;
   my $n_layers = $in->{layers} ;

#2. iterate over ROI and peel layers
   my $peeled_roi = alnmnr::roiops::peel_roi({
      roidef => $in->{roidef},
      n_layers => $n_layers,
   }) ;

#3. print out all ROI to file
   if (exists $in->{out_fn}) {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         out_fn => $in->{out_fn},
         roidef => $peeled_roi
      }) ;
   } else {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         roidef => $peeled_roi }) ;
   }

}


=head2 _run_roi_subtract()

   Title:       _run_roi_subtract()
   Function:    Computes the difference between two ROIs
   Args:        $_->{roidef} - region of interestd definitions
                $_->{roi1} - name of first ROI
                $_->{roi2} - name of second ROI
                $_->{new_name} - new ROI name

=cut

sub _run_roi_subtract {

   my $in = shift ; #should come in with roidefs
   my $roidef = $in->{roidef} ;
   my $roi1 = $in->{roi1} ;
   my $roi2 = $in->{roi2} ;
   my $new_name = $in->{new_name} ;


#2. iterate over ROI and peel layers
   my $new_roi = alnmnr::roiops::subtract_roi({
      roidef => $in->{roidef},
      roi1 => $roi1,
      roi2 => $roi2,
      new_name => $new_name,
   }) ;

#3. print out all ROI to file
   if (exists $in->{out_fn}) {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         out_fn => $in->{out_fn},
         roidef => $new_roi
      }) ;
   } else {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         roidef => $new_roi }) ;
   }

}


=head2 _run_roi_rename()

   Title:       _run_roi_rename()
   Function:    Renames an ROI
   Args:        $_->{roidef} - region of interest definitions
                $_->{old_name} - old ROI name
                $_->{new_name} - new ROI name

=cut

sub _run_roi_rename {

   my $in = shift ; #should come in with roidefs
   my $roidef = $in->{roidef} ;
   my $old_name = $in->{old_name} ;
   my $new_name = $in->{new_name} ;

#2. iterate over ROI and peel layers
   my $new_roi = alnmnr::roiops::rename_roi({
      roidef => $in->{roidef},
      old_name => $old_name,
      new_name => $new_name,
   }) ;

#3. print out all ROI to file
   if (exists $in->{out_fn}) {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         out_fn => $in->{out_fn},
         roidef => $new_roi
      }) ;
   } else {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         roidef => $new_roi }) ;
   }

}



=head2 _run_roi_separate_adjacent()

   Title:       _run_roi_separate_adjacent()
   Function:    Removes ROI edge layers that border neighboring ROIs
   Args:        $_->{roidef} - region of interest definitions
                $_->{layers} - number of edge layers to remove

=cut

sub _run_roi_separate_adjacent {

   my $in = shift ; #should come in with roidefs
   my $roidef = $in->{roidef} ;
   my $n_layers = $in->{layers} ;

#1. build new ROIdefs by removing adjacent points
   my $shaven_roi = separate_adjacent_roi({
      roidef => $in->{roidef},
      n_layers => $n_layers,
   }) ;

#2. print out ROI to file
   if (exists $in->{out_fn}) {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         out_fn => $in->{out_fn},
         roidef => $shaven_roi
      }) ;
   } else {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         roidef => $shaven_roi }) ;
   }

}



=head2 _run_roi_partition()

   Title:       _run_roi_partition()
   Function:    computes axis-aligned sub-rois for entropy or gradient search
   Args:        $_->{roidef} - region of interest definition
                $_->{axes} - AP (x), DV (y), or ML (z)

=cut

sub _run_roi_partition {

   my $in = shift ;
   my $fitted_cuts = 0;
   if (exists $in->{fitted_cuts} && $in->{fitted_cuts} == 1) {
      $fitted_cuts = 1;}

   my $roi_basename = '' ;
   if (exists $in->{roi_basename}) {
      $roi_basename = $in->{roi_basename}; }
   my $roidef_partitions = alnmnr::roiops::calc_roi_partitions({
      roi_basename => $roi_basename,
      fitted_cuts => $fitted_cuts,
      roidef => $in->{roidef},
      numbins => $in->{numbins},
      axes => $in->{axes},
   }) ;

   if (exists $in->{out_fn}) {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         out_fn => $in->{out_fn},
         roidef => $roidef_partitions
      }) ;
   } else {
      alnmnr::io::display_roidef({
         age    => $in->{age},
         roidef => $roidef_partitions
      }) ;
   }

   if (exists $in->{pdb_fn}) {
      alnmnr::io::display_roidef_pdb({
         roidef => $roidef_partitions,
         out_fn => $in->{pdb_fn}
      }) ;
   }

}


=head2 _run_roi2pdb()

   Title:       _run_roi2pdb()
   Function:    computes axis-aligned sub-rois for entropy or gradient search
   Args:        $_->{roidef} - region of interest definition
                $_->{axes} - AP (x), DV (y), or ML (z)
                $_->{

=cut

sub _run_roi2pdb {

   my $in = shift ;

   alnmnr::io::display_roidef_pdb({
      age    => $in->{age},
      roidef => $in->{roidef},
      out_fn => $in->{pdb_fn}
   }) ;

}



sub _make_report {
# Goal: take -search or -simsearch results and format to html for easier review
# NOTE: DOESNT WORK FOR ALL SEARCH MODES YET 

   my $in = shift ;
   my $specs ;
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   my $roidef = $in->{roidef} ;
   my $age = $in->{age} ;

# um = view-volume * $scaling
   my $scaling_voxelnum_2_um =
      $specs->{view_volume}->{$in->{age}}->{spacing} ;

# 1. figure out search and displaly specs ----------------------------------------
# * figure out salient lines from search results: record_type
# * figure out search field - either user-specified or default for this search type

   #file handles to store splits of search results.
   my ($perbatch_fh, $perbatch_fn) ;

   my $search_specs = {} ; # Search type, sorting prefs info
   my $seen_rectypes = {} ; my $rec2f2i = {};
   open(RESF, $in->{search_results_fn}) ;
   while (my $line = <RESF>) {
      chomp $line;
      if ($line =~ /^\#/) {
         $line =~ s/^#// ;
         my @t = split(/\t/, $line) ;
         $seen_rectypes->{$t[0]}++ ;
         map {$rec2f2i->{$t[0]}->{$t[$_]} = $_; } (0 .. $#t) ;
      }
   }
   close(RESF) ;

   my $numrectypes = keys %{$seen_rectypes} ;
   if ($numrectypes == 1) {
      my ($rectype) = keys %{$seen_rectypes} ;
      $search_specs->{search} = $specs->{report}->{rec2specs}->{
                                   $rectype}->{search} ;
      $search_specs->{record} = $rectype ;
   } else {
      $search_specs->{search} = 'search_patterned' ;
   }


# Determine sort by field
   if (exists $in->{sortby}) {
      $search_specs->{sortby} = $in->{sortby} ;

# Set record type if not already there - which means search_patterned
      if (!exists $search_specs->{record}) {
         foreach my $rectype (keys %{$seen_rectypes}) {
            if (exists $rec2f2i->{$rectype}->{$in->{sortby}}) {
               $search_specs->{record} = $rectype;
               last;
            }
         }
         if (!exists $search_specs->{record})  {
            die "ERROR: unrecognized sortby" ; }
      }
   } else {
      if (!exists $search_specs->{record}) {
         $search_specs->{record} =
            $specs->{report}->{search2rec}->{$search_specs->{search}} ; }

      $search_specs->{sortby} =
         $specs->{report}->{rec2specs}->{$search_specs->{record}}->{sortby} ;
   }
   $search_specs->{splitby} = $specs->{report}->{rec2specs}->{
                                 $search_specs->{record}}->{splitby};
   $search_specs->{f2i} = $rec2f2i->{$search_specs->{record}} ;
   my $f2i = $search_specs->{f2i} ;

#   print STDERR "SEARCH: ".$search_specs->{search}."\n";
#   print STDERR "RECORD: ".$search_specs->{record}."\n";
#   print STDERR "SPLIT BY: ".join(", ", @{$search_specs->{splitby}})."\n";

   $search_specs->{splitby_ind} = [] ;
   foreach my $f ( @{$search_specs->{splitby}} ) {
      push @{$search_specs->{splitby_ind}},
         $search_specs->{f2i}->{$f} ; }

# BY NOW WE HAVE: search_specs->{search, record, sortby, splitby, f2i}

# 2. store results for each (splitby field) into separate files --------------

#   print STDERR "SPLITTING RESULTS FILE ".$in->{search_results_fn}.": " ;
   open(RESF, $in->{search_results_fn}) ;
   my $total_numhits = {} ;
   while (my $line = <RESF>) {
      chomp $line;

      my @t = split(/\t/, $line) ;

# If record_type is 'P', also need ROI names from 'S' lines
      if ($search_specs->{record} eq 'P' && $t[0] eq 'S') {
         my $roi = $t[$rec2f2i->{S}->{ROI}] ;
         my ($roi_basename) = ($roi =~ /(.+)_part[0-9]+$/) ;
         $search_specs->{roi_basename2roi}->{$roi_basename}->{$roi}++;
      }

# If this is the wrong record type, go to next line
      if ($t[0] ne $search_specs->{record}) {next;}

#NOTE: if record_type is 'E', ignore ROI=internal_background lines
      if ($search_specs->{record} eq 'E' &&
          $line =~ /internal_background/) {next; }


      if ((exists $in->{query_plane} &&
           $t[$f2i->{query_plane}] ne $in->{query_plane}) ||

          (exists $in->{hit_plane} &&
           $t[$f2i->{hit_plane}] ne $in->{hit_plane}) ||

          (exists $in->{plane} &&
           $t[$f2i->{plane}] ne $in->{plane})

         ) {next;}

# join splitby fields to get fh
      my $cursplit = join("\t", @t[@{$search_specs->{splitby_ind}}]);

      if (!exists $perbatch_fh->{$cursplit}) {
         ($perbatch_fh->{$cursplit}, $perbatch_fn->{$cursplit}) = tempfile() ; }
      print {$perbatch_fh->{$cursplit}} $line."\n";
      $total_numhits->{$cursplit}++ ;
   }
   close(RESF) ;
   map {close($perbatch_fh->{$_}); } (keys %{$perbatch_fh}) ;
#   print STDERR "X\n";



#3. Sort and format individual batches ----------------------------

   my $roi_slice_bounds = {} ;
   alnmnr::roiops::calc_roi_slicebounds({
      roi_slice_bounds => $roi_slice_bounds,
      roidef           => $roidef }) ;

# if simsearch, header rows = ref atlas, query slices
# if search,    header rows = ref atlas slices

   my $html_header_rows = '' ;
   my $html_header = "<html>
<head>
<title>ALLENMINER results</title>
<style>
body {
   font-family: \"trebuchet ms\", Helvetica, Verdana, Geneva, Arial;
   font-size: smaller;
}
td.headerrow {
   background-color: lightgrey
}
th {
   border-bottom: navy solid thin;
   border-top: navy solid thin;
   border-left: none;
   border-right: none;
   text-align: left;
}
</style>
</head>
<body>\n" ;


# sort by user-specified field -> tempfile2
   my $sort_field = $f2i->{$search_specs->{sortby}} + 1 ;
   my $hit_seen = {} ;
   my $navbar = '';
   foreach my $cursplit (keys %{$perbatch_fh}) {

      my $total_numpages = POSIX::ceil($total_numhits->{$cursplit} /
                                        $specs->{report}->{hits_per_page}) ;

      my @rois = $cursplit ;
      if ($search_specs->{record} eq 'FE') {
         @rois = split(/\t/, $cursplit) ;
      } elsif ($search_specs->{record} eq 'P') {
         @rois = sort keys %{$search_specs->{roi_basename2roi}->{$cursplit}} ;
      }
#      print STDERR "ROIs are: @rois\n";

      $html_header .= "<h1>ALLENMINER ".$search_specs->{search}.
                      " over $cursplit</h1>\n";

      my $html_page_num = 1 ;
      my $html_fh ;
      open($html_fh, ">".$in->{out_prefix}."_p".$html_page_num.".html") ;
      my $nextpage_url = $in->{out_prefix}."_p".($html_page_num + 1).".html" ;

      print {$html_fh} $html_header ;
      $navbar = "<p><h3>Page $html_page_num (of $total_numpages)";
#      $navbar = "<p><h3>Page $html_page_num";
      if ($html_page_num < $total_numpages) {
         $navbar .= " | <a href=\"$nextpage_url\">Go to next page</a>"; }
      $navbar .= "</h3></p>\n";
      print {$html_fh} $navbar ;
      print {$html_fh} "<table>\n" ;

      my $tcom = "sort -k$sort_field,$sort_field"."nr ".
                 $perbatch_fn->{$cursplit}." | " ;


# Get relevant list of ROIs
# if rectype not P, is the splitfield
# if rectype=P, use search_specs->{roi_basename2roi} list.


      my $numimages = $in->{num_images} || $specs->{report}->{num_images} ;

      my $refrows = []; my $total_numslices = 0 ;
      foreach my $roi (@rois) { # Get reference atlas images
         my $refatlas_planes = {} ;
         if (exists $in->{query_plane} && exists $in->{hit_plane} &&
             $in->{query_plane} eq $in->{hit_plane}) {
            $refatlas_planes->{$in->{query_plane}}++ ;
         } elsif (exists $in->{plane}) {
            $refatlas_planes->{$in->{plane}}++ ;
         } elsif ($age eq 'P56') {
            $refatlas_planes->{'coronal'}++ ;
            $refatlas_planes->{'sagittal'}++ ;
         } else {
            $refatlas_planes->{'sagittal'}++ ;
         }

         my $planenum = 0 ;
         foreach my $refatlas_plane (sort keys %{$refatlas_planes}) {

            my $planeinc = $planenum * 3 ; ; #if second plane to show, new rows

            my $images = _run_roi3d_to_2d({
               roidef   => $roidef,
               roi_name => $roi,
               specs    => $specs,
               age      => $age,
               oneimage_per_z3d => 1,
               refatlas_fl => 1,
               refatlas_plane => $refatlas_plane,
            }) ;
            my $plane_type = 'x2yz' ;
            my $width_var = 'z' ;
            my $plane_var = 'x' ;
            if ($refatlas_plane eq 'sagittal') {
               $plane_var = 'z' ;
               $width_var = 'x' ;
               $plane_type = 'z2xy' ; }

            my $numslices = keys %{$roi_slice_bounds->{$roi}->{$plane_type}} ;

            push @{$refrows->[(0 + $planeinc)]}, "<td colspan=\"$numslices\">".
               "Allen reference atlas ($age, $refatlas_plane): $roi</td>\n" ;

            my @td_headers ; my @td_images ; my @td_info ;
            foreach my $z3d_voxels (sort {$a<=>$b}
               keys %{$roi_slice_bounds->{$roi}->{$plane_type}}) {
   
               my $z3d_mm = sprintf("%.2f", $z3d_voxels *
                                    $scaling_voxelnum_2_um/ 1000) ;
   
               if (exists $images->{z3d_to_maximage}->{$z3d_mm}) {
                  my $z_image = $images->{z3d_to_maximage}->{$z3d_mm}->{z} ;
                  my $image_id = $images->{info}->{$z_image}->{id} ;
                  push @td_headers, "<th>$plane_var=$z3d_mm mm</th>" ;
   
                  push @td_images, "<td><a href=\"".
                       $images->{urls}->{$z_image}->{interact}."\">".
                     "<img src=\"".$images->{urls}->{$z_image}->{ish}.
                     "\"></a></td>\n" ;
                  push @td_info, "<td>image ".$images->{info}->{$z_image}->{id}.
                                 '</td>' ;
   
               } else {
                  push @td_headers, '<th></th>' ;
                  push @td_images, "<td></td>\n";
                  push @td_info, '' ;
               }
            }
            push @{$refrows->[(1 + $planeinc)]}, @td_headers ;
            push @{$refrows->[(2 + $planeinc)]}, @td_images ;
            $total_numslices += $numslices ;
            $planenum++ ;
         }
      }
      foreach my $i ( 0 .. $#{$refrows}) {
         $html_header_rows .= "<tr>".join(" ", @{$refrows->[$i]})."</tr>\n" ; }
      $html_header_rows .= "<tr><td colspan=\"$total_numslices\" ".
                            "height=20></td></tr>\n" ;

      print {$html_fh} $html_header_rows ;

      my $hitnum = 0 ;
      open(SPLITRESF, $tcom) ;
      while (my $line = <SPLITRESF>) {
         chomp $line;
         my @t = split(/\t/, $line) ;

# 1. New page if necesary

         $hitnum++ ;
         if ($hitnum % $specs->{report}->{hits_per_page} == 0) {
            print {$html_fh} "</table>$navbar</body></html>\n" ;
            close($html_fh) ;

            $html_page_num++ ;
            open($html_fh, ">".$in->{out_prefix}."_p".$html_page_num.".html") ;

            my $prevpage_url = $in->{out_prefix}."_p".($html_page_num - 1).
                               ".html" ;
            my $nextpage_url = $in->{out_prefix}."_p".($html_page_num + 1).
                               ".html" ;

            $navbar = "<p><h3>Page $html_page_num (of $total_numpages) | ".
                      "<a href=\"$prevpage_url\">Go to previous page</a>";
            if ($html_page_num < $total_numpages) {
               $navbar .= " | <a href=\"$nextpage_url\">next page</a>"; }
            $navbar .= "</h3></p>\n";

            print {$html_fh} $html_header ;
            print {$html_fh} $navbar ;
            print {$html_fh} "<table>\n" ;
            print {$html_fh} $html_header_rows ;
         }


# 1. Info about current hit for header row.
# HERENOW - NOTE: THIS IS RECORD SPECIFIC, currently only for simsearch
         my $curentry ;
         my $hit_header ;
         my $numslices ;
         if ($search_specs->{record} eq 'SIM') {
   
            $curentry->{hit_gene} = $t[$f2i->{hit_gene}] ;
            $curentry->{query_gene} = $t[$f2i->{query_gene}] ;
   
            ($curentry->{hit_experiment_id}) =
               ($t[$f2i->{hit_file}] =~ /_([0-9]+)\.xpz/) ;
   
            ($curentry->{query_experiment_id}) =
               ($t[$f2i->{query_file}] =~ /_([0-9]+)\.xpz/) ;
   
            $curentry->{hit_url} =
               $specs->{download}->{URL}->{interactive_imageseries}.
                  $curentry->{hit_experiment_id} ;
   
            $curentry->{query_url} =
               $specs->{download}->{URL}->{interactive_imageseries}.
                  $curentry->{query_experiment_id} ;
   
            my $query_url = $specs->{download}->{URL}->{interactive_imageseries}.
               $curentry->{query_experiment_id} ;
   
            $curentry->{hit_plane} = $t[$f2i->{hit_plane}] ;
            $curentry->{plane_type} = 'x2yz' ;
            $curentry->{plane_var} = 'x' ;
            if ($curentry->{hit_plane} eq 'sagittal') {
               $curentry->{plane_var} = 'z' ;
               $curentry->{plane_type} = 'z2xy' ; }

            my $roi = $t[$f2i->{ROI}] ;
   
            $numslices = keys %{$roi_slice_bounds->{$roi}->{
                                $curentry->{plane_type}}} ;

            $hit_header = "$hitnum. ".
               "<b>".$curentry->{hit_gene}." </b> (<a href=\"".
               $curentry->{hit_url}."\">".$curentry->{hit_experiment_id}.", ".
               $t[$f2i->{hit_plane}]."</a>) ".  
               "similarity score=".$t[$f2i->{similarity_score}].
               ", mean expression = ".$t[$f2i->{hit_mean_cube_level}].
               " over ".$t[$f2i->{num_shared_cubes}].
               " ROI voxels shared with query ".
               $curentry->{query_gene}." (<a href=\"".$curentry->{query_url}.
               "\">".$curentry->{query_experiment_id}.", ".
               $t[$f2i->{query_plane}]."</a>) ".  
               "(mean expression ".$t[$f2i->{query_mean_cube_level}].")" ;
   
            if ($curentry->{query_experiment_id} == 
                $curentry->{hit_experiment_id}) {
               $hit_header = "<tr><td colspan=\"$numslices\" class=\"headerrow\">".
                             $hit_header."</td></tr>\n" ;
               $html_header_rows .= $hit_header;
            } else {
               $hit_header = "<tr><td colspan=\"$numslices\">$hit_header</td></tr>\n";
            }

            $curentry->{hit_file} = $t[$f2i->{hit_file}] ;

         } elsif ($search_specs->{record} eq 'P') {

            $curentry->{hit_gene} = $t[$f2i->{gene_name}] ;
            $curentry->{hit_plane} = $t[$f2i->{plane}] ;
            ($curentry->{hit_experiment_id}) =
               ($t[$f2i->{xpz_name}] =~ /_([0-9]+)\.xpz/) ;
            $curentry->{hit_url} =
               $specs->{download}->{URL}->{interactive_imageseries}.
               $curentry->{hit_experiment_id} ;

            $curentry->{plane_type} = 'x2yz' ;
            $curentry->{plane_var} = 'x' ;
            if ($curentry->{hit_plane} eq 'sagittal') {
               $curentry->{plane_var} = 'z' ;
               $curentry->{plane_type} = 'z2xy' ; }

            $numslices = 0 ;
            foreach my $roi (keys %{$search_specs->{roi_basename2roi}->{
                                    $cursplit}}) {
               $numslices += keys %{$roi_slice_bounds->{$roi}->{
                                $curentry->{plane_type}}} ;
            }

            $hit_header = "$hitnum. ".
               "<b>".$curentry->{hit_gene}." </b> (<a href=\"".
               $curentry->{hit_url}."\">".$curentry->{hit_experiment_id}.", ".
               $curentry->{hit_plane}."</a>) " ;

            foreach my $measure (qw/numpoints density intensity energy/) {
               foreach my $score (qw/infocontent infocontent_nl gradient_score/) {
                  $hit_header .= $measure."_".$score."=".
                                 $t[$f2i->{$measure."_".$score}]." " ;
               }
            }
            $hit_header =
               "<tr><td colspan=\"$numslices\" class=\"headerrow\">".
               $hit_header."</td></tr>\n" ;
            $curentry->{hit_file} = $t[$f2i->{xpz_name}] ;
         }
         print {$html_fh} $hit_header ;


         my $xpz_fn = $specs->{allen_xpz_dir}."/".$in->{age}."/".
                      substr($curentry->{hit_file}, 0, 1)."/".
                      substr($curentry->{hit_file}, 0, 3)."/".
                      $curentry->{hit_file} ;

         if (!exists $hit_seen->{$curentry->{hit_file}} && $numimages > 0) {
            $numimages-- ;

            my @td_headers ; my @td_images ;

            foreach my $roi (@rois) {

               my $images = _run_roi3d_to_2d({
                  xpz_fn   => $xpz_fn,
                  roidef   => $roidef,
                  roi_name => $roi,
                  specs    => $specs,
                  age      => $age,
                  oneimage_per_z3d => 1,
               }) ;
   
               foreach my $z3d_voxels (sort {$a<=>$b}
                  keys %{$roi_slice_bounds->{$roi}->{$curentry->{plane_type}}}) {
   
                  my $z3d_mm = sprintf("%.2f", $z3d_voxels *
                                       $scaling_voxelnum_2_um/ 1000) ;
   
                  if (exists $images->{z3d_to_maximage}->{$z3d_mm}) {
                     my $z_image = $images->{z3d_to_maximage}->{$z3d_mm}->{z} ;
                     my $image_id = $images->{info}->{$z_image}->{id} ;
                     push @td_headers, "<th>".
                                       $curentry->{plane_var}."=$z3d_mm mm, ".
                                       $curentry->{hit_gene}.
                                       "</th>" ;
   
                     my $inturl = $images->{urls}->{$z_image}->{interact} ;
                     push @td_images, "<td><a href=\"$inturl\">".
                        "<img src=\"".$images->{urls}->{$z_image}->{ish}.
                        "\"></a></td>\n" ;
   
                  } else {
                     push @td_headers, '<th></th>' ;
                     push @td_images, "<td></td>\n";
                  }
               }
            }

            print {$html_fh} "<tr>".join(' ', @td_headers)."</tr>
<tr>".join(' ', @td_images)."</tr>\n" ;

            $hit_seen->{$curentry->{hit_file}}++ ;
         }
         print {$html_fh} "<tr><td colspan=\"$numslices\" height=5>".
                          "</td></tr>\n";
      }

      print {$html_fh} "</table>$navbar</body></html>\n" ;
      close($html_fh) ;
      close(SPLITRESF) ;
   }

}



=head2 _run_roi3d_to_2d()

   Title:       _run_roi3d_to_2d()
   Function:    returns 2D points corresponding to 3D ROI
   Args:        $_->{xpz_fn}
                $_->{roi}

   COORDINATE SYSTEMS:
      roi def = voxel counts
      - convert to um: using $specs->{view_volume}->{$in->{age}}->{spacing} ;
      - convert to mm: / 1000.

=cut

sub _run_roi3d_to_2d {

   my $in = shift ;
   my $specs ;
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   my $xpz ;
   if (exists $in->{refatlas_fl}) {
      my $refatlas_fn ;
      
      if ($in->{age} eq 'P56') {
         $refatlas_fn = $specs->{allen_adult_atlas_dir}.
            '/Spaces/P56/'.$in->{refatlas_plane}.'_images.xml' ;
      } else {
         $refatlas_fn = $specs->{allen_devel_atlas_dir}.
            '/Spaces/'.$in->{age}.'/'.$in->{refatlas_plane}.'_images.xml' ;
      }

      $xpz = {};
      alnmnr::io::parse_image_series_xml({fn => $refatlas_fn, data => $xpz});

   } else {
      $xpz = alnmnr::io::aba_parse_xpz({fn => $in->{xpz_fn}}) ;
   }
   my $roidef = $in->{roidef};

# um = view-volume * $scaling
   my $scaling = $specs->{view_volume}->{$in->{age}}->{spacing} ;

   my $max = 20 ;
   my $bounds2d = {}; #->{z}->{min|max x|y}
   my $image_urls = {} ; my $image_info = {} ;
   foreach my $xyz (keys %{$roidef->{point2roi}}) {

      if (exists $in->{roi} &&
          !exists $roidef->{point2roi}->{$xyz}->{$in->{roi_name}}) { next; }

      my ($x, $y, $z) = split(',', $xyz) ;
      my $x_3d = [ $x * $scaling,
                   $y * $scaling,
                   $z * $scaling ] ;
#      print STDERR "PLANE: ".$xpz->{meta}->{plane}."\n";
#      print STDERR "   CONVERTING (".join(", ", @{$x_3d}).") to 2D\n";
      my $x_2d = alnmnr::core::convert_3d_to_2d_coords({
         xpz    => $xpz,
         x_3d   => $x_3d,
         specs  => $specs,
      }) ;

      if ($xpz->{meta}->{plane} eq 'coronal') {
         $image_urls->{$x_2d->{z}}->{z3d} = sprintf("%.2f", ($x_3d->[0]/1000)) ;
      } else {
         $image_urls->{$x_2d->{z}}->{z3d} = sprintf("%.2f", ($x_3d->[2]/1000)) ;
      }

      my $numkeys = keys %{$bounds2d} ;

      if (!exists $bounds2d->{$x_2d->{z}}) {
         $bounds2d->{$x_2d->{z}}->{min_x} = $x_2d->{x}->[0] ;
         $bounds2d->{$x_2d->{z}}->{max_x} = $x_2d->{x}->[0] ;
         $bounds2d->{$x_2d->{z}}->{min_y} = $x_2d->{x}->[1] ;
         $bounds2d->{$x_2d->{z}}->{max_y} = $x_2d->{x}->[1] ;

      } else {
         if ($x_2d->{x}->[0] < $bounds2d->{$x_2d->{z}}->{min_x}) {
            $bounds2d->{$x_2d->{z}}->{min_x} = $x_2d->{x}->[0] ; }
         if ($x_2d->{x}->[0] > $bounds2d->{$x_2d->{z}}->{max_x}) {
            $bounds2d->{$x_2d->{z}}->{max_x} = $x_2d->{x}->[0] ; }

         if ($x_2d->{x}->[1] < $bounds2d->{$x_2d->{z}}->{min_y}) {
            $bounds2d->{$x_2d->{z}}->{min_y} = $x_2d->{x}->[1] ; }
         if ($x_2d->{x}->[1] > $bounds2d->{$x_2d->{z}}->{max_y}) {
            $bounds2d->{$x_2d->{z}}->{max_y} = $x_2d->{x}->[1] ; }
      }
   }


   my $z3d_to_maximage = {};
   foreach my $z (sort {$a <=> $b} keys %{$bounds2d}) {
      my $z3d = $image_urls->{$z}->{z3d} ;

      my $cur_width = $bounds2d->{$z}->{max_x} - $bounds2d->{$z}->{min_x} + 1;
      my $cur_height = $bounds2d->{$z}->{max_y} - $bounds2d->{$z}->{min_y} + 1;

# Determine # zoom levels ("tiers") available for this image, store in $xpz
      if (!exists $xpz->{meta}->{images}->{$z}->{numtiers}) {
         my $info_url = $specs->{download}->{URL}->{image_info}.
                        $xpz->{meta}->{images}->{$z}->{imagepath};
# URL: http://mouse.brain-map.org/aba/api/image/info?path=/external/aibssan/production4/Pcp4_Baylor_8338/zoomify/primary/null/Pcp4_66_null_A.aff

         my ($t_fh, $t_fn) = tempfile() ;
         my $wget_com = "wget --quiet \"$info_url\" -O $t_fn" ;
#         print STDERR "\n$info_url\n";
         my $wget_status = system($wget_com) ;
         open(INFOF, $t_fn); my $info = <INFOF>;
         chomp $info; close(INFOF) ;
         unlink $t_fn ;
# string: <IMAGE_PROPERTIES WIDTH="6160" HEIGHT="4520" NUMTILES="619" NUMTIERS="6" NUMIMAGES="1" VERSION="1.8" TILESIZE="256"/>

         my ($numtiers) = ($info =~ /NUMTIERS="([0-9]+)"/) ;
         $xpz->{meta}->{images}->{$z}->{numtiers} = $numtiers ;
      }

# Image zoom logic
      my $zoomout ;
      if ($xpz->{meta}->{images}->{$z}->{numtiers} >= 6) {
         $zoomout = 6 - $specs->{report}->{image_zoom} ;
      } else {
         $zoomout = $xpz->{meta}->{images}->{$z}->{numtiers} -
                    $specs->{report}->{image_zoom} ;
      }

# Note: Older image series (id <= 2478 - Rfn160 (P56 sagittal)) were produced in Baylor and have a weird resolution - can never match the newer images exactly - either ~0.8 or 1.6x size...

      my $url_zoom_scale = 2**$zoomout ;
      my $url_zoom_scale_inv = 1 / $url_zoom_scale ;
      my $url_zoom_level =
            ($xpz->{meta}->{images}->{$z}->{numtiers} - 1) - $zoomout ;

# calc top/left/width/height to request.
# NOTE! if devel - need to consider image-top/image-left (mult images/slide)
# HERENOW 120208_1734 

#  every zoom step is a 2x reduction.

      my $rawwidth  = $bounds2d->{$z}->{max_x} - $bounds2d->{$z}->{min_x} + 1;
      my $rawheight = $bounds2d->{$z}->{max_y} - $bounds2d->{$z}->{min_y} + 1;

      my $margin_x = $rawwidth * $specs->{report}->{image_margins_perc} ;
      my $margin_y = $rawheight * $specs->{report}->{image_margins_perc} ;

      my $interact_top   = POSIX::floor($bounds2d->{$z}->{min_y} - $margin_x) ;
      my $img_top        = $interact_top + $xpz->{meta}->{images}->{$z}->{top} ;

      my $interact_left  = POSIX::floor($bounds2d->{$z}->{min_x} - $margin_x) ;
      my $img_left       = $interact_left + $xpz->{meta}->{images}->{$z}->{left} ;

#      print STDERR "IMAGESERIES: ".$xpz->{meta}->{imageseries_id}."\n";
#      print STDERR " NAME: ".$xpz->{meta}->{gene_name}."\n";
#      print STDERR " image: ".$xpz->{meta}->{images}->{$z}->{id}."\n";
#      print STDERR " IMAGE.XML TOP: ".$xpz->{meta}->{images}->{$z}->{top}."\n";
#      print STDERR " IMAGE.XML LEFT: ".$xpz->{meta}->{images}->{$z}->{left}."\n";

      my $width = POSIX::floor($url_zoom_scale_inv * 
                               ($rawwidth + 2 * $margin_x)) ;
      my $height= POSIX::floor($url_zoom_scale_inv * 
                               ($rawheight + 2 * $margin_y)) ;

      $image_urls->{$z}->{ish}  = $specs->{download}->{URL}->{image_part1}.
                $url_zoom_level.
                $specs->{download}->{URL}->{image_part2}.
                $xpz->{meta}->{images}->{$z}->{imagepath}.
                "&top=$img_top&left=$img_left&width=$width&height=$height";

      if (!exists $in->{refatlas_fl}) {
         $image_urls->{$z}->{expr} = $specs->{download}->{URL}->{image_part1}.
                $url_zoom_level.
                $specs->{download}->{URL}->{image_part2}.
                $xpz->{meta}->{images}->{$z}->{"expression-imagepath"}.
                "&top=$img_top&left=$img_left&width=$width&height=$height";
      }

# NOTE interactive coordinates for multi-image differs from the ISH/EXPR images
#      - it is only relative to the individual image
# eg, developmental series, some of the adult atlas slides
      $image_urls->{$z}->{interact} =
                $specs->{download}->{URL}->{interactive_image_part1}.
                $xpz->{meta}->{images}->{$z}->{id}.
                $specs->{download}->{URL}->{interactive_image_part2}.
               'x='.int($interact_left + ($rawwidth / 2)).
               '&y='.int($interact_top + ($rawheight / 2)).
               '&z=50' ;

#      print STDERR "NOTE!! image z=$z. URL: ".$image_urls->{$z}->{ish} ;

#NOTE 120204_1228 - interactive URL z should be proportional to image zoom/width
# but for now, stick with 50, default for Brain Explorer click-through...

      my $area = $width * $height ;
      if (!exists $z3d_to_maximage->{$z3d} ||
          $area > $z3d_to_maximage->{$z3d}->{area}) {
         $z3d_to_maximage->{$z3d}->{area} = $area ;
         $z3d_to_maximage->{$z3d}->{z} = $z ;
         $image_info->{$z}->{id} = $xpz->{meta}->{images}->{$z}->{id} ;
      }
   }

   if (exists $in->{oneimage_per_z3d} && $in->{oneimage_per_z3d} == 1) {
      my @allz = keys %{$bounds2d} ;
      my $keepz = {};
      map {$keepz->{$z3d_to_maximage->{$_}->{z}}++;} (keys %{$z3d_to_maximage});
      my @deletez ;
      foreach my $z (keys %{$bounds2d}) {
         if (!exists $keepz->{$z}) {push @deletez, $z;} }
      map {delete $image_urls->{$_};} @deletez ;
   }

   return {
      z3d_to_maximage   => $z3d_to_maximage,
      info              => $image_info,
      urls              => $image_urls,
   } ;

}


=head2 _run_xpr2pdb()

   Title:       _run_xpr2pdb()
   Function:    converts XPR file to PDB format
   Args:        $_->{xpr_fn} - input xpr filename
                $_->{pdb_fn} - output pdb filename

=cut

sub _run_xpr2pdb {

   my $in = shift ;

   alnmnr::io::aba_parse_xpr({fn        => $in->{xpr_fn},
                              out_pdb   => $in->{pdb_fn}}) ;

}


=head2 _run_xpz2pdb()

   Title:       _run_xpz2pdb()
   Function:    converts XPZ file to PDB format
   Args:        $_->{xpz_fn} - input xpz filename
                $_->{pdb_fn} - output pdb filename

=cut

sub _run_xpz2pdb {

   my $in = shift ;

   alnmnr::io::aba_parse_xpz({fn        => $in->{xpz_fn},
                              out_pdb   => $in->{pdb_fn}}) ;

}


=head2 _run_sva2pdb()

   Title:       _run_sva2pdb()
   Function:    converts SVA file to PDB format
   Args:        $_->{sva_fn} - input SVA filename
                $_->{pdb_fn} - output pdb filename

=cut

sub _run_sva2pdb {

   my $in = shift ;

   alnmnr::io::aba_parse_sva({fn        => $in->{sva_fn},
                              out_pdb   => $in->{pdb_fn}}) ;

}


=head2 _run_xprarchive2pdb()

   Title:       _run_xprarchive2pdb()
   Function:    converts all local XPR files to PDB format
   Args:        $_->{out_dir} - directory to place PDB files
                $_->{compress_fl} - will gzip output PDB files

=cut

sub _run_xprarchive2pdb {

   my $in = shift ;
   if (!exists $in->{out_dir}) {die "out_dir not specified";}

   my $compress_fl = 0 ;
   if (exists $in->{compress_fl} && $in->{compress_fl} == 1) {
      $compress_fl = 1 ; }

   my $out_dir = $in->{out_dir} ;

   my $specs = alnmnr::getspecs() ;
   my $xpr_list = [] ;
   @{$xpr_list} = glob($specs->{allen_xpr_dir}.'/*/*xpr') ;

   print STDERR "now on:  0";
   foreach my $j (0 .. $#{$xpr_list}) {
      print STDERR "\b"x(length($j)).($j + 1) ;
      my $xpr_fn = $xpr_list->[$j] ;
      if (!-s $xpr_fn) {
         print STDERR "\tWARNING: $xpr_fn not found\nnow on:  ".($j + 1)."\n";
         next;
      }
      my $basename = basename($xpr_fn) ;
      $basename =~ s/\.xpr$// ;

      my $out_pdb_dir = $out_dir.'/'.substr($basename,0,1).'/' ;
      if (!-s $out_pdb_dir) { mkpath($out_pdb_dir);}

      my $out_pdb_fn = $out_pdb_dir.'/'.$basename.'.pdb' ;
      if ($compress_fl) { $out_pdb_fn .= '.gz' ; }
      if (-s $out_pdb_fn) {
         print STDERR "\tSKIPPING: $xpr_fn already converted\nnow on:  ".($j + 1)."\n";
         next;
      }

      my $tcom = $0." -mode xpr2pdb -xpr_fn $xpr_fn -pdb_fn $out_pdb_fn";
      system("perl $tcom") ;
      if (!-s $out_pdb_fn) {
         print STDERR "\tWARNING: $xpr_fn not converted properly\nnow on:  ".
                      ($j + 1)."\n";
         next;
      }

   }
   print STDERR "\n" ;
}


=head2 _run_xpr2txt()

   Title:       _run_xpr2txt()
   Function:    converts XPR file to PDB format
   Args:        $_->{xpr_fn} - input xpr filename
                $_->{out_fn} - output pdb filename

=cut

sub _run_xpr2txt {

   my $in = shift ;

   alnmnr::io::aba_parse_xpr({fn        => $in->{xpr_fn},
                              out_xyz   => $in->{txt_fn}}) ;

}


=head2 _run_xpz2txt()

   Title:       _run_xpz2txt()
   Function:    converts XPZ file to PDB format
   Args:        $_->{xpz_fn} - input xpz filename
                $_->{out_fn} - output pdb filename

=cut

sub _run_xpz2txt {

   my $in = shift ;

   alnmnr::io::aba_parse_xpz({fn        => $in->{xpz_fn},
                              out_xyz   => $in->{txt_fn}}) ;

}


=head2 _run_sva2txt()

   Title:       _run_sva2txt()
   Function:    converts SVA file to PDB format
   Args:        $_->{sva_fn} - input SVA filename
                $_->{out_fn} - output pdb filename

=cut

sub _run_sva2txt {

   my $in = shift ;

   alnmnr::io::aba_parse_sva({fn        => $in->{sva_fn},
                              out_xyz   => $in->{txt_fn}}) ;

}


sub _run_fastsearch {
#Purpose: performs fast ABA boolean searches using precomputed 
#         search() results for all ABA reference atlas regions.

# Allows 2 types of query:
# 1. enrichment: roi1 vs roi2
# 2. expression level: roi1

   my $in = shift;
   my $age = $in->{age} ;
   $in->{headers} = alnmnr::core::calc_expr({ return_headers_only => 1}) ;

   my $specs ;
   if (exists $in->{specs}) {
      $specs = $in->{specs} ;
   } else {
      $specs = alnmnr::getspecs() ;
   }

   my $f2i ;
   my $raw_headers = $in->{headers}->{E} ;
   foreach my $j ( 0 .. $#{$in->{headers}->{E}}) {
      $f2i->{$in->{headers}->{E}->[$j]} = $j ; }

   my $queries ;
   my $get_rois ;
   open(QUERYSPECS, $in->{query_specs_fn}) ;
   while (my $line = <QUERYSPECS>) {
      chomp $line;
      if ($line =~ /^\#/ || $line =~ /^$/) {next;}
      my ($name, $type, @roi)= split(/\t/, $line) ;
      if ($type eq 'enrichment') {
         $queries->{$name}->{type} = $type ;
         $queries->{$name}->{roi_1} = $roi[0] ;
         $queries->{$name}->{roi_2} = $roi[1] ;
         $get_rois->{$roi[0]}++ ;
         $get_rois->{$roi[1]}++ ;
      } elsif ($type eq 'level') {
         $queries->{$name}->{type} = $type ;
         $queries->{$name}->{roi} = $roi[0] ;
         $get_rois->{$roi[0]}++ ;
      }
   }
   close(QUERYSPECS) ;

   my $roi_data ;

   foreach my $slice (qw/sagittal coronal/) {
      $roi_data->{xpzlist}->{$slice} = [] ;
      $roi_data->{xpz2ind}->{$slice} = {};
      $roi_data->{roi_data}->{$slice} = {};
   }
   my $roidef_info = {};

# iterate over all roi expressions and come up with an equation to
#  describe the final value.

   my $roi_eqn ;
   my $roidef ;
   foreach my $roi ( keys %{$get_rois}) {

      $roidef->{$roi}->{points} = alnmnr::io::get_brainstrx_voxels({
         specs => $specs,
         brainstrx => $roi,
         onpoints_only => 1,
         age => $age,
      }) ;
      $roidef->{$roi}->{num_roidef_points} = keys %{$roidef->{$roi}->{points}};

      if ($roi =~ /\,/) {
         my @subroi = split(/\,/, $roi) ;
         my $subroi_details ;

         my $roi_kids ;
         foreach my $j ( 0 .. $#subroi) {
            $subroi[$j] =~ s/ //g ;
            my ($name, $status) ;
            if ($subroi[$j] =~ /\=/) {
               ($name, $status) = ($subroi[$j] =~ /(.+)\=(.+)/) ;
            } else {
               $name = $subroi[$j] ;
               $status = 'on' ;
            }
            $subroi_details->{$name}->{status} = $status ;
            $roi_kids->{$name} = get_brainstrx_info({brainstrx => $name,
                                                     age => $age}) ;
         }

# figure out relation and equation for scoring expression
         my @all_rois = sort keys %{$subroi_details} ;
         my $roi_kid2parent ;
         foreach my $roi_1 ( @all_rois) { #candidate parent
            foreach my $roi_2 ( @all_rois) { #candidate kid
               if ($roi_1 eq $roi_2) {next;}
               if ($subroi_details->{$roi_2}->{status} eq 'disappear'){
                  next;}

               if (val_is_in_arr({array => $roi_kids->{$roi_1}->{kids},
                      value => $roi_kids->{$roi_2}->{brainstrx_id}})) {
                  $roi_kid2parent->{$roi_2} = $roi_1 ;
                  if ($subroi_details->{$roi_2}->{status} eq 'off') {
                     $subroi_details->{$roi_1}->{minus}->{$roi_2}++; } 
                  $subroi_details->{$roi_2}->{status} = 'disappear' ;
               }
            }
         }

# all offs by now should either be subtacted from their parents if there is one
#  or disappeared if no parent is found; just warn about those
         foreach my $subroi_name (keys %{$subroi_details}) {
            if ($subroi_details->{$subroi_name}->{status} eq 'off' &&
                !exists $roi_kid2parent->{$subroi_name}) {
               $subroi_details->{$subroi_name}->{status}  = 'disappear' ;
               print STDERR "WARNING: parent structure of $subroi_name (specified as off) was not found in the roi definition - deleting from specs\n" ;
            }
         }

# bring pieces together
         my $finaleqn_parts ;
         foreach my $subroi_name (keys %{$subroi_details}) {
            if ($subroi_details->{$subroi_name}->{status} eq 'disappear'){
               next;}
            my $cur_part = $subroi_name ;
            if (exists $subroi_details->{$subroi_name}->{minus}) {
               foreach my $minus_roi (keys %{$subroi_details->{$subroi_name}->{minus}}) {
                  $cur_part .= " - $minus_roi" ;
               }
            }
            push @{$finaleqn_parts}, $cur_part ;
         }

         $roi_eqn->{$roi} = join(" + ", @{$finaleqn_parts}) ;

      } else { # if nothing specifeid about the roi, just add it
         $roi_eqn->{$roi} = $roi ;
      }
#      print STDERR "EQUATION for $roi: ".$roi_eqn->{$roi}."\n" ;
   }

# now iterate over the input, and for each gene read, output the
#   roi info for each query type.

   my ($last_xpzname, $last_slices) ;
   my $cur_xpz_info ;
   my $entry ;

# default is to have individual queries go to individual files.
# if out_fn is defined, then everything gets put there.

   my @headers_enrich = (qw/FE xpz slice query_name query_type_enrichment/,
      qw/ROI1 roi1_density roi1_intensity roi1_energy/,
      qw/ROI2 roi2_density roi2_intensity roi2_energy/,
      qw/enrichment_density enrichment_intensity enrichment_energy/) ;
   my @headers_level = (qw/F xpz slice query_name query_type_level/,
                        qw/ROI roi_density roi_intensity roi_energy/);

   my ($out_fh, $out_fn) ;
   if (defined $in->{out_fn}) {
      my $fh; open($fh, ">".$in->{out_fn}) ;
      print {$fh} '#'.join("\t", @headers_enrich)."\n" ;
      print {$fh} '#'.join("\t", @headers_level)."\n" ;
      foreach my $name (keys %{$queries}) {
         $out_fn->{$name} = $in->{out_fn} ;
         $out_fh->{$name} = $fh ;
      }
   } elsif (defined $in->{outfile_prefix}) {
      foreach my $name (keys %{$queries}) {
         $out_fn->{$name} = $in->{outfile_prefix}."_$name.allenminer.out" ;
         open($out_fh->{$name}, ">".$out_fn->{$name}) ;
         print {$out_fh->{$name}} '#'.join("\t", @headers_enrich)."\n" ;
         print {$out_fh->{$name}} '#'.join("\t", @headers_level)."\n" ;
      }
   } else {
      foreach my $name (keys %{$queries}) {
         $out_fn->{$name} = "aba_fast_query.$name.$$.allenminer.out" ;
         open($out_fh->{$name}, ">".$out_fn->{$name}) ;
         print {$out_fh->{$name}} '#'.join("\t", @headers_enrich)."\n" ;
         print {$out_fh->{$name}} '#'.join("\t", @headers_level)."\n" ;
      }
   }

   my $fastsearch_data_fn = $specs->{fastsearch_fn}->{$age} ;
   if (!defined $fastsearch_data_fn || ! -s $fastsearch_data_fn) {
      die "ERROR: can't find FASTSEARCH expression data file: ".
          $fastsearch_data_fn."\n";
   }
   if ($fastsearch_data_fn =~ /gz$/) {
      open(FASTSEARCHF, 'zcat '.$fastsearch_data_fn.' |');
   } else {
      open(FASTSEARCHF, $fastsearch_data_fn.' |');
   }

   while (my $line = <FASTSEARCHF>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      $entry = {};

      my @t = split(/\t/, $line) ;
      if ($t[0] ne 'E') {next;}
      foreach my $j (0 .. $#t) {
         $entry->{$raw_headers->[$j]} = $t[$j] ; }

      if ($entry->{plane} ne 'sagittal' &&
          $entry->{plane} ne 'coronal') { next; }
      my $slice = $entry->{plane} ;


# if new xpz: calc roi enrichment scores, IC metric and specific for last one
# and store in results

      my $lastline_fl = 0 ;
      if (eof(FASTSEARCHF)) {
         foreach my $f (sort keys %{$f2i}) {
            if ($entry->{$f} eq 'NA') {$entry->{$f} = 0;}
            $cur_xpz_info->{$entry->{ROI}}->{$f} = $entry->{$f} ; }
         $lastline_fl = 1 ;
      }

      if (defined $last_xpzname && $entry->{xpz_name} ne $last_xpzname ||
          $lastline_fl == 1) { # one xpz data complete: score query rois

# push onto master xpz list
         push @{$roi_data->{xpzlist}->{$last_slices}}, $last_xpzname ;
         $roi_data->{xpz2ind}->{$last_slices}->{$last_xpzname} =
            $#{$roi_data->{xpzlist}->{$last_slices}} ;

# iterate over roi_eqn and score
         my $curxpz_queryroi ;
         foreach my $query_roi (keys %{$roi_eqn}) {
            foreach my $data_type (qw/density intensity energy/) {
               my $cur_expr = $roi_eqn->{$query_roi} ;
#               print STDERR "$query_roi $data_type: $cur_expr turning to " ;
               $cur_expr =~ s/([a-zA-Z]+)/$cur_xpz_info->{$1}->{$data_type}/g ;
               $curxpz_queryroi->{$query_roi}->{$data_type}= eval($cur_expr) ;
#               print STDERR "   = $cur_expr = ".
#                  $curxpz_queryroi->{$query_roi}->{$data_type}."\n" ;
            }
         }

# iterate over queries and display results from roi
         foreach my $query_name (keys %{$queries}) {
# output format: xpz roi density intensity energy
# if enrichment: followed by second roi info, followed by enrichment info
            my $query = $queries->{$query_name} ;
            my @outvals = ('F', $last_xpzname, $last_slices,
                           $query_name, $query->{type}) ;

            if ($query->{type} eq 'enrichment') {
               $outvals[0] = 'FE' ;
               my $roi_1 = $query->{roi_1} ;
               my $roi_2 = $query->{roi_2} ;
               foreach my $cur_roi ($roi_1, $roi_2) {
                  push @outvals, $cur_roi ;
                  push @outvals,
                     $curxpz_queryroi->{$cur_roi}->{density},
                     $curxpz_queryroi->{$cur_roi}->{intensity},
                     $curxpz_queryroi->{$cur_roi}->{energy};
               }
               my $escore_denom =
                  $roidef->{$roi_1}->{num_roidef_points} /
                  $roidef->{$roi_2}->{num_roidef_points} ;
               foreach my $data_type (qw/density intensity energy/) {
                  if ($curxpz_queryroi->{$roi_1}->{$data_type} == 0 ||
                      $curxpz_queryroi->{$roi_2}->{$data_type} == 0) {
                      push @outvals, 0 ;
                  } else {
                     push @outvals, (($curxpz_queryroi->{$roi_1}->{$data_type} /
                                   $curxpz_queryroi->{$roi_2}->{$data_type})/
                                  $escore_denom) ;
                  }
               }
            } else {
               push @outvals, $query->{roi},
                  $curxpz_queryroi->{$query->{roi}}->{density},
                  $curxpz_queryroi->{$query->{roi}}->{intensity},
                  $curxpz_queryroi->{$query->{roi}}->{energy};
            }

            print {$out_fh->{$query_name}} join("\t", @outvals)."\n" ;
         }
         $cur_xpz_info = {} ;
      }

      $last_xpzname = $entry->{xpz_name} ;
      $last_slices = $entry->{plane} ;
      foreach my $f (sort keys %{$f2i}) {
         if ($entry->{$f} eq 'NA') {$entry->{$f} = 0;}
         $cur_xpz_info->{$entry->{ROI}}->{$f} = $entry->{$f} ; }

   }
   close(FASTSEARCHF) ;

   foreach my $query_name (keys %{$queries}) {
      close($out_fh->{$query_name}) ; }

   return $roi_data ;

}


sub create_fastsearch_data {
# Creates an ROI definition file that includes all defined regions for each age, run in roi-list mode, and compress the output into the curdir()

   my $in = shift ;
   my $specs = $in->{specs} ;
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   foreach my $age (keys %{$specs->{dims}}) {
      my ($roidef_fh, $roidef_fn) =
         tempfile("roi_allregions_$age.XXXXX", SUFFIX => ".roi") ;

      my $brainstrx = alnmnr::io::load_aba_brainstructures({age => $age}) ;
      foreach my $strx (keys %{$brainstrx->{abbreviation2ID}}) {
         print {$roidef_fh} join("\t", $strx, '', 'brainstrx', $strx)."\n"; }
      close($roidef_fh) ;

      my $out_fn = "fastsearch_$age.txt.gz" ;
      my $tcom =
         "allenminer.pl -mode search -age $age -roidef_fn $roidef_fn ".
         "-out_fn $out_fn" ;
      print $tcom."\n" ;
   }
   
}

1 ;
