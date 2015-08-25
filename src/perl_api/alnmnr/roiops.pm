=head1 NAME

alnmnr::roiops.pm - routines to manipulate ROI

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

package alnmnr::roiops;
use strict;
use warnings;
use Cwd;
use File::Path qw/mkpath/ ;
use File::Temp qw/tempfile tempdir/ ;
use POSIX qw/floor ceil/ ;

=head2 peel_roi()

   Title:       peel_roi()
   Function:    Remove the outer n layers of an ROI
   Args:        $_->{roidef} = ROIdef structure: $_->{x,y,z}->{roiname_subname}
                $_->{n} = number of layers to peel off - optional, defaults to 1
                          (positive integer)
   Returns:     $_ = new peeled ROIdef structure

=cut

sub peel_roi {
#purpose remove outer n layers of an ROI.

   my $in = shift ;
   my $points = $in->{points}  ;
   my $n_layers = 1;
   if (exists $in->{n_layers}) {
      $n_layers = $in->{n_layers};}
   my $roidef = $in->{roidef} ;

   my $peeled_roidef ;
   if ($n_layers > 0) {
   $peeled_roidef = $roidef ;
   foreach my $j ( 1 .. $n_layers) {
      $peeled_roidef = get_edge_points({
         return_shaven_roidef => 1,
         roidef => $peeled_roidef }) ;
   }

   } elsif ($n_layers < 0) {

# build n_layer shell ROIdef
      $peeled_roidef = build_roi_shell({
         roidef => $roidef,
         n_layers => $n_layers * -1
      }) ;

# supplement shell ROIdef with the original ROIdef
      foreach my $point (keys %{$roidef->{point2roi}}) {
         map {$peeled_roidef->{point2roi}->{$point}->{$_}++;}
            (keys %{$roidef->{point2roi}->{$point}}) ; }

   }

   return $peeled_roidef ;

}


=head2 subtract_roi()

   Title:       subtract_roi()
   Function:    Subtract roi2 from roi1
   Args:        $_->{roidef} = ROIdef structure: $_->{x,y,z}->{roiname_subname}
                $_->{roi1} = name of first roi
                $_->{roi2} = name of second roi
                $_->{new_name} = name of new roi
   Returns:     $_ = new ROIdef structure

=cut

sub subtract_roi {
#purpose: calculate new_roi = roi1 - roi2

   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $roi1 = $in->{roi1} ;
   my $roi2 = $in->{roi2} ;
   my $new_name = $in->{new_name} ;

   my $new_roidef ;
   foreach my $point (keys %{$roidef->{point2roi}}) {
      if (!exists $roidef->{point2roi}->{$point}->{$roi1} ||
          exists $roidef->{point2roi}->{$point}->{$roi2}) {next;}
      $new_roidef->{point2roi}->{$point}->{$new_name}++ ;
   }

   return $new_roidef ;
}


=head2 rename_roi()

   Title:       rename_roi()
   Function:    Rename ROI
   Args:        $_->{roidef} = ROIdef structure: $_->{x,y,z}->{roiname_subname}
                $_->{old_name} = old ROI name
                $_->{new_name} = new ROI name
   Returns:     $_ = new ROIdef structure

=cut

sub rename_roi {
#purpose: calculate new_roi = roi1 - roi2

   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $old_name = $in->{old_name} ;
   my $new_name = $in->{new_name} ;

   my $new_roidef ;
   foreach my $point (keys %{$roidef->{point2roi}}) {
      if (!exists $roidef->{point2roi}->{$point}->{$old_name}) {next;}
      $new_roidef->{point2roi}->{$point}->{$new_name}++ ;
   }

   return $new_roidef ;
}


=head2 flip_roi()

   Title:       flip_roi()
   Function:    Flips a ROI across the mid-sagittal plane to convert it to
                the equivalent definition in the other hemisphere
   Args:        $_->{points}->{x,y,z} = hash of points in an ROI
                $_->{n} = number of adjacent edge layers to remove
                          (positive integer)
   Returns:     $_->{x,y,z} = hash of points in peeled ROI

=cut

sub flip_roi {
#NOTE: if age=adult, ASSUMES original 100sva format ROIdef (133x81x115)
#HERENOW111206_1031  double check that proper dims are used for P56

   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $age = $in->{age} ;
   my $specs = $in->{specs};
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   my $roi_hash = {} ;
   if (exists $in->{roi_list}) {
      map {$roi_hash->{$_}++;} @{$in->{roi_list}} ; }

   print STDERR "AGE = $age\n";
   my $z_midpoint = $specs->{dims}->{$age}->{z} / 2;

   my $flipped_roi ;
   foreach my $xyz (keys %{$roidef->{point2roi}}) {
      my ($x,$y,$z) = split(',',$xyz) ;
      my $newz = $z - 2 * ($z - $z_midpoint) ;
      print STDERR "converting oldz=$z to newz=$newz\n";
      map {$flipped_roi->{point2roi}->{$x.','.$y.','.$newz}->{"flip_".$_}++; }
         (keys %{$roidef->{point2roi}->{$xyz}}) ;
   }

   return $flipped_roi ;

}


=head2 separate_adjacent_roi()

   Title:       separate_adjacent_roi()
   Function:    Removes n ROI edge layers between adjacnet ROIs 
                - effectively imposes a minimum separation.
   Args:        $_->{points}->{x,y,z} = hash of points in an ROI
                $_->{n} = number of adjacent edge layers to remove
                          (positive integer)
   Returns:     $_->{x,y,z} = hash of points in peeled ROI

=cut

sub separate_adjacent_roi {

   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $n_layers = $in->{n_layers} ;

   my $shave_points ;
   my $shell_roidef = build_roi_shell({ roidef => $roidef,
                                        n_layers => $n_layers }) ;
   foreach my $point (keys %{$shell_roidef->{point2roi}}) {
      if (!exists $roidef->{point2roi}->{$point}) {next;}

      my @curshells = sort keys %{$shell_roidef->{point2roi}->{$point}} ;
      foreach my $roi1 (@curshells) {
         my @currois = sort keys %{$roidef->{point2roi}->{$point}} ;
         foreach my $roi2 (@currois) {
            if ($roi1 eq $roi2) {next;}
            $shave_points->{$point}->{$roi2}++ ;
         }
      }
   }

   my $shaven_roidef = {};
   foreach my $point (keys %{$roidef->{point2roi}}) {
      foreach my $roi (keys %{$roidef->{point2roi}->{$point}}) {
         if (exists $shave_points->{$point} &&
             exists $shave_points->{$point}->{$roi}) {next;}
         $shaven_roidef->{$point}->{$roi}++ ;
      }
   }

   return {point2roi => $shaven_roidef} ;

}


sub get_edge_points {
# purpose: given ROIdef sends back edge points

   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $point2roi = $roidef->{point2roi} ;

   my $return_shaven_roidef_fl = 0 ;
   if (exists $in->{return_shaven_roidef} &&
       $in->{return_shaven_roidef} == 1) {$return_shaven_roidef_fl = 1;}

   my $roi2edge_points = {};
   my $peeled_roidef ;

   foreach my $point (keys %{$point2roi}) {
      my ($x,$y,$z) = split(/\,/,$point) ;
      foreach my $roi (keys %{$point2roi->{$point}}) {
         my $on_edge = 0 ;
# check existence of 26 neighbors
         foreach my $t_x ( ($x - 1) .. ($x + 1)) {
            if ($t_x < 0) {next;}
            foreach my $t_y ( ($y - 1) .. ($y + 1)) {
               if ($t_y < 0) {next;}
               foreach my $t_z ( ($z - 1) .. ($z + 1)) {
                  if ($t_z < 0) {next;}
                  if ($x == $t_x && $y == $t_y && $z == $t_z) {next;}
                  if (!exists $point2roi->{$t_x.','.$t_y.','.$t_z} ||
                      !exists $point2roi->{$t_x.','.$t_y.','.$t_z}->{$roi}) {
                        $roi2edge_points->{$roi}->{$point}++ ;
                        $on_edge++ ;
                        last;
                  }
               }
               if ($on_edge) {last;}
            }
            if ($on_edge) {last;}
         }
         if ($return_shaven_roidef_fl && $on_edge == 0) {
            $peeled_roidef->{$point}->{$roi}++ ; }
      }
   }

   if ($return_shaven_roidef_fl) {
      return {point2roi => $peeled_roidef} ;
   } else {
      return $roi2edge_points ;
   }

}


sub build_roi_shell {
# given ROIdef, return ROIdef of shell points?

   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $n_layers = $in->{n_layers} ;

   my $edge_points = get_edge_points({roidef => $roidef}) ;
   my $shell_roidef ;
   foreach my $roi (keys %{$edge_points}) {
      foreach my $point (keys %{$edge_points->{$roi}}) {
         my ($x,$y,$z) = split(/\,/,$point) ;
         foreach my $t_x ( ($x - $n_layers) .. ($x + $n_layers)) {
            if ($t_x < 0) {next;}
            foreach my $t_y ( ($y - $n_layers) .. ($y + $n_layers)) {
               if ($t_y < 0) {next;}
               foreach my $t_z ( ($z - $n_layers) .. ($z + $n_layers)) {
                  if ($t_z < 0) {next;}
                  if ($x == $t_x && $y == $t_y && $z == $t_z) {next;}
                  if (exists $roidef->{$t_x.','.$t_y.','.$t_z} &&
                      exists $roidef->{$t_x.','.$t_y.','.$t_z}->{$roi}) {next;}
                  $shell_roidef->{$t_x.','.$t_y.','.$t_z}->{$roi}++ ; }}}}}

   return {point2roi => $shell_roidef} ;

}

# given an roidef (boxbounds), enumerates all voxels inside
sub calc_roi_boxbounds_pointlist {

   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $specs ;
   if (exists $in->{specs}) {
      $specs = $in->{specs}
   } else {
      $specs = alnmnr::getspecs() ;
   }
   my $bindefs = $specs->{bindefs} ;

   my $search_bounds ;
   foreach my $bound_axes (qw/x y z/) {
    foreach my $side (qw/min max/) {
      if ($roidef->{$bound_axes.$side}->{type} eq 'constant') {
         $search_bounds->{$bound_axes.$side} =
            $roidef->{$bound_axes.$side}->{value} ;
      } else {
# HERENOW111123_1840 CHANGE THIS TO use $specs->{dims} (or view volume??)
         $search_bounds->{$bound_axes.$side} =
            $bindefs->{$bound_axes."_coord_$side"} ;
      }
    }
   }

   my $roidef_points ;
   my @t = ($search_bounds->{xmin} .. $search_bounds->{xmax}) ;
   foreach my $cur_x ($search_bounds->{xmin} .. $search_bounds->{xmax}) {
    foreach my $cur_y ($search_bounds->{ymin} .. $search_bounds->{ymax}) {
     foreach my $cur_z ($search_bounds->{zmin} .. $search_bounds->{zmax}) {
         my $cur_point = {
            x => $cur_x,
            y => $cur_y,
            z => $cur_z,
         } ;
         my $outofbounds = 0 ;
         foreach my $bound_axes (qw/x y z/) {
            if ($roidef->{$bound_axes.'min'}->{type} ne 'constant') {
               my $bound_val = $roidef->{$bound_axes.'min'}->{value} ;
               $bound_val =~ s/([xyz])/$cur_point->{$1}/g ;
               $bound_val = eval($bound_val) ;
               if ($cur_point->{$bound_axes} < $bound_val) {
                  $outofbounds = 1;
                  last;
               }
            }

            if ($roidef->{$bound_axes.'max'}->{type} ne 'constant') {
               my $bound_val = $roidef->{$bound_axes.'max'}->{value} ;
               $bound_val =~ s/([xyz])/$cur_point->{$1}/g ;
               $bound_val = eval($bound_val) ;
               if ($cur_point->{$bound_axes} > $bound_val) {
                  $outofbounds = 1;
                  last;
               }
            }
         }
         if (!$outofbounds) {
            $roidef_points->{$cur_x.",".$cur_y.",".$cur_z}++ ; }
     }
    }
   }

   return $roidef_points ;

}


sub calc_roi_slicebounds {
   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $out    = $in->{roi_slice_bounds} ;

   foreach my $xyz (keys %{$roidef->{point2roi}}) {
      foreach my $roi (keys %{$roidef->{point2roi}->{$xyz}}) {
         my ($x, $y, $z) = split(',', $xyz) ;
         if (!exists $out->{$roi} || !exists $out->{$roi}->{x2yz}->{$x}) {
            $out->{$roi}->{x2yz}->{$x}->{y_min} = $y ;
            $out->{$roi}->{x2yz}->{$x}->{y_max} = $y ;
            $out->{$roi}->{x2yz}->{$x}->{z_min} = $z ;
            $out->{$roi}->{x2yz}->{$x}->{z_max} = $z ;
         } else {
            if ($y < $out->{$roi}->{x2yz}->{$x}->{y_min}) {
               $out->{$roi}->{x2yz}->{$x}->{y_min} = $y; }
            if ($y > $out->{$roi}->{x2yz}->{$x}->{y_max}) {
               $out->{$roi}->{x2yz}->{$x}->{y_max} = $y; }
            if ($z < $out->{$roi}->{x2yz}->{$x}->{z_min}) {
               $out->{$roi}->{x2yz}->{$x}->{z_min} = $z; }
            if ($z > $out->{$roi}->{x2yz}->{$x}->{z_max}) {
               $out->{$roi}->{x2yz}->{$x}->{z_max} = $z; }
         }

         if (!exists $out->{$roi}->{z2xy}->{$z}) {
            $out->{$roi}->{z2xy}->{$z}->{x_min} = $x ;
            $out->{$roi}->{z2xy}->{$z}->{x_max} = $x ;
            $out->{$roi}->{z2xy}->{$z}->{y_min} = $y ;
            $out->{$roi}->{z2xy}->{$z}->{y_max} = $y ;
         } else {
            if ($x < $out->{$roi}->{z2xy}->{$z}->{x_min}) {
               $out->{$roi}->{z2xy}->{$z}->{x_min} = $x; }
            if ($x > $out->{$roi}->{z2xy}->{$z}->{x_max}) {
               $out->{$roi}->{z2xy}->{$z}->{x_max} = $x; }
            if ($y < $out->{$roi}->{z2xy}->{$z}->{y_min}) {
               $out->{$roi}->{z2xy}->{$z}->{y_min} = $y; }
            if ($y > $out->{$roi}->{z2xy}->{$z}->{y_max}) {
               $out->{$roi}->{z2xy}->{$z}->{y_max} = $y; }

         }

      }
   }

#   foreach my $roi (keys %{$out}) {
#      foreach my $x (keys %{$out->{$roi}->{x2yz}}) {
#         print STDERR "ROI = $roi\n";
#         print STDERR "   ROI slice $x:\n";
#         print STDERR "      y_min:".$out->{$roi}->{x2yz}->{$x}->{y_min}."\n";
#         print STDERR "      y_max:".$out->{$roi}->{x2yz}->{$x}->{y_max}."\n";
#         print STDERR "      z_min:".$out->{$roi}->{x2yz}->{$x}->{z_min}."\n";
#         print STDERR "      z_max:".$out->{$roi}->{x2yz}->{$x}->{z_max}."\n";
#      }
#   }

   return ;
}


sub calc_roi_partitions {

   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $axes_key = {
      RC => 'x',
      AP => 'x',
      DV => 'y',
      ML => 'z',
   };

   my $axes_num = {
      'x' => 0,
      'y' => 1,
      'z' => 2,
   };

   if (!exists $in->{axes}||
       !exists $axes_key->{$in->{axes}}) {
      die "partition axes not recognized" ; }
   if (!exists $in->{numbins}) {die "numbins not specified";}

   my $cut_axes = $axes_key->{$in->{axes}} ;

# 1. iterate over points: figure out x,y,z box bounds
   my $roi_bounds;
   foreach my $xyz (keys %{$roidef->{point2roi}}) {
      my @t = split(/\,/, $xyz) ;
      my $point = { x => $t[0], y => $t[1], z => $t[2] } ;
      my @rangecheck = @t ;
      $rangecheck[$axes_num->{$cut_axes}] = '*';
      my $rangecheck = join(',', @rangecheck);
      foreach my $roi (keys %{$roidef->{point2roi}->{$xyz}}) {
         if (exists $in->{fitted_cuts} && $in->{fitted_cuts} == 1) {
            if (!exists $roi_bounds->{$roi}->{$rangecheck}->{min} ||
                $point->{$cut_axes} < $roi_bounds->{$roi}->{$rangecheck}->{min}) {
               $roi_bounds->{$roi}->{$rangecheck}->{min} = $point->{$cut_axes} ; }

            if (!exists $roi_bounds->{$roi}->{$rangecheck}->{max} ||
                $point->{$cut_axes} > $roi_bounds->{$roi}->{$rangecheck}->{max}) {
               $roi_bounds->{$roi}->{$rangecheck}->{max} = $point->{$cut_axes} ; }

         } else {
            if (!exists $roi_bounds->{$roi}->{min} ||
                $point->{$cut_axes} < $roi_bounds->{$roi}->{min}) {
               $roi_bounds->{$roi}->{min} = $point->{$cut_axes} ; }

            if (!exists $roi_bounds->{$roi}->{max} ||
                $point->{$cut_axes} > $roi_bounds->{$roi}->{max}) {
               $roi_bounds->{$roi}->{max} = $point->{$cut_axes} ; }

         }
      }
   }

# 2. divvy up specified axes into specified number of bins.
   if (exists $in->{fitted_cuts} && $in->{fitted_cuts} == 1) {
      foreach my $roi (keys %{$roi_bounds}) {
         foreach my $rangecheck (keys %{$roi_bounds->{$roi}}) {
            $roi_bounds->{$roi}->{$rangecheck}->{range} =
               $roi_bounds->{$roi}->{$rangecheck}->{max} -
               $roi_bounds->{$roi}->{$rangecheck}->{min} ;
            if ($roi_bounds->{$roi}->{$rangecheck}->{range} == 0) {
               $roi_bounds->{$roi}->{$rangecheck}->{range} = 1; }
         }
      }
   } else {
      foreach my $roi (keys %{$roi_bounds}) {
         $roi_bounds->{$roi}->{range} = $roi_bounds->{$roi}->{max} -
                                       $roi_bounds->{$roi}->{min} ;
      }
   }

# 3. iterate over points: assign each point to its bin
   my $newroidef ;
   foreach my $xyz (keys %{$roidef->{point2roi}}) {
      my @t = split(/\,/, $xyz) ;
      my $point = { x => $t[0], y => $t[1], z => $t[2] } ;
      my @rangecheck = @t ;
      $rangecheck[$axes_num->{$cut_axes}] = '*';
      my $rangecheck = join(',', @rangecheck);
      foreach my $roi (keys %{$roidef->{point2roi}->{$xyz}}) {
         my $roibin_no ;
         if (exists $in->{fitted_cuts} && $in->{fitted_cuts} == 1) {
            $roibin_no =floor((($point->{$cut_axes} -
                               $roi_bounds->{$roi}->{$rangecheck}->{min})/
                               $roi_bounds->{$roi}->{$rangecheck}->{range})
                               * $in->{numbins}) ;
         } else {
            $roibin_no = floor((($point->{$cut_axes} -
                                 $roi_bounds->{$roi}->{min}) /
                                $roi_bounds->{$roi}->{range}) * $in->{numbins});
         }

         if ($roibin_no == $in->{numbins}) {
            $roibin_no  = $roibin_no - 1;}
         my $newroi = '';
         if (exists $in->{roi_basename} &&
             defined $in->{roi_basename}) {
            $newroi = $in->{roi_basename}."_";
         }
         $newroi .= $roi."_part$roibin_no" ;
#         print STDERR "newroi is $newroi e\n" ;
         $newroidef->{point2roi}->{$xyz}->{$newroi}++ ;
      }
   }

# 4. return new point-specified roi specs 
   return $newroidef ;

}


sub convert_roidef_100_to_200 {
#purpose: convert ROIdef from 100 um to 200 um grid for SVA parsing
# ROI 100 grid (default ABA XPR format):   133 x 81 x 115
# ROI 200 grid (SVA smoothed energy files): 67 x 41 x  58

   my $in = shift ;
   my $roidef_100 = $in->{roidef} ;

   my $roidef_200 ;
   foreach my $point (keys %{$roidef_100->{point2roi}}) {
      my ($x,$y,$z) = split(/\,/, $point) ;
      my $new_x = POSIX::floor($x / 2) ;
      my $new_y = POSIX::floor($y / 2) ;
      my $new_z = POSIX::floor($z / 2) ;

      map {$roidef_200->{point2roi}->{$new_x.','.$new_y.','.$new_z}->{$_}++;}
         (keys %{$roidef_100->{point2roi}->{$point}}) ;
   }

# Count up number of ROI points
   foreach my $xyz (keys %{$roidef_200->{point2roi}}) {
      map {$roidef_200->{num_roidef_points}->{$_}++ ;}
         (keys %{$roidef_200->{point2roi}->{$xyz}}) ; }

# Also copy over meta info for completeness sake
   $roidef_200->{roi_meta} = {} ;
   %{$roidef_200->{roi_meta}} = %{$roidef_100->{roi_meta}} ;
   return $roidef_200 ;

}


sub create_lrcheck_roi {
# given ROIDEF input, for each ROI, union a flip about sagittal center; split into L/R halves
# and append new L/R ROIs to the input ->{roidef}

# By default splits all ROIs, optionally uses ->{roi_list}.

# if a -lr-check is given to calc_expr or calc_expr_sim, independent operations for L/R sides are done, and a diff computed.

   my $in       = shift ;
   my $roidef   = $in->{roidef} ;
   my $age      = $in->{age} ;
   my $specs = $in->{specs};
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   my $roi_list = [];
   if (exists $in->{roi_list}) {
      $roi_list = $in->{roi_list} ;
   } else {
      @{$roi_list} = keys %{$roidef->{roi_meta}} ;
   }

   foreach my $orig_roi (@{$roi_list}) {
      print STDERR "* creating lrcheck ROI for $orig_roi\n";

# 1. flip across mid-sagittal plane
      my $flipped_roi = flip_roi({
         age            => $in->{age},
         roidef         => $in->{roidef},
         roi_list       => [$orig_roi],
      }) ;

      my $new_name = $orig_roi."-LR_joint" ;
# 2. join with orig ROI
      my $lr_joint_roi = join_roi({
         new_name       => $new_name,
         roidef1        => $roidef,
         roi_name1      => $orig_roi,
         roidef2        => $flipped_roi,
         roi_name2      => "flip_".$orig_roi,
      }) ;

# 3. split across mid-sagittal plane; if midplane included, give to both sides.
      my $lname = $orig_roi."-LR_part0" ;
      my $rname = $orig_roi."-LR_part1" ;
      my $lr_split_roi = midsagittal_split_roi({
         roidef         => $lr_joint_roi,
         roi_name       => $new_name,
         age            => $age,
         specs          => $specs,
         l_name         => $lname,
         r_name         => $rname
      }) ;

      foreach my $newroi ($lname, $rname) {
         $roidef->{roi_meta}->{$newroi} = {} ;
         $roidef->{num_roidef_points}->{$newroi} = 0;
      }

# 4. add new rois 'ORIG-LR_part0' and 'ORIG-LR_part1' to orig $roidef ;
      foreach my $xyz (keys %{$lr_split_roi->{point2roi}}) {
         foreach my $roi (keys %{$lr_split_roi->{point2roi}->{$xyz}}) {
            $roidef->{point2roi}->{$xyz}->{$roi}++;
            $roidef->{num_roidef_points}->{$roi}++; } }
      foreach my $newroi ($lname, $rname) {
         print STDERR "$newroi: ".$roidef->{num_roidef_points}->{$newroi}." points\n"; }
   }

}


sub join_roi {
# joins ->{roi1} and ->{roi2} into an roi named ->{new_name}
# if roi_name1 or roi_name2 are specified, then only those ROI
#  out of the original roidef will be merged.

   my $in       = shift ;

   my $new_roidef = {};
   foreach my $i (1, 2) {
      foreach my $xyz (keys %{$in->{"roidef$i"}->{point2roi}}) {
         if (exists $in->{"roi_name$i"} &&
             !exists $in->{"roidef$i"}->{point2roi}->{
                 $xyz}->{$in->{"roi_name$i"}}) {next;}
         $new_roidef->{point2roi}->{$xyz}->{$in->{new_name}}++; } }

   return $new_roidef ;
}


sub midsagittal_split_roi {
# splits ->{roi_name} about mid-sagittal plane into ->{l_name} and ->{r_name}

   my $in = shift ;
   my $roidef = $in->{roidef} ;
   my $age = $in->{age} ;
   my $specs = $in->{specs};
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }
   my $l_name   = $in->{l_name} ;
   my $r_name   = $in->{r_name} ;
   my $roi_name = $in->{roi_name} ;

   my $l_max = POSIX::ceil($specs->{dims}->{$age}->{z} / 2);
   my $r_min = POSIX::floor($specs->{dims}->{$age}->{z} / 2);

   my $new_roidef = {};
   foreach my $xyz (keys %{$roidef->{point2roi}}) {
      if (!exists $roidef->{point2roi}->{$xyz}->{$roi_name}) {next;}
      my ($x,$y,$z) = split(',',$xyz) ;
      if ($z < $l_max) {
#         print STDERR "assigned $xyz to LEFT ($z < $l_max)\n";
         $new_roidef->{point2roi}->{$xyz}->{$l_name}++ ;
      }
      if ($z > $r_min) {
#         print STDERR "assigned $xyz to RIGHT ($z > $r_min)\n";
         $new_roidef->{point2roi}->{$xyz}->{$r_name}++ ;
      }
   }

   return $new_roidef ;
}

1 ;
