=head1 NAME

alnmnr::io.pm - routines to interact with ABA and ALLENMINER files

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

package alnmnr::io ;
use strict;
use warnings;
use Cwd;
use File::Path qw/mkpath/ ;
use File::Temp qw/tempfile tempdir/ ;
use POSIX qw/floor ceil/ ;


=head2 aba_parse_xpz()

   Title:       aba_parse_xpz()
   Function:    Parses binary ABA XPZ files and converts it to a tab-delimited
                  text format or a PDB file (XPZ = Developing mouse atlas)
   Args:        $_->{xpz} - XPZ file path
                $_->{out_xyz} - file name for XYZ output [optional]
                $_->{out_pdb} - file name for PDB format output [optional]
                $_->{out_data} - flag to return hash holding data
                $_->{roidef} - optional; if specified, only stores points in ROI
                $_->{energy_fl} - optional; if specified, sends back xyz hash
                                  pointing to `energy' surrogate.

                     Energy = Intensity * sqrt(Density)

                     We don't have the raw pixel/cube mapping to get the 'real'
                     energy per ABA White Paper. this should do for now.

   Returns:     $->{} (if $_->{out_data} is specified}
   File out:    XYZ output (in $_->{out_xyz}, if specified}
                  1. x
                  2. y
                  3. z
                  4. density
                  5. intensity

  NOTE: Density values in the Brain Explorer are floor(VALUE/65)
  NOTE: Density values also stored as 'num_expressors' and
        intensity values as 'expression_level' so XPR routines can
        be reused on XPZ files without significant modification.

=cut

sub aba_parse_xpz {
   require File::Basename ;

# XPZ file are ZIP archives of:
# 1. image_series.xml - contains meta info (ascii)
# 2. intensity.mhd - header file describing intensity file format (ascii)
# 3. intensity.raw - raw intensity data (binary format)
# 4. density.mhd - header file describing density file format (ascii)
# 5. density.raw - raw density data (binary format)

   my $in = shift ;
   my ($fn, $fh) ;
   if (!exists $in->{fn}) {
      die "XPZ file not specified" ; }
   if (!-s $in->{fn})     {
      print STDERR "XPZ file $in->{fn} not found " ;
      return {error_fl => "file not found "} ;
   }
   $fn->{xpz} = $in->{fn} ;


   my $data = {};
   $data->{meta}->{xpz_name} = File::Basename::basename($fn->{xpz}) ;

# 0. Make a temporary directory and unzip the file.
   
   my $pwd = getcwd ;
   my $unzip_dir = tempdir("unzip_xpz_XXXXX") ;
   alnmnr::core::safe_copy($fn->{xpz}, $unzip_dir) ;
   chdir($unzip_dir);

   system("unzip ".$data->{meta}->{xpz_name}." > /dev/null") ;
   chdir($pwd) ;

   my @expected_files = ('image_series.xml', 'intensity.mhd', 'intensity.raw',
                         'density.mhd', 'density.raw') ;
   foreach my $cur_fn (@expected_files) {
      if (!-s $unzip_dir.'/'.$cur_fn) {
         return {error_fl =>
                 "Fatal error: $fn->{xpz} XPZ file did not contain $cur_fn"};
      }
   }

# energy_only flag set if only vector of energy values wanted for simsearch
   my $energy_only = 0 ;
   if (exists $in->{energy_only} && $in->{energy_only} == 1) {
      $energy_only = 1 ; }

# 1. Parse image-series meta info
#   $data->{meta}->{imageseries} = read_xpz_imageseries_file("image_series.xml");

   parse_image_series_xml({
      data      => $data,
      fn        => $unzip_dir.'/image_series.xml'
   }) ;

# 2. Parse intensity header file
#-------------------------------
# ObjectType = Image
# NDims = 3
# BinaryData = True
# BinaryDataByteOrderMSB = False
# CompressedData = False
# TransformMatrix = 1 0 0 0 1 0 0 0 1
# Offset = 0 0 -160
# CenterOfRotation = 0 0 0
# AnatomicalOrientation = RAI
# ElementSpacing = 80 80 80
# DimSize = 70 75 40
# ElementType = MET_UCHAR
# ElementDataFile = intensity.raw

   my $raw_data = {};
   my $max_val = {} ;
   foreach my $file_type (qw/intensity density/) {
      open(MHD, $unzip_dir.'/'.$file_type.'.mhd') ;
      while(my $line = <MHD>) {
         chomp $line ;
         my ($key, $val) = split(' = ', $line) ;
         $data->{meta}->{$file_type}->{$key} = $val ;
      }
      close(MHD) ;

      ($data->{meta}->{dims_x},
       $data->{meta}->{dims_y},
       $data->{meta}->{dims_z}) = split(' ',
         $data->{meta}->{$file_type}->{DimSize}) ;
      my $element_type = $data->{meta}->{$file_type}->{ElementType} ;

#MET_UCHAR = 1 byte; MET_USHORT = 2 byte
      my ($read_size, $unpack_tmpl) ;
      if ($element_type eq 'MET_UCHAR') {
         $read_size = 1 ;
         $unpack_tmpl = 'C8' ;
      } elsif ($element_type eq 'MET_USHORT') {
         $read_size = 2 ;
         $unpack_tmpl = 'S' ;
      } else {
         return {error_fl => "FATAL ERROR: $fn->{xpz}/$file_type.raw".
                             " ElementType $element_type not recognized\n" };
      }

      open(RAW, $unzip_dir.'/'.$file_type.'.raw') ;
      my $buffer ;
      my ($cur_x, $cur_z, $cur_y) = (0,0,0) ;
      while (read(RAW, $buffer, $read_size)) {
         my $val = unpack($unpack_tmpl, $buffer) ;
         push @{$raw_data->{$file_type}}, $val ;
         if (!exists $max_val->{$file_type} ||
             $val > $max_val->{$file_type}) {
            $max_val->{$file_type} = $val ;
         }
      }
      close(RAW) ;
#      print STDERR "MAX $file_type: ".$max_val->{$file_type}."\n";
   }


   my $xy_plane_size = $data->{meta}->{dims_x} * $data->{meta}->{dims_y} ;
   $data->{xyz} = [] ;
   foreach my $j ( 0 .. $#{$raw_data->{intensity}}) {

#NOTE: fpd111224_0715  - what is the conversion code for density?
# per e-notebook 2010-Oct-22 Brain Explorer 2 displays density as floor(VAL/65)
      $raw_data->{density}->[$j] /= 65535 ;

      if ($raw_data->{density}->[$j]   == 0 ||
          $raw_data->{intensity}->[$j] == 0) {next;}

      my $z = POSIX::floor($j / $xy_plane_size) ;
      my $y = POSIX::floor(($j - $z * $xy_plane_size) /
                              $data->{meta}->{dims_x}) ;

      my $x = $j - $z * $xy_plane_size - $y * $data->{meta}->{dims_x} ;

# if roidef specified, only store points within ROI
      if (defined $in->{roidef} &&
          !exists $in->{roidef}->{point2roi}->{$x.','.$y.','.$z}) {next;}

#HERENOW 111220_1524 - double check `energy' equation
         $raw_data->{energy}->[$j] = sprintf("%.3f",
            $raw_data->{intensity}->[$j] * $raw_data->{density}->[$j]) ;

      if ($energy_only) {
         $data->{energy}->{"$x,$y,$z"} = $raw_data->{energy}->[$j] ;
      } else {
         push @{$data->{xyz}}, {
            x              => $x,
            y              => $y,
            z              => $z,
            intensity      => $raw_data->{intensity}->[$j],
            density        => $raw_data->{density}->[$j],
            energy         => $raw_data->{energy}->[$j],
         } ;
      }
   }
   $raw_data = {} ;

   if (exists $in->{out_xyz}) {
      my $fh_xyz ;
      open($fh_xyz, "> $in->{out_xyz}" );

      foreach my $metainfo (@{$data->{meta}->{gene}}) {
         print {$fh_xyz} "#meta\t$metainfo\n" ;
      }
      print {$fh_xyz} "#imageseries_id\t".$data->{meta}->{imageseries_id}."\n" ;
      print {$fh_xyz} "#gene\t".$data->{meta}->{gene_name}."\n" ;
      print {$fh_xyz} "#stage\t".$data->{meta}->{stage}."\n" ;
      print {$fh_xyz} "#slice\t".$data->{meta}->{plane}."\n" ;
      my @fields = qw/x y z intensity density energy/;
      print {$fh_xyz} '#'.join("\t", @fields)."\n" ;
      my $xy_plane_size = $data->{meta}->{dims_x} *
                          $data->{meta}->{dims_y} ;
      foreach my $point (@{$data->{xyz}}) {
         my @outvals = ($point->{x},
                        $point->{y},
                        $point->{z},
                        $point->{intensity},
                        sprintf("%.5f", $point->{density}),
                        $point->{energy}) ;
         print {$fh_xyz} join("\t", @outvals)."\n";
      }
      close($fh_xyz) ;
   }

   if (exists $in->{out_pdb}) {
      my $out_pdb_fn = $in->{out_pdb} ;
      my $compress_fl = 0 ;
      if ($out_pdb_fn =~ /\.gz$/) {
         $compress_fl = 1 ; $out_pdb_fn =~ s/\.gz// ;}

      aba_displayxpr_pdbformat({
         out_pdb => $out_pdb_fn,
         data => $data,
      }) ;

      if ($compress_fl) {
         system("gzip $out_pdb_fn") ; }
   }

# Clean up temporary unzip directory
   foreach my $file (@expected_files) {
      unlink $unzip_dir.'/'.$file ; }
   unlink $unzip_dir.'/'.$data->{meta}->{xpz_name} ;
   chdir $pwd ;
   rmdir $unzip_dir ;

   return $data ;

}


=head2 aba_parse_xpr()

   Title:       aba_parse_xpr()
   Function:    Parses binary ABA XPR files and converts it to a tab-delimited
                  text format or a PDB file
   Args:        $_->{xpr} - XPR file path
                $_->{out_xyz} - file name for XYZ output [optional]
                $_->{out_pdb} - file name for PDB format output [optional]
                $_->{out_data} - flag to return hash holding data

   Returns:     $->{} (if $_->{out_data} is specified}
   File out:    XYZ output (in $_->{out_xyz}, if specified}
                ^# comment lines holding meta data, header_id, riboprobe info
                  1. grid_id
                  2. x
                  3. y
                  4. z
                  5. cell_diameter
                  6. grid_area
                  7. expression_level
                  8. section
                  9. location_1
                 10. location_2
                 11. num_expressors

                PDB output (in $_->{out_pdb}, if specified}
                  PDB ATOM RECORD format - per aba_displayxpr_pdbformat()

=cut


# NOTE: what is missing:
# * header info (first 128 bytes after uncompression)
#   - only extracts (i) identifier number, and (ii) number of aff entries
#
# * also dont know what is in the 32 bytes betwen aff entries
#  (AFF entries = primary and express zoomify image locations - probably used
#     by their software to link to in situ images)
#
# inspiration from: 
# http://www.openrce.org/articles/full_view/16
# got python zlib uncompression script (and strategy) from web

sub aba_parse_xpr {

   require File::Basename ;

#cut off first 16 bytes, uncompress rest with zlib

   my $in = shift ;
   my ($fn, $fh) ;
   if (!exists $in->{fn}) {
      die "XPR file not specified" ; }
   if (!-s $in->{fn}) {
      die "XPR file $in->{fn} not found " ; }
   $fn->{xpr} = $in->{fn} ;

   ($fh->{xpr_frag}, $fn->{xpr_frag}) =
      tempfile("xpr_fragment.XXXXX", SUFFIX => ".frag") ;
   close($fh->{xpr_frag}) ;

   _aba_parse_xpr_uncompress_xpr({
      fn_in => $fn->{xpr},
      fn_out => $fn->{xpr_frag}
   }) ;


   my $data ;
   $data->{meta} = {} ;
   if ($fn->{xpr} =~ /sagittal/) {
      $data->{meta}->{plane} = 'sagittal' ;
   } else {
      $data->{meta}->{plane} = 'coronal' ;
   }

   {
      my $a  = File::Basename::basename($fn->{xpr}) ;
      $data->{meta}->{xpr_name} = $a ;
      $a =~ s/_coronal.*// ;
      $a =~ s/_sagittal.*// ;
      $data->{meta}->{gene_name} = $a ;
   }



#read file header

   my $buffer ;
   open( XPRDATA, $fn->{xpr_frag}) ;
   read( XPRDATA,$data->{header},64);
   $data->{meta}->{header_id} = unpack("I",substr($data->{header},0,4));
   $data->{hexstring_identifier} = unpack("H4",substr($data->{header},0,4));
#   print STDERR "id $data->{identifier} from hexstring $data->{hexstring_identifier}\n" ;

   $data->{hexstring_num_aff_records} = unpack("H4",substr($data->{header},56,4));
   $data->{num_aff_records} = unpack("I",substr($data->{header},56,4));
#   print STDERR "number of aff records: $data->{num_aff_records} from $data->{hexstring_num_aff_records}\n" ;

#read gene meta data
   while (!exists $data->{meta}->{rp}) {
      read(XPRDATA,$buffer, 1) ;
      my $cur_reclength = unpack("c",$buffer);

      if ($cur_reclength <= 0) {
         close(XPRDATA) ;
         my $error = "ERROR: xpr parsing error (".__LINE__.
                      "): negative record length" ;
         print STDERR $error."\n" ;
         return {error_fl => $error};
      }

#      my $hex_cur_reclength = unpack("H2",$buffer);
#      print STDERR "length of $cur_reclength (from $hex_cur_reclength)\n" ;

      read(XPRDATA,$buffer,$cur_reclength) ;
      my $hex_cur_rec = unpack("H".$cur_reclength, $buffer) ;
      my $cur_rec = unpack("A".$cur_reclength, $buffer) ;

#      print STDERR "$hex_cur_rec giving $cur_rec\n";

      if ($cur_rec =~ /^RP/) {
         $data->{meta}->{rp} = $cur_rec ;
      } else {
         push @{$data->{meta}->{gene}}, $cur_rec ;
      }
   }

   my $j = 0 ;
   while ($j < $data->{num_aff_records}) {
      read(XPRDATA,$buffer, 16) ;
      my $hex_pre_aff = unpack("H16", $buffer) ;
      push @{$data->{aff_header}}, $hex_pre_aff ;

      my $aff_id = unpack("c", substr($buffer,0,2)) ;
      push @{$data->{aff_id}}, $aff_id ;
#      print STDERR "hex is $hex_pre_aff, affid = $aff_id\n" ;

      read(XPRDATA,$buffer, 1) ;
      my $hex_primary_aff_length = unpack("c", $buffer) ;

# fpd101205_2132 - catches corrupted files like Cks2_coronal.
      if ($hex_primary_aff_length <= 0) {
         close(XPRDATA) ;
         my $error = "ERROR: xpr parsing error (".__LINE__.
                      "): negative record length" ;
         print STDERR $error."\n" ;
         return {error_fl => $error};
      }

      read(XPRDATA,$buffer,$hex_primary_aff_length) ;
      my $hex_primary_aff = unpack("H".$hex_primary_aff_length,$buffer) ;
      my $primary_aff = unpack("A".$hex_primary_aff_length,$buffer) ;
      push @{$data->{primary_aff}}, $primary_aff ;
#      print STDERR "$hex_primary_aff to\n$primary_aff\n" ;

      read(XPRDATA,$buffer, 1) ;
      my $hex_express_aff_length = unpack("c", $buffer) ;

      if ($hex_express_aff_length <= 0) {
         close(XPRDATA) ;
         my $error = "ERROR: xpr parsing error (".__LINE__.
                      "): negative record length" ;
         print STDERR $error."\n" ;
         return {error_fl => $error};
      }

      read(XPRDATA,$buffer,$hex_express_aff_length) ;
      my $hex_express_aff = unpack("H".$hex_express_aff_length,$buffer) ;
      my $express_aff = unpack("A".$hex_express_aff_length,$buffer) ;
      push @{$data->{express_aff}}, $express_aff ;
#      print STDERR "$hex_express_aff to\n$express_aff\n" ;

      $j++ ;
   }

   $data->{read_something} = 0 ;
   while (read(XPRDATA,$buffer,24)) {
      my ($hexes, $vals) ;
      $hexes->{grid_id} = unpack("H8", substr($buffer,0,4)) ;
      $vals->{grid_id} = unpack("I8", substr($buffer,0,4)) ;
#      print STDERR "grid is $hexes->{grid_id} to $vals->{grid_id}\n" ;
      my $z = floor(($vals->{grid_id} -1 )/ 10773) ;
      my $y = floor(($vals->{grid_id} - 1 - ($z * 10773)) / 133) ;
      my $x = $vals->{grid_id} - 1 - $z * 10773 - $y * 133 ;
      ($vals->{x}, $vals->{y}, $vals->{z}) = ($x, $y, $z) ;
      $vals->{x} = sprintf("%d", $vals->{x}) ;
      $vals->{y} = sprintf("%d", $vals->{y}) ;
      $vals->{z} = sprintf("%d", $vals->{z}) ;

      $hexes->{cell_diameter} = unpack("H8", substr($buffer,4,4)) ;
      $vals->{cell_diameter} = sprintf("%.3f", unpack("f", substr($buffer,4,4))) ;
#      print STDERR "cell diameter is $hexes->{cell_diameter} to $vals->{cell_diameter}\n" ;

      $hexes->{grid_area} = unpack("H8", substr($buffer,8,4)) ;
      $vals->{grid_area} = sprintf("%.2f", unpack("f", substr($buffer,8,4))) ;
#      print STDERR "grid_area is $hexes->{grid_area} to $vals->{grid_area}\n" ;

      $hexes->{expression_level} = unpack("H8", substr($buffer,12,4)) ;
      $vals->{expression_level} = sprintf("%.3f", unpack("f", substr($buffer,12,4))) ;
#      print STDERR "expression_level is $hexes->{expression_level} to $vals->{expression_level}\n" ;

      $hexes->{section} = unpack("H4", substr($buffer,16,2)) ;
      $vals->{section} = unpack("s", substr($buffer,16,2)) ;
#      print STDERR "section is $hexes->{section} to $vals->{section}\n" ;

      $hexes->{location_1} = unpack("H4", substr($buffer,18,2)) ;
      $vals->{location_1} = unpack("s", substr($buffer,18,2)) ;
#      print STDERR "location_1 is $hexes->{location_1} to $vals->{location_1}\n" ;

      $hexes->{location_2} = unpack("H4", substr($buffer,20,2)) ;
      $vals->{location_2} = unpack("s", substr($buffer,20,2)) ;
#      print STDERR "location_2 is $hexes->{location_2} to $vals->{location_2}\n" ;

      $hexes->{num_expressors} = unpack("H4", substr($buffer,22,2)) ;
      $vals->{num_expressors} = unpack("s", substr($buffer,22,2)) ;
#      print STDERR "num_expressors is $hexes->{num_expressors} to $vals->{num_expressors}\n" ;

      push @{$data->{xyz}}, $vals ;
      $data->{read_something} = 1 ;
#      print STDERR "\n\n\n" ;
   }
   close( XPRDATA) ;
   unlink $fn->{xpr_frag} ;

   if (exists $in->{out_xyz}) {
      my $fh_xyz ;
      open($fh_xyz, "> $in->{out_xyz}" );

      foreach my $metainfo (@{$data->{meta}->{gene}}) {
         print {$fh_xyz} "#meta\t$metainfo\n" ;
      }
      print {$fh_xyz} "#header_id\t$data->{meta}->{header_id}\n" ;
      print {$fh_xyz} "#riboprobe\t$data->{meta}->{rp}\n" ;

      my @fields = qw/grid_id x y z cell_diameter grid_area expression_level section location_1 location_2 num_expressors/ ;
      print {$fh_xyz} '#'.join("\t", @fields)."\n" ;
      foreach my $point (@{$data->{xyz}}) {
         print {$fh_xyz} join("\t", @{$point}{@fields})."\n" ;
      }
      close($fh_xyz) ;
   }

   if (exists $in->{out_pdb}) {

      my $out_pdb_fn = $in->{out_pdb} ;
      my $compress_fl = 0 ;
      if ($out_pdb_fn =~ /\.gz$/) {
         $compress_fl = 1 ; $out_pdb_fn =~ s/\.gz// ;}

      aba_displayxpr_pdbformat({
         out_pdb => $out_pdb_fn,
         data => $data,
      }) ;

      if ($compress_fl) {
         system("gzip $out_pdb_fn") ; }
   }

   return $data ;

}



=head2 aba_displayxpr_pdbformat()

   Title:       aba_displayxpr_pdbformat
   Function:    Prints out parsed ABA XPR (or XPZ) data as a PDB file
   Args:        $_->{data} = parsed ABA XPR data - from aba_parse_xpr()
                                      (or XPZ file from aba_parse_xpz())
                $_->{out_pdb} - output PDB filename
   Returns:     nothing
   Files out:   Output PDB file ($_->{out_pdb})
                 GLY CA atoms - chains from a-z,A-Z, updated when residue #
                 reaches 99999
                 occupancy field = density
                 B-factor field = energy
                 charge = intensity

=cut

sub aba_displayxpr_pdbformat {

   my $in = shift ;
   my $data = $in->{data} ;

   my @chains = (1 .. 9);
   push @chains, ('a' .. 'z'), ('A' .. 'Z');
   push @chains, ('a' .. 'z'), ('A' .. 'Z');
   push @chains, ('a' .. 'z'), ('A' .. 'Z');

   my $fh_pdb ;
   open($fh_pdb, "> $in->{out_pdb}") ;

   my $atomno = 1;
   my $chainid = shift @chains ;
   foreach my $point (@{$data->{xyz}}) {
#      my $resno = floor($atomno / 10) ;
      my $resno = $atomno ;

      print {$fh_pdb} "ATOM  ".
            sprintf("%5d", $atomno),' ',
            " CA ",
            " ",
            "GLY",
            ' ',
            $chainid,
            sprintf("%4d", $resno),
            " ",
            '   ',
            sprintf("%8.3f",$point->{x}),
            sprintf("%8.3f",$point->{y}),
            sprintf("%8.3f",$point->{z}),
            sprintf("%6.2f",$point->{density}),
            sprintf("%6.2f",$point->{intensity}),
            '           ',
            'CA',
            $point->{energy}."\n" ;

      $atomno++ ;
      if ($atomno > 9999) {
         $atomno = 1 ;
         $resno  = 1 ;
         $chainid = shift @chains ;
      }
   }
   close($fh_pdb) ;

}


=head2 _aba_parse_xpr_uncompress_xpr {

   Title:       _aba_parse_xpr_uncompress_xpr()
   Function:    creates and runs a python zlib uncompression script to
                 decompress the zlibbed XPR block
   Args:        $_->{fn_out} - output file name
   Files in:    XPR file ($_->{fn_in}) ;
   Files out:   uncompressed XPR block ($_->{fn_out})

=cut

sub _aba_parse_xpr_uncompress_xpr {

   my $in = shift ;

   my $xpr_uncompress = "
import zlib, struct, sys

try:
    import psyco
    psyco.full()
except ImportError:
    pass

def main(args):
    fp = file('".$in->{fn_in}."', 'rb')
    data = fp.read(16) #starts 16 in
    file('".$in->{fn_out}."', 'wb').write(zlib.decompress(fp.read()))

if __name__=='__main__':
    sys.exit(main(sys.argv[1:]))
" ;

   my ($fh, $fn) = tempfile() ;
   print $fh ($xpr_uncompress) ;
   close($fh) ;

   system("python $fn") ;

   if (!-s $in->{fn_out}) {
      die "XPR uncompression didn't work" ; }
   return ;

}


=head2 aba_parse_sva()

   Title:       aba_parse_sva()
   Function:    Parses binary ABA SVA files

   Args:        $_->{sva} - XPR file path
                $_->{roidef} - optionally specify ROIdef structure, if 
                               only those ROI expression is desired
   Returns:     $->{x,y,z} = expression energy

=cut

sub aba_parse_sva {

   require File::Basename ;
   my $in = shift ;

   my $fn ;
   $fn->{sva} = $in->{fn} ;

#Comment:Smoothed energy volume for imageseriesId 71919026
#Dimensions:67,41,58
#63,9,8,1.43691e-06
#64,9,8,1.69753e-05

   open(SVAF, $fn->{sva}) ;

   my $comment_line = <SVAF> ;
   my $dim_line = <SVAF> ; chomp $dim_line; 

   my $xpr_like = 0 ;
   if (exists $in->{xpr_like} ||
       exists $in->{out_xyz}  ||
       exists $in->{out_pdb}    ) { $xpr_like = 1 ; }

   my $xpr_data ;
   if ($xpr_like) {
      $xpr_data->{meta} = {} ;
      if ($fn->{sva} =~ /sagittal/) {
         $xpr_data->{meta}->{plane} = 'sagittal' ;
      } else {
         $xpr_data->{meta}->{plane} = 'coronal' ;
      }

      $xpr_data->{meta}->{stage} = 'P56' ;
      $xpr_data->{meta}->{xpr_name} = File::Basename::basename($fn->{sva}) ;
      ($xpr_data->{meta}->{gene_name},
       $xpr_data->{meta}->{plane},
       $xpr_data->{meta}->{imageseries_id}) =
       split('_', $xpr_data->{meta}->{xpr_name});
      $xpr_data->{meta}->{imageseries_id} =~ s/\.sva$// ;

      ($xpr_data->{meta}->{dims_x},
       $xpr_data->{meta}->{dims_y},
       $xpr_data->{meta}->{dims_z}) =
         ($dim_line =~ /Dimensions:([0-9]+),([0-9]+),([0-9]+)/) ;
      $xpr_data->{xyz} = [] ;
   }

   my $data ;
   while (my $line = <SVAF>) {
      chomp $line;
      my ($x,$y,$z,$val) = split(/\,/, $line) ;
      my $point = $x.','.$y.','.$z ;
      if (defined $in->{roidef_200} && 
          !exists $in->{roidef_200}->{point2roi}->{$point}) {next;}
      if ($xpr_like) {
         push @{$xpr_data->{xyz}}, {
            x => $x,
            y => $y,
            z => $z,
            energy => $val,
            expression_level => $val,
            num_expressors => $val,
            grid_area => $val,
            cell_diameter => $val,
         } ;
      }
      $data->{$point} = $val ;
   }
   close(SVAF) ;

   if (exists $in->{out_xyz}) {
      my $fh_xyz ;
      open($fh_xyz, "> $in->{out_xyz}" );

      foreach my $metainfo (@{$xpr_data->{meta}->{gene}}) {
         print {$fh_xyz} "#meta\t$metainfo\n" ;
      }
      print {$fh_xyz} "#imageseries_id\t".
                      $xpr_data->{meta}->{imageseries_id}."\n" ;
      print {$fh_xyz} "#gene\t".$xpr_data->{meta}->{gene_name}."\n" ;
      print {$fh_xyz} "#stage\t".$xpr_data->{meta}->{stage}."\n" ;
      print {$fh_xyz} "#slice\t".$xpr_data->{meta}->{plane}."\n" ;
      my @fields = qw/x y z energy/;
      print {$fh_xyz} '#'.join("\t", @fields)."\n" ;
      my $xy_plane_size = $xpr_data->{meta}->{dims_x} *
                          $xpr_data->{meta}->{dims_y} ;
      foreach my $point (@{$xpr_data->{xyz}}) {
         print {$fh_xyz} join("\t", $point->{x},
                                    $point->{y},
                                    $point->{z},
                                    $point->{energy})."\n";
      }
      close($fh_xyz) ;
   }

   if (exists $in->{out_pdb}) {
      my $out_pdb_fn = $in->{out_pdb} ;
      my $compress_fl = 0 ;
      if ($out_pdb_fn =~ /\.gz$/) {
         $compress_fl = 1 ; $out_pdb_fn =~ s/\.gz// ;}

      aba_displayxpr_pdbformat({
         out_pdb => $out_pdb_fn,
         data => $xpr_data,
      }) ;

      if ($compress_fl) {
         system("gzip $out_pdb_fn") ; }
   }

   if (exists $in->{xpr_like}) {
      return $xpr_data ;
   } else {
      return $data ;
   }

}



#sub readin_roi_list_results {
#
#   my $in = shift;
#   my $f2i ;
#   foreach my $j ( 0 .. $#{$in->{headers}->{headers}}) {
#      $f2i->{$in->{headers}->{headers}->[$j]} = $j ; }
#
#   open(RESULTSF, $in->{fn}) ;
#   my $roi_data ;
#   while (my $line = <RESULTSF>) {
#      chomp $line;
#      if ($line =~ /^#/) {next;}
#      my $entry ;
#      my @t = split(/\t/, $line) ;
#      foreach my $j (0 .. $#t) {
#         $entry->{$in->{headers}->{headers}->[$j]} = $t[$j] ;
#      }
#      if ($entry->{plane} ne 'sagittal' &&
#          $entry->{plane} ne 'coronal') {
#         next; }
#      foreach my $f (sort keys %{$f2i}) {
#         push @{$roi_data->{$entry->{roi_name}}->{$entry->{plane}}->{$f}},
#            $entry->{$f} ; }
#   }
#   close(RESULTSF) ;
#   return $roi_data ;
#
#}


=head2 readin_roi_file()

   Title:       readin_roi_file()
   Function:    Read in the ROI definition
   Args:        $_->{fn} [optional - if not will read in from STDIN]
   Returns:     $_->{}->{point2strx}->{x,y,z} = brain_strx_id ;
                $_->{}->{strx2points}->{brain_strx_id}->[i] = "x_i,y_i,z_i"
   Files in:    AtlasAnnotation100.sva

=cut

sub readin_roi_file {

   my $in = shift ;
   my $age = $in->{age} ;

   my $options = {
      roi_axes => {
         xmin => 1,
         xmax => 1,
         ymin => 1,
         ymax => 1,
         zmin => 1,
         zmax => 1,
      }
   } ;


   my $specs ;
   if (exists $in->{specs}) {
      $specs = $in->{specs} ;
   } else {
      $specs = alnmnr::getspecs() ;
   }

   if (!exists $specs->{dims}->{$age}) {
      die "FATAL ERROR: age $age is not one of the ABA timepoints: ".
          join(", ", sort keys %{$specs->{dims}})."\n" ;
   }

   if (!-s $in->{fn}) {
      die "FATAL ERROR: ROI definition file ".$in->{fn}." not found\n"; }

# view-volume = mm / $scaling
   my $scaling = $specs->{view_volume}->{$in->{age}}->{spacing} / 1000;
   my $scaling_inv = 1000 / $specs->{view_volume}->{$in->{age}}->{spacing} ;

# don't divide by $scaling, just multiply by $scaling_inv if needed
# because of precision issue.
# eg, a = 11.6 / (1 / 5 ) = 57.999999999999992894572642398998, not 58.
# so if we build an array of (a .. n), will start from 57, not 58.
# alt, call ceil() to be safe...

   print STDERR "Reading ROI file $in->{fn}\n" ;
   open(ROIDEF, $in->{fn}) ;
   my $roidef ;
   my $results = {};

#120104_1824 - convert all coordinates as they are read from mm to view volume
   while (my $line = <ROIDEF>) {
      if ($line =~ /^\#/ || $line =~ /^$/) {next;}
      chomp $line;
      my ($name, $subname, $roi_type, @t) = split(/\t/, $line) ;
      if ($subname eq '') {$subname = $name;}
      if (!exists $roidef->{$name}->{$subname}) {
         $results->{roi_meta}->{$name}->{$subname}->{roi_spec} = $line ;
         $results->{roi_meta}->{$name}->{$subname}->{roi_type} = $roi_type ;
         $roidef->{$name}->{$subname}->{roi_type} = $roi_type ;
      }

      if ($roi_type =~ /^boxbounds/) {
         my ($axes_type, $value) = @t ;
         if (!exists $options->{roi_axes}->{$axes_type}) {
            die "ERROR: roi axes $axes_type not recognized\n".
                "must be one of: ".join(", ",keys %{$options->{roi_axes}})."\n";
         }

         $value =~ s/x/(x \* $scaling)/g ;
         $value =~ s/y/(y \* $scaling)/g ;
         $value =~ s/z/(z \* $scaling)/g ;

         if ($value =~ /[xyz]/)  {
            $value = '( '.$value." ) * $scaling_inv" ;
            $roidef->{$name}->{$subname}->{$axes_type}->{value} = $value ;
            $roidef->{$name}->{$subname}->{$axes_type}->{type} = 'variable' ;
         } else {
            $roidef->{$name}->{$subname}->{$axes_type}->{type} = 'constant' ;

# even though constant, have to apply scaling factor
            $roidef->{$name}->{$subname}->{$axes_type}->{value} =
               $value * $scaling_inv ;
         }

      } elsif ($roi_type eq 'point') {

         my ($x,$y,$z) ;
         if ($t[0] =~ /\,/) {
            ($x,$y,$z) = split(',', $t[0]) ;
         } else {
            ($x,$y,$z) = ($t[0], $t[1], $t[2]) ;
         }
#fpd120405_1825 - internal coordinates must be integer (as in xpz files)
         my $x2 = sprintf("%.0f", $x * $scaling_inv) ;
         my $y2 = sprintf("%.0f", $y * $scaling_inv) ;
         my $z2 = sprintf("%.0f", $z * $scaling_inv) ;

         $roidef->{$name}->{$subname}->{points}->{"$x2,$y2,$z2"}++ ;

      } elsif ($roi_type eq 'brainstrx') {
         $roidef->{$name}->{$subname}->{brainstrx} = $t[0];
      }
   }
   close(ROIDEF) ;

   foreach my $roi_name (keys %{$roidef}) {
      my $a = {};
#      print STDERR "   read in roi $roi_name\n" ;

      my $off_points ; #enables use of strx=off AND boxbounds for same ROI

      foreach my $roi_subname (keys %{$roidef->{$roi_name}}) {
         my $cur_roidef = $roidef->{$roi_name}->{$roi_subname} ;

         if ($cur_roidef->{roi_type} eq 'brainstrx') {
            if (!exists $specs->{parsed_annotation}) {
               $specs->{parsed_annotation}= load_aba_Annotation({age => $age});}
            $cur_roidef->{points} = get_brainstrx_voxels({
               age => $age,
               specs => $specs,
               brainstrx => $cur_roidef->{brainstrx}
            })  ;
         }

         if ($cur_roidef->{roi_type} eq 'boxbounds') {
            $cur_roidef->{points} = alnmnr::roiops::calc_roi_boxbounds_pointlist({
               roidef => $cur_roidef,
               pointlist => 1,
            }) ;
         }
         $cur_roidef->{num_roidef_points} = keys %{$cur_roidef->{points}} ;
         foreach my $point (keys %{$cur_roidef->{points}}) {
            if ($cur_roidef->{points}->{$point} > 0) { #if positive
               $results->{point2roi}->{$point}->{$roi_name}++ ;
               $a->{$point}++ ;
            } else { #otherwise, an off point
               $off_points->{$point}++ ;
            }
         }
      }

      foreach my $off_point (keys %{$off_points}) {
         if (exists $results->{point2roi}->{$off_point}->{$roi_name}) {
            delete $results->{point2roi}->{$off_point}->{$roi_name} ; }
         if (exists $a->{$off_point}) {
            delete $a->{$off_point} ; }
      }

      $results->{num_roidef_points}->{$roi_name} = keys %{$a} ;
   }

   foreach my $roi_name (keys %{$roidef}) {
      print STDERR "  $roi_name: ";
      if ($results->{num_roidef_points}->{$roi_name} > 0) {
         print STDERR $results->{num_roidef_points}->{$roi_name}." points\n";
         next;
      }
      print STDERR "ignored because no voxels are defined\n";
      delete $results->{roi_meta}->{$roi_name} ;
      delete $results->{num_roidef_points}->{$roi_name} ;
   }

   return $results;

}


sub PERLZLIB_aba_parse_xpz {
# fpd111124_1821 - slower by ~5-6x than external unzip
#   use IO::Uncompress::Unzip qw/unzip $UnzipError/ ;
   require File::Basename ;

# XPZ file are ZIP archives of:
# 1. image_series.xml - contains meta info (ascii)
# 2. intensity.mhd - header file describing intensity file format (ascii)
# 3. intensity.raw - raw intensity data (binary format)
# 4. density.mhd - header file describing density file format (ascii)
# 5. density.raw - raw density data (binary format)

   my $in = shift ;
   my ($fn, $fh) ;
   if (!exists $in->{fn}) {
      die "XPZ file not specified" ; }
   if (!-s $in->{fn})     {
      print STDERR "XPZ file $in->{fn} not found " ;
      return {error_fl => "file not found "} ;
   }
   $fn->{xpz} = $in->{fn} ;


   my $data = {};
   $data->{meta}->{xpr_name} = File::Basename::basename($fn->{xpz}) ;

# Check that the XPZ is complete
   {
      my $xpz_u = new IO::Uncompress::Unzip $fn->{xpz} || 
         return {error_fl => "Fatal error: couldn't open ".$fn->{xpz}} ;

      my $expected_files = {} ;
      map {$expected_files->{$_}++; } ('image_series.xml', 'intensity.mhd',
                           'intensity.raw', 'density.mhd', 'density.raw') ;
      for (my $status = 1; !$xpz_u->eof(); $status = $xpz_u->nextStream()) {
         my $internal_fn = $xpz_u->getHeaderInfo()->{Name} ;
         delete $expected_files->{$internal_fn} ;
      }
      $xpz_u->close();

      if (keys %{$expected_files} > 0) {
         return {error_fl =>
                 "Fatal error: $fn->{xpz} XPZ file did not contain: ".
                 join(", ", sort keys %{$expected_files})};
      }
   }

   my $raw_data = {};
   my $xpz_u = new IO::Uncompress::Unzip $fn->{xpz} ||
      return {error_fl => "Fatal error: couldn't open ".$fn->{xpz}} ;
   my $raw_read_size = {} ; my $raw_unpack_tmpl = {} ;

   for (my $status = 1; !$xpz_u->eof(); $status = $xpz_u->nextStream()) {
      my $internal_fn = $xpz_u->getHeaderInfo()->{Name} ;

      if ($internal_fn eq 'image_series.xml') {
# 1. Parse image-series meta info

# <section-thickness> only in adult series...
         while (! $xpz_u->eof()) {
            my $line = $xpz_u->getline() ;
            chomp $line;
            if ($line =~ /^  <id/) {
               ($data->{meta}->{imageseries_id}) =
                  ($line =~ /<id>([0-9]+)<\/id>/) ;
            } elsif ($line =~ /<plane-of-section>/) {
               ($data->{meta}->{plane}) =
                  ($line =~ /<plane-of-section>(.+)<\/plane-of-section>/) ;
            } elsif ($line =~ /<section-thickness>/) {
               ($data->{meta}->{"section-thickness"}) =
                  ($line =~ /<section-thickness>(.+)<\/section-thickness>/) ;
            } elsif ($line =~ /<age>/) {
               ($data->{meta}->{stage}) =
                  ($line =~ /<age>(.+)<\/age>/) ;
            } elsif ($line =~ /<gene-symbol>/) {
               ($data->{meta}->{gene_name}) =
                  ($line =~ /<gene-symbol>(.+)<\/gene-symbol>/) ;
            } elsif ($line =~ /<\/gene>$/) {
               last;
            }
         }

      } elsif ($internal_fn =~ /mhd$/) {

         my ($file_type) = ($internal_fn =~ /(.+)\.mhd$/) ;

         while (! $xpz_u->eof()) {
            my $line = $xpz_u->getline() ;
            chomp $line;
            my ($key, $val) = split(' = ', $line) ;
            $data->{meta}->{$file_type}->{$key} = $val ;
         }

         ($data->{meta}->{dims_x},
          $data->{meta}->{dims_y},
          $data->{meta}->{dims_z}) = split(' ',
            $data->{meta}->{$file_type}->{DimSize}) ;
         my $element_type = $data->{meta}->{$file_type}->{ElementType} ;

#MET_UCHAR = 1 byte; MET_USHORT = 2 byte
         my ($read_size, $unpack_tmpl) ;
         if ($element_type eq 'MET_UCHAR') {
            $raw_read_size->{$file_type} = 1 ;
            $raw_unpack_tmpl->{$file_type} = 'C8' ;
         } elsif ($element_type eq 'MET_USHORT') {
            $raw_read_size->{$file_type} = 2 ;
            $raw_unpack_tmpl->{$file_type} = 'S' ;
         } else {
            return {error_fl => "FATAL ERROR: $fn->{xpz}/$file_type.raw".
                             " ElementType $element_type not recognized\n" };
         }

      } elsif ($internal_fn =~ /raw$/) {

         my ($file_type) = ($internal_fn =~ /(.+)\.raw$/) ;

         my $buffer ;
         my ($cur_x, $cur_z, $cur_y) = (0,0,0) ;
         while (!$xpz_u->eof()) {
            read($xpz_u, $buffer, $raw_read_size->{$file_type}) ;
            my $val = unpack($raw_unpack_tmpl->{$file_type}, $buffer) ;
            push @{$raw_data->{$file_type}}, $val ;
         }

      }
   }
   $xpz_u->close();


# 2. Parse intensity header file
#-------------------------------
# ObjectType = Image
# NDims = 3
# BinaryData = True
# BinaryDataByteOrderMSB = False
# CompressedData = False
# TransformMatrix = 1 0 0 0 1 0 0 0 1
# Offset = 0 0 -160
# CenterOfRotation = 0 0 0
# AnatomicalOrientation = RAI
# ElementSpacing = 80 80 80
# DimSize = 70 75 40
# ElementType = MET_UCHAR
# ElementDataFile = intensity.raw

   my $energy = {};
   my $xy_plane_size = $data->{meta}->{dims_x} * $data->{meta}->{dims_y} ;
   $data->{xyz} = [] ;
   foreach my $j ( 0 .. $#{$raw_data->{intensity}}) {
      if ($raw_data->{density}->[$j] == 0 &&
          $raw_data->{intensity}->[$j] == 0) {next;}

      my $z = POSIX::floor($j / $xy_plane_size) ;
      my $y = POSIX::floor(($j - $z * $xy_plane_size) /
                              $data->{meta}->{dims_x}) ;
      my $x = $j - $z * $xy_plane_size - $y * $data->{meta}->{dims_x} - 1 ;

#fpd101109_1354 - if roidef specified, only store points within ROI
      if (defined $in->{roidef} &&
          !exists $in->{roidef}->{point2roi}->{$x.','.$y.','.$z}) {next;}

      if (exists $in->{energy_fl} &&
          $raw_data->{density}->[$j] > 0 &&
          $raw_data->{intensity}->[$j] > 0) {
         $energy->{$x.','.$y.','.$z} = $raw_data->{intensity}->[$j]  *
                                          sqrt($raw_data->{density}->[$j]) ; }

      push @{$data->{xyz}}, {
         x => $x,
         y => $y,
         z => $z,
         expression_level => $raw_data->{intensity}->[$j],
         num_expressors => $raw_data->{density}->[$j],
         grid_area => $raw_data->{density}->[$j],
         cell_diameter => $raw_data->{intensity}->[$j],
      } ;
   }
   $raw_data = {} ;

   if (exists $in->{out_xyz}) {
      my $fh_xyz ;
      open($fh_xyz, "> $in->{out_xyz}" );

      foreach my $metainfo (@{$data->{meta}->{gene}}) {
         print {$fh_xyz} "#meta\t$metainfo\n" ;
      }
      print {$fh_xyz} "#imageseries_id\t".$data->{meta}->{imageseries_id}."\n" ;
      print {$fh_xyz} "#gene\t".$data->{meta}->{gene_name}."\n" ;
      print {$fh_xyz} "#stage\t".$data->{meta}->{stage}."\n" ;
      print {$fh_xyz} "#plane\t".$data->{meta}->{plane}."\n" ;
      my @fields = qw/x y z intensity density/;
      print {$fh_xyz} '#'.join("\t", @fields)."\n" ;
      my $xy_plane_size = $data->{meta}->{dims_x} *
                          $data->{meta}->{dims_y} ;
      foreach my $point (@{$data->{xyz}}) {
         print {$fh_xyz} join("\t", $point->{x},
                                    $point->{y},
                                    $point->{z},
                                    $point->{expression_level},
                                    $point->{num_expressors})."\n";
      }
      close($fh_xyz) ;
   }

   if (exists $in->{out_pdb}) {
      my $out_pdb_fn = $in->{out_pdb} ;
      my $compress_fl = 0 ;
      if ($out_pdb_fn =~ /\.gz$/) {
         $compress_fl = 1 ; $out_pdb_fn =~ s/\.gz// ;}

      aba_displayxpr_pdbformat({
         out_pdb => $out_pdb_fn,
         data => $data,
      }) ;

      if ($compress_fl) {
         system("gzip $out_pdb_fn") ; }
   }

## Clean up temporary unzip directory
#   foreach my $file (@expected_files) {
#      unlink $unzip_dir.'/'.$file ; }
#   unlink $unzip_dir.'/'.$data->{meta}->{xpr_name} ;
#   chdir $pwd ;
#   rmdir $unzip_dir ;

# For exprsim_query runs: 111122_1811 - dont just return energy here...
   if (exists $in->{energy_fl}) {
      $data->{energy} = $energy; }

# So the XPR routines can be easily adapted, level=intensity; expressors=density
# fpd101206_1135 - actually, already took care of this above in the {xyz} set.
#                  these two lines are most likely unnecessary, though harmless
   $data->{expression_level} = $data->{intensity} ;
   $data->{num_expressors}   = $data->{density} ;

   return $data ;

}


sub display_roidef {

   my $in = shift ;

   my $specs ;
   if (exists $in->{specs}) {
      $specs = $in->{specs} ;
   } else {
      $specs = alnmnr::getspecs() ;
   }

# mm = view-volume * $scaling
   my $scaling = $specs->{view_volume}->{$in->{age}}->{spacing} / 1000;

   my $out_fh ;
   if (exists $in->{out_fn}) {
      open($out_fh, ">$in->{out_fn}") ;
   } else {
      open($out_fh, ">-") ;
   }

   foreach my $xyz (keys %{$in->{roidef}->{point2roi}}) {
      my ($x, $y, $z) = split(',', $xyz) ;

# Convert to mm coordinate system
      my $x2 = sprintf("%.2f", $x * $scaling) ;
      my $y2 = sprintf("%.2f", $y * $scaling) ;
      my $z2 = sprintf("%.2f", $z * $scaling) ;
      foreach my $roi (keys %{$in->{roidef}->{point2roi}->{$xyz}}) {
         print {$out_fh} join("\t", $roi, "", "point", "$x2,$y2,$z2")."\n" ; }
   }

   close($out_fh) ;

}


sub display_roidef_pdb {

   my $in = shift ;
   my $out_fh ;
   if (exists $in->{out_fn}) {
      open($out_fh, ">$in->{out_fn}") ;
   } else {
      open($out_fh, ">-") ;
   }

   my @chains = (1 .. 9);
   push @chains, ('a' .. 'z'), ('A' .. 'Z');
   push @chains, ('a' .. 'z'), ('A' .. 'Z');
   push @chains, ('a' .. 'z'), ('A' .. 'Z');
   my $atomno = 1; my $chainid = shift @chains ;
   my $roi2num = {} ; my @roi_names ;
   foreach my $xyz (keys %{$in->{roidef}->{point2roi}}) {
      my $resno = floor($atomno / 10) ;
      foreach my $roi (keys %{$in->{roidef}->{point2roi}->{$xyz}}) {
         if (!exists $roi2num->{$roi}) {
            push @roi_names, $roi ;
            $roi2num->{$roi} = $#roi_names ;
         }
         my ($x, $y, $z) = split(/\,/, $xyz) ;
         print {$out_fh} "ATOM  ".
            sprintf("%5d", $atomno),
            " CA ",
            " ",
            "GLY",
            ' ',
            $chainid,
            sprintf("%4d", $resno),
            " ",
            '   ',
            sprintf("%8.3f",$x),
            sprintf("%8.3f",$y),
            sprintf("%8.3f",$z),
            sprintf("%6.2f",$roi2num->{$roi}),
            sprintf("%6.2f",$roi2num->{$roi}),
            '           ',
            'CA',
            $roi2num->{$roi}."\n" ;

         if ($atomno > 99999) {
               $atomno = 1 ;
               $chainid = shift @chains ;
         }
      }
      $atomno++ ;
   }

   close($out_fh) ;

}


=head2 downsample_aba_Annotation()

   Title:       load_aba_Annotation()
   Function:    loads the ABA voxel - brain structure annotation mapping
                 from the developing or adult reference atlases
   Args:        $_->{fn} [optional - if not will call set_locations()]
   Returns:     $_->{}->{point2strx}->{x,y,z} = brain_strx_id ;
                $_->{}->{strx2points}->{brain_strx_id}->[i] = "x_i,y_i,z_i"
   Files in:    AtlasAnnotation100.sva (adult) or Annotation (developing)

=cut

sub downsample_aba_Annotation {

   my $in = shift ;
   my $out_fn = $in->{out_fn};
   my $specs ;
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }
   my $age = $in->{age} ;

   my $data = {};
      my $fn ;
      if ($age eq 'P56') {
         $fn = $specs->{allen_adult_atlas_dir}."/Spaces/P56/Annotation" ;
      } else {
         $fn = $specs->{allen_devel_atlas_dir}."/Spaces/$age/Annotation" ;
      }
      if (!-s $fn) {
         die "Can't find Annotation file for $age, expecting: $fn" ;}

      my $xdim = $specs->{atlas_dims}->{$age}->{x} ;
      my $ydim = $specs->{atlas_dims}->{$age}->{y} ;
      my $zdim = $specs->{atlas_dims}->{$age}->{z} ;
      my $xy_plane_size = $xdim * $ydim ;

      my $downsample = 1 ;
      if (exists $specs->{atlas_downsample}->{$age}) {
         $downsample = $specs->{atlas_downsample}->{$age}->{from} /
                       $specs->{atlas_downsample}->{$age}->{to} ; }

      my $unpacktype = "C" ; #developing atlases have max 232 entries.
      my $read_length = 1 ;
      if ($age eq 'P56') { #adult atlas has more entries, max Id range to 10737
         $unpacktype = "S";
         $read_length = 2 ;
      }

      my $buffer ;
      my $entry = 0 ;
      open(ANNF, $fn);
      while(read(ANNF, $buffer, $read_length)) {
         $entry++ ;
         my $point_num = $entry ;
         my $z = POSIX::floor($point_num / $xy_plane_size) ;
         my $y = POSIX::floor(($point_num - $z * $xy_plane_size) / $xdim) ;
         my $x = $point_num - $z * $xy_plane_size - $y * $xdim - 1 ;

         $x = POSIX::floor($downsample * $x) ;
         $y = POSIX::floor($downsample * $y) ;
         $z = POSIX::floor($downsample * $z) ;

         my $strx = unpack($unpacktype, $buffer) ;

         if ($strx == 0) {next;}

         $data->{point2strx}->{"$x,$y,$z"} = $strx ;
         $data->{strx2points}->{$strx}->{"$x,$y,$z"}++ ;
      }
      close(ANNF) ;

# NOTE 111122_2005 - danger now is if down-sampling results in a single new-coord-system point being associated with multiple sructures... though luckily {point2strx} doesn't seem to be used...

   open(OUTF, ">".$out_fn) ;
   foreach my $strx (keys %{$data->{strx2points}}) {
      foreach my $point (keys %{$data->{strx2points}->{$strx}}) {
         print OUTF $strx."\t".$point."\n"; } }
   close(OUTF) ;

   return $data ;
}


=head2 load_aba_brainstructures()

   Title:       load_aba_brainstructures()
   Function:    loads the ABA's brainstructure.csv file
   Args:        $_->{fn} [optional - if not will call set_locations()]
   Returns:     $_->{abbreviation2ID}->{abbreviation} = informaticsId ;
                $_->{ID2info}->{informaticsId} = {
                  StructureName =>,
                  Abbreviation =>,
                  ParentStruct =>,
                  red => ,
                  green =>,
                  blue =>,
                  informaticsId =>,
                  StructureId =>,
                }
   Files in:    brainstructures.csv

=cut

sub load_aba_brainstructures {

   my $in = shift ;
   my $specs ;
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   my $fn ;
   if (exists $in->{fn}) {
      $fn = $in->{fn} ;
   } else {
      if ($in->{age} eq 'P56') {
#fpd111122_0644          $fn = $specs->{fn}->{brainstructures} ;
         $fn = $specs->{allen_adult_atlas_dir}.'/ontology.csv' ;
      } else {
         $fn = $specs->{allen_devel_atlas_dir}.'/ontology.csv' ;
      }
   }

   if (!-s $fn) {
      die "Can't find brainstructures.csv, expecting: ".$fn ; }

   open(INF, $fn) ;

   my @fields ;
   my $f2i ;
#   if ($in->{age} eq 'adult') { #brainstructures.csv has a header line
#      my $header = <INF> ;
#      @fields = split(/\,/ , $header) ;
#   } else { #These are my inferred headers for ontology.csv
      @fields = qw/StructureName Abbreviation ParentStruct red green blue informaticsId unk1 unk2 ontology_level/ ;
#   }

   foreach my $j ( 0 .. $#fields) {
      $f2i->{$fields[$j]} = $j; }

   my $data = {};
   while (my $line = <INF>) {
      chomp $line;
      $line =~ s/\, /\; /g ;
      $line =~ s/"//g ;
      my @t = split(/\,/, $line) ; 

      $data->{abbreviation2ID}->{$t[$f2i->{Abbreviation}]} =
         $t[$f2i->{informaticsId}] ;

      foreach my $f (keys %{$f2i}) {
         $data->{ID2info}->{$t[$f2i->{informaticsId}]}->{$f}= $t[$f2i->{$f}] ;}
   }
   close(INF) ;

   $data->{tree} = {} ;
   foreach my $node (keys %{$data->{ID2info}}) {
      if ($data->{ID2info}->{$node}->{ParentStruct} eq '') {next;}
      my $parent_name = $data->{ID2info}->{$node}->{ParentStruct} ;
      my $parent_id = $data->{abbreviation2ID}->{$parent_name} ;
      $data->{tree}->{kids}->{$parent_id}->{$node}++ ;
      $data->{tree}->{parent}->{$node}->{$parent_id}++ ;
   }

   return $data ;

}


=head2 get_brainstrx_voxels()

   Title:       get_brainstrx_voxels()
   Function:    Returns list of voxels assigned to a specific ABA brain region
   Args:        $_->{brainstrx}  either the informatics_ID or abbreviation for a
                  region in the ABA reference atlas

                $_-{onpoints_only} if set to 1, only returns on points,
                                    otherwise, by default returns off points

   Returns:     $_->{x,y,z} = 1  if point (x,y,z) is in the specified ROI
                              -1 if point is specified OFF in the ROI

=cut

sub get_brainstrx_voxels {

   my $in = shift ;
   my $brainstrx = $in->{brainstrx} ;
   my $age = $in->{age} ;
   my $specs = $in->{specs} ;

   my $strx2status ;
   if ($brainstrx =~ /\,/ || $brainstrx =~ /\=on/ || $brainstrx =~ /\=off/) {
      my @parts = split(/\,/, $brainstrx) ;
      foreach my $j ( 0.. $#parts) {
         my ($strx, $status) = ($parts[$j] =~ /(.+)=(.+)/) ;
         $strx2status->{$strx} = $status ;
      }
   } else {
      $strx2status->{$brainstrx} = 'on' ;
   }

   my $atlas_annotations = {};
   if (!exists $specs->{parsed_annotation}) {
      $atlas_annotations = load_aba_Annotation({age => $age}) ;
   } else {
      $atlas_annotations = $specs->{parsed_annotation} ;
   }

   my $strx2points ;
   foreach my $brainstrx (keys %{$strx2status}) {
      my $brainstrx_info = get_brainstrx_info({brainstrx => $brainstrx,
                                               age => $age });
      my $brainstrx_id = $brainstrx_info->{brainstrx_id} ;

      my $brainstrx_points ;
      push @{$brainstrx_points},
         keys %{$atlas_annotations->{strx2points}->{$brainstrx_id}} ;
      foreach my $node (@{$brainstrx_info->{kids}}) {
         push @{$brainstrx_points},
            keys %{$atlas_annotations->{strx2points}->{$node}};}

      $strx2points->{$brainstrx} = $brainstrx_points ;
   }

   my $point_light ;
   foreach my $brainstrx (keys %{$strx2status}) {
      foreach my $point (@{$strx2points->{$brainstrx}}) {
         if (!exists $point_light->{$point}) {
            if ($strx2status->{$brainstrx} eq 'on') {
               $point_light->{$point} = 1 ;
            } else {
               $point_light->{$point} = 0 ;
            }
         } else {
            if ($strx2status->{$brainstrx} eq 'off') {
               $point_light->{$point} = 0 ;
            }
         }
      }
   }

# NOTE: If developing atlas, have to convert between Annotation (tissue-volume)
#       and XPZ (view-volume) coordinate systems
   my $conversion_factor = 1 ;
   if ($age ne 'P56') {
      $conversion_factor = $specs->{dims}->{$age}->{spacing} /
                           $specs->{view_volume}->{$age}->{spacing} ;
   }

   my $results ;
   foreach my $point (keys %{$point_light}) {

      my ($x, $y, $z) = split(',', $point) ;
      my $newpoint = $point ;
      if ($age ne 'P56') {
         $x = POSIX::floor($x * $conversion_factor) ;
         $y = POSIX::floor($y * $conversion_factor) ;
         $z = POSIX::floor($z * $conversion_factor) ;
         $newpoint = $x.','.$y.','.$z ;
      }

      if ($point_light->{$point} == 1) {
         $results->{$newpoint}++ ;
      } else {
         if (!exists $in->{onpoints_only} ||
             $in->{onpoints_only} == 0) {
            $results->{$newpoint} = -1 ;
         }
      }
   }

   return $results ;
}


=head2 get_brainstrx_info()

   Title:       get_brainstrx_info()
   Function:    loads annotation for voxel - brain structure mapping
   Args:        $_->{brainstrx}  either the informatics_ID or abbreviation for a
                  region in the ABA reference atlas
   Returns:     $_->{brainstrx_id} - ABA informatics ID of the brain region
                $_->{kids} = [kid1, kid2,...] - arrayref of ABA informatics IDs
                   for structures that are children of the specified brainstrx
   Files in:    AtlasAnnotation100.sva

=cut

sub get_brainstrx_info {

   my $in = shift ;
   my $brainstrx = load_aba_brainstructures({age => $in->{age}}) ;

   my $brainstrx_id ;
   if ($in->{brainstrx} =~ /[A-Za-z]/) {
      $brainstrx_id = $brainstrx->{abbreviation2ID}->{$in->{brainstrx}} ;
   } else {
      $brainstrx_id = $in->{brainstrx} ;
   }

#fpd101107_1612  - not used.
#   my $atlas_annotations = load_aba_AtlasAnnotation100() ;

   my $brainstrx_points ;

   my @kids = get_kids({ # also need to get points from kids.
      node => $brainstrx_id,
      tree => $brainstrx->{tree}
   }) ;

   return  {
      brainstrx_id => $brainstrx_id,
      kids => \@kids
   } ;

}


=head2 load_aba_Annotation()

   Title:       load_aba_Annotation()
   Function:    loads the ABA voxel - brain structure annotation mapping
                 from the developing or adult reference atlases
   Args:        $_->{fn} [optional - if not will call set_locations()]
   Returns:     $_->{}->{point2strx}->{x,y,z} = brain_strx_id ;
                $_->{}->{strx2points}->{brain_strx_id}->[i] = "x_i,y_i,z_i"
   Files in:    AtlasAnnotation100.sva (adult) or Annotation (developing)

=cut

sub load_aba_Annotation {

   my $in = shift ;
   my $specs ;
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }
   my $age = $in->{age} ;

   my $data = {};
   my $fn ;
   if ($age eq 'P56') {
      $fn = $specs->{allen_adult_atlas_dir}."/Spaces/P56/Annotation" ;
   } else {
      $fn = $specs->{allen_devel_atlas_dir}."/Spaces/$age/Annotation" ;
   }
   if (!-s $fn) {
      die "Can't find Annotation file for $age, expecting: $fn" ;}

   if (exists $specs->{atlas_downsample}->{$age}) {
      my $downsample_fn = $fn.'.allenminer2-downsample.txt';
      if (!-s $downsample_fn) {
         print STDERR "Note: Must downsample atlas Annotation file, ".
                      " one-time ~10min process: ";
         downsample_aba_Annotation({
            age         => $age,
            specs       => $specs,
            orig_fn     => $fn,
            out_fn      => $downsample_fn,
         }) ;
         print STDERR " Done\n";
      }

      open(DOWNSAMPLEF, $downsample_fn) ;
      while (my $line = <DOWNSAMPLEF>) {
         chomp $line;
         my ($strx, $point) = split(/\t/, $line);
         $data->{point2strx}->{$point}->{$strx}++ ;
         $data->{strx2points}->{$strx}->{$point}++ ;
      }
      close(DOWNSAMPLEF) ;

   } else {

      my $xdim = $specs->{dims}->{$age}->{x} ;
      my $ydim = $specs->{dims}->{$age}->{y} ;
      my $zdim = $specs->{dims}->{$age}->{z} ;
      my $xy_plane_size = $xdim * $ydim ;

      my $buffer ;
      my $entry = 0 ;
      open(ANNF, $fn);
      while(read(ANNF, $buffer, 1)) {
         $entry++ ;
         my $point_num = $entry ;
         my $z = POSIX::floor($point_num / $xy_plane_size) ;
         my $y = POSIX::floor(($point_num - $z * $xy_plane_size) / $xdim) ;
         my $x = $point_num - $z * $xy_plane_size - $y * $xdim - 1 ;

         my $strx = unpack("C", $buffer) ;
         if ($strx == 0) {next;}

         $data->{point2strx}->{"$x,$y,$z"} = $strx ;
         $data->{strx2points}->{$strx}->{"$x,$y,$z"}++ ;
      }
      close(ANNF) ;
   }
   return $data ;
}



sub get_kids {

   my $in = shift ;
   my @kids ;

   if (!exists $in->{tree}->{kids}->{$in->{node}}) {
      return () ;
   } ;
   foreach my $kid (keys %{$in->{tree}->{kids}->{$in->{node}}}) {
      push @kids, $kid ;
      push @kids, get_kids({node => $kid, tree => $in->{tree}}) ;
   }

   return @kids ;
}


=head2 readin_1col_arr()

   Title:       readin_1col_arr()
   Function:    Reads in a single column and returns an arrayref of the values
   Args:        $_->{fn} = input filename
   Returns:     $_->[i] = value of the ith lines in the file

=cut

sub readin_1col_arr {

   my $in = shift;
   my $fn = $in->{fn} ;

   open (FH, $in->{fn}) ;
   my $results; 
   while (my $line = <FH>) {
      chomp $line;
      push @{$results}, $line ;
   }

   return $results;
}


sub parse_image_series_xml {

   my $in       = shift;
   my $data     = $in->{data} ;
   my $fn       = $in->{fn} ;


   my $image_keys = {} ;
   map {$image_keys->{$_}++;}
      (qw/id imagepath expression-imagepath height width left top/,
       qw/parentimageheight parentimagewidth/,
       qw/specimen-tissue-index section-index/,
       qw/tsv tvs z-value/) ;


   open(XMLINFO, $fn) ;
# <section-thickness> only in adult series...
   while (my $line = <XMLINFO>) {
      chomp $line;
      if ($line =~ /^  <id/) {
         ($data->{meta}->{imageseries_id}) = ($line =~ /<id>([0-9]+)<\/id>/) ;
      } elsif ($line =~ /<plane-of-section>/) {
         ($data->{meta}->{plane}) = ($line =~ /<plane-of-section>(.+)<\/plane-of-section>/) ;
      } elsif ($line =~ /<section-thickness>/) {
         ($data->{meta}->{"section-thickness"}) =
            ($line =~ /<section-thickness>(.+)<\/section-thickness>/) ;
      } elsif ($line =~ /<age>/) {
         ($data->{meta}->{stage}) = ($line =~ /<age>(.+)<\/age>/) ;
      } elsif ($line =~ /<trv>/) {
         ($data->{meta}->{trv}) = ($line =~ /<trv>(.+)<\/trv>/) ;
      } elsif ($line =~ /<tvr>/) {
         ($data->{meta}->{tvr}) = ($line =~ /<tvr>(.+)<\/tvr>/) ;
      } elsif ($line =~ /<gene-symbol>/) {
         ($data->{meta}->{gene_name}) = ($line =~ /<gene-symbol>(.+)<\/gene-symbol>/) ;
#      } elsif ($line =~ /<\/gene>$/) {

      } elsif ($line =~ /^    <image>/) {

         my $curimage ;
         while ($line = <XMLINFO>) {
            if ($line =~ /^    <\/image>/) {last;}
            chomp $line;
            my ($key, $val) = ($line =~ /<([^>]+)>(.+)<\//) ;
            if (!defined $key || !exists $image_keys->{$key}) {next;}
            $curimage->{$key} = $val ;
         }
         $data->{meta}->{images}->{$curimage->{"z-value"}} = $curimage ;

      } elsif ($line =~ /<\/images>$/) {
         last;
      }
   }
# now read in image data
   close(XMLINFO) ;

   return ;

}



#fpd111219_1114 
#sub readbitvec_roi_list_results {
#
#   my $in = shift;
#   my $f2i ;
#   foreach my $j ( 0 .. $#{$in->{headers}->{headers}}) {
#      $f2i->{$in->{headers}->{headers}->[$j]} = $j ; }
#
## 1. read through first gene to get roi names, keep levels in buffer
## 2. convert first gene to bit vector
## 3. read ahead and assign bit vectors on the fly
## 4. send back absence/presence bit vectors
#
#   open(RESULTSF, $in->{fn}) ;
#   my $roi_data ;
#   while (my $line = <RESULTSF>) {
#      chomp $line;
#      if ($line =~ /^#/) {next;}
#      my $entry ;
#      my @t = split(/\t/, $line) ;
#      foreach my $j (0 .. $#t) {
#         $entry->{$in->{headers}->{headers}->[$j]} = $t[$j] ;
#      }
#      if ($entry->{plane} ne 'sagittal' &&
#          $entry->{plane} ne 'coronal') {
#         next; }
#      foreach my $f (sort keys %{$f2i}) {
#         push @{$roi_data->{$entry->{roi_name}}->{$entry->{plane}}->{$f}},
#            $entry->{$f} ; }
#   }
#   close(RESULTSF) ;
#   return $roi_data ;
#
#}


1;
