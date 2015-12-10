=head1 NAME

alnmnr::download.pm - routines to download ABA data

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

package alnmnr::download ;
use strict;
use warnings;
use Cwd;
use File::Path qw/mkpath/ ;
use File::Temp qw/tempfile tempdir/ ;
use POSIX qw/floor ceil/ ;

sub _run_download_data {

   my $in = shift ;
   my $age = $in->{age} ;
   my $roidef = $in->{roidef} ;
   my $specs = $in->{specs} ;
   if (!exists $in->{specs}) {
      $specs = getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   print STDERR "Your data directory is currently set to: ".
      $specs->{allen_data_dir}."\n\n".
      "To change this, exit and type the full path to your desired ".
      "allen_data_dir on line 65 of alnmnr.pm\n\n".
      "To continue with the current path, type Y: " ;

   my $answer = <STDIN> ; chomp $answer;
   if ($answer ne 'Y' && $answer ne 'y') {die "\n";}

   if ($in->{atlas} eq 'adultbrain') {
      download_atlas({specs => $specs, atlas_type => 'adult'}) ;

      make_adult_xpz_list({specs => $specs}) ;

      download_adultbrain_data({specs => $specs}) ;

   } elsif ($in->{atlas} eq 'develbrain') {

      download_atlas({specs => $specs, atlas_type => 'devel'}) ;
      download_develbrain_data({specs => $specs}) ;

   } elsif ($in->{atlas} eq 'spinalcord') {

      alnmnr::spinalcord::download_spinalcord_data({specs => $specs}) ;

   } elsif ($in->{atlas} eq 'fastsearch') {

      download_fastsearch_data({specs => $specs}) ;

   } else {

      die "must specify -atlas adultbrain|develbrain|spinalcord|fastsearch" ;

   }

}


sub download_atlas {

   my $in = shift ;
   my $type = $in->{atlas_type} ;

   my $specs = $in->{specs} ;
   if (!exists $in->{specs}) {
      $specs = getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   my $atlas_dir = $specs->{"allen_".$type."_atlas_dir"} ;
   if (! -s $atlas_dir) { File::Path::mkpath($atlas_dir); }
   my $url = $specs->{download}->{URL}->{$type."_atlas"} ;
   my $atlas_fn = $url; $atlas_fn =~ s/.*\/// ;

   if (-s "$atlas_dir/$atlas_fn") {
      print STDERR "NOTE: $type atlas already installed; ".
         "To reinstall, remove $atlas_dir/$atlas_fn and rerun\n";
      return;
   }

   print STDERR "Installing $type atlas:\n";
   print STDERR "  Downloading: ";
   my $atlas_com = "wget -q \'$url\' -O $atlas_dir/$atlas_fn" ;
   system($atlas_com) ;
   print STDERR "X\n";

   my $curdir = getcwd() ;

   chdir $atlas_dir ;
   print STDERR "  Uncompressing: ";
   system("unzip -q $atlas_fn") ;
   print STDERR "X\n" ;
   chdir $curdir ;

}

sub download_fastsearch_data {
# download adult and developing expression summaries

   my $in = shift ;
   my $specs = $in->{specs} ;
   if (!exists $in->{specs}) {
      $specs = getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   if (! -s $specs->{allen_fastsearch_dir}) {
      File::Path::mkpath($specs->{allen_fastsearch_dir}); }

   print STDERR "Retrieving fastsearch files:\n";
   foreach my $age (keys %{$specs->{atlas_dims}}) {
      print STDERR "   age $age: \n";
      if (-s $specs->{fastsearch_fn}->{$age}) {
         print STDERR "Fastquery file for $age already installed; ".
            "To reinstall, remove ".$specs->{fastsearch_fn}->{$age}.
            "and rerun\n";
         next;
      }
      my $url = $specs->{download}->{URL}->{fastsearch_prefix}."$age.txt.gz" ;
      my $wget_com = "wget -q \'$url\' -O ".$specs->{fastsearch_fn}->{$age} ;
      system($wget_com) ;
      print STDERR "X\n";
   }

}


sub download_develbrain_data {

   my $in = shift ;
   my $specs = $in->{specs} ;
   if (!exists $in->{specs}) {
      $specs = getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }


   my $curdir = getcwd() ;
   my $tempdir = tempdir() ;
   chdir $tempdir ;

   my $ages = {}; map {$ages->{$_} = 1;} qw/E11.5 E13.5 E15.5 E18.5 P4 P14 P28/;

#need to search for terms explicitly
   my @terms ; push @terms, 0 .. 9; push @terms, 'A' .. 'Z' ;

   print STDERR "Retrieving developing brain expression files:\n";
   print STDERR "   Getting gene list: ";
   my $gene_info = {};
   foreach my $term (@terms) {
      print STDERR "$term" ;
      my ($out_fh, $out_fn) ;
      my $num_hits;

# Have to switch to adult atlas-like per-term search, grep number of rows, then
#  grab them chunk by chunk.

      my $initsearch_url =
         $specs->{download}->{URL}->{devel_search_part1}."*$term*".
         $specs->{download}->{URL}->{devel_search_part2}."0".
         $specs->{download}->{URL}->{devel_search_part3}."5".
         $specs->{download}->{URL}->{devel_search_part4} ;

      {
         ($out_fh->{xml}, $out_fn->{xml}) = tempfile() ; close($out_fh->{xml}) ;
         my $tcom = "wget -q '$initsearch_url' -O ".$out_fn->{xml} ;
         system($tcom) ;
         open(P1, $out_fn->{xml}) ;
         while (my $line = <P1>) {
            if ($line =~ /total_rows/) {
               chomp $line;
               ($num_hits) = ($line =~ /total_rows='([0-9]+)'/) ;
               last;
            }
         }
         close(P1) ;

         unlink $out_fn->{xml} ;
      }



      if ($num_hits == 0) { next;}
      my $num_chunks = POSIX::ceil($num_hits /
                        $specs->{download}->{devel_xpz_search_chunk}) ;

      foreach my $j ( 0 .. $num_chunks) {
         my $startrow = $j * $specs->{download}->{devel_xpz_search_chunk} ;
         my $search_url =
               $specs->{download}->{URL}->{devel_search_part1}."*$term*".
               $specs->{download}->{URL}->{devel_search_part2}.$startrow.
               $specs->{download}->{URL}->{devel_search_part3}.
                  $specs->{download}->{devel_xpz_search_chunk}.
               $specs->{download}->{URL}->{devel_search_part4} ;

         my ($t_fh, $t_fn) = tempfile() ; close($t_fh) ;
         my $tcom = "wget -q \'$search_url\' -O $t_fn" ;
         system($tcom) ;

         if (!-s $t_fn) {
            print STDERR "ERROR: couldn't download developmental gene ".
                         "list XML file: $tcom\n" ;
            die ;
         }

         aba_parse_devel_xml({ xml_fn    => $t_fn,
                               gene_info => $gene_info}) ;
         unlink $t_fn ;
      }
   }
   chdir $curdir ;
   rmdir $tempdir ;
   print STDERR ". X\n";


   my ($out_fh, $out_fn) ;
   ($out_fh->{devel_wget}, $out_fn->{devel_wget}) =
      tempfile("allen_devel_wget.XXXXX", SUFFIX => '.out') ;

   print STDERR "   Building list of expression file URLs: ";
   foreach my $gene (keys %{$gene_info}) {
      foreach my $imageseries_id ( keys %{$gene_info->{$gene}->{imageseries}}) {
         my $imageseries =
            $gene_info->{$gene}->{imageseries}->{$imageseries_id} ;
         if (!exists $ages->{$imageseries->{age}}) {next;}
         my $xpz_dir = $specs->{allen_xpz_dir}."/".$imageseries->{age}.'/'.
                       substr($gene,0,1).'/'.substr($gene,0,3) ;
         if (! -s $xpz_dir) { File::Path::mkpath($xpz_dir); }

         my $url = $specs->{download}->{URL}->{devel_xpz}.$imageseries->{id} ;
         my $xpz_fn = $gene.'_'.
                      $imageseries->{age}.'_'.
                      $imageseries->{"plane-of-section"}.'_'.
                      $imageseries_id ;
         if ($imageseries->{failed} ne 'false') { $xpz_fn .= '_failed'; }
         $xpz_fn .= '.xpz' ;
                      
         my $xpz_command = "wget -q \'$url\' -O \'$xpz_dir/$xpz_fn\'" ;
         print {$out_fh->{devel_wget}} $xpz_command."\n" ;
      }
   }
   close($out_fh->{devel_wget}) ;
   print STDERR "X\n";

   my ($wget_err_fh, $wget_err_fn) =
      tempfile("devel_wget.XXXXX", SUFFIX => ".err") ;
      close($wget_err_fh) ;
   my ($wget_out_fh, $wget_out_fn) =
      tempfile("devel_wget.XXXXX", SUFFIX => ".out") ;
      close($wget_out_fh) ;

   print STDERR "  About to start download (commands in ".
                 $out_fn->{devel_wget}.") - OK? (Y/n) " ;
   my $answer = <STDIN> ; chomp $answer;
   if ($answer ne '' && $answer ne 'Y' && $answer ne 'y') {die "\n";}

   print STDERR "  Downloading files: ";
   my $download_com = "bash ".$out_fn->{devel_wget}.
                      " 2> $wget_err_fn >$wget_out_fn" ;
   system($download_com) ;
   print STDERR "X\n";

   ($out_fh->{missing_files}, $out_fn->{missing_files}) =
      tempfile("allen_devel_wget_MISSING_FILES.XXXXX", SUFFIX => '.err') ;

# Check that the expected files were downloaded
   open(WGETF, $out_fn->{devel_wget}) ;
   my $num_missing = 0 ;
   while (my $line = <WGETF>) {
      chomp $line;
      my ($filepath) = ($line =~ /.*-O '(.+)'$/);
      if (-s $filepath) {next;}
      print {$out_fh->{missing_files}} $line."\n" ;
      $num_missing++ ;
   }
   close(WGETF) ;
   close($out_fh->{missing_files}) ;

   if ($num_missing == 0) {
      print STDERR "All files were downloaded\n";
      unlink $out_fn->{missing_files} ;
   } else {
      print STDERR "Warning: $num_missing files did not download\n";
      "Retry by running: bash ".$out_fn->{missing_files}."\n" ;
   }

}


=head2 aba_parse_devel_xml()

   Title:       aba_parse_devel_xml()
   Function:    Parses developmental ABA XML files to get gene/stage/imageseries
   Args:        ->{xml_fn} = XML filenmae 
                ->{gene_info} = hash to hold output
   Returns:     Nothing
   Output:      $in->{gene_info}->{gene-symbol}->{imageseries}->[i] =
                  { imageseries_id =>
                    age            =>
                    slices         => }

=cut

sub aba_parse_devel_xml {

   my $in = shift;
   my $xml_fn = $in->{xml_fn} ;         #XML file to parse
   my $gene_info = $in->{gene_info} ;   #->{gene_symbol}->{imageseries}->[i]=
                                        # {imageseries_id|age|slices}
                                        #
   my $ageid2text = {
      1 => 'E18.5',
      2 => 'E15.5',
      4 => 'E13.5',
      5 => 'E11.5',
      7 => 'P4',
      11 => 'P14',
      14 => 'P28',
      15 => 'P56'
   } ;

   my $planeid2text = {
      1 => 'coronal',
      2 => 'sagittal'
   } ;

   my $cur_gene = {};
   my $cur_imageseries = {};

   open(SEARCHRESF, $in->{xml_fn}) ;
   while (my $line = <SEARCHRESF>) {
      chomp $line;

      if ($line =~ /^ {4}<id>/) {                                       #GOOD
         ($cur_gene->{id}) = ($line =~ /^    <id>(.+)<\/id>$/) ;

      } elsif ($line =~ /^ {4}<name>/) {                                #GOOD
         ($cur_gene->{name}) = ($line =~ /<name>(.+)<\/name>$/) ;

      } elsif ($line =~ /^ {4}<entrez-id>/) {                           #GOOD
         ($cur_gene->{"entrez-id"}) =
            ($line =~ /<entrez-id>(.+)<\/entrez-id>$/) ;

      } elsif ($line =~ /^ {4}<acronym>/) {                             #GOOD
         ($cur_gene->{"gene-symbol"}) =
            ($line =~ /<acronym>(.+)<\/acronym>$/) ;

      } elsif ($line =~ /^ {8}<id type="integer">/) {                   #GOOD
         ($cur_imageseries->{id}) =
            ($line =~ /<id type="integer">(.+)<\/id>$/) ;

      } elsif ($line =~ /^ {8}<plane-of-section-id type="integer">/) { #GOOD
         my ($plane_id) =
            ($line =~ /<plane-of-section-id type="integer">(.+)<\/plane-of-section-id>$/) ;

         ($cur_imageseries->{"plane-of-section"}) = $planeid2text->{$plane_id};

      } elsif ($line =~ /^ {12}<age-id type="integer">/) {              #GOOD
         my ($age_id) =
            ($line =~ /<age-id type="integer">(.+)<\/age-id>$/) ;
         ($cur_imageseries->{"age"}) = $ageid2text->{$age_id};

      } elsif ($line =~ /^ {8}<failed type="boolean">/) {               #GOOD
         ($cur_imageseries->{"failed"}) =
            ($line =~ /<failed type="boolean">(.+)<\/failed>$/) ;

      } elsif ($line =~ /^      <\/data-set>$/) {                       #GOOD
         if (!exists $gene_info->{$cur_gene->{"gene-symbol"}}) {
            $gene_info->{$cur_gene->{"gene-symbol"}}->{imageseries} = {} ;
            foreach my $gene_feat (keys %{$cur_gene}) {
               $gene_info->{$cur_gene->{"gene-symbol"}}->{$gene_feat} =
                  $cur_gene->{$gene_feat} ; } }

         if (!exists $gene_info->{$cur_gene->{"gene-symbol"}}->{imageseries}->{$cur_imageseries->{id}}) {
            foreach my $image_series_feat (keys %{$cur_imageseries}) {
$gene_info->{$cur_gene->{"gene-symbol"}}->{imageseries}->{$cur_imageseries->{id}}->{$image_series_feat}=
                  $cur_imageseries->{$image_series_feat} ;
            }
         }

         $cur_imageseries= {} ;

      } elsif ($line =~ /  <\/gene>$/) {
         $cur_gene = {} ;
      }

   }
   close(SEARCHRESF) ;

   return ;

}


=head2 PRE20150514_aba_parse_devel_xml()

   Title:       PRE20150514_aba_parse_devel_xml() -- DEPRECATED
   Function:    Parses developmental ABA XML files to get gene/stage/imageseries
   Args:        ->{xml_fn} = XML filenmae 
                ->{gene_info} = hash to hold output
   Returns:     Nothing
   Output:      $in->{gene_info}->{gene-symbol}->{imageseries}->[i] =
                  { imageseries_id =>
                    age            =>
                    slices         => }

=cut

sub PRE20150514_aba_parse_devel_xml {

   my $in = shift;
   my $xml_fn = $in->{xml_fn} ;         #XML file to parse
   my $gene_info = $in->{gene_info} ;   #->{gene_symbol}->{imageseries}->[i]=
                                        # {imageseries_id|age|slices}
   my $cur_gene = {};
   my $cur_imageseries = {};

   open(SEARCHRESF, $in->{xml_fn}) ;
   while (my $line = <SEARCHRESF>) {
      chomp $line;

      if ($line =~ /^    <id>/) {
         ($cur_gene->{id}) = ($line =~ /^    <id>(.+)<\/id>$/) ;
      } elsif ($line =~ /^    <name>/) {
         ($cur_gene->{name}) = ($line =~ /^    <name>(.+)<\/name>$/) ;
      } elsif ($line =~ /^    <entrez-id>/) {
         ($cur_gene->{"entrez-id"}) =
            ($line =~ /^    <entrez-id>(.+)<\/entrez-id>$/) ;
      } elsif ($line =~ /^    <gene-symbol>/) {
         ($cur_gene->{"gene-symbol"}) =
            ($line =~ /^    <gene-symbol>(.+)<\/gene-symbol>$/) ;
      } elsif ($line =~ /^        <id>/) {
         ($cur_imageseries->{id}) =
            ($line =~ /<id>(.+)<\/id>$/) ;
      } elsif ($line =~ /^        <plane-of-section>/) {
         ($cur_imageseries->{"plane-of-section"}) =
            ($line =~ /<plane-of-section>(.+)<\/plane-of-section>$/) ;
      } elsif ($line =~ /^        <age>/) {
         ($cur_imageseries->{"age"}) =
            ($line =~ /<age>(.+)<\/age>$/) ;
      } elsif ($line =~ /^        <failed>/) {
         ($cur_imageseries->{"failed"}) =
            ($line =~ /<failed>(.+)<\/failed>$/) ;

      } elsif ($line =~ /^      <\/image-series>$/) {
         if (!exists $gene_info->{$cur_gene->{"gene-symbol"}}) {
            $gene_info->{$cur_gene->{"gene-symbol"}}->{imageseries} = {} ;
            foreach my $gene_feat (keys %{$cur_gene}) {
               $gene_info->{$cur_gene->{"gene-symbol"}}->{$gene_feat} =
                  $cur_gene->{$gene_feat} ; } }

         if (!exists $gene_info->{$cur_gene->{"gene-symbol"}}->{imageseries}->{$cur_imageseries->{id}}) {
            foreach my $image_series_feat (keys %{$cur_imageseries}) {
$gene_info->{$cur_gene->{"gene-symbol"}}->{imageseries}->{$cur_imageseries->{id}}->{$image_series_feat}=
                  $cur_imageseries->{$image_series_feat} ;
            }
         }

         $cur_imageseries= {} ;

      } elsif ($line =~ /  <\/gene>$/) {
         $cur_gene = {} ;
      }

   }
   close(SEARCHRESF) ;

   return ;

}


sub make_adult_xpz_list {

   my $in = shift ;
   my $specs = $in->{specs} ;
   if (!exists $in->{specs}) {
      $specs = getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   my $xpz_list = {};
   print STDERR "Retrieving list of expression files: ";
   print STDERR "." ;
   my $initsearch_url =
      $specs->{download}->{URL}->{adult_search_part1}.
      $specs->{download}->{URL}->{adult_search_part2}.'0'.
      $specs->{download}->{URL}->{adult_search_part3}.'50' ;

   my $num_hits ;
   {
   my ($t_fh, $t_fn) = tempfile() ; close($t_fh) ;
   my $tcom = "wget -q \'$initsearch_url\' -O $t_fn"  ;
   system($tcom) ;

   open(SEARCHRESF, $t_fn) ;
   while (my $line = <SEARCHRESF>) {
      if ($line =~ /total_rows/) {
         chomp $line;
         ($num_hits) = ($line =~ /total_rows='([0-9]+)'/) ;
         last;
      }
   }
   close(SEARCHRESF) ;
   unlink $t_fn ;
   }

   if ($num_hits == 0) { die "ERROR: No expression files found on server";}

   my $num_chunks = POSIX::ceil($num_hits /
                     $specs->{download}->{adult_xpz_search_chunk}) ;
   foreach my $j ( 0 .. $num_chunks) {
      my $startrow = $j * $specs->{download}->{adult_xpz_search_chunk} ;
      my $search_url =
            $specs->{download}->{URL}->{adult_search_part1}.
            $specs->{download}->{URL}->{adult_search_part2}.$startrow.
            $specs->{download}->{URL}->{adult_search_part3}.
            $specs->{download}->{adult_xpz_search_chunk} ;
      my ($t_fh, $t_fn) = tempfile() ; close($t_fh) ;
      my $tcom = "wget -q \'$search_url\' -O $t_fn" ;
      system($tcom) ;

      aba_parse_adult_search_results({ fn => $t_fn, hits => $xpz_list }) ;
      unlink $t_fn ;
   }

   print STDERR "X\n";
   print STDERR "X\n";

   if (!-s $specs->{allen_xpz_dir}) {
      File::Path::mkpath($specs->{allen_xpz_dir}) ; }

   open(OUTF, ">".$specs->{download}->{fn}->{adult_xpz_list}) ;
   my @headers = qw/gene gene-id name entrez-id imageseries-id/;
   push @headers, qw/plane age failed/ ;
   print OUTF '#'.join("\t", @headers)."\n";
   foreach my $gene (sort keys %{$xpz_list}) {
      foreach my $imageseries (sort {$a <=> $b}
            keys %{$xpz_list->{$gene}->{imageseries}}) {

         if ($xpz_list->{$gene}->{imageseries}->{$imageseries}->{failed} eq
             'true') { next; }

         if (!exists $xpz_list->{$gene}->{"entrez-id"}) {
            $xpz_list->{$gene}->{"entrez-id"} = '' ; }

         if (!exists $xpz_list->{$gene}->{imageseries}->{$imageseries}->{age}){
            $xpz_list->{$gene}->{imageseries}->{$imageseries}->{age} = 'P56';}

         my @outvals = ( $gene,
            $xpz_list->{$gene}->{"id"},
            $xpz_list->{$gene}->{"name"},
            $xpz_list->{$gene}->{"entrez-id"},
            $xpz_list->{$gene}->{imageseries}->{$imageseries}->{id},
$xpz_list->{$gene}->{imageseries}->{$imageseries}->{"plane-of-section"},
$xpz_list->{$gene}->{imageseries}->{$imageseries}->{age},
$xpz_list->{$gene}->{imageseries}->{$imageseries}->{failed} ) ;
         print OUTF join("\t", @outvals)."\n" ;
      }
   }
   close(OUTF) ;

}


sub aba_parse_adult_search_results {

   my $in = shift ;

   my $hits ;
   if (exists $in->{hits}) { # if hits hash provided, append.
      $hits = $in->{hits} ;
   } else {
      $hits = {};
   }

   my $planeid2text = {
      1 => 'coronal',
      2 => 'sagittal'
   } ;

   my $cur_gene = {};
   my $cur_imageseries = {};

### example of query response
#<Response success='true' start_row='0' num_rows='50' total_rows='19991'><objects>
#  <object>
#    <acronym>0610005G16Rik*</acronym>
#    <alias-tags nil="true"/>
#    <chromosome-id nil="true"/>
#    <ensembl-id nil="true"/>
#    <entrez-id nil="true"/>
#    <homologene-id nil="true"/>
#    <id>149630</id>
#    <legacy-ensembl-gene-id nil="true"/>
#    <name>RIKEN cDNA 0610005G16 gene (non-RefSeq)</name>
#    <organism-id>2</organism-id>
#    <sphinx-id>158073</sphinx-id>
#    <data-sets>
#      <data-set>
#        <blue-channel nil="true"/>
#        <delegate type="boolean">false</delegate>
#        <expression type="boolean">true</expression>
#        <failed type="boolean">false</failed>
#        <failed-facet type="integer">734881840</failed-facet>
#        <green-channel nil="true"/>
#        <id type="integer">74357755</id>
#        <name nil="true"/>
#        <plane-of-section-id type="integer">2</plane-of-section-id>
#        <qc-date type="datetime">2009-05-02T23:03:44Z</qc-date>
#        <red-channel nil="true"/>
#        <reference-space-id type="integer">10</reference-space-id>
#        <rnaseq-design-id type="integer" nil="true"/>
#        <section-thickness type="float">25.0</section-thickness>
#        <specimen-id type="integer">74302711</specimen-id>
#        <sphinx-id type="integer">121773</sphinx-id>
#        <storage-directory>/external/aibssan/production32/prod341/image_series_74357755/</storage-directory>
#        <weight type="integer">5270</weight>
#        <products type="array">
#          <product>
#            <abbreviation>Mouse</abbreviation>
#            <description>Genome-wide, high-resolution ISH data detailing gene expression throughout the adult mouse brain.</description>
#            <id type="integer">1</id>
#            <name>Mouse Brain</name>
#            <product-name-facet type="integer">2757464828</product-name-facet>
#            <resource>Allen Mouse Brain Atlas</resource>
#            <species>Mouse</species>
#            <species-name-facet type="integer">1861523945</species-name-facet>
#            <tags nil="true"/>
#          </product>
#        </products>
#      </data-set>
#    </data-sets>
#  </object>
###



   open(SEARCHRESF, $in->{fn}) ;
   while (my $line = <SEARCHRESF>) {
      chomp $line;

      if ($line =~ /^ {4}<id>/) {
         ($cur_gene->{id}) = ($line =~ /^    <id>(.+)<\/id>$/) ;

      } elsif ($line =~ /^ {4}<name>/) {
         ($cur_gene->{name}) = ($line =~ /<name>(.+)<\/name>$/) ;

      } elsif ($line =~ /^ {4}<entrez-id>/) {
         ($cur_gene->{"entrez-id"}) =
            ($line =~ /<entrez-id>(.+)<\/entrez-id>$/) ;

      } elsif ($line =~ /^ {4}<acronym>/) {
         ($cur_gene->{"gene-symbol"}) =
            ($line =~ /<acronym>(.+)<\/acronym>$/) ;

      } elsif ($line =~ /^ {8}<id type="integer">/) {
         ($cur_imageseries->{id}) =
            ($line =~ /<id type="integer">(.+)<\/id>$/) ;

      } elsif ($line =~ /^ {8}<plane-of-section-id type="integer">/) {
         my ($plane_id) =
            ($line =~ /<plane-of-section-id type="integer">(.+)<\/plane-of-section-id>$/) ;

         ($cur_imageseries->{"plane-of-section"}) = $planeid2text->{$plane_id};

      } elsif ($line =~ /^ {8}<failed type="boolean">/) {
         ($cur_imageseries->{"failed"}) =
            ($line =~ /<failed type="boolean">(.+)<\/failed>$/) ;


      } elsif ($line =~ /^      <\/data-set>$/) {
         $cur_imageseries->{"age"} = "P56" ;

         if (!exists $hits->{$cur_gene->{"gene-symbol"}}) {
            $hits->{$cur_gene->{"gene-symbol"}}->{imageseries} = {} ;
            foreach my $gene_feat (keys %{$cur_gene}) {
               $hits->{$cur_gene->{"gene-symbol"}}->{$gene_feat} =
                  $cur_gene->{$gene_feat} ; } }

         if (!exists $hits->{$cur_gene->{"gene-symbol"}}->{imageseries}->{$cur_imageseries->{id}}) {
            foreach my $image_series_feat (keys %{$cur_imageseries}) {
$hits->{$cur_gene->{"gene-symbol"}}->{imageseries}->{$cur_imageseries->{id}}->{$image_series_feat} =
                  $cur_imageseries->{$image_series_feat} ;
            }
         }

         $cur_imageseries= {} ;


      } elsif ($line =~ /  <\/gene>$/) {
         $cur_gene = {} ;
      }

   }
   close(SEARCHRESF) ;

   return $hits ;

}


sub download_adultbrain_data {

   my $in = shift ;
   my $specs = $in->{specs} ;
   if (!exists $in->{specs}) {
      $specs = getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   my ($out_fh, $out_fn) ;
   ($out_fh->{adult_wget}, $out_fn->{adult_wget}) =
      tempfile("allen_adult_wget.XXXXX", SUFFIX => '.out') ;

   open(LISTF, $specs->{download}->{fn}->{adult_xpz_list}) ;

   my $header = <LISTF>; $header =~ s/^#// ; chomp $header;
   my @headers = split(/\t/, $header) ;

   print STDERR "Retrieving adult brain expression files:\n";
   my $f2i ; map {$f2i->{$headers[$_]} = $_;} (0 .. $#headers) ;
   my $num_alreadyexist = 0 ;
   while (my $line = <LISTF>) {
      chomp $line;
      my @t = split(/\t/, $line) ;

      my $gene = $t[$f2i->{gene}] ;
      my $imageseries_id = $t[$f2i->{"imageseries-id"}] ;
      my $plane = $t[$f2i->{plane}] ;
      my $age = $t[$f2i->{age}] ;
      my $failed = $t[$f2i->{failed}] ;

      my $xpz_fn = $gene.'_'.$age.'_'.$plane.'_'.$imageseries_id ;
      if ($failed ne 'false') {$xpz_fn .= '_failed';}
      $xpz_fn .= '.xpz' ;

      my $xpz_dir = $specs->{allen_xpz_dir}."/P56/".
                    substr($gene,0,1).'/'.substr($gene,0,3).'/' ;
      if (! -s $xpz_dir) { File::Path::mkpath($xpz_dir); }

      my $url = $specs->{download}->{URL}->{adult_xpz}.$imageseries_id ;
      my $xpz_command = "wget -q \'$url\' -O \'$xpz_dir/$xpz_fn\'" ;

      if (-s "$xpz_dir/$xpz_fn") {
         $num_alreadyexist++ ;
      } else {
         print {$out_fh->{adult_wget}} $xpz_command."\n" ;
      }
   }
   close($out_fh->{adult_wget}) ;
   print STDERR "Skipping $num_alreadyexist expression files that already exist\n";
   close(LISTF) ;

   my ($wget_err_fh, $wget_err_fn) =
      tempfile("adult_wget.XXXXX", SUFFIX => ".err") ;
      close($wget_err_fh) ;
   my ($wget_out_fh, $wget_out_fn) =
      tempfile("adult_wget.XXXXX", SUFFIX => ".out") ;
      close($wget_out_fh) ;

   print STDERR " running wget: " ;
   my $download_com = "bash ".$out_fn->{adult_wget}.
                      " 2> $wget_err_fn >$wget_out_fn" ;
   system($download_com) ;
   print STDERR "X\n";


   ($out_fh->{missing_files}, $out_fn->{missing_files}) =
      tempfile("allen_adult_wget_MISSING_FILES.XXXXX", SUFFIX => '.err') ;

# Check that the expected files were downloaded
   open(WGETF, $out_fn->{adult_wget}) ;
   my $num_missing = 0 ;
   while (my $line = <WGETF>) {
      chomp $line;
      my ($filepath) = ($line =~ /.*-O '(.+)'$/);
      if (-s $filepath) {next;}
      print {$out_fh->{missing_files}} $line."\n" ;
      $num_missing++ ;
   }
   close(WGETF) ;
   close($out_fh->{missing_files}) ;

   if ($num_missing == 0) {
      print STDERR "All files were downloaded\n";
      unlink $out_fn->{missing_files} ;
   } else {
      print STDERR "Warning: $num_missing files did not download\n";
      "Retry by running: bash ".$out_fn->{missing_files}."\n" ;
   }
   
   return ;
}

1 ;
