=head1 NAME

alnmnr::spinalcord.pm - retrieve and search ABA mouse spinal cord data

=head1 AUTHOR

Fred P. Davis (fredpdavis@gmail.com)

=head1 SYNOPSIS

SEARCHSPINE offers customized searches of the rostrocaudal and dorso-ventral
expression patterns of juvenile (P4) and adult (P56) spinal cord available from
the Allen Brain Atlas (http://mousespinal.brain-map.org)

SEARCHSPINE searches the ABA spinal cord data (in the text file from
download_spinalcord_data() ) with user-specified criteria and produces HTML
output.

Each search takes under a second.

=head1 DESCRIPTION

=cut


package alnmnr::spinalcord ;
use strict;
use warnings;
use File::Temp qw/tempdir tempfile/ ;

use Exporter ;
our @ISA = qw/Exporter/ ;
our @EXPORT = qw/_run_spinesearch/ ;

sub _run_spinesearch {
# purpose: searches text format aba spinal cord data for query criteria

   my $in = shift ;
   my $specs = $in->{specs} ;
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }

   my $usage = "allenminer.pl [-rcbits RCPATTERN] [-dvbits DVPATTERN] [-min_rc INTEGER] [-maxrc INTEGER] [-min_dv INTEGER] [-max_dv INTEGER] [-ignore_missing 1] [-color 0|1]

 -color 0|1: 0 = don't use color, default 1 = use colors in output

 -rcbits RCPATTERN: 16-character long string of 0,1,or .(wildcard)
 -dvbits DVPATTERN: 11-character long string of 0,1,or .(wildcard)

 -min_rc (1-16): Minimum number of RC bits with expression
 -min_dv (1-11): Minimum number of DV bits with expression

 -max_rc (1-16): Maximum number of RC bits with expression
 -max_dv (1-11): Maximum number of DV bits with expression

 -ignore_missing 0: will treat missing bits as mis-matches (default)
                 1: will allow missing bits for '0' and '1' query bits
                 2: will allow missing bits for '0' query bits
                 3: will allow missing bits for '1' query bits

  Bit reminder: 0=no expression, 1=expression, .=wildcard, x=missing data
  * RC: ~4 bits each in cervical, thoracic, lumbar, sacral/coccygeal
  * DV: 1. Laminae 1-3
        2. Laminae 4-6
        3. Laminae 7-8
        4. Laminae 9
        5. Intermediolateral Column
        6. Gray Matter
        7. White Matter
        8. Central Canal
        9. Ventral-dorsal Midline in Gray Matter
        10. Radially Arrayed in White Matter
        11. Vascular-like in Gray and White Matter
";

   if (!exists $in->{ignore_missing}) { $in->{ignore_missing} = '0';}
   if (!exists $in->{color})          { $in->{color} = '1';}


   if (!-s $specs->{allen_spinalcord_fn}) {
      print STDERR "ERROR: spinal cord data file not found: ".
                   $specs->{allen_spinalcord_fn}."\n";
      die $usage ;
   }

   if (exists $in->{h} || exists $in->{help}) {
      die $usage ; }


   my $out_fh ;
   if (exists $in->{out_fn}) {
      open($out_fh, ">".$in->{out_fn}) ;
   } else {
      ($out_fh, undef) = tempfile("spinesearch_results.XXXXX",SUFFIX=>".html");
   }


# search options:
#  rcbits regex: .=wildcard, 1=on, 0=off, x=missing
#  dvbits regex: .=wildcard, 1=on, 0=off
#  later: calc-RC gradient (bins)

   my $gene_url_strx = {
      prefix => 'http://mousespinal.brain-map.org/imageseries/detail/',
      suffix => '.html'
   };
   my $image_url_strx = {
      prefix => 'http://mousespinal.brain-map.org/imageseries/show.html?id=',
      suffix => ''
   };


# Double check search criteria correctness
   my $exp_bit_length = {
      rcbits => 16,
      dvbits => 11,
   } ;
   foreach my $type (qw/rcbits dvbits/) {
      if (!exists $in->{$type}) {next;}
      if (length($in->{$type}) != $exp_bit_length->{$type} ||
          $in->{$type} !~ /^[01\.\*]+$/) {
            die 'ERROR: -'.$type.' must be a '.$exp_bit_length->{$type}.
                '-character string containing only 0, 1, . (or *)' ;
      }
      $in->{$type} =~ s/\*/\./g ;

      if ($in->{ignore_missing} == 1) {
         $in->{$type} =~ s/1/\[1x\]/g ;
         $in->{$type} =~ s/0/\[0x\]/g ;
      } elsif ($in->{ignore_missing} == 2) {
         $in->{$type} =~ s/0/\[0x\]/g ;
      } elsif ($in->{ignore_missing} == 3) {
         $in->{$type} =~ s/1/\[1x\]/g ;
      }
   }

   print {$out_fh} "<html>\n" ;
   print {$out_fh} "<body>\n" ;
   print {$out_fh} "<h2>ALLENMINER results (".localtime().")</h2>\n" ;
   my @options ;
   print {$out_fh} "<b>Search options:</b><br>\n";
   print {$out_fh} '<table style="border-width: 1px; border-color:#000000; border-style: solid;">'."\n";
   foreach my $key (sort keys %{$in}) {
      if ($key eq 'ARGV' || $key eq 'specs') {next;}
      my $val = $in->{$key} ;
      if ($key =~ /bits/) {
         $val = '<tt>'.$val.'</tt>' ;
      }
      print {$out_fh} '<tr><th align="left">'.$key.'</th><td>'.$val.'</td></tr>'."\n";
   }
   print {$out_fh} "</table>\n";
   print {$out_fh} "<p><b>Command line:</b><br>\n";
   print {$out_fh} "<tt> allenminer.pl ".join(" ", @ARGV)."</tt><br>\n";
   print {$out_fh} "<p>\n";
   print {$out_fh} "<b>Search results:</b>\n";
   print {$out_fh} "<table align=center cellspacing=5>\n";

   {
      my @outheaders = qw/No. Gene Age R-C_pattern D-V_pattern Imageseries/ ;
      foreach my $j ( 0 .. $#outheaders) {
         $outheaders[$j] = '<th>'.$outheaders[$j].'</th>' ; }
      print {$out_fh} '<tr>'.join(" ", @outheaders).'</tr>'."\n" ;
   }

   if ($specs->{allen_spinalcord_fn} =~ /gz$/) {
      open(DATAF, "zcat ".$specs->{allen_spinalcord_fn}." |") ;
   } else {
      open(DATAF, $specs->{allen_spinalcord_fn}) ;
   }

   my $f2i = {};
   {
      my $headers = <DATAF> ; chomp $headers; $headers =~ s/^\#// ;
      my @t = split(/\t/, $headers) ;
      map {$f2i->{$t[$_]} = $_; } (0 .. $#t) ;
   }
   my $hitnum = 0 ;
   while (my $line = <DATAF>) {
      chomp $line;
      if ($line =~ /^#/) {next;}
      my @t = split(/\t/, $line) ;

      if (exists $in->{rcbits} &&
          $t[$f2i->{rcbits}] !~ /$in->{rcbits}/) {next;}

      if (exists $in->{dvbits} &&
          $t[$f2i->{dvbits}] !~ /$in->{dvbits}/) {next;}

      if (exists $in->{min_rc}) {
         my $num_rc_on = $t[$f2i->{rcbits}] =~ tr/1//;
         if ($num_rc_on < $in->{min_rc}) {next;}
      }

      if (exists $in->{max_rc}) {
         my $num_rc_on = $t[$f2i->{rcbits}] =~ tr/1//;
         if ($num_rc_on > $in->{max_rc}) {next;}
      }

      if (exists $in->{min_dv}) {
         my $num_dv_on = $t[$f2i->{dvbits}] ; $num_dv_on = tr/1//;
         if ($num_dv_on < $in->{min_dv}) {next;}
      }

      if (exists $in->{max_dv}) {
         my $num_dv_on = $t[$f2i->{dvbits}] ; $num_dv_on = tr/1//;
         if ($num_dv_on > $in->{max_dv}) {next;}
      }

      $hitnum++ ;
      my $gene_url = $gene_url_strx->{prefix}.$t[$f2i->{imageseries}].
                     $gene_url_strx->{suffix} ;

      my $image_url = $image_url_strx->{prefix}.$t[$f2i->{imageseries}].
                     $image_url_strx->{suffix} ;
      my $rcbits = $t[$f2i->{rcbits}];
      my $dvbits = $t[$f2i->{dvbits}];
      if ($in->{color}) {
         $rcbits =~ s/x/<font color="grey">x<\/font>/g ;
         $rcbits =~ s/1/<font color="green">1<\/font>/g ;
         $rcbits =~ s/0/-/g ;
         $dvbits =~ s/x/<font color="grey">x<\/font>/g ;
         $dvbits =~ s/1/<font color="blue">1<\/font>/g ;
         $dvbits =~ s/0/-/g ;
      }

      my @outfields ;
      push @outfields, $hitnum ;
      push @outfields, '<a href="'.$gene_url.'">'.$t[$f2i->{gene_id}].'</a>' ;
      push @outfields, $t[$f2i->{age}] ;
      push @outfields, '<tt>'.$rcbits.'</tt>' ;
      push @outfields, '<tt>'.$dvbits.'</tt>' ;
      push @outfields, '<a href="'.$image_url.'">'.
                       $t[$f2i->{imageseries}].'</a>' ;
      foreach my $j ( 0 .. $#outfields) {
         $outfields[$j] = '<td>'.$outfields[$j].'</td>' ; }
      if ($in->{color} && ($hitnum % 2 == 1)) {
         print {$out_fh} '<tr bgcolor="lightgray">'.join(" ", @outfields).'</tr>'."\n" ;
      } else {
         print {$out_fh} '<tr>'.join(" ", @outfields).'</tr>'."\n" ;
      }
   }
   print {$out_fh} "</table><hr>\n" ;
   print {$out_fh} "</body>" ;
   print {$out_fh} "</html>" ;

}


sub download_spinalcord_data {
# purpose: retrieves ABA spinal cord web pages and extracts expression patterns

   my $in = shift ;
   my $type = $in->{atlas_type} ;

   my $specs = $in->{specs} ;
   if (!exists $in->{specs}) {
      $specs = alnmnr::getspecs() ;
   } else {
      $specs = $in->{specs} ;
   }


   my ($temp_fh, $temp_fn) = tempfile() ; close($temp_fh) ;
   my $url = $specs->{download}->{URL}->{spinalcord_part1}.'1'.
             $specs->{download}->{URL}->{spinalcord_part2};
   system("wget --quiet \'$url\' -O $temp_fn") ;
   if (!-s $temp_fn) { die "ERROR retrieving first ABA webpage: $url\n"; }

   my $numpages ;
   open(PAGEONE, $temp_fn) ;
   while (my $line = <PAGEONE>) {
      chomp $line;
      if ($line !~ /"pagination"/) {next;}
      while ($line =~ /([0-9]+)\<\/a/g) {$numpages = $1;}
      print STDERR "$numpages total pages\n";
      last;
   }
   close(PAGEONE) ;

   my $out_fh; open($out_fh, ">".$specs->{allen_spinalcord_fn}) ;

   print {$out_fh} '#'.join("\t", qw/imageseries gene_id age rcbits dvbits/)."\n";

   extract_expression_patterns({in_fn => $temp_fn, out_fh => $out_fh}) ;
   unlink $temp_fn ;

   print STDERR "now on page: 1";
   foreach my $pagenum (2 .. $numpages) {
      print STDERR "\b"x(length($pagenum - 1)).$pagenum ;
      $url = $specs->{download}->{URL}->{spinalcord_part1}.$pagenum.
             $specs->{download}->{URL}->{spinalcord_part2};
      system("wget --quiet \'$url\' -O $temp_fn") ;
      if (!-s $temp_fn) { die "ERROR retrieving page $url\n";}
      extract_expression_patterns({in_fn => $temp_fn, out_fh => $out_fh}) ;
      unlink $temp_fn ;
   }
   close($out_fh) ;
   print STDERR "\nDone!\n" ;

}


sub extract_expression_patterns {

   my $in = shift ;
   my $in_fn = $in->{in_fn} ;
   my $out_fh = $in->{out_fh} ;

   open(DATA, $in_fn) ;
   my $in_table = 0 ;
   my $cur_record = {};
   while(my $line = <DATA>) {
      chomp $line;
      if ($in_table) {
         if ($line =~ /"checkbox"/) {

            if (exists $cur_record->{imageseries}) {
               print {$out_fh} join("\t",
                @$cur_record{qw/imageseries gene_id age rcbits dvbits/})."\n"; }

            my ($imageseries) = ($line =~ /id="([0-9]+)"/) ;

            my @keys = keys %{$cur_record} ;
            map {delete $cur_record->{$_};} @keys ;

            $cur_record->{imageseries} = $imageseries;
            $cur_record->{rcbits} = '' ;
            $cur_record->{dvbits} = '' ;

         } elsif ($line =~ /^\S4/) {
            $cur_record->{age} = 4 ;
          
         } elsif ($line =~ /"tableLink"/) {

            my ($gene_id) = ($line =~ /tableLink"\>(.+)\<\/a/) ;
            $cur_record->{gene_id} = $gene_id ;
            my $line2 = <DATA> ; #</td>
            $line2 = <DATA> ;    #<td >
            $line2 = <DATA> ;    #gene name
            $line2 = <DATA> ;    #</td>
            $line2 = <DATA> ;    #<td >
            $line2 = <DATA> ;    #<td > 
            my ($age) = ($line2 =~ /([0-9]+)/) ;
            $cur_record->{age} = $age ;

         } elsif ($line =~ /ection/) {

            if ($line =~ /missingSection/) {
               $cur_record->{rcbits} .= 'x' ;
            } else {
               my $exprtype = '';
               if ($line =~ /sectionB/)    { $exprtype = 'rcbits'; }
               elsif ($line =~ /sectionC/) { $exprtype = 'dvbits'; }

               if ($line =~ /100\%/) {$cur_record->{$exprtype} .= '1';}
               else                  {$cur_record->{$exprtype} .= '0';}
            }

         } elsif ($line =~ /Compare selected/) {
            last;
         }
      } elsif ($line =~ /Images/) {
         $in_table = 1;
      }
   }
   print {$out_fh} join("\t",
      @$cur_record{qw/imageseries gene_id age rcbits dvbits/})."\n";

}
