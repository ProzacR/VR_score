#!/usr/bin/perl
#
# try to get Kd from pdb.org using id
#
# VR

use warnings;

die('usage: ./get_Kd.pl id') if (@ARGV != 1);
my $id = $ARGV[0];


my $url = 'http://pdb.org/pdb/explore/explore.do?structureId=' . $id;

use LWP::Simple;
my $content = get $url;
die "Couldn't get $url" unless defined $content;
@parts = split(/"External Ligand Annotations"/, $content);
@lines = split(/\n/, $parts[0]);
foreach(@lines) {
  if ($_ =~ m/LigandSummary/i) {
   @data1 = split(/\t/, $_);
  }
  if($_ =~ m/Kd\:/i) {
    @data = split(/>|</, $_);
    last;
  }
}


#foreach (@data1) {
#print "$_\n";
#}


die "no data!" unless defined @data;
print $id,",",$data1[4],",",$data[4],",",$data[10],"\n";
