#!/usr/bin/perl
#
# scoring funtion implemented with Perl
#
# VR

use warnings;
#use Data::Dumper;

use read_mol2;
use get_atom_parameter;
use move;

#Scoring weights:
%Weight  = (
          #'Contact' => -1.7e-1,
          'Repulsion' => -2.7e-1,
          #'Gauss1' => -5.9e-2,
          #'Hydrophobic' => 4.3e-2,
          #'Hydrophobic1' => 1,
          #'Hydrophobic2' => 1,
          'Hydrophobic3' => 5.1e-2,
          #'Hydrogen1' => 3.6,
          #'Hydrogen11' => 1,
          'Hydrogen12' => 5.1,
          #'Hydrogen13' => 1,
          'Hydrogen2' => 1.2,
          #'Hydrogen21' => 1,
          'Hydrogen22' => 1.8,
          #'Hydrogen23' => 1,
          #'Gap' => 1.8e-2,
          'Clash' => 1,
          'Charge' => 1.2e+2, #negative means good
          'Combined' => 1
           );


#take file names
if ((@ARGV == 3) || (@ARGV == 2)) {
$protein = $ARGV[0];
$ligand = $ARGV[1];
@move = (0, 0, 0);
@move = split /:/, $ARGV[2] if ($ARGV[2]);
} else {
die("usage: VR_score.pl protein.mol2 ligand.mol2");
}


#parse files
print STDERR "reading protein...\n";
($head, $atom, $atom_matrix, $foot) = read_mol2::read_mol2($protein);
#@protein_head = @$head;
@protein_atom = @$atom;
@protein_atom_matrix = @$atom_matrix;
#@protein_foot = @$foot;
print STDERR "reading ligand...\n";
($head, $atom, $atom_matrix, $foot) = read_mol2::read_mol2($ligand);
@ligand_head = @$head;
@ligand_atom = @$atom;
@ligand_atom_matrix = @$atom_matrix;
@ligand_foot = @$foot;
#that case probably mol2 files mixed places:
die("Ligand has 200+ atoms. Usage: VR_score.pl protein.mol2 ligand.mol2") if (@ligand_atom > 200);
#print STDERR Dumper \@ligand_head, \@ligand_foot;


#main function
%toppoints = score();
while (($key, $value) = each %toppoints)
{
  print "$key", " ", $value, " ";
}
$main = 0;
#move away check score
while ($main < $move[2]) {
#score ligand pose in protein and write if better
%points = score();
while (($key, $value) = each %points)
{
  print "$key", " ", $value, "\n";
}
@ligand_atom_matrix = move::random_move(@ligand_atom_matrix);
if (($points{'Combined'} < $toppoints{'Combined'}) && ($points{'Repulsion'} < 1)) {
 print STDERR "\n writing ligand...\n";
 write_ligand();
 %toppoints = %points;
}
#print STDERR Dumper \@ligand_atom_matrix[0];
$main++;
}
#write_ligand();


##############################
#
# subroutines here
#
##############################


#distance between atoms
#this one must be called for for each atom pair
#FIXME has to be faster
sub distance_sqared {
my $dxs = $_[0][0]-$_[1][0];
my $dys = $_[0][1]-$_[1][1];
my $dzs = $_[0][2]-$_[1][2];
return $dxs*$dxs+$dys*$dys+$dzs*$dzs;
}


#scoring function
sub score {
#initial score
my $all = 1.1;
%score  = (
          #'Contact' => 0,
          'Repulsion' => 0,
          #'Gap' => 0,
          #'Hydrophobic' => 0,
          #'Hydrophobic1' => 0,
          #'Hydrophobic2' => 0,
          'Hydrophobic3' => 0,
          #'Hydrogen1' => 0,
          #'Hydrogen11' => 0,
          'Hydrogen12' => 0,
          #'Hydrogen13' => 0,
          'Hydrogen2' => 0,
          #'Hydrogen21' => 0,
          'Hydrogen22' => 0,
          #'Hydrogen23' => 0,
          #'Gauss1' => 0,
          'Charge' => 0,
          'Clash' => 0,
          'Combined' => 1.5 #initial value
           );


#electrostatic force, d matrix, Gauss and repulsion
#d=distance-R_ligand_atom-R_protein_atom
$x=0;
while($ligand_atom[$x]{'atom_type'}[0]) {
 $y=0;
 my $lig_radius = get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'radius');
 die("$ligand_atom[$x]{'atom_type'}[0] radius undefined!") if (!($lig_radius));
 while($protein_atom[$y]) {
 $d[$x][$y] = distance_sqared($ligand_atom_matrix[$x],$protein_atom_matrix[$y]);
 if ($d[$x][$y] < 100) { #skip very distant atom pairs
  $score{'Charge'} += $ligand_atom[$x]{'charge'}*$protein_atom[$y]{'charge'}/$d[$x][$y];
  $d[$x][$y] = sqrt($d[$x][$y]) 
  - $lig_radius
  - get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'radius');
  #calculate Gauss1 and Gauss2
  if ($d[$x][$y] < 2) {
   #$score{'Gauss1'} += exp(-8*($d[$x][$y]**2)); #means if abs distance 0 then +1 else +less
   #$score{'Gap'}++ if ($d[$x][$y] > 0); # count just gap
   #$score{'Contact'}++ if (abs($d[$x][$y]) < 0.25);
    if ($d[$x][$y] < 0) {
     #calculate repulsion:
     $score{'Repulsion'}++;
     $score{'Clash'}++ if ($d[$x][$y] < -2);
     #print STDERR "\n $ligand_atom[$x]{'atom_id'} repeals $protein_atom[$y]{'atom_id'}\n";
    }
  }
 }
 $y++;
 }
$x++;
}
#print STDERR "d size ", $x, " ", $y, "\n";


#calculate hydrophobic
$x = 0;
my $hydrophobic = 0;
while($d[$x]) {
 $y = 0;
 if (get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'hydrophobic')) {
  while($d[$x][$y]) {
   if ((0.25 < $d[$x][$y]) && ($d[$x][$y] < 2)) {
   if (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'hydrophobic')) {
    #if ($d[$x][$y] < 0.5) {
      #$score{'Hydrophobic'}++;
      #$score{'Hydrophobic1'}++ if (abs($d[$x][$y]) < 0.25);
      #$score{'Hydrophobic2'}++ if ($d[$x][$y] < -0.25);
      $score{'Hydrophobic3'}++;
   }
   }
   $y++;
  }
 }
$x++;
}


#ligand as hydrogen bond donor
$x = 0;
my $hydrogenbd = 0;
while($d[$x]) {
 $y = 0;
 if (($ligand_atom[$x]{'charge'} > 0.1) && ($ligand_atom[$x]{'atom_type'}[0] eq 'H')) {
  while($d[$x][$y]) {
   if ((-2 < $d[$x][$y]) && ($d[$x][$y] < -0.25)) {
   if (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'H_acceptor')) {
      #$score{'Hydrogen1'}++ if ($d[$x][$y] < 0);
      #$score{'Hydrogen11'}++ if (abs($d[$x][$y]) < 0.25);
      $score{'Hydrogen12'}++;
      #$score{'Hydrogen13'}++ if ($d[$x][$y] > 0.25);
   }
   }
   $y++;
  }
 }
$x++;
}


#ligand as hydrogen bond acceptor
$x = 0;
my $hydrogenba = 0;
while($d[$x]) {
 $y = 0;
 if (get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'H_acceptor')) {
  while($d[$x][$y]) {
   if ((-2 < $d[$x][$y]) && ($d[$x][$y] < 0)) {
   if (($protein_atom[$y]{'charge'} > 0.1) && ($protein_atom[$y]{'atom_type'}[0] eq 'H')) {
     $score{'Hydrogen2'}++;
     #$score{'Hydrogen21'}++ if (abs($d[$x][$y]) < 0.25);
     $score{'Hydrogen22'}++ if ($d[$x][$y] < -0.25);
     #$score{'Hydrogen23'}++ if ($d[$x][$y] > 0.25);
   }
   }
   $y++;
  }
 }
$x++;
}


#* by Weight and calculate combined score
foreach my $key ( keys %score )
{
   if ($key ne 'Combined') {
    $score{$key} *= $Weight{$key};
    $score{'Combined'} += $score{$key};
   }
}
return %score;
}


#write edited ligand as *.mol2 file
sub write_ligand {
#gather lines
my @lines;
foreach my $line (@ligand_head) {
 push @lines, $line;
}
$x = 0;
while ($ligand_atom[$x]{'atom_type'}[0]) {
$ligand_atom[$x]{'atom_type'}[1]=" " if (!$ligand_atom[$x]{'atom_type'}[1]);
 $line =    ' '.
            $ligand_atom[$x]{'atom_id'}.' '.
            $ligand_atom[$x]{'atom_name'}.' '.
            $ligand_atom_matrix[$x][0].' '.
            $ligand_atom_matrix[$x][1].' '.
            $ligand_atom_matrix[$x][2].' '.
#            $ligand_atom[$x]{'status_bit'}.' '.
            $ligand_atom[$x]{'atom_type'}[0].'.'.
            $ligand_atom[$x]{'atom_type'}[1].' '.
            $ligand_atom[$x]{'subst_id'}.' '.
            $ligand_atom[$x]{'subst_name'}.' '.
            $ligand_atom[$x]{'charge'};
 #print STDERR "$line\n";
 push @lines, $line;
 $x++;
}
foreach my $line (@ligand_foot) {
 push @lines, $line;
}

#ligand file name
$name = ">out_$toppoints{'Combined'}_$ligand";


#print them
open (LIGAND, $name) or die('can not write ligand output file');
foreach my $line (@lines) {
 print LIGAND $line;
}
close(LIGAND);
}
