#!/usr/bin/perl
#
# scoring funtion implemented with Perl
#
# VR

use warnings;
use Data::Dumper;

use read_mol2;
use get_atom_parameter;
use move;


#take file names
if ((@ARGV == 3) || (@ARGV == 2)) {
$protein = $ARGV[0];
$ligand = $ARGV[1];
@move = (0, 1, 1);
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
print STDERR "readling ligand...\n";
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
  print "$key", " ", $value, "\n";
}
$main = 0;
#move away check score
while ($main < $move[2]) {
#score ligand pose in protein and write if better
%points = score();
#while (($key, $value) = each %points)
#{
#  print "$key", " ", $value, "\n";
#}
@ligand_atom_matrix = move::random_move(@ligand_atom_matrix);
if ($points{'Combined'} < $toppoints{'Combined'}) {
 print STDERR "writing ligand...\n";
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
sub distance_sqared {
my $dxs = $_[0][0]-$_[1][0];
my $dys = $_[0][1]-$_[1][1];
my $dzs = $_[0][2]-$_[1][2];
return $dxs*$dxs+$dys*$dys+$dzs*$dzs;
}


#scoring function
sub score {
#lower means better:
%Weight  = (
          'Repulsion' => 1,
          'Gauss1' => -1e-2,
          'Hydrophobic' => -1e-1,
          'hydrogen1' => -1,
          'hydrogen2' => -1,
          'Gauss2' => -1e-3,
          'Electrostatic' => 1 #negative means good
           );


#electrostatic force, d matrix, Gauss and repulsion
#d=distance-R_ligand_atom-R_protein_atom
my $F = 0;
my $Gauss1 = 0;
my $Gauss2 = 0;
my $repulsion = 0;
$x=0;
while($ligand_atom[$x]{'atom_type'}[0]) {
 $y=0;
 my $lig_radius = get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'radius');
 while($protein_atom[$y]) {
 $d[$x][$y] = distance_sqared($ligand_atom_matrix[$x],$protein_atom_matrix[$y]);
 if ($d[$x][$y] < 100) { #skip very distant atom pairs
  $F += $ligand_atom[$x]{'charge'}*$protein_atom[$y]{'charge'}/$d[$x][$y];
  $d[$x][$y] = sqrt($d[$x][$y]) 
  - $lig_radius
  - get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'radius');
  #calculate Gauss1 and Gauss2
  if ($d[$x][$y] < 10) {
   $Gauss1 += exp(-(($d[$x][$y]*2)**2));
   $Gauss2 += exp(-((($d[$x][$y]-3)/2)**2));
    if ($d[$x][$y]<0) {
     #calculate repulsion
     $repulsion += $d[$x][$y]**2;
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
   if ($d[$x][$y] < 1.5) {
   if (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'hydrophobic')) {
    if ($d[$x][$y] < 0.5) {
      $hydrophobic++;
     } else {
      $hydrophobic += -$d[$x][$y] + 1.5; #so linearly interpolated
     }
   }
   }
   $y++;
  }
 }
$x++;
}
#print STDERR "Hydrophobic: ", $hydrophobic, "\n";


#ligand as hydrogen bond donor
$x = 0;
my $hydrogenbd = 0;
while($d[$x]) {
 $y = 0;
 if (($ligand_atom[$x]{'charge'} > 0.1) && ($ligand_atom[$x]{'atom_type'}[0] eq 'H')) {
  while($d[$x][$y]) {
   if ($d[$x][$y] < 0) {
   if (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'H_acceptor')) {
    if ($d[$x][$y] < -0.7) {
      $hydrogenbd++;
     } else {
      $hydrogenbd += -1.45 * $d[$x][$y]; #so linearly interpolated
     }
   }
   }
   $y++;
  }
 }
$x++;
}
#print STDERR "Ligand as hydrogen bond donor: ", $hydrogenbd, "\n";


#ligand as hydrogen bond acceptor
$x = 0;
my $hydrogenba = 0;
while($d[$x]) {
 $y = 0;
 if (get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'H_acceptor')) {
  while($d[$x][$y]) {
   if ($d[$x][$y] < 0) {
   if (($protein_atom[$y]{'charge'} > 0.1) && ($protein_atom[$y]{'atom_type'}[0] eq 'H')) {
    if ($d[$x][$y] < -0.7) {
      $hydrogenba++;
     } else {
      $hydrogenba += -1.45 * $d[$x][$y]; #so linearly interpolated
     }
   }
   }
   $y++;
  }
 }
$x++;
}
#print STDERR "Ligand as hydrogen bond acceptor: ", $hydrogenba, "\n";


#return sore
my %score = ("Electrostatic" => $F, "Repulsion" => $repulsion, "Gauss1" => $Gauss1, "Gauss2" => $Gauss2,
"Hydrophobic" => $hydrophobic, "hydrogen1" => $hydrogenbd, "hydrogen2" => $hydrogenba);
#* by Weight and combined score
my $all;
foreach my $key ( keys %score )
{
   $score{$key} *= $Weight{$key};
   $all += $score{$key};
}
#print STDERR "Combined: ", $all, "\n";
$score{'Combined'} = $all;
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
