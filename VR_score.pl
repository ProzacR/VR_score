#!/usr/bin/perl
#
# scoring funtion implemented with Perl
#
# VR

use warnings;
use Data::Dumper;

use read_mol2;
use get_atom_parameter;

use constant e_math => 2.71828;


#take file names
if (@ARGV == 2) {
$protein = $ARGV[0];
$ligand = $ARGV[1];
} else {
die("usage: VR_score.pl protein.mol2 ligand.mol2");
}


#parse files
print STDERR "reading protein...\n";
@protein_atom = read_mol2::read_mol2($protein);
print STDERR "readling ligand...\n";
@ligand_atom = read_mol2::read_mol2($ligand);
#that case probably mol2 files mixed places:
die("Ligand has 200+ atoms. Usage: VR_score.pl protein.mol2 ligand.mol2") if (@ligand_atom > 200);


$main = 0;
#move away check score
while ($main < 1) {
#score ligand pose in protein:
%points = score();
while (($key, $value) = each %points)
{
  print "$key", " ", $value, "\n";
}
move_ligand('x', 'p');
$main++;
}


##############################
#
# subroutines here
#
##############################


#move ligand
#use ex. move_ligand('x', 'p')
sub move_ligand {
my $p;
my $x = 0;
while ($ligand_atom[$x]{$_[0]}) {
 # + or - direction?
 $p = 1 if ($_[1] eq 'p');
 $p = -1 if ($_[1] eq 'n');
 $ligand_atom[$x]{$_[0]} += $p;
 $x++;
}
return 1;
}


#distance between atoms
sub distance_sqared {
my $dxs = ($_[0]{'x'}-$_[1]{'x'})**2;
my $dys = ($_[0]{'y'}-$_[1]{'y'})**2;
my $dzs = ($_[0]{'z'}-$_[1]{'z'})**2;
return $dxs+$dys+$dzs;
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


#electrostatic force
my $F;
$x=0;
while($ligand_atom[$x]{'atom_type'}[0]) {
 $y=0;
 while($protein_atom[$y]) {
 $d[$x][$y] = distance_sqared($ligand_atom[$x],$protein_atom[$y]);
 $F += $ligand_atom[$x]{'charge'}*$protein_atom[$y]{'charge'}/$d[$x][$y];
 $y++;
 }
$x++;
}
#print STDERR "Electrostatic force: ", $F, "\n";


#calculate d matrix (d values)
#d=distance-R_ligand_atom-R_protein_atom
print STDERR "calculating d";
$x=0;
while($ligand_atom[$x]{'atom_type'}[0]) {
 print STDERR $ligand_atom[$x]{'atom_type'}[0];
 $y=0;
 while($protein_atom[$y]) {
 $d[$x][$y] = sqrt($d[$x][$y]) 
 - get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'radius')
 - get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'radius') if (abs($d[$x][$y]) < 100); 
 $y++;
 }
$x++;
print STDERR ".";
}
print STDERR "OK\n";
#print STDERR Dumper \@d


#calculate repulsion
$x = 0;
my $repulsion;
while($d[$x]) {
 $y = 0;
 while($d[$x][$y]) {
  if ($d[$x][$y]<0) {
  $repulsion += $d[$x][$y]**2;
  }
 $y++;
 }
$x++;
}
#print STDERR "Repulsion: ", $repulsion, "\n";


#calculate Gauss1
$x = 0;
my $Gauss1;
while($d[$x]) {
 $y = 0;
 while($d[$x][$y]) {
 $Gauss1 += e_math ** (-(($d[$x][$y]*2)**2));
 $y++;
 }
$x++;
}
#print STDERR "Gauss1: ", $Gauss1, "\n";


#calculate Gauss2
$x = 0;
my $Gauss2;
while($d[$x]) {
 $y = 0;
 while($d[$x][$y]) {
 $Gauss2 += e_math ** (-((($d[$x][$y]-3)/2)**2));
 $y++;
 }
$x++;
}
#print STDERR "Gauss2: ", $Gauss2, "\n";


#calculate hydrophobic
$x = 0;
my $hydrophobic;
while($d[$x]) {
 $y = 0;
 if (get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'hydrophobic')) {
  while($d[$x][$y]) {
   if (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'hydrophobic')) {
    if ($d[$x][$y] < 0.5) {
      $hydrophobic++;
     } elsif ($d[$x][$y] < 1.5) {
      $hydrophobic += -$d[$x][$y] + 1.5; #so linearly interpolated
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
my $hydrogenbd;
while($d[$x]) {
 $y = 0;
 if (($ligand_atom[$x]{'charge'} > 0.1) &&($ligand_atom[$x]{'atom_type'}[0] eq 'H')) {
  while($d[$x][$y]) {
   if (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'H_acceptor')) {
    if ($d[$x][$y] < -0.7) {
      $hydrogenbd++;
     } elsif ($d[$x][$y] < 0) {
      $hydrogenbd += -1.45 * $d[$x][$y]; #so linearly interpolated
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
my $hydrogenba;
while($d[$x]) {
 $y = 0;
 if (get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'H_acceptor')) {
  while($d[$x][$y]) {
   if (($protein_atom[$y]{'charge'} > 0.1) &&($protein_atom[$y]{'atom_type'}[0] eq 'H')) {
    if ($d[$x][$y] < -0.7) {
      $hydrogenba++;
     } elsif ($d[$x][$y] < 0) {
      $hydrogenba += -1.45 * $d[$x][$y]; #so linearly interpolated
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
