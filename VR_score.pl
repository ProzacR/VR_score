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
print STDERR "readling ligand....\n";
@ligand_atom = read_mol2::read_mol2($ligand);
#that case probably mol2 files mixed places:
die("Ligand has 200+ atoms. Usage: VR_score.pl protein.mol2 ligand.mol2") if (@ligand_atom > 200);


#electrostatic force
$x=0;
while($ligand_atom[$x]) {
 $y=0;
 while($protein_atom[$y]) {
 $F += $ligand_atom[$x]{'charge'}*$protein_atom[$y]{'charge'}/distance_sqared($ligand_atom[$x],$protein_atom[$y]);
 $y++;
 }
$x++;
}
print "Electrostatic force: ", $F, "\n";


#distance between atoms
sub distance_sqared {
my $dxs = ($_[0]{'x'}-$_[1]{'x'})**2;
my $dys = ($_[0]{'y'}-$_[1]{'y'})**2;
my $dzs = ($_[0]{'z'}-$_[1]{'z'})**2;
return $dxs+$dys+$dzs;
}


#calculate d matrix (d values)
#d=distance-R_ligand_atom-R_protein_atom
print STDERR "calculating d";
$x=0;
while($ligand_atom[$x]) {
 #print STDERR $ligand_atom[$x]{'atom_type'}[0];
 $y=0;
 while($protein_atom[$y]) {
 $d[$x][$y] = sqrt(distance_sqared($ligand_atom[$x],$protein_atom[$y])) 
 - get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'radius')
 - get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'radius');
 $y++;
 }
$x++;
print STDERR ".";
}
print STDERR "OK\n";
#print STDERR Dumper \@d


#calculate repulsion
$x=0;
while($d[$x]) {
 $y=0;
 while($d[$x][$y]) {
  if ($d[$x][$y]<0) {
  $repulsion += $d[$x][$y]**2;
  }
 $y++;
 }
$x++;
}
print "Repulsion: ", $repulsion, "\n";


#calculate Gauss1
$x=0;
while($d[$x]) {
 $y=0;
 while($d[$x][$y]) {
 $Gauss1 += e_math ** (-(($d[$x][$y]*2)**2));
 $y++;
 }
$x++;
}
print "Gauss1: ", $Gauss1, "\n";


#calculate Gauss2
$x=0;
while($d[$x]) {
 $y=0;
 while($d[$x][$y]) {
 $Gauss2 += e_math ** (-((($d[$x][$y]-3)/2)**2));
 $y++;
 }
$x++;
}
print "Gauss2: ", $Gauss2, "\n";


#calculate hydrophobic
$x=0;
while($d[$x]) {
 $y=0;
 if (get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'hydrophobic')) {
  while($d[$x][$y]) {
   if (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'hydrophobic')) {
    if ($d[$x][$y]<0.5) {
      $hydrophobic++;
     } elsif ($d[$x][$y]<1.5) {
      $hydrophobic += -$d[$x][$y]+1.5; #so linearly interpolated
     }
   }
   $y++;
  }
 }
$x++;
}
print "Hydrophobic: ", $hydrophobic, "\n";


#ligand as hydrogen bond donor
$x=0;
while($d[$x]) {
 $y=0;
 if (($ligand_atom[$x]{'charge'}>0.1) &&($ligand_atom[$x]{'atom_type'}[0] eq 'H')) {
  while($d[$x][$y]) {
   if (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'H_acceptor')) {
    if ($d[$x][$y]<-0.7) {
      $hydrogenbd++;
     } elsif ($d[$x][$y]<0) {
      $hydrogenbd += -1.45*$d[$x][$y]; #so linearly interpolated
     }
   }
   $y++;
  }
 }
$x++;
}
print "Ligand as hydrogen bond donor: ", $hydrogenbd, "\n";


#ligand as hydrogen bond acceptor
$x=0;
while($d[$x]) {
 $y=0;
 if (get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'H_acceptor')) {
  while($d[$x][$y]) {
   if (($protein_atom[$y]{'charge'}>0.1) &&($protein_atom[$y]{'atom_type'}[0] eq 'H')) {
    if ($d[$x][$y]<-0.7) {
      $hydrogenba++;
     } elsif ($d[$x][$y]<0) {
      $hydrogenba += -1.45*$d[$x][$y]; #so linearly interpolated
     }
   }
   $y++;
  }
 }
$x++;
}
print "Ligand as hydrogen bond acceptor: ", $hydrogenba, "\n";



