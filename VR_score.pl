#!/usr/bin/perl
#
# scoring funtion implemented with Perl
#
# VR

use warnings;
use Data::Dumper;

use read_mol2;
use get_atom_parameter;


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

#simple parameter get:
$par = get_atom_parameter::get_atom_parameter('C', 'radius');
print $par, "\n";
