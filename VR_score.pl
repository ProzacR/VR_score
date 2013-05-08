#!/usr/bin/perl
#
# scoring funtion implemented with Perl
#
# VR

use warnings;
use Data::Dumper;

use read_mol2;

#take file names
if (@ARGV == 2) {
$protein = $ARGV[0];
$ligand = $ARGV[1];
} else {
die("usage: VR_score.pl protein.mol2 ligand.mol2");
}

#parse files
@protein_atom = read_mol2::read_mol2($protein);
@ligand_atom = read_mol2::read_mol2($ligand);
print Dumper \$protein_atom[0];

