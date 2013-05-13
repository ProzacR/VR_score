package read_mol2;

#
# simple *.mol2 file reader
# http://www.tripos.com/data/support/mol2.pdf
#
# VR

use warnings;
use Data::Dumper;


##take file name
#if (@ARGV == 1) {
#$file = $ARGV[0];
#} else {
#die("usage: read_mol2.pl file.mol2");
#}

sub read_mol2 {
#read file into lines
@lines = '';
open(INFO, $_[0]) or die("Could not open file.");
while(<INFO>) {
 push (@lines, $_);
}
close(INFO);
#print Dumper \@lines;


$line_number=0;
#skip to first @<TRIPOS>MOLECULE
while(!($lines[$line_number] =~ /MOLECULE/)) {
#just ignore comments and empty lines
#print STDERR "comment: ", @lines[$line_number];
$line_number++;
}
$line_number++;
#read MOLECULE part (dump for now)
while(!($lines[$line_number] =~ /ATOM/)) {
#print STDERR "molecule: ", $lines[$line_number];
$line_number++;
}
$line_number++;
$x=0;
#read ATOM part to data structure
my @atoms;
my @atom;
while(!($lines[$line_number] =~ /BOND/)) {
#print STDERR "atom: ", @lines[$line_number];
#split lines into symbols
push @{$atoms[$x]}, split(/ +/, $lines[$line_number]);
$x++;
$line_number++;
}
$line_number++;
#print Dumper \@atoms[0];
#ignore BOND part for now....


#make formated hash:
$x = 0;
while($atoms[$x]) {
#COLUMNS DATA TYPE CONTENTS
#--------------------------------------------------------------------------------
# atom_id (integer) = the ID number of the atom at the time the file was
# created. This is provided for reference only and is not used when the
# .mol2
# file is read into SYBYL.
$atom[$x]{'atom_id'} = @{$atoms[$x]}[1];

# atom_name (string) = the name of the atom.
$atom[$x]{'atom_name'} = @{$atoms[$x]}[2];

#
$atom[$x]{'x'} = @{$atoms[$x]}[3];

#
$atom[$x]{'y'} = @{$atoms[$x]}[4];

#
$atom[$x]{'z'} = @{$atoms[$x]}[5];

# atom_type (string) = the SYBYL atom type for the atom
push @{$atom[$x]{'atom_type'}}, split('\.', @{$atoms[$x]}[6]);

# subst_id (integer) = the ID number of the substructure containing the
# atom.
$atom[$x]{'subst_id'} = @{$atoms[$x]}[7];

# subst_name (string) = the name of the substructure containing the atom
$atom[$x]{'subst_name'} = @{$atoms[$x]}[8];

# charge (real) = the charge associated with the atom
$atom[$x]{'charge'} = @{$atoms[$x]}[9];

# status_bit (string) = the internal SYBYL status bits associated with the
# atom
$atom[$x]{'status_bit'} = @{$atoms[$x]}[10];

$x++;
}

#print STDERR Dumper \$atom[0];
return @atom;
}

1;
