#!/usr/bin/perl
#
# simple *.mol2 file reader
# http://www.tripos.com/data/support/mol2.pdf
#
# VR

#use warnings;
use Data::Dumper;


#take file name
if (@ARGV == 1) {
$file = $ARGV[0];
} else {
die("usage: read_mol2.pl file.mol2");
}


#read file into lines
open(INFO, $file) or die("Could not open file.");
while(<INFO>) {
 push (@lines, $_);
}
close(INFO);
#print Dumper \@lines;


$line_number=0;
#skip to first @<TRIPOS>MOLECULE
while(!(@lines[$line_number] =~ /MOLECULE/)) {
#just ignore comments and empty lines
#print STDERR "comment: ", @lines[$line_number];
$line_number++;
}
$line_number++;
#read MOLECULE part to STDERR now
while(!(@lines[$line_number] =~ /ATOM/)) {
#print STDERR "molecule: ", @lines[$line_number];
$line_number++;
}
$line_number++;
$x=0;
#read ATOM part to data structure
while(!(@lines[$line_number] =~ /BOND/)) {
#print STDERR "atom: ", @lines[$line_number];
#split lines into symbols
push @{$atoms[$x]}, split(//, $lines[$line_number]);
$x++;
$line_number++;
}
$line_number++;
print Dumper \@atoms[0];
#ignore BOND part for now....


#make formated hash: FIXME
$x = 0;
while($atoms[$x]) {
#COLUMNS DATA TYPE CONTENTS
#--------------------------------------------------------------------------------
# 1 - 6 Record name "ATOM "
$atom[$x]{'record'} = join "", @{$atoms[$x]}[0..5];

# 7 - 11 Integer Atom serial number.
$atom[$x]{'number'} = join "", @{$atoms[$x]}[6..10];

#13 - 16 Atom Atom name.
$atom[$x]{'name'} = join "", @{$atoms[$x]}[12..15];

#17 Character Alternate location indicator.
$atom[$x]{'location'} = join "", @{$atoms[$x]}[16];

#18 - 20 Residue name Residue name.
$atom[$x]{'resname'} = join "", @{$atoms[$x]}[17..19];

#22 Character Chain identifier.
$atom[$x]{'chain'} = join "", @{$atoms[$x]}[21];

#23 - 26 Integer Residue sequence number.
$atom[$x]{'res_nr'} = join "", @{$atoms[$x]}[22..25];

#27 AChar Code for insertion of residues.
$atom[$x]{'ins_code'} = join "", @{$atoms[$x]}[26];

#31 - 38 Real(8.3) Orthogonal coordinates for X in Angstroms.
$atom[$x]{'x'} = join "", @{$atoms[$x]}[30..37];

#39 - 46 Real(8.3) Orthogonal coordinates for Y in Angstroms.
$atom[$x]{'y'} = join "", @{$atoms[$x]}[38..45];

#47 - 54 Real(8.3) Orthogonal coordinates for Z in Angstroms.
$atom[$x]{'z'} = join "", @{$atoms[$x]}[46..53];

#55 - 60 Real(6.2) Occupancy.
$atom[$x]{'occupancy'} = join "", @{$atoms[$x]}[54..59];

#61 - 66 Real(6.2) Temperature factor (Default = 0.0).
$atom[$x]{'temp_factor'} = join "", @{$atoms[$x]}[60..65];

#73 - 76 LString(4) Segment identifier, left-justified.
$atom[$x]{'segment'} = join "", @{$atoms[$x]}[72..75];

#77 - 78 LString(2) Element symbol, right-justified.
$atom[$x]{'element'} = join "", @{$atoms[$x]}[76..77];

#79 - 80 LString(2) Charge on the atom.
$atom[$x]{'charge'} = join "", @{$atoms[$x]}[78..79];

$x++;
}

#print Dumper \@atom;
