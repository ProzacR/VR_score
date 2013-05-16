package move;

use warnings;


#move ligand
#use ex. move_ligand('0', 'p')
sub move_ligand {
my $p;
my $x = 0;
while ($ligand_atom_matrix[$x][$_[0]]) {
 # + or - direction?
 $p = 1 if ($_[1] eq 'p');
 $p = -1 if ($_[1] eq 'n');
 $ligand_atom[$x]{$_[0]} += $p;
 $x++;
}
return 1;
}


1;
