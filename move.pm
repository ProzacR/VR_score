package move;

use warnings;
use Data::Dumper;


#make random move
#use random_move(ligand_atom_matrix)
sub random_move {
#print STDERR Dumper \@_;
$move[0] = int(rand(3));
if(int(rand(2))) { #so 50/50 chance
 $move[1] = 1;
} else {
 $move[1] = -1;
}
#print STDERR Dumper \@move;
@ligand_atom_matrix = move::move_ligand(\@_, $move[0], $move[1]);
return @ligand_atom_matrix;
}


#move ligand
#use ex. move_ligand(ligand_ref, x y or z (0 1 or 2), p or n (-1 or 1))
sub move_ligand {
my $ref = $_[0];
my $x = 0;
while ($$ref[$x][$_[1]]) {
 # + or - direction?
 $$ref[$x][$_[1]] += $_[2];
 $x++;
}
return @$ref;
}


1;
