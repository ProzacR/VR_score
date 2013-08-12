package move;

use warnings;
use Data::Dumper;
use Math::Trig;


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
if(int(rand(2))) { #so 50/50 chance rotate or move
 print STDERR "done\n";
 @ligand_atom_matrix = move::move_ligand(\@_, $move[0], $move[1]);
} else {
 print STDERR "rotation\n";
 @ligand_atom_matrix = move::rotate_ligand(\@_, $move[0]);
}

return @ligand_atom_matrix;
}


#move ligand
#use ex. move_ligand(ligand_ref, x y or z (0 1 or 2), p or n (-1 or 1))
sub move_ligand {
my $ref = $_[0];
my $x = 0;
while ($$ref[$x][$_[1]]) {
 $$ref[$x][$_[1]] += $_[2];
 $x++;
}
#print STDERR Dumper \@$ref;
return @$ref;
}


#rotate ligand
#use ex. rotate(ligand_ref, x y or z (0 1 or 2))
sub rotate_ligand {
my $ref = $_[0];
my $x = 0;
my $step = pi/8;

#print STDERR Dumper \@$ref;
@rotation_point = ($$ref[$x][0], $$ref[$x][1], $$ref[$x][2]);
#move to 0,0,0:
while (defined($$ref[$x][0])) {
$$ref[$x][0] -= $rotation_point[0];
$$ref[$x][1] -= $rotation_point[1];
$$ref[$x][2] -= $rotation_point[2];
$x++;
}
print STDERR Dumper \@$ref;


$x = 0;
#if rotx
if ($_[1] == 0) {
while (defined($$ref[$x][0])) {
 $$ref[$x][1] = ($$ref[$x][1]*cos($step) - $$ref[$x][2]*sin($step));
 $$ref[$x][2] = ($$ref[$x][1]*sin($step) + $$ref[$x][2]*cos($step));
 $x++;
}
print STDERR "rotx\n";
}


#if roty
if ($_[1] == 1) {
while (defined($$ref[$x][0])) {
 $$ref[$x][0] = ($$ref[$x][0]*cos($step) + $$ref[$x][2]*sin($step));
 $$ref[$x][2] = (-$$ref[$x][0]*sin($step) + $$ref[$x][2]*cos($step)); 
 $x++;
}
print STDERR "roty\n";
}


#if rotz
if ($_[1] == 2) {
while (defined($$ref[$x][0])) {
 $$ref[$x][0] = ($$ref[$x][0]*cos($step) - $$ref[$x][1]*sin($step));
 $$ref[$x][1] = ($$ref[$x][0]*sin($step) + $$ref[$x][1]*cos($step)); 
 $x++;
}
print STDERR "rotz\n";
}


print STDERR Dumper \@$ref;
#move to initial location:
$x = 0;
while (defined($$ref[$x][0])) {
$$ref[$x][0] += $rotation_point[0];
$$ref[$x][1] += $rotation_point[1];
$$ref[$x][2] += $rotation_point[2];
$x++;
}


print STDERR Dumper \@$ref;
return @$ref;
}

1;
