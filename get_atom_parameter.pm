package get_atom_parameter;

use warnings;
use Data::Dumper;


sub get_atom_parameter {
#use:
# get_atom_parameter(name, parameter);

my @atom;
my $x;

#take correct value from data table
$x = 0 if($_[1] eq "radius");
$x = 1 if($_[1] eq "depth");
$x = 2 if($_[1] eq "solvation_parameter");
$x = 3 if($_[1] eq "volume");
$x = 4 if($_[1] eq "covalent_radius");
$x = 5 if($_[1] eq "hydrophobic");
$x = 6 if($_[1] eq "H_acceptor");
$x = 7 if($_[1] eq "MW");

 
# {'name'}= radius, depth, solvation_parameter, volume, covalent_radius, hydrophobic, H_acceptor, MW
$atom{'C'}= [2.00000, 0.15000, -0.00143, 33.51030, 0.77, 1, 0, 12];
$atom{'A'}= [2.00000, 0.15000, -0.00052, 33.51030, 0.77, 0, 0, ''];
$atom{'N'}= [1.75000, 0.16000, -0.00162, 22.44930, 0.75, 0, 1, 14]; 
$atom{'O'}= [1.60000, 0.20000, -0.00251, 17.15730, 0.73, 0, 1, 16];
$atom{'P'}= [2.10000, 0.20000, -0.00110, 38.79240, 1.06, 0, 0, 30];
$atom{'S'}= [2.00000, 0.20000, -0.00214, 33.51030, 1.02, 0, 1, 32];
$atom{'H'}= [1.00000, 0.02000, 0.00051, 0.00000, 0.37, 0, 0, 1];
$atom{'F'}= [1.54500, 0.08000, -0.00110, 15.44800, 0.71, 1, 0, 19];
$atom{'I'}= [2.36000, 0.55000, -0.00110, 55.05850, 1.33, 1, 0, 127];
$atom{'NA'}= [1.75000, 0.16000, -0.00162, 22.44930, 0.75, 0, 0, 14];
$atom{'OA'}= [1.60000, 0.20000, -0.00251, 17.15730, 0.73, 0, 0, 16];
$atom{'SA'}= [2.00000, 0.20000, -0.00214, 33.51030, 1.02, 0, 0, 32];
$atom{'HD'}= [1.00000, 0.02000, 0.00051, 0.00000, 0.37, 0, 0, 1];
$atom{'Mg'}= [0.65000, 0.87500, -0.00110, 1.56000, 1.30, 0, 0, 24];
$atom{'Mn'}= [0.65000, 0.87500, -0.00110, 2.14000, 1.39, 0, 0, 55];
$atom{'Zn'}= [0.74000, 0.55000, -0.00110, 1.70000, 1.31, 0, 1, 65];
$atom{'Ca'}= [0.99000, 0.55000, -0.00110, 2.77000, 1.74, 0, 0, 40];
$atom{'Fe'}= [0.65000, 0.01000, -0.00110, 1.84000, 1.25, 0, 0, 56];
$atom{'Cl'}= [2.04500, 0.27600, -0.00110, 35.82350, 0.99, 1, 0, 35.5];
$atom{'Br'}= [2.16500, 0.38900, -0.00110, 42.56610, 1.14, 1, 0, 80];


return $atom{$_[0]}[$x];
}

1;
