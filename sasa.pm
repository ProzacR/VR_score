package sasa;

# SASA using Shrake-Rupley algorithm
# calculate vapor to water free energy of transfer
# using Solvent Accessible Surface Area
# SASA
# VR 2013

use Math::Trig;
#use Data::Dumper;
use warnings;


use get_atom_parameter;

#FIXME nesutampa su get_atom_parameter?!
##solvation parameter
#$solv_par{'H'}=0;
#$solv_par{'C'}=27;
#$solv_par{'N'}=-116;
#$solv_par{'O'}=-116;


#H2O radius:
$prad=1.4;
#tasku sk. aplink atomus (daugiau tiksliau: ~simtai):
$M=500;

sub sasa {
$E = 0;
#print STDERR Dumper \@_;
($atom_matrix, $atom_type) = @_;
#print STDERR $$atom_matrix[0][0];
#atoms example 0,1,...  (x,y,z,radius,solvation parameter):
#@N[0] = [1,2,3,4,5];
#set atom new way, find radius
$x = 0;
while ($$atom_type[$x]) {
  #print STDERR $$atom_type[$x], "\n";
  #$N[$x] = [x, y, z, radii, solv_par];
  $N[$x] = [$$atom_matrix[$x][0], $$atom_matrix[$x][1], $$atom_matrix[$x][2],
  get_atom_parameter::get_atom_parameter($$atom_type[$x], 'radius'),
  get_atom_parameter::get_atom_parameter($$atom_type[$x], 'solvation_parameter')];
  $x++;
}
#print STDERR Dumper \@N;
#print STDERR $N[0][4];


$i = 0;
$irad = 0;
$k = 0;


while($N[$i]) {
 #print "skersmuo: ", $N[$i][3];
 $irad=$N[$i][3] + $prad;
 while ($k<$M) {
  $u=rand();
  $v=rand();
  $theta=2*pi*$u;
  $phi=acos(2*$v-1);
  $x=cos($theta)*sin($phi);
  $y=sin($theta)*sin($phi);
  $z=cos($phi);
  #sukuriamas taskas:
  $point[0]=$N[$i][0]+$x*$irad;
  $point[1]=$N[$i][1]+$y*$irad;
  $point[2]=$N[$i][2]+$z*$irad;
  push @{$pts[$i][$k]}, @point;
 $k++;
 }
$i++;
$k = 0;
}

#tasku visuma:
#print STDERR Dumper \@pts;

#eit per atomus ir istrinti taskus kurie kito atomo kelyje:
$i = 0;
$k = 0;


while($N[$i]) {
 $irad=$N[$i][3] + $prad;
 $Mp=$M;
 while ($k<$M) {
  $fail=0;
  $j=$i+1; #strange here...sometimes jumps over next while at all!!?
  while($N[$j]) {
   $jrad=$N[$j][3] + $prad;
   $r=sqrt(( ($pts[$i][$k][0]-$N[$j][0])+($pts[$i][$k][1]-$N[$j][1])+($pts[$i][$k][2]-$N[$j][2]) )**2);
   if($r <= $jrad) {
    $fail=1;
   }
   $j++;
  }
  if($fail) {
   $Mp--;
  }
 $k++;
 }
 #SASA atomu i:
 $sasa[$i]=4*pi*$irad*$irad*$Mp/$M;
$i++;
$k = 0;
}
#print STDERR "sasa: \n";
#print STDERR Dumper \@sasa;


#solvation energy:
$x=0;
while (@sasa >= $x) {
$E=$E+$sasa[$x]*$N[$x][4];
$x++;
}
#print STDERR "x: ", $x, "\n";
#print STDERR "solvation energy: ", $E/1000, " kcal/mol\n";
return $E;
}


1;
