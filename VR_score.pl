#!/usr/bin/perl
#
# scoring funtion implemented with Perl
#
# VR


use warnings;
#use Data::Dumper;


use read_mol2;
use get_atom_parameter;
use move;
#use sasa;


##Scoring weights:
#%Weight  = (
#          #'Contact' => -1.7e-1,
#          'Repulsion' => -1.5e-1,
#          'Gauss1' => -6.3e-2,
#          #'Hydrophobic' => 4.3e-2,
#          #'Hydrophobic1' => 1,
#          #'Hydrophobic2' => 1,
#          'Hydrophobic3' => 5.5e-2,
#          #'Hydrogen1' => 3.6,
#          #'Hydrogen11' => 1,
#          'Hydrogen12' => 1.3,
#          #'Hydrogen13' => 1,
#          'Hydrogen2' => 1.3,
#          #'Hydrogen21' => 1,
#          #'Hydrogen22' => 1.8,
#          'Hydrogen3' => 2.6e-1,
#          #'Gap' => 1.8e-2,
#          'Clash' => 1,
#          'Charge' => 9.8e+1, #negative means good
#          'MWs' => 1,
#          #'Combined' => 1
#           );


#minimum alowed atom distance (overlaping if -)
$colision = -1.28;


#take file names
if ((@ARGV == 3) || (@ARGV == 2)) {
$protein = $ARGV[0];
$ligand = $ARGV[1];
@move = (0, 0, 0);
@move = split /:/, $ARGV[2] if ($ARGV[2]);
} else {
die("usage: VR_score.pl protein.mol2 ligand.mol2");
}


#parse files
print STDERR "reading protein $protein ...\n";
($head, $atom, $atom_matrix, $foot) = read_mol2::read_mol2($protein);
#@protein_head = @$head;
@protein_atom = @$atom;
@protein_atom_matrix = @$atom_matrix;
#@protein_foot = @$foot;
print STDERR "reading ligand $ligand ...\n";
($head, $atom, $atom_matrix, $foot) = read_mol2::read_mol2($ligand);
@ligand_head = @$head;
@ligand_atom = @$atom;
@ligand_atom_matrix = @$atom_matrix;
@ligand_foot = @$foot;
#that case probably mol2 files mixed places:
die("Ligand has 300+ atoms. Usage: VR_score.pl protein.mol2 ligand.mol2") if (@ligand_atom > 300);
#print STDERR Dumper \@ligand_head, \@ligand_foot;


#main function
%toppoints = score();
while (($key, $value) = each %toppoints)
{
  print "$key", ",";
}
print "\n";
while (($key, $value) = each %toppoints)
{
  print $value, ",";
}
#print STDERR "\n";
#print just one term for debug:
print STDERR " MWs: $toppoints{'MWs'} \n";


$main = 0;
#move away check score
while ($main < $move[2]) {
#score ligand pose in protein and write if better
%points = score();
while (($key, $value) = each %points)
{
  print "$key", " ", "\n";
}
while (($key, $value) = each %points)
{
  print $value, "\n";
}
@ligand_atom_matrix = move::random_move(@ligand_atom_matrix);
if (($points{'Combined'} < $toppoints{'Combined'}) && ($points{'Clash'} < 1) && ($points{'Repulsion'} > -30)) {
 print STDERR "\n writing ligand...\n\n";
 write_ligand();
 %toppoints = %points;
}
#print STDERR Dumper \@ligand_atom_matrix[0];
$main++;
}
#write_ligand();


##############################
#
# subroutines here
#
##############################


#distance between atoms
#this one must be called for for each atom pair
#FIXME has to be faster
sub distance_sqared {
my $dxs = $_[0][0]-$_[1][0];
my $dys = $_[0][1]-$_[1][1];
my $dzs = $_[0][2]-$_[1][2];
return $dxs*$dxs+$dys*$dys+$dzs*$dzs;
}


#scoring function
sub score {
#initial score
my @d = ();
%score  = (
          'Contact' => 0,
#          'Repulsion' => 0,
#          'Gap' => 0,
          'Hydrophobic1' => 0,
#          'Hydrophobic2' => 0,
          'Hydrophobic3' => 0,
#          'Hydrophobic4' => 0,
          'Hydrophobic5' => 0,
#          'Hydrophobic6' => 0,
#          'Hydrophobic7' => 0,
#          'Hydrophobic8' => 0,
          'Hydrophobic9' => 0,
#          'Hydrophobic10' => 0,
#          'Hydrogen11' => 0,
#          'Hydrogen12' => 0,
#          'Hydrogen13' => 0,
#          'Hydrogen14' => 0,
#          'Hydrogen15' => 0,
#          'Hydrogen16' => 0,
          'Hydrogen21' => 0,
#          'Hydrogen22' => 0,
#          'Hydrogen23' => 0,
          'Hydrogen24' => 0,
          'Hydrogen25' => 0,
          'Hydrogen26' => 0,
#          'Hydrogen31' => 0,
#          'Hydrogen32' => 0,
#          'Hydrogen33' => 0,
#          'Hydrogen34' => 0,
#          'Hydrogen35' => 0,
#          'Hydrogen36' => 0,
#          'HH_31' => 0,
          'HH_32' => 0,
#          'HH_33' => 0,
          'HH_34' => 0,
#          'HH_35' => 0,
#          'HH_36' => 0,
#          'HH2_31' => 0,
#          'HH2_32' => 0,
#          'HH2_33' => 0,
#          'HH2_34' => 0,
#          'HH2_35' => 0,
#          'HH2_36' => 0,
#          'lHH1' => 0,
          #'lHH2' => 0,
          #'lHH2_33' => 0,
          #'lHH2_34' => 0,
          #'lHH2_35' => 0,
          #'lHH2_36' => 0,
          'Gauss1' => 0,
          'Gauss2' => 0,
          'Gauss3' => 0,
#          'Charge' => 0,
#          'Clash' => 0,
          'MWs' => 0,
#          'noContact1' => 0,
#          'noContact2' => 0,
#          'noContact3' => 0,
#          'ZH_31' => 0,
          'ZH_32' => 0,
          'ZH_33' => 0,
#          'ZH_34' => 0,
#          'ZH_35' => 0,
#          'ZH_36' => 0,
          'SO1' => 0,
#          'SO2' => 0,
#          'SO3' => 0,
#          'SO4' => 0,
#          'SO5' => 0,
#          'SO6' => 0,
#          'SS' => 0,
          'AA1' => 0,
#          'AA2' => 0,
#          'AA3' => 0,
#          'AA4' => 0,
          'AA5' => 0,
#          'AA6' => 0,
#          'AN1' => 0,
#          'AN2' => 0,
#          'AN3' => 0,
          'AN4' => 0,
          'AN5' => 0,
          'AN6' => 0,
          'C3' => 0,
          'X3' => 0,
          'N3' => 0,
          'O3' => 0,
          #'SASA1' => 0,
          #'SASA2' => 0,
          #'SASA3' => 0,
          #'Combined' => 1.8 #initial value
           );


#MW of ligand
#respect substructures
$x = 0;
%ids = ();
while ($ligand_atom[$x]{'atom_type'}[0]) {
  #print STDERR $ligand_atom[$x]{'subst_id'};
  $score{'MWs'} += get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'MW');
  $score{'C3'} ++ if (($ligand_atom[$x]{'atom_type'}[0] eq 'C') && ($ligand_atom[$x]{'atom_type'}[1] eq '3'));
  $score{'N3'} ++ if (($ligand_atom[$x]{'atom_type'}[0] eq 'N') && ($ligand_atom[$x]{'atom_type'}[1] eq '3'));
  $score{'O3'} ++ if (($ligand_atom[$x]{'atom_type'}[0] eq 'O') && ($ligand_atom[$x]{'atom_type'}[1] eq '3'));
  $score{'X3'} ++ if ($ligand_atom[$x]{'atom_type'}[1] eq '3');
  $id = $ligand_atom[$x]{'subst_id'};
  $ids{$id} = 1;
  $x++;
}
#print STDERR Dumper \%ids;
$SS = keys %ids;
#$score{'SS'} = $SS;
$score{'MWs'} = $score{'MWs'}/$SS;
die ("multilig!") if ($SS != 1);


#electrostatic force, d matrix, Gauss and repulsion
#d=distance-R_ligand_atom-R_protein_atom
$x=0;
while($ligand_atom[$x]{'atom_type'}[0]) {
 $y=0;
 my $lig_radius = get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'radius');
 die("$ligand_atom[$x]{'atom_type'}[0] radius undefined!") if (!($lig_radius));
 die ("ligand_atom_matrix element $x undefined!") if !($ligand_atom_matrix[$x]);
 while($protein_atom[$y]) {
 $d[$x][$y] = distance_sqared($ligand_atom_matrix[$x],$protein_atom_matrix[$y]);
 if ($d[$x][$y] < 100) { #skip very distant atom pairs
#  $score{'Charge'} += $ligand_atom[$x]{'charge'}*$protein_atom[$y]{'charge'}/$d[$x][$y];
  $d[$x][$y] = sqrt($d[$x][$y]) 
  - $lig_radius
  - get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'radius');
  #calculate Gauss1 and Gauss2
  if ($d[$x][$y] < 2) {
   $score{'Gauss1'} += exp(-2*($d[$x][$y]**2)); #means if abs distance 0 then +1 else +less
   $score{'Gauss2'} += exp(-4*($d[$x][$y]**2));
   $score{'Gauss3'} += exp(-8*($d[$x][$y]**2));
#   $score{'Gap'}++ if ($d[$x][$y] > 0); # count just gap
   $score{'Contact'}++ if (abs($d[$x][$y]) < 0.25);
    if ($d[$x][$y] < 0) {
     #calculate repulsion:
#     $score{'Repulsion'}++;
     if ($d[$x][$y] < $colision) {
#          $score{'Clash'}++;
          #print STDERR "colision detected! $d[$x][$y]\n";
     }
     #print STDERR "\n $ligand_atom[$x]{'atom_id'} repeals $protein_atom[$y]{'atom_id'}\n";
    }
  }
 }
 $y++;
 }
$x++;
}
#print STDERR "d size ", $x, " ", $y, "\n";
#print STDERR "x: ", $ligand_atom[2]{'atom_type'}[0], " y: ", $protein_atom[3147]{'atom_type'}[0], " d:", $d[2][3147], "\n";


#calculate hydrophobic
$x = 0;
while($d[$x]) {
 $y = 0;
 if (get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'hydrophobic')) {
  while($d[$x][$y]) {
   if (($colision < $d[$x][$y]) && ($d[$x][$y] < 2)) {
   if (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'hydrophobic')) {
    #if ($d[$x][$y] < 0.5) {
      $score{'Hydrophobic1'}++;
#      $score{'Hydrophobic2'}++ if ($d[$x][$y] < 0.25);
      $score{'Hydrophobic3'}++ if (abs($d[$x][$y]) < 0.25);
#      $score{'Hydrophobic4'}++ if ($d[$x][$y] < -0.25);
      $score{'Hydrophobic5'}++ if ($d[$x][$y] < 0);
#      $score{'Hydrophobic6'}++ if ($d[$x][$y] > 0.25);
#      $score{'Hydrophobic7'}++ if ($d[$x][$y] > 0.5);
#      $score{'Hydrophobic8'}++ if ($d[$x][$y] > 0.75);
      $score{'Hydrophobic9'}++ if ($d[$x][$y] > 1);
#      $score{'Hydrophobic10'}++ if ($d[$x][$y] > 1.25);
   }
   }
   $y++;
  }
 }
$x++;
}


#ligand as hydrogen bond donor
$x = 0;
while($d[$x]) {
 $y = 0;
 $count = 0;
 if (($ligand_atom[$x]{'charge'} > 0.1) && ($ligand_atom[$x]{'atom_type'}[0] eq 'H')) {
  while($d[$x][$y]) {
   if (($colision < $d[$x][$y]) && ($d[$x][$y] < 2)) {
   if (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'H_acceptor')) {
      $count++ if ($d[$x][$y] < 0.25);
#      $score{'Hydrogen11'}++;
#      $score{'Hydrogen12'}++ if ($d[$x][$y] < 0.25);
#      $score{'Hydrogen13'}++ if (abs($d[$x][$y]) < 0.25);
#      $score{'Hydrogen14'}++ if ($d[$x][$y] < -0.25);
#      $score{'Hydrogen15'}++ if ($d[$x][$y] < 0);
#      $score{'Hydrogen16'}++ if ($d[$x][$y] > 0.25);
   }
   }
   $y++;
  }
 }
# $score{'lHH1'}++ if ($count > 1);
 $x++;
}


#ligand as hydrogen bond acceptor
$x = 0;
while($d[$x]) {
 $y = 0;
 if (get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'H_acceptor')) {
  while($d[$x][$y]) {
   if (($colision < $d[$x][$y]) && ($d[$x][$y] < 2)) {
   if (($protein_atom[$y]{'charge'} > 0.1) && ($protein_atom[$y]{'atom_type'}[0] eq 'H')) {
     $score{'Hydrogen21'}++;
#     $score{'Hydrogen22'}++ if ($d[$x][$y] < 0.25);
#     $score{'Hydrogen23'}++ if (abs($d[$x][$y]) < 0.25);
     $score{'Hydrogen24'}++ if ($d[$x][$y] < -0.25);
     $score{'Hydrogen25'}++ if ($d[$x][$y] < 0);
     $score{'Hydrogen26'}++ if ($d[$x][$y] > 0.25);
   }
   }
   $y++;
  }
 }
$x++;
}


##non polar H
#$x = 0;
#while($d[$x]) {
# $y = 0;
# if ((abs($ligand_atom[$x]{'charge'}) < 0.1) && ($ligand_atom[$x]{'atom_type'}[0] eq 'H')) {
#  while($d[$x][$y]) {
#   if (($colision < $d[$x][$y]) && ($d[$x][$y] < 2)) {
#   if ((abs($protein_atom[$y]{'charge'}) < 0.1) && ($protein_atom[$y]{'atom_type'}[0] eq 'H')) {
#     $score{'Hydrogen31'}++;
#     $score{'Hydrogen32'}++ if ($d[$x][$y] < 0.25);
#     $score{'Hydrogen33'}++ if (abs($d[$x][$y]) < 0.25);
#     $score{'Hydrogen34'}++ if ($d[$x][$y] < -0.25);
#     $score{'Hydrogen35'}++ if ($d[$x][$y] < 0);
#     $score{'Hydrogen36'}++ if ($d[$x][$y] > 0.25);
#   }
#   }
#   $y++;
#  }
# }
#$x++;
#}


#hydrophobic - hydrophilic
$x = 0;
while($d[$x]) {
 $y = 0;
 if ((abs($ligand_atom[$x]{'charge'}) < 0.1) || (get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'hydrophobic'))) {
  while($d[$x][$y]) {
   if (($colision < $d[$x][$y]) && ($d[$x][$y] < 2)) {
   if ((abs($protein_atom[$y]{'charge'}) > 0.1) || !(get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'hydrophobic'))) {
#     $score{'HH_31'}++;
     $score{'HH_32'}++ if ($d[$x][$y] < 0.25);
#     $score{'HH_33'}++ if (abs($d[$x][$y]) < 0.25);
     $score{'HH_34'}++ if ($d[$x][$y] < -0.25);
#     $score{'HH_35'}++ if ($d[$x][$y] < 0);
#     $score{'HH_36'}++ if ($d[$x][$y] > 0.25);
   }
   }
   $y++;
  }
 }
$x++;
}


#hydrophilic - hydrophobic
$x = 0;
while($d[$x]) {
 $y = 0;
 if ((abs($ligand_atom[$x]{'charge'}) > 0.1) || !(get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'hydrophobic'))) {
  while($d[$x][$y]) {
   if (($colision < $d[$x][$y]) && ($d[$x][$y] < 2)) {
   if ((abs($protein_atom[$y]{'charge'}) < 0.1) || (get_atom_parameter::get_atom_parameter($protein_atom[$y]{'atom_type'}[0], 'hydrophobic'))) {
#     $score{'HH2_31'}++;
#     $score{'HH2_32'}++ if ($d[$x][$y] < 0.25);
#     $score{'HH2_33'}++ if (abs($d[$x][$y]) < 0.25);
#     $score{'HH2_34'}++ if ($d[$x][$y] < -0.25);
#     $score{'HH2_35'}++ if ($d[$x][$y] < 0);
#     $score{'HH2_36'}++ if ($d[$x][$y] > 0.25);
   }
   }
   $y++;
  }
 }
$x++;
}


#Zn
$x = 0;
while($d[$x]) {
 $y = 0;
 if ((abs($ligand_atom[$x]{'charge'}) > 0.1) || !(get_atom_parameter::get_atom_parameter($ligand_atom[$x]{'atom_type'}[0], 'hydrophobic'))) {
  while($d[$x][$y]) {
   if (($colision < $d[$x][$y]) && ($d[$x][$y] < 2)) {
   if ($protein_atom[$y]{'atom_type'}[0] eq 'Zn') {
#     $score{'ZH_31'}++;
#very good score for HCA:
     $score{'ZH_32'}++ if ($d[$x][$y] < 0.25);
     #print STDERR $ligand_atom[$x]{'atom_type'}[0];
     $score{'ZH_33'}++ if (abs($d[$x][$y]) < 0.25);
#     $score{'ZH_34'}++ if ($d[$x][$y] < 0.2);
#     $score{'ZH_35'}++ if ($d[$x][$y] < 0.1);
#     $score{'ZH_36'}++ if ($d[$x][$y] < 0.05);
   }
   }
   $y++;
  }
 }
$x++;
}


#S.O2
$x = 0;
while($d[$x]) {
 $y = 0;
 if (($ligand_atom[$x]{'atom_type'}[0] eq 'S') && (($ligand_atom[$x]{'atom_type'}[1] eq 'O2') || ($ligand_atom[$x]{'atom_type'}[1] == 3))) {
  while($d[$x][$y]) {
   if (($colision < $d[$x][$y]) && ($d[$x][$y] < 2)) {
   if ((abs($protein_atom[$y]{'charge'}) > 0.1) && ($protein_atom[$y]{'atom_type'}[0] eq 'H')) {
     $score{'SO1'}++;
#     $score{'SO2'}++ if ($d[$x][$y] < 0.25);
#     $score{'SO3'}++ if (abs($d[$x][$y]) < 0.25);
#     $score{'SO4'}++ if ($d[$x][$y] < -0.25);
#     $score{'SO5'}++ if ($d[$x][$y] < 0);
#     $score{'SO6'} = $score{'SO5'}/$SS;
   }
   }
   $y++;
  }
 }
$x++;
}


#ar.-ar.
$x = 0;
while($d[$x]) {
 $y = 0;
 if ($ligand_atom[$x]{'atom_type'}[1] eq 'ar') {
  while($d[$x][$y]) {
   if (($colision < $d[$x][$y]) && ($d[$x][$y] < 2)) {
   if ($protein_atom[$y]{'atom_type'}[1] eq 'ar') {
     $score{'AA1'}++;
#     $score{'AA2'}++ if ($d[$x][$y] < 0.25);
#     $score{'AA3'}++ if (abs($d[$x][$y]) < 0.25);
#     $score{'AA4'}++ if ($d[$x][$y] < -0.25);
     $score{'AA5'}++ if ($d[$x][$y] < 0);
#     $score{'AA6'}++ if ($d[$x][$y] > 0.25);
   }
   }
   $y++;
  }
 }
$x++;
}


# N - ar.
$x = 0;
while($d[$x]) {
 $y = 0;
 if (($ligand_atom[$x]{'atom_type'}[0] eq 'N')) {
  while($d[$x][$y]) {
   if (($colision < $d[$x][$y]) && ($d[$x][$y] < 2)) {
   if ($protein_atom[$y]{'atom_type'}[1] eq 'ar') {
#     $score{'AN1'}++;
#     $score{'AN2'}++ if ($d[$x][$y] < 0.25);
#     $score{'AN3'}++ if (abs($d[$x][$y]) < 0.25);
     $score{'AN4'}++ if ($d[$x][$y] < -0.25);
     $score{'AN5'}++ if ($d[$x][$y] < 0);
     $score{'AN6'}++ if ($d[$x][$y] > 0.25);
   }
   }
   $y++;
  }
 }
$x++;
}


##ligand atom no contact
#$x = 0;
#while($d[$x]) {
# $y = 0;
# $contact1 = 0;
# $contact2 = 0;
# $contact3 = 0;
#   while($d[$x][$y]) {
#    if (($colision < $d[$x][$y]) && ($d[$x][$y] < 0)) {
#     $contact1++;
#    }
#    $contact2++ if ($d[$x][$y] < -0.25);
#    $contact3++ if ($d[$x][$y] < -0.5);
#  $y++;
#  }
# $score{'noContact1'}++ if !($contact1);
# $score{'noContact2'}++ if !($contact2);
# $score{'noContact3'}++ if !($contact3);
#$x++;
#}


##dSASA
#$x = 0;
#while ($ligand_atom[$x]{'atom_type'}[0]) {
#  push @ligand_atom_type, $ligand_atom[$x]{'atom_type'}[0];
#  $x++;
#}
#$x = 0;
#$maxd = 6;
##supaprastinimas: surinkti tik tuos balt atom kurie arti ligando
#while ($protein_atom[$x]{'atom_type'}[0]) {
#  $y = 0;
#    while ($ligand_atom[$y]{'atom_type'}[0]) {
#      if ($d[$y][$x] < $maxd) {
#        push (@protein_atom_type, $protein_atom[$x]{'atom_type'}[0]);
#        last;
#      } 
#      $y++;
#      push (@protein_atom_type, 0); #signal to skip it
#    }
#  $x++;
#}
#@atom_matrix = (@ligand_atom_matrix, @protein_atom_matrix);
#@atom_type = (@ligand_atom_type, @protein_atom_type);
#$score{'SASA1'} = sasa::sasa(\@atom_matrix, \@atom_type) -  sasa::sasa(\@ligand_atom_matrix, \@ligand_atom_type)
#- sasa::sasa(\@protein_atom_matrix, \@protein_atom_type);
##VERY slow:- sasa::sasa(\@protein_atom_matrix, \@protein_atom_type);


##* by Weight and calculate combined score
foreach my $key ( keys %score )
{
   #$score{$key} = $score{$key}/$SS;
   #if ($key ne 'Combined') {
   # $score{$key} *= $Weight{$key};
   # $score{'Combined'} += $score{$key};
   #}
}


return %score;
}


#write edited ligand as *.mol2 file
sub write_ligand {
#gather lines
my @lines;
foreach my $line (@ligand_head) {
 push @lines, $line;
}
$x = 0;
while ($ligand_atom[$x]{'atom_type'}[0]) {
$ligand_atom[$x]{'atom_type'}[1]=" " if (!$ligand_atom[$x]{'atom_type'}[1]);
 $line =    ' '.
            $ligand_atom[$x]{'atom_id'}.' '.
            $ligand_atom[$x]{'atom_name'}.' '.
            $ligand_atom_matrix[$x][0].' '.
            $ligand_atom_matrix[$x][1].' '.
            $ligand_atom_matrix[$x][2].' '.
#            $ligand_atom[$x]{'status_bit'}.' '.
            $ligand_atom[$x]{'atom_type'}[0].'.'.
            $ligand_atom[$x]{'atom_type'}[1].' '.
            $ligand_atom[$x]{'subst_id'}.' '.
            $ligand_atom[$x]{'subst_name'}.' '.
            $ligand_atom[$x]{'charge'};
 #print STDERR "$line\n";
 push @lines, $line;
 $x++;
}
foreach my $line (@ligand_foot) {
 push @lines, $line;
}

#ligand file name
$name = ">out_$points{'Combined'}_$ligand";


#print them
open (LIGAND, $name) or die('can not write ligand output file');
foreach my $line (@lines) {
 print LIGAND $line;
}
close(LIGAND);
}
