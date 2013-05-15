VR_score
========

scoring function implemented with Perl

File preparation
========

Ligand:
Add hydrogens with Avogadro save as *.mol2

Protein:
Chimera  do Dock prep and save as *.mol2

Test use
========

````bash
 ./VR_score.pl 2HOC_.mol2 ligH.mol2
````

Moving ligand
========
````bash
 ./VR_score.pl 2HOC_.mol2 ligH.mol2 x:p
````
Where:
x direction (x, y, z)
p sign (p, n)

