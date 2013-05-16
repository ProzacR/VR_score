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
 ./VR_score.pl 2HOC_.mol2 ligH.mol2 x:s:n
````
Where:
x direction (x is 0, y is 1, z is 2)
s sign (-1 or 1)
n number of 1A moves
(now moves random so just ::n)


