# Explanation for files

## Json files

these were the first storage files for the first results most follow the simple concept: name: results_inp.pk contains Bravais matrix, orbit states and occupation and the first iteration distance, mostly used for unaries and tests

## pickle

arrays or dicts of calcnodes with the corresponding first iteration distance. Used mainy for oxides.

S charge transfer is the charge transfer tests from s1/2 to the p of oxygen for the metals in the s group, p charge transfer for the metals with 1 electron in the p orb (Ga, Al , In etc...) , s charge d contains the nodes for the d metals with half filled s orbitals (Cu , Ag etc)

Copper: s charge transfer s 0.9: set the s1/2 valence orb to 0.9 and transfered charge from d3/2 to oxygen 

s charge transfer d 1.6 set the d3/2 to 1.6 (transfered charge to s1/2) then transfered charge to s1/2


## EXTRAS

I switched to the tag like extras system (from saving the node pk in json/pickle). Exisiting Extras:

for inp calcs: 

Z atom number of the oxides 

element name

configuration what oxide configurations was used eg. XO , XO2 etc


for starting density/charge transfer calculations:

element: element used in oxides

type  what type of structrue eg oxides

calc type what was done with the calculation

transfered charge, amount of charge transfered from X to oxygen( in the case of oxides)

oxide struct what oxide was used eg. XO XO2 etc