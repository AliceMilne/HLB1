
This code was written to describe the epidemiology of the serious disease in citus -- Huanglongbing disease (HLB).  
This is an acute bacterial disease that threatens the sustainability of citrus production across the world.  
For the disease to be controlled sufficient numbers of growers must adopt the control - therefore behaviour 
becomes important. To this end we linked our epidemiological model with a model that uses opinion dynamics to 
simulate the behaviour of individual growers. This work is published in

Alice E. Milne, Tim Gottwald, Stephen R. Parnell, Vasthi Alonso Chavez, and Frank van den Bosch, 2020 
What makes or breaks a campaign to stop an invading plant pathogen? PLOS Comp Biol 

The code was produced using Microsoft visual studio but it should be straightforward to compile with other systems.

We have used NAG libraries to randomly generate nummbers according to specified distributions and perform integration (see www.nag.co.uk).
Calls to nag are of forms similar to "G05KGF" and instructions on what this part of the code does can be found on Nags website 
(e.g. https://www.nag.co.uk/numeric/fl/nagdoc_fl24/). If you do not have access to NAG libraries then it should be straightforward 
to replace these calls with a different random number generators.  

Please contact us with any quieries or problems 

- the authors

