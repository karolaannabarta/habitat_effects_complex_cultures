# habitat_effects_complex_cultures
Repository for simulation scripts on the effects of habitat characteristics on the diversity of complex cultural systems.

We argue that changes in habitat structure can influence the diversity of animal cultures by affecting social distancing of individuals. To explore this, we proposed a spatial modelling framework that focuses on complex cultures. It incorporates parameters such as hearing distance, edge effect controlled by habitat shape, and carrying capacity of the habitat represented by population density. Two models were built based on the sample of tutors (the whole population or neighbours only) whom the young individuals acquire their repertoire representing two species-specific learning strategies.

This repository was created to store and share scripts of the simulation exploring the above mentioned scientific questions. The two different model setups relating to how young individuals acquire their initial repertoire have been written separately. The file "complex_culture_spatial_model_allpop.py" contains model I, where young individuals acquire their initial repertoire from the whole population. Model II, in which young individuals acquire their initial repertoire from their neighbours, can be found in "complex_culture_spatial_model_allpop_younginds_rep_from_neighbours.py".

The scripts use multiprocessing.
