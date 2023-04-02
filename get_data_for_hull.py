import requests
import numpy as np

def get_formation_energies(species_list,model=None):
    """
    Query OpenKIM repository for test results or reference data for constructing a thermodynamic hull.
    Args:
        species_list:
            List of all species to query for. Every result or data containing any combination of one or more of these will be returned
        model:
            OpenKIM model to query for. If this is None, reference data are queried for instead

    Returns:
        * (Nxlen(species_list)) ndarray representing N points to buld the hull out of. The last column is the energy, the others are the molar fractions of each element in species_list except the first
        * N-length list of prototype labels
    """
    assert len(species_list>1)
    if model==None:
        model_query=''
        meta_type='rd'
    else:
        model_query='"meta.subject.extended-id":"%s"'%model
        meta_type='tr'
    result = requests.post("https://query.openkim.org/api",data={'query': '{"meta.type":"%s",%s,"property-id":"tag:staff@noreply.openkim.org,2023-02-21:property/binding-energy-crystal","stoichiometric-species.source-value":{"$in":%s}}'%(meta_type,model_query,str(species_list).replace("'",'"')), 'fields': '{"prototype-label.source-value":1,"stoichiometric-species.source-value":1,"binding-potential-energy-per-formula.source-value":1}', 'database': 'data'}).json()
    elemental_references = {}
    # initialize dict for elemental references
    for species in species_list:
        elemental_references[species]=0.
    for entry in result:
        species=entry["stoichiometric-species"]["source-value"]
        if len(species)==1:
            energy = entry["binding-potential-energy-per-formula"]["source-value"]
            species = species[0]
            if energy < elemental_references[species]:
                elemental_references[species] = energy
    print(elemental_references)
    #now build the hull list of dicts
    #each entry:
    # array of fractions of each element
    # formation energy per atom
    # prototype label
    for entry in result:
    

            