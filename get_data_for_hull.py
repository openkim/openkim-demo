import requests

def get_trs_for_hull(model,species_list):
    result = requests.post("https://query.openkim.org/api",data={'query': '{"meta.type":"tr","meta.subject.extended-id":"%s","property-id":"tag:staff@noreply.openkim.org,2023-02-21:property/binding-energy-crystal","stoichiometric-species.source-value":{"$in":%s}}'%(model,str(species_list).replace("'",'"')), 'fields': '{"prototype-label.source-value":1,"stoichiometric-species.source-value":1,"binding-potential-energy-per-formula.source-value":1}', 'database': 'data'}).json()
    elemental_references = {}
    # initialize dict for elemental reference
    for species in species_list:
        elemental_references[species]=0.
    for entry in result:
        species=entry["stoichiometric-species"]["source-value"]
        if len(species)==1:
            energy = entry["binding-potential-energy-per-formula"]["source-value"]
            species = species[0]
            if energy < elemental_references[species]:
                elemental_references[species] = energy
    # format of a 
    print(elemental_references)
            
            