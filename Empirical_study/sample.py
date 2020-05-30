from collections import OrderedDict


#### Primate Scandentia Glires ####
psgInfo = OrderedDict()
psgInfo["Primate"] = [# (Common name, Latin name)
    ("Human", "Homo_sapiens"),
    ("Chimpanzee", "Pan_troglodytes"),
    ("Gorilla", "Gorilla_gorilla"),
    ("Rhesus_monkey", "Macaca_mulatta"),
    ("Orangutan", "Pongo_abelii"),
    ("Marmoset", "Callithrix_jacchus"),
    ("Tarsier", "Carlito_syrichta"),
    ("Mouse_lemur","Microcebus_murinus"),
]
psgInfo["Scandentia"] = [
    ("Northern_treeshrew", "Tupaia_belangeri"),
    ("Chinese_treeshrew", "Tupaia_chinensis"),
]
psgInfo["Glires"] = [
    ("Mouse", "Mus_musculus"),
    ("Rat", "Rattus_norvegicus"),
    ("Squirrel", "Ictidomys_tridecemlineatus"),
    ("Guinea_pig", "Cavia_porcellus"),
    ("Rabbit", "Oryctolagus_cuniculus"),
    ("Pika", "Ochotona_princeps"),
]    
psgInfo["Outgroup"] = [
    ("Dog", "Canis_familiaris"),
    ("Pig", "Sus_scrofa"),
    ("Little_brown_bat", "Myotis_lucifugus"),
]

