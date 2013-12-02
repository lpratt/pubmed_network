# Queries the Entrez/Medline databases and creates a JSON file with the results

import sys
import json
from Bio import Entrez
from Bio import Medline

# Tell NCBI who I am
Entrez.email = "lwrpratt@gmail.com"

# Search term taken from command line input
searchterm = sys.argv[1]

# Searching the db for a term, creating a list of paper IDs
handle = Entrez.esearch(db="pubmed", term=searchterm)
record_for_id = Entrez.read(handle)
idlist = record_for_id["IdList"]
handle_medline = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
records = Medline.parse(handle_medline)
records = list(records)

# Creating output JSON file, currently overwrites file, need to append to file
data = {}
with open("{}.json".format(searchterm.replace(' ','_')), "w") as f:
    for record in records:
        data[record.get("PMID", "?")] = {"title": record.get("TI", "?"), "author(s)": record.get("FAU", "?"), "keywords": record.get("MH", "?")}
    json.dump(data, f, sort_keys=True, indent=4, separators=(',', ': '))
f.closed