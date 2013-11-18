import sys
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

# Just as a check
for record in records:
     print "title:", record.get("TI", "?")
     print "authors:", record.get("FAU", "?")
     print "keywords:", record.get("MH", "?")
     print