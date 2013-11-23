#!/usr/local/bin/python python

"""
Structure and Dynamics of Complex Networks Final Project

Lily Pratt & Erik Islo
25 November 2013

Begins by querying Entrez/Medline databases through the Bio python module.
From this data, a JSON file is created. Subsequent runs of the program
append to the existing JSON file, excluding duplicates.
"""

import sys, json
from Bio import Entrez, Medline

def gatherData(search_term, email="lwrpratt@gmail.com"):
  # Tell NCBI who we are (they track daily usage)
  Entrez.email = email

  # What are we looking for? Def'd by terminal argument.
  search_term = sys.argv[1]

  # Search for the term, create list of matching paper's ID's
  handle = Entrez.esearch(db="pubmed", term=searchterm)
  record_for_id = Entrez.read(handle)
  idlist = record_for_id["IdList"]
  handle_medline = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
  records = Medline.parse(handle_medline)
  records = list(records)

  data = {}
  with open("{}.json".format(searchterm.replace(' ', '_')), "w") as f:
    for record in records:
      data[record.get("PMID", "?")] = {"title" : record.get("MH", "?")}, "author(s)" : record.get("FAU", "?"), "keywords" : record.get("MH", "?")}
      json.dump(data, f, sort_keys=True, indent=4, separators=(',', ': '))
    f.closed


def readData(fname="data.json"):
  json_data = open(fname, "r")
  read_data = json.load(json_data)
  json_data.close()
  print "json file read in successfully"

  return read_data


def appendToJSON(fname="data.json"):
  searchterm =  sys.argv[1]

  readData = readData(fname)

  json_termfile = open(searchterm, "r")
  read_term_data = json.load(json_termfile)
  json_termfile.close()
  print "file read in"

  # appends new data together
  read_data.update(read_term_data)

  with open(fname, "w") as f:
    json.dump(read_data, f, sort_keys=True, indent=4, separators=(', ', ': '))
  f.close

# Script starts here
print "arguments: 1-search term, 2-filename," # 3-email"
search_term = sys.argv[1]
fname = sys.argv[2]
# email = sys.argv[3] not implemented yet

print "gathering data"
getData(search_term)
print "writing to file"
appendToJSON(fname)
