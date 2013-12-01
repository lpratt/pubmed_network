#!/usr/local/bin/python python

"""
Structure and Dynamics of Complex Networks Final Project

Lily Pratt & Erik Islo
25 November 2013

Begins by querying Entrez/Medline databases through the Bio python module.
From this data, a JSON file is created. Subsequent runs of the program
append to the existing JSON file, excluding duplicates.

Authors are stored as a dictionary with keys as Author names and values
as a Counter object with all linked authors and their respective connection
weights.
"""

import sys, json
from collections import Counter
from Bio import Entrez, Medline

authors = {}

def authorWeight(authors)
  """
  Defines an edge weight for coauthors of a given paper.
  Determined so that the more authors on a paper the weaker their weight.
  """

  # returns a dictionary of linked authors and respective edge weights
  return 


def addAuthor(author, coauthors)
  if author in authors:
    authors.update({author: authorWeight(coauthors)})

def updateAuthor(author, coauthors)

def writeToJSON()


def gatherData(search_term, email="lwrpratt@gmail.com"):

  Entrez.email = email       # NCBI identification
  search_term = sys.argv[1]  # what are we looking for?

  # See how many papers match this query
  handle = Entrez.egquery(term=search_term)
  record = Entrez.read(handle)
  for row in record["eGQueryResult"]:
    if row["DbName"]=="pubmed":
      num_papers = row["Count"]
  print num_papers + " papers found with egquery."

  handle = Entrez.esearch(db="pubmed", term=search_term)
  record = Entrez.read(handle)
  idlist = record["IdList"]
  print "Search for '%s' resulted in %s results.\n" % (search_term, str(len(idlist)))

  handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode ="text")
  records = Medline.parse(handle)
  print "efetch complete"


  for record in records:
    au = record.get('AU', '?')
    for a in au:
      if au in authors:
        authors.update{}
      else:

      

"""
  # Put Paper Information in a JSON file
  data = {}
  with open("data.json", "wb") as f:
    for record in records:
      data[record.get("PMID", "?")] = {"title": record.get("TI", "?"), "author(s)": record.get("FAU", "?"), "keywords": record.get("MH", "?")}
      json.dump(data, f, sort_keys=True, indent=4, separators=(',', ': '))
    f.closed
"""

#-------------------------------------------------------------------------#
#                          Script starts here                             #
#-------------------------------------------------------------------------#
print "ARGUMENTS: 1-search term, 2-filename\n"
search_term = sys.argv[1]
# fname = sys.argv[2]

print "gathering data..."
gatherData(search_term)
