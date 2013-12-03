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

def authorWeight(coauthors):
  """
  Defines an edge weight for coauthors of a given paper.
  Determined so that the more authors on a paper the weaker their weight.
  """
  coauth = {}
  num_authors = len(coauthors) + 1        # doesn't include current author being updated
  weight = 1.0/(float(num_authors) - 1.0) # each author on paper given same link weight

  for a in coauthors:
    coauth.update({a: weight})

  # returns a dictionary of linked authors and respective edge weights
  return coauth


def addAuthor(authors, author, coauthors):
  if coauthors == []:
    authors.update({author : {}})
  else:
    authors.update({author : authorWeight(coauthors)})

  return authors


def updateAuthor(authors, author, coauthors):
  if coauthors == []:
    pass
  else:
    authors[author] = dict(Counter(authors[author]) + Counter(authorWeight(coauthors)))

  return authors 

#-------------------------------------------------------------------------#
#                END author network updating functions                    #
#-------------------------------------------------------------------------#

def writeToJSON(authors, file_name="authors.json"):
  # Put Paper Information in a JSON file
  with open(file_name, 'w') as f:
     json.dump(authors, f, sort_keys=True, indent=4, separators=(',', ': '))
     f.closed

def makeAdjList(authors, file_name="pubmed_authors.txt"):
  with open(file_name, 'w') as f:
    for a in authors:
      au = a.replace(' ','_')
      for co in authors[a]:
        coau = co.replace(' ','_')
        f.write(au + ' ' + coau + ' ' + str(authors[a][co])+'\n')
  f.close()

# TODO check about adding weights for multiple occurence of same paper
# TODO handle cycles from authors connecting to each other
  # though perhaps this is handled by cytoscape/networkx
def gatherData(search_term, email="lwrpratt@gmail.com"):
  """
  Queries pubmed for papers with the given search_term using the BioPython
  toolkit's Entrez and Medline modules.

  Saves each paper's authors in a dict with all connected authors and
  their repsective connection strength
  """

  try:
    with open('authors.json', 'r') as auths:
      authors = json.load(auths)
      print authors
      print "successfully imported"
    auths.closed
  except IOError:
    authors = {}
    print "no old data to import"

  Entrez.email = email       # NCBI identification

  # See how many papers match this query
  handle = Entrez.esearch(db='pubmed', term=search_term)
  result = Entrez.read(handle)
  num_papers = result['Count']
  print num_papers + " papers found with esearch."

  #TODO Why does this give a huge number when ids (below) is so short?

  # Get the info on these papers
  ids = result['IdList']
  handle = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
  records = Medline.parse(handle)
  print "efetch complete"


  for record in records:
    au = record.get('AU', '?')
    au.sort()
    for a in au:
      #sys.stdout.write("Author: %s, " % str(a)) # debug
      if a in authors:
        #print "existing"                        # debug
        au.remove(a)
        authors = updateAuthor(authors, a, au)
        au.append(a) # TODO this is gross and should not stay this way
        au.sort()
      else:
        #print "new"                             # debug
        au.remove(a)
        authors = addAuthor(authors, a, au)
        au.append(a)
        au.sort()

  return authors

      

#-------------------------------------------------------------------------#
#                          Script starts here                             #
#-------------------------------------------------------------------------#
if __name__ == '__main__':

  print "ARGUMENTS: 1-search term, 2-filename\n"
  search_term = sys.argv[1]
  # fname = sys.argv[2]
  print "gathering data..."
  authors = gatherData(search_term)

  writeToJSON(authors)

    authors[author] = dict(Counter(authors[author]) + Counter(authorWeight(coauthors)))

  makeAdjList(authors)
