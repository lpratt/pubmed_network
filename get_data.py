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

import sys, json, re
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

def writeToJSON(obj_to_write, file_name):
  # Put Paper Information in a JSON file
  with open(file_name, 'w') as f:
     json.dump(obj_to_write, f, sort_keys=True, indent=4, separators=(',', ': '))
     f.closed

def makeAdjList(authors, file_name="pubmed_authors.txt"):
  with open(file_name, 'w') as f:
    for a in authors:
      au = a.replace(' ','_')
      for co in authors[a]:
        coau = co.replace(' ','_')
        f.write(au + ' ' + coau + ' ' + str(authors[a][co]) + '\n')
  f.close()

def writeForInfomap(authors, file_name="pubmed_authors_infomap.txt"):
  author_ids = {}
  id_count = 1
  with open(file_name, 'w') as f:
    for a in authors:
      if a not in author_ids:
        author_ids[a] = id_count
        id_count += 1
      au = author_ids[a]
      for co in authors[a]:
        if co not in author_ids:
          author_ids[co] = id_count
          id_count += 1
        coau = author_ids[co]
        f.write(str(au)+' '+str(coau)+' '+str(authors[a][co])+'\n')
  f.close()

# TODO check about adding weights for multiple occurence of same paper
# TODO handle cycles between authors?  Might be cytoscape/nx parsing

def gatherData(search_term, email="lwrpratt@gmail.com"):
  """
  Queries pubmed for papers with the given search_term using the BioPython
  toolkit's Entrez and Medline modules.

  Saves each paper's authors in a dict with all connected authors and
  their repsective connection strength
  """

  try:               # import existing author data
    with open('authors.json', 'r') as auths:
      authors = json.load(auths)
      auths.closed
  except IOError:    # none to import
    print "no old data to import, creating new authors dict"
    authors = {}

  try:               # import existing paper ids
    with open('papers.json', 'r') as paps:
      papers = json.load(paps)
      paps.closed
  except IOError:    # none to import
    print "no old data to import, creating new papers dict"
    papers = {}

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
  print "efetch complete. %s results." % str(len(ids))


  for record in records:
    idnum = record.get('PMID', '?') 

    if idnum in papers:
      print "Already have %s" % str(idnum)
      pass # Already have this paper

    else:
      papers.update({idnum:{"Title" : record.get('TI', '?'), "Authors" : record.get('AU','?'),"Keywords":record.get('MH','?')}})

      au = record.get('AU', '?')
      au.sort()
      for a in au:
        if a in authors:
          au.remove(a)
          authors = updateAuthor(authors, a, au)
          au.append(a) # TODO this is gross and should not stay this way
          au.sort()
        else:
          au.remove(a)
          authors = addAuthor(authors, a, au)
          au.append(a)
          au.sort()

  return (authors, papers)

def recurseKeywords(papers, keyword_queue):
  for paper in papers:
#    print paper
    for keyword in papers[paper]['Keywords']:
      # remove *, separate / ,
      keyword = keyword.replace('*','')
      keyword_list = re.split('[/|,|&]+', keyword)
#      print keyword_list
      for add_keyword in keyword_list:
        if add_keyword not in keyword_queue and add_keyword != '?':
          keyword_queue[add_keyword] = 0
#  print keyword_queue
  # part of function that recurses
  for search_term in keyword_queue:
    print search_term
    if keyword_queue[search_term] == 0 and (search_term != 'Insulin-Secreting Cells' or search_term != 'Ocean and Seas':
      (Ignore, newPapers) = gatherData(search_term)
      # save data structures for future runs of gatherData()
      print "writing JSON files..."
      writeToJSON(Ignore, 'authors.json')
      writeToJSON(newPapers, 'papers.json')
      # creating maps
      makeAdjList(Ignore)
      writeForInfomap(Ignore)
      keyword_queue[search_term] = 1
      recurseKeywords(newPapers, keyword_queue)
  

#-------------------------------------------------------------------------#
#                          Script starts here                             #
#-------------------------------------------------------------------------#
if __name__ == '__main__':

  print "ARGUMENTS: 1-search term, 2-filename\n"
  search_term = sys.argv[1]
  # fname = sys.argv[2]
  print "gathering data..."
  (Authors, Papers)  = gatherData(search_term)

  # testing purposes
  print "keyword queue:"
  keyword_queue = {}
  recurseKeywords(Papers, keyword_queue)

  # save data structures for future runs of gatherData()
  print "writing JSON files..."
  writeToJSON(Authors, 'authors.json')
  writeToJSON(Papers, 'papers.json')

  # create file for networkx or cytoscape to read in as network
  print "creating author adjacency list..."
  makeAdjList(Authors)

  # create file for infomap to read in as network
  print "creating numbered author adjacency list..."
  writeForInfomap(Authors)
