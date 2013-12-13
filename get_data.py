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

# TODO - cytoscape input file so that only one edge between nodes appears

def writeToJSON(obj_to_write, file_name):
  # Put Paper Information in a JSON file
  with open(file_name, 'w') as f:
     json.dump(obj_to_write, f, sort_keys=True, indent=4, separators=(',', ': '))
     f.closed

def makeAdjList(authors, file_name="pubmed_authors.txt"):
  added_authors = {}
  with open(file_name, 'w') as f:
    for a in authors:
      au = a.replace(' ','_')
      for co in authors[a]:
        coau = co.replace(' ','_')
        if (au,coau) not in added_authors and (coau,au) not in added_authors: # add to file and to added_authors
          added_authors[(au,coau)] = authors[a][co]
          added_authors[(coau,au)] = authors[a][co]
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

def authorKeyword(target_file="author_keyword.json", f="papers.json"):
#  with open(target_file, 'r+') as f:
#    authkey = json.load(target_file)

  with open(f, 'r') as pap:
    papers = json.load(pap)
    pap.closed

  try:
    with open(target_file, 'r') as target:
      authkey = json.load(target)
      target.closed
  except IOError:
    authkey = {}

  for id in papers:
    for author in papers[id]['Authors']:
      id_set = [id]
      keyword_set = []
      for key in papers[id]['Keywords']:
        keyword_set.append(key)
      if author in authkey:
        authkey[author]['Keywords'].append(keyword_set)
        authkey[author]['Paper ID'].append(id)
      else:
        authkey.update({author: {'Keywords': keyword_set, 'Paper ID': id_set}})

  writeToJSON(authkey, 'author_keyword.json')
#    json.dump()
#    f.closed
  return authkey


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
      # add search terms
      if search_term not in papers[idnum]["Search Terms"]:
        papers[idnum]["Search Terms"].append(search_term)
    else:
      search_term_set = [search_term]
      papers.update({idnum:{"Search Terms": search_term_set, "Title" : record.get('TI', '?'), "Authors" : record.get('AU','?'),"Keywords":record.get('MH','?')}})

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

# TODO: breaks on certain words, also does not deal with empty keyword_queue
def recurseKeywords(papers, keyword_queue):
  for paper in papers:
#    print paper
    for keyword in papers[paper]['Keywords']:
      # remove *, separate / ,
      keyword = keyword.replace('*','')
      keyword_list = re.split('[/|,|&]+', keyword)
#      print keyword_list
      for add_keyword in keyword_list:
        if add_keyword != '?' and add_keyword not in keyword_queue:
          keyword_queue[add_keyword] = 0
#  print keyword_queue
  # part of function that recurses
  for search_term in keyword_queue:
    if keyword_queue[search_term] == 0 and search_term != 'Insulin-Secreting Cells' and search_term != 'Oceans and Seas' and search_term != "Time Factors" and search_term != "HIV Infections" and search_term != "Lymphopoiesis":
      print search_term
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
  

# TODO : function that takes in a community and labels the keywords associated with it

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
