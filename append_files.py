# Updates data.json to include the file given
# Would it be useful for it to create a file with names of already appended files?

import sys
import json

# Search term taken from command line input
searchterm = sys.argv[1]

# Reads in json files
json_data = open("data.json", "r")
read_data = json.load(json_data)
json_data.close()

json_termfile = open(searchterm, "r")
read_term_data = json.load(json_termfile)
json_termfile.close()
#print read_term_data

# Appends new data together
read_data.update(read_term_data)
#print read_data["24194598"]

# Puts appended data back into data.json file
with open("data.json", "w") as f:
    json.dump(read_data, f, sort_keys=True, indent=4, separators=(',', ': '))
f.close