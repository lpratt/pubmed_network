# Reads data.json into a python dictionary

import json

# Reads in json file
json_data = open("data.json", "r")
read_data = json.load(json_data)
json_data.close()
print read_data