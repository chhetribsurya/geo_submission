#!/usr/bin/env python2
# -*- coding: latin-1 -*-
'''GET the results of a search from an ENCODE server'''

import requests, json, re
from collections import defaultdict, OrderedDict


# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}

# This URL locates the ENCODE biosample with accession number ENCBS000AAA
fastq_acc ="ENCFF000PPW"
fastq_acc = "ENCFF000POC"
fastq_acc = "ENCFF000PQC"
fastq_acc = "ENCFF698KPI"
fastq_acc = "ENCFF088SMS"
fastq_acc = "ENCFF000POQ"
fastq_acc = "ENCFF000PMN"
fastq_acc = "ENCFF000PMN"
fastq_acc = "ENCFF002EEF"
fastq_acc_list = ["ENCFF854PST"]
fastq_acc_list = ["ENCFF000POQ", "ENCFF000POQ", "ENCFF002EEF", "ENCFF002EEF", "ENCFF000PMN", "ENCFF000PMN"]



GEO_ID_dict = defaultdict(list)
for fastq_acc in fastq_acc_list:
	URL = "https://www.encodeproject.org/biosample/" + fastq_acc + "/?frame=embedded" # frame=object; for short, brief and robust description
	response = requests.get(URL, headers=HEADERS) # GET the object

	# Extract the JSON response as a python dict
	response_json_dict = response.json() # Print the object; print json.dumps(response_json_dict, indent=4, separators=(',', ': '))
	expt_acc = response_json_dict["dataset"].split("/")[-2] # extracts the experiment accession from the fastq accession:

	expURL = "https://www.encodeproject.org/biosample/" + expt_acc + "/?frame=object"
	exp_response = requests.get(expURL, headers=HEADERS)
	exp_response_json_dict = response.json()

	extdb_accession_list = exp_response_json_dict["replicate"]["experiment"]["dbxrefs"]  # use; exp_response_json_dict.items() to see the data in tuple format:	
	mystring = "___".join(extdb_accession_list)
	regexes = [
			re.compile("GEO*"), 
			re.compile("..."), 
			re.compile("..."),
			]

	if any(regex.match(mystring) for regex in regexes):
		for i,each in enumerate(extdb_accession_list):
		    if each.startswith("GEO"):
		    	GEO_ID = each
		    	GEO_ID_dict[fastq_acc].append(GEO_ID)
		    	#GEO_ID_list.append(GEO_ID)
		        print GEO_ID, ": found at position", i #GEO:GSM803530 : found at position 0				
	else:
		print "GEO_ID not available, so using curl -L"
		cmd = """curl -L -H "Accept: application/json" https://www.encodeproject.org/biosamples/"""
		cmd_call = join(cmd, fastq_acc + "/")
		ouput_file = fastq_acc + ".json"
		#os.system(cmd > "tmp.json" )
		os.system('cmd_call')
		with open(ouput_file) as data_file:
			data = json.load(data_file)

		#"assigning N/A"
		#GEO_ID = "NA"
		GEO_ID_dict[fastq_acc].append(GEO_ID)

	print "Here's the GEO_ID for %s fastq file - %s\n\n" %(fastq_acc, GEO_ID)

geo_id_df = pd.DataFrame.from_dict(GEO_ID_dict, orient="index") # orient="columns"; if keys to be columns 

# import re

# regexes = [
#     # your regexes here
#     re.compile('hi'),
# #    re.compile(...),
# #    re.compile(...),
# #    re.compile(...),
# ]

# mystring = 'hi'

# if any(regex.match(mystring) for regex in regexes):
#     print 'Some regex matched!'


