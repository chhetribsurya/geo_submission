#!/usr/bin/env python2
# -*- coding: latin-1 -*-
'''GET the results of a search from an ENCODE server'''

import requests, json, re, os, subprocess
from collections import defaultdict, OrderedDict
import pandas as pd
from os.path import expanduser, join

# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}

output_dir = os.path.expanduser("~/geo_submission/")
GEO_ID_outdir = join(output_dir, "GEO_ID_output")
if not os.path.exists(GEO_ID_outdir):
    os.makedirs(GEO_ID_outdir)

#xls_df = pd.read_excel(expanduser("~/Dropbox/encode_3/final_figures/supplementary/Encode_full_hepg2_datasets.xlsx"), sheetname=None) #All sheets as a dictionary of DataFrames
xls_df = pd.read_excel(expanduser("~/geo_submission/Encode_full_hepg2_datasets.xlsx"), sheetname=None) #All sheets as a dictionary of DataFrames
xls_df.keys() #to see all the sheetnames from dictionary of dataframes

### choose the columns containing (samplename, rep1_acc, rep1_control_acc, rep2_control_acc) and the set the index using sample name:
xls_df_1 = xls_df["myers_lab_acc"].iloc[:,[4,1,3,7,9,10]]
xls_df_1.columns = ["sample_name", "rep1", "control1", "rep2", "control2", "lab"]
xls_df_2 = xls_df["other_labs_acc"].iloc[:,[4,0,1,5,6,10]]
xls_df_2.columns = ["sample_name", "rep1", "control1", "rep2", "control2", "lab"]
sample_fastq_acc_xls_df = pd.concat([xls_df_1, xls_df_2], ignore_index=True) 
sample_fastq_acc_df = sample_fastq_acc_xls_df.set_index(["sample_name", "lab"])
#sample_fastq_acc_df = sample_fastq_acc_df.sample(30)

sample_fastq_acc_dict = sample_fastq_acc_df.to_dict(orient="split")
sample_fastq_acc_ziplist = zip(sample_fastq_acc_dict["index"], sample_fastq_acc_dict["data"])

GEO_ID_dict = defaultdict(list)
for sample_name,fastq_acc_list in sample_fastq_acc_ziplist:
    for fastq_acc in fastq_acc_list:
        if not fastq_acc.startswith("ENC"):
            GEO_ID = fastq_acc
            GEO_ID_dict[sample_name].append(GEO_ID)
            print "ENC### not available, so Upload the fastq for", sample_name 
            continue

        else:
            # This URL locates the ENCODE biosample with accession number ENCBS000AAA...(eg)
            URL = "https://www.encodeproject.org/biosample/" + fastq_acc + "/?frame=embedded" # frame=object; for short, brief and robust description
            response = requests.get(URL, headers=HEADERS) # GET the object

            # Extract the JSON response as a python dict
            response_json_dict = response.json() # Print the object; print json.dumps(response_json_dict, indent=4, separators=(',', ': '))     
            try:
                expt_acc = response_json_dict["dataset"].split("/")[-2] # extracts the experiment accession from the fastq accession:       
            except KeyError:
                print "Data not availabe to public for %s \n\n\n" %(fastq_acc)
                GEO_ID = "NA_public"
                GEO_ID_dict[sample_name].append(GEO_ID)
                continue

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
                        GEO_ID_dict[sample_name].append(GEO_ID)
                        #GEO_ID_list.append(GEO_ID)
                        print GEO_ID, ": found at position", i #GEO:GSM803530 : found at position 0             

            else:
                print "GEO_ID not available, so using curl -L for deep search"
                json_out = subprocess.check_output(["curl", "-L", "-H", "Accept: application/json", join("https://www.encodeproject.org/biosamples", fastq_acc + "/" )])
                #str_args = [ str(x) for x in args ]
                json_response = json.loads(json_out)

                extdb_accession_list = json_response["replicate"]["experiment"]["dbxrefs"]  # use; json_response.items() to see the data in tuple format:   
                mystring = "___".join(extdb_accession_list)
                regexes = [
                        re.compile("GEO*"), 
                        re.compile("..."), 
                        re.compile("..."),
                        ]

                if any(regex.match(mystring) for regex in regexes):
                    for i,each in enumerate(extdb_accession_list):
                        if each.startswith("GEO"):
                            GEO_ID = "NA_movedTo " + each
                            GEO_ID_dict[fastq_acc].append(GEO_ID)
                            #GEO_ID_list.append(GEO_ID)
                            print GEO_ID, ": found at position after deep search", i #GEO:GSM803530 : found at position 0

                else: 
                    print "GEO_ID not available, so assigning NA_nywhere"
                    GEO_ID = "NA_nywhere"
                    GEO_ID_dict[sample_name].append(GEO_ID)

        print "Here's the GEO_ID for %s %s fastq(ENC### Acc)  - %s\n\n" %(sample_name, fastq_acc, GEO_ID)


GEO_ID_df = pd.DataFrame.from_dict(GEO_ID_dict, orient="index") # orient="columns"; if keys to be columns 
GEO_ID_df.columns = ["Rep1", "Control1", "Rep2", "Control2"]

GEO_ID_df_reset = GEO_ID_df.reset_index()
GEO_ID_df_reset[["TF_name", "Lab"]] = GEO_ID_df_reset["index"].apply(pd.Series)
GEO_ID_df_final = GEO_ID_df_reset.loc[:, ["TF_name", "Rep1", "Control1", "Rep2", "Control2", "Lab"]]

GEO_ID_df_final.to_excel(join(GEO_ID_outdir, "GEO_ID_final_details_v2.xlsx"))
### Or, print this way to see the clean list of tf_name and GEO_ID:
for key,values in GEO_ID_dict.iteritems():
    print key, ", ".join(values)

#GEO_ID_df.reset_index()[GEO_ID_df.reset_index()["index"].str.contains("human")].shape

# import json
# from pprint import pprint

# with open('json_test.json') as data_file:    
#     data = json.load(data_file)

# pprint(data)

