#!/usr/bin/env python2
# -*- coding: latin-1 -*-
'''GET the results of a search from an ENCODE server'''

import requests, json, re, os, subprocess, pickle
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
xls_df_1 = xls_df["myers_lab_acc"].iloc[:,[4,0,1,2,3,6,7,8,9,10]]
xls_df_1.columns = ["sample_name", "sl_rep1", "rep1", "sl_control1", "control1", "sl_rep2", "rep2", "sl_control2", "control2", "lab"]
xls_df_1["SL_rep1_comb"] = xls_df_1["sl_rep1"].astype(str) + "," + xls_df_1["rep1"].astype(str)
xls_df_1["SL_control1_comb"] = xls_df_1["sl_control1"].astype(str) + "," + xls_df_1["control1"].astype(str)
xls_df_1["SL_rep2_comb"] = xls_df_1["sl_rep2"].astype(str) + "," + xls_df_1["rep2"].astype(str)
xls_df_1["SL_control2_comb"] = xls_df_1["sl_control2"].astype(str) + "," + xls_df_1["control2"].astype(str)

xls_df_2 = xls_df["other_labs_acc"].iloc[:,[4,2,0,3,1,7,5,8,6,10]]
xls_df_2.columns = ["sample_name", "sl_rep1", "rep1", "sl_control1", "control1", "sl_rep2", "rep2", "sl_control2", "control2", "lab"]
xls_df_2["SL_rep1_comb"] = xls_df_2["sl_rep1"].astype(str) + "," + xls_df_2["rep1"].astype(str)
xls_df_2["SL_control1_comb"] = xls_df_2["sl_control1"].astype(str) + "," + xls_df_2["control1"].astype(str)
xls_df_2["SL_rep2_comb"] = xls_df_2["sl_rep2"].astype(str) + "," + xls_df_2["rep2"].astype(str)
xls_df_2["SL_control2_comb"] = xls_df_2["sl_control2"].astype(str) + "," + xls_df_2["control2"].astype(str)

sample_fastq_acc_xls_df = pd.concat([xls_df_1, xls_df_2], ignore_index=True) 
select_col = ["sample_name", "SL_rep1_comb", "SL_control1_comb", "SL_rep2_comb", "SL_control2_comb", "lab"]
sample_fastq_acc_xls_df = sample_fastq_acc_xls_df.loc[:, select_col]
sample_fastq_acc_df = sample_fastq_acc_xls_df.set_index(["sample_name", "lab"])
#sample_fastq_acc_df = sample_fastq_acc_df.sample(15)

sample_fastq_acc_dict = sample_fastq_acc_df.to_dict(orient="split")
sample_fastq_acc_ziplist = zip(sample_fastq_acc_dict["index"], sample_fastq_acc_dict["data"])

SL_num_list = []
GEO_ID_dict = defaultdict(list)
for sample_name,fastq_acc_list in sample_fastq_acc_ziplist:
    for fastq_acc in fastq_acc_list:
        SL_num = fastq_acc.split(",")[0]
        fastq_acc_new = fastq_acc.split(",")[1]
        if not fastq_acc_new.startswith("ENC"):
            GEO_ID = SL_num+".fastq.gz"
            GEO_ID_dict[sample_name].append(GEO_ID)
            SL_num_list.append(SL_num)
            print "ENC### not available, so Upload the fastq for", sample_name 
            continue

        else:
            # This URL locates the ENCODE biosample with accession number ENCBS000AAA...(eg)
            URL = "https://www.encodeproject.org/biosample/" + fastq_acc_new + "/?frame=embedded" # frame=object; for short, brief and robust description
            response = requests.get(URL, headers=HEADERS) # GET the object

            # Extract the JSON response as a python dict
            response_json_dict = response.json() # Print the object; print json.dumps(response_json_dict, indent=4, separators=(',', ': '))     
            try:
                expt_acc = response_json_dict["dataset"].split("/")[-2] # extracts the experiment accession from the fastq accession:       
            except KeyError:
                print "Data not availabe to public for %s \n\n\n" %(fastq_acc_new)
                GEO_ID = fastq_acc_new + "_in_Progress"
                GEO_ID_dict[sample_name].append(GEO_ID)
                #SL_num_list.append(SL_num) # don't upload these fastq
                continue

            expURL = "https://www.encodeproject.org/biosample/" + expt_acc + "/?frame=object"
            exp_response = requests.get(expURL, headers=HEADERS)
            exp_response_json_dict = response.json()

            extdb_accession_list = exp_response_json_dict["replicate"]["experiment"]["dbxrefs"]  # use; exp_response_json_dict.items() to see the data in tuple format: 
            # platform_list = response_json_dict_expt["files"][0]["platform"]["aliases"]
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
                json_out = subprocess.check_output(["curl", "-L", "-H", "Accept: application/json", join("https://www.encodeproject.org/biosamples", fastq_acc_new + "/" )])
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
                            GEO_ID_dict[fastq_acc_new].append(GEO_ID)
                            #GEO_ID_list.append(GEO_ID)
                            print GEO_ID, ": found at position after deep search", i #GEO:GSM803530 : found at position 0

                else: 
                    print "GEO_ID not available, so assigning NA_nywhere"
                    #GEO_ID = SL_num + ".fastq.gz"
                    GEO_ID = fastq_acc_new
                    GEO_ID_dict[sample_name].append(GEO_ID)
                    #SL_num_list.append(SL_num)

        print "Here's the GEO_ID for %s %s fastq(ENC### Acc)  - %s\n\n" %(sample_name, fastq_acc_new, GEO_ID)


GEO_ID_df = pd.DataFrame.from_dict(GEO_ID_dict, orient="index") # orient="columns"; if keys to be columns 
GEO_ID_df.columns = ["Rep1", "Control1", "Rep2", "Control2"]

GEO_ID_df_reset = GEO_ID_df.reset_index()
GEO_ID_df_reset[["TF_name", "Lab"]] = GEO_ID_df_reset["index"].apply(pd.Series)
GEO_ID_df_final = GEO_ID_df_reset.loc[:, ["TF_name", "Rep1", "Control1", "Rep2", "Control2", "Lab"]]
GEO_ID_df_final.to_excel(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v4.xlsx"))

SL_num_df = pd.Series(SL_num_list)
SL_num_df.to_csv(join(GEO_ID_outdir, "SL_num_list_v4.txt"), header=False, index=False, sep="\t")

### Or, save as pickle, load(when needed), and print this way to see the clean list of tf_name and GEO_ID:
GEO_ID_df_final.to_pickle(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v4.pkl"))
GEO_ID_df_final = pd.read_pickle(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v4.pkl"))
pickle.dump(GEO_ID_dict, open(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v4_dict.pkl"), "wb"))  # save it into a file named open('save.p', "wb")
GEO_ID_dict = pickle.load(open(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v4_dict.pkl"), "rb"))

for key,values in GEO_ID_dict.iteritems():
    print key, ", ".join(values)

for key,values in GEO_ID_dict.iteritems():
    print key[0], "\t\t\t\t\t\t\t", key[1]


# import json
# from pprint import pprint

# with open('json_test.json') as data_file:    
#     data = json.load(data_file)
# pprint(data)


df = pd.read_excel("/Users/suryachhetri/Dropbox/encode_3/geo_submission/GEO_ID_output/GEO_sample_details_w_fastq_v4.xlsx")
df["Rep1"] = df["Rep1"].str.replace(r"GEO:\w+|ENC\w+", "np.nan")
df["Rep2"] = df["Rep2"].str.replace(r"GEO:\w+|ENC\w+", "np.nan")
df["Control1"] = df["Control1"].str.replace(r"GEO:\w+|ENC\w+", "np.nan")
df["Control2"] = df["Control2"].str.replace(r"GEO:\w+|ENC\w+", "np.nan")

df["combined_col"] = df["Rep1"].astype(str) + "_" + df["Control1"].astype(str) + "_" + df["Rep2"].astype(str) + "_" + df["Control2"].astype(str)
#df["combined_col"].str.replace("_np.nan", "").str.split("_").apply(pd.Series)
df[["upload1", "upload2", "upload3", "upload4"]] = df["combined_col"].str.replace("_np.nan", "").str.split("_").apply(pd.Series)
df.to_excel("/Users/suryachhetri/Dropbox/encode_3/geo_submission/GEO_ID_output/upload_data/GEO_ID_details_final_processed_for_uploads.xlsx")

ls /gpfs/gpfs1/home/mmackiewicz/TAG-CHIP/*/*merged*


#!/bin/bash

### File to be appended on:
append_file="/gpfs/gpfs1/home/schhetri/geo_submission/GEO_ID_output/SL_num_list_v4_path_for_gzip.txt"

### This is important since some of the dir could be just linked and find func won't find the linked dir:
#for each in $(ls $(pwd)/*batch_libraries); do
input_dir="/gpfs/gpfs1/home/schhetri/for_encode/spp/Libraries"
find_dir=$(readlink -f $input_dir)
for each in $(grep -v "SL161" "SL_num_list_v4.txt"|sort |uniq); do find $find_dir -name $each -type d; done >> ${append_file}

input_dir="/gpfs/gpfs1/home/schhetri/for_encode/spp/II_batch_Libraries/links"
find_dir=$(readlink -f $input_dir)
for each in $(grep -v "SL161" "SL_num_list_v4.txt" |sort |uniq); do find $find_dir -name $each -type d; done >> $append_file

input_dir="/gpfs/gpfs1/home/schhetri/for_encode/spp/VI_batch_libraries"
find_dir=$(readlink -f $input_dir)
for each in $(grep -v "SL161" "SL_num_list_v4.txt" |sort |uniq); do find $find_dir -name $each -type d; done >> ${append_file}

input_dir="/gpfs/gpfs1/home/schhetri/for_encode/spp/VII_batch_libraries"
find_dir=$(readlink -f $input_dir)
for each in $(grep -v "SL161" "SL_num_list_v4.txt" |sort |uniq); do find $find_dir -name $each -type d; done >> ${append_file}

### Submit the job for all the collection of files:
for each_SL in $(cat ${append_file}); do file_name=$(basename ${each_SL}); bsub -J "geo submission" -o "geo_submit_gzip.out" "zcat ${each_SL}/* | gzip > $output_dir/${file_name}.fastq.gz"; done

### Submit job for Mark dir containing SL161 series:
output_dir="/gpfs/gpfs1/home/schhetri/geo_submission/surya_chhetri/surya_chhetri_ChIPSeq"
input_dir="/gpfs/gpfs1/home/mmackiewicz/TAG-CHIP/*/*"
find_dir=$(readlink -f $input_dir)
for each in $(grep "SL161" SL_num_list_v4.txt); do file_loc=$(find $find_dir -name ${each}*merged.fastq.gz -type f); bsub -J "geo submission" -o "geo_submit.out" "cp $file_loc $output_dir/${each}.fastq.gz"; done 



import GEOparse
gse = GEOparse.get_GEO(geo="GSE96509", destdir="./")
sorted(gse.gsms.keys())

#!zgrep "GSM" ./GSE96509_family.soft.gz | head
history

# Data structure:
test_dict = {"tf1":{"chip" : ["ENCSS","ENCSS","ENCS2"], 
                    "control" : ["ENCFF", "ENCFF"]
                    }, 
            "tf2" : {"chip" : ["GSM1", "GSM2"], 
                    "control" : ["GSM", "GSM", "GSM3"]
                    }
            }








