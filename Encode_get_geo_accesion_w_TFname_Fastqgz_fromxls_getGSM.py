#!/usr/bin/env python2
# -*- coding: latin-1 -*-
'''GET the results of a search from an ENCODE server'''

import requests, json, GEOparse, re, os, subprocess, pickle
from collections import defaultdict, OrderedDict
import pandas as pd
from os.path import expanduser, join

# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}

output_dir = os.path.expanduser("~/geo_submission/")
GEO_ID_outdir = join(output_dir, "GEO_ID_output")
if not os.path.exists(GEO_ID_outdir):
    os.makedirs(GEO_ID_outdir)

GSE_parse_dir = join(output_dir, "GEO_parsed_dir")
if not os.path.exists(GSE_parse_dir):
    os.makedirs(GSE_parse_dir)

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
#sample_fastq_acc_df = sample_fastq_acc_df.sample(30)

sample_fastq_acc_dict = sample_fastq_acc_df.to_dict(orient="split")
sample_fastq_acc_ziplist = zip(sample_fastq_acc_dict["index"], sample_fastq_acc_dict["data"])

SL_num_list = []
#GEO_ID_all_list = []
GEO_ID_dict = defaultdict(list)
for sample_name,fastq_acc_list in sample_fastq_acc_ziplist:
    for fastq_acc in fastq_acc_list:
        SL_num = fastq_acc.split(",")[0]
        fastq_acc_new = fastq_acc.split(",")[1]
        if not fastq_acc_new.startswith("ENC"):
            #GEO_ID = SL_num+".fastq.gz"
            GEO_ID = [SL_num]
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
                #GEO_ID = [expt_acc + " in_progress"]
                GEO_ID = [fastq_acc_new + "_in_progress"]
                GEO_ID_dict[sample_name].append(GEO_ID)
                #SL_num_list.append(SL_num) # don't upload these fastq
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
                        if each.startswith("GEO:GSE"):
                            GEO_ID = each.split(":")[1] # as format is GEO:GSE1234, function only recognises "GSE...." and not "GEO"
                            gse = GEOparse.get_GEO(geo=GEO_ID, destdir=GSE_parse_dir)
                            # gse.download_SRA() ## download_SRA(email, directory=’series’, filterby=None, **kwargs) #Download SRA files for each GSM in series
                            GSM_ID_list = gse.gsms.keys()
                            GEO_ID_dict[sample_name].append(GSM_ID_list)

                        else: #assuming it to be GEO:GSM(which is what we need)
                            GEO_ID = [each]
                            GEO_ID_dict[sample_name].append(GEO_ID)
                            #GEO_ID_list.append(GEO_ID)
                            print GEO_ID, ": found at position", i #GEO:GSM803530 : found at position 0             

            else:
                print "GEO_ID not available, so using curl -L for deep search"
                json_out = subprocess.check_output(["curl", "-L", "-H", "Accept: application/json", join("https://www.encodeproject.org/biosamples", fastq_acc_new + "/" )])
                #str_args = [ str(x) for x in args ]
                json_response = json.loads(json_out)

                extdb_accession_list = json_response["replicate"]["experiment"]["dbxrefs"]  # use; json_response.items() to see the data in tuple format:   
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
                            GEO_ID = ["NA_movedTo " + each]
                            GEO_ID_dict[sample_name].append(GEO_ID)
                            #GEO_ID_list.append(GEO_ID)
                            print GEO_ID, ": found at position after deep search", i #GEO:GSM803530 : found at position 0

                else: 
                    print "GEO_ID not available, so assigning NA_nywhere"
                    #GEO_ID = SL_num + ".fastq.gz"
                    #GEO_ID = [expt_acc]
                    GEO_ID = [fastq_acc_new]
                    GEO_ID_dict[sample_name].append(GEO_ID)
                    #SL_num_list.append(SL_num)

        print "Here's the GEO_ID for %s %s fastq(ENC### Acc)  - %s\n\n" %(sample_name, fastq_acc_new, GEO_ID)


GEO_ID_df = pd.DataFrame.from_dict(GEO_ID_dict, orient="index") # orient="columns"; if keys to be columns 
GEO_ID_df.columns = ["Rep1", "Control1", "Rep2", "Control2"]

GEO_ID_df_reset = GEO_ID_df.reset_index()
GEO_ID_df_reset[["TF_name", "Lab"]] = GEO_ID_df_reset["index"].apply(pd.Series)
GEO_ID_df_final = GEO_ID_df_reset.loc[:, ["TF_name", "Rep1", "Control1", "Rep2", "Control2", "Lab"]]

### Arrange the TF's with reference to standard TF order submitted to GEO:
ref_df = pd.read_excel(expanduser("~/geo_submission/Reference_order_of_TFs.xlsx"), header=None)
ref_df = ref_df.set_index([0])
GEO_ID_df_final_reset = GEO_ID_df_final.set_index(["TF_name"])
final_xls_df = GEO_ID_df_final_reset.reindex(ref_df.index)
final_xls_df["chip"] = final_xls_df["Rep1"] + final_xls_df["Rep2"]
final_xls_df["control"] = final_xls_df["Control1"] + final_xls_df["Control2"]
final_xls_df.to_excel(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v5_final.xlsx"))

final_xls_df_select = final_xls_df.loc[:, ["chip", "control"]]
final_xls_df_uniq = final_xls_df_select.applymap(lambda x: " ".join(list(set(x))))
final_xls_df_uniq["chip"] = final_xls_df_uniq["chip"].str.replace("GEO:", "")
final_xls_df_uniq["control"] = final_xls_df_uniq["control"].str.replace("GEO:", "")

## For Controls and ChIP antibody column:
final_xls_df_uniq["ChIP_antibody"] = final_xls_df_uniq["TF_name"].map(lambda x: "-".join(x.split("[")[::-1]).replace("FLAG]", "FLAG"))
final_xls_df_uniq.to_excel(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v5_final_chip_control.xlsx"))

#### For the merging of the datasets w.r.t ot :
#pd.Series(df["control"].unique()).str.split()[::1].sort_values() ## order of control sample to be concatenated at the end:
df_new = final_xls_df_uniq.loc[:, ["TF_name", "control"]]
df_dict = df_new.set_index(["control"]).to_dict(orient="split").keys()
zipped_dict = zip(df_dict["index"], df_dict["data"])

dup_elim_dict = defaultdict(list)
for gsm_id, values in zipped_dict:
    dup_elim_dict[gsm_id].append("Input for " + values[0])
gsm_id_df = pd.DataFrame.from_dict(dup_elim_dict, orient="index")
gsm_id_df_sorted = gsm_id_df.sort_index()

### Print the data in bedtools merge style:
for key, values in dup_elim_dict.iteritems():
    print "%s\t%s"%(key, ", ".join(values))

### Or, save as pickle, load(when needed), and print this way to see the clean list of tf_name and GEO_ID:
GEO_ID_df_final_reset.to_pickle(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v5_final.pkl"))
GEO_ID_df_final_reset = pd.read_pickle(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v5_final.pkl"))
pickle.dump(GEO_ID_dict, open(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v5_dict.pkl"), "wb"))  # save it into a file named open('save.p', "wb")
GEO_ID_dict = pickle.load(open(join(GEO_ID_outdir, "GEO_ID_final_details_w_fastqgz_v5_dict.pkl"), "rb"))

for key,values in GEO_ID_dict.iteritems():
    print key, ", ".join(values)

for key,values in GEO_ID_dict.iteritems():
    print key[0], "\t\t\t\t\t\t\t", key[1]


#### FOR hudsonalphs SLs##
SL_dict = defaultdict(list)
SL_num_df = pd.Series(SL_num_list)
for each_SL in SL_num_df.unique():
    URL="http://hts.hudsonalpha.org/api/htwos/library/simple/" + each_SL
    response = requests.get(URL, headers=HEADERS)
    response_json_dict = response.json()
    read_info = response_json_dict["run_condition"]
    SL_dict[each_SL].append(read_info)

SL_read_len_df = pd.DataFrame.from_dict(SL_dict, orient="index")
SL_read_len_df.to_excel(join(GEO_ID_outdir, "SL_read_len_v5_final.xlsx"))
SL_num_df.to_csv(join(GEO_ID_outdir, "SL_num_list_v5.txt"), header=False, index=False, sep="\t")


#### Finding the sequencer or platform name from ENCODE accession   
#response_json_dict_expt["files"][0]["platform"]["aliases"]

gsm_id_list = final_xls_df_select["chip"] + final_xls_df_select["control"]
gsm_list = []
for each_gsm_list in gsm_id_list:
    for each_gsm in each_gsm_list:
        if each_gsm.startswith("GSM"):

            gsm_list.append(each_gsm)

gsm_df = pd.Series(gsm_list)
gsm_df.to_csv(join(GEO_ID_outdir, "GSM_ID_list_v5.txt"), header=False, index=False, sep="\t")


library(GEOquery)
library(dplyr)
library(data.table)
library(xlsx)
## use : names(Meta(gsm)) or names(Meta(gpl)) to see the keys or column names:
gsm_list_df <- fread("~/geo_submission/GEO_ID_output/GSM_ID_list_v5.txt", header=FALSE)
names(gsm_list_df) <- "gsm_num"
gsm_list_df$GPL_ID <- lapply(gsm_list_df$gsm_num, function(x) {gsm_id <- getGEO(x); Meta(gsm_id)$platform})
gsm_list_df$PLATFORM_ID <- lapply(gsm_list_df$GPL_ID, function(x) {platform_id <- getGEO(x); Meta(platform_id)$title})
### Remove the duplicates either by whole dataframe or by duplicated(x) on series:
gsm_gpl_platform_df <- unique(gsm_list_df %>% as.data.frame)
rownames(gsm_gpl_platform_df) <-  1:nrow(gsm_gpl_platform_df)
colnames(gsm_gpl_platform_df) <- c("GSM_ID", "GPL_ID", "INSTRUMENT_MODEL")
write.xlsx(gsm_gpl_platform_df, "/gpfs/gpfs1/home/schhetri/geo_submission/GEO_ID_output/HepG2_ChIPSeq_GSM_GPLplatform_metadata_final.xlsx", sheetName="Sheet1", row.names=FALSE)

### FOR WGBS, individual check
#gsm_list_df_wgbs <- data.frame("gsm_num" : c("GSM2308630", "GSM2308631"))
gsm_id <- getGEO("GSM2308630")
gsm_id <- getGEO("GSM2308631")
gpl_platform <- Meta(gsm_id)$platform
gpl <- getGEO(gpl_platform) #getGEO('GPL11154')
sequencer_name = Meta(gpl)$title ##
# gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform})
# head(gsmplatforms)

#u'title': u'TOWARD A COMPREHENSIVE FUNCTIONAL ANNOTATION OF THE HUMAN GENOME',
#  u'url': u'http://projectreporter.nih.gov/project_info_details.cfm?aid=8402461',




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



### Extract the read length from the SLs:
LIB_LIST='SL81123'
for LIB in ${LIB_LIST}; do /opt/python2.7/bin/python ~snewberry/bin/eat_json.py hts.morgan.haib.org/api/htwos/library/simple/$LIB/ 
    'print a["library_id"],a["cell_line"],a["library_name"],a["ref_genome"]'|sort -k2 | awk '{print $1"\t"$3}' >> ./"lib_info.txt"
done


import json
import subprocess
import sys

#TODO Boring error checking
url = sys.argv[1]
cmd = sys.argv[2]

# Could use a python library to make the call, check for 200 return, etc
#  but this should work on any system and not require anyhting but default
#  packages
a = json.loads(subprocess.check_output(['curl', '-s', url]))
exec cmd



### Data structure, to make dictionary containing dictionary with list:
In [49]: d = defaultdict(lambda: defaultdict(list))

In [50]: d
Out[50]: defaultdict(<function __main__.<lambda>>, {})

In [51]: d[tf_name_key]["chip"].append(fastq_acc) ## d[tf_name]["chip"].append(fastq_acc)

In [52]: d
Out[52]: 
defaultdict(<function __main__.<lambda>>,
            {'test': defaultdict(list, {'chip': 0})})

In [53]: d[tf_name_key]["control"].append(fastq_acc)














