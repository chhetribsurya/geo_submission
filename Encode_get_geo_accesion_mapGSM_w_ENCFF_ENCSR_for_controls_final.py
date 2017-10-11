import pandas as pd
import os
from collections import defaultdict
from os.path import join
from os.path import basename
from os.path import splitext
from os.path import expanduser

input_dir = expanduser("~/Dropbox/encode_3/geo_submission/GEO_ID_output")
output_dir = expanduser("~/Dropbox/encode_3/geo_submission/GEO_ID_output")

### Read GSM_ID xlsx file:
gsm_id_df = pd.read_excel(join(input_dir, "GSM_ID_w_TF_name_final.xlsx"))
gsm_id_df.columns = ["TF_name", "gsm_chip", "gsm_control"]

### Read ENCFF_ID xlsx file:
encff_id_df = pd.read_excel(join(input_dir, "ENCFF_ID_w_TF_name_final.xlsx"))
encff_id_df.columns = ["TF_name", "encff_chip", "encff_control"]

### Read ENCSR_ID xlsx file:
encsr_id_df = pd.read_excel(join(input_dir, "ENCSR_ID_w_TF_name_final.xlsx"))
encsr_id_df.columns = ["TF_name", "encsr_chip", "encsr_control"]

### Merge all the files on "TF_name" column:
df_list = [gsm_id_df, encff_id_df, encsr_id_df]
merged_id_df = reduce(lambda left, right: pd.merge(left, right, on=["TF_name"]), df_list)

### Extract GSM chip, ENCFF chip and ENCSR chip:
select_cols = ["TF_name", "gsm_chip", "encff_chip", "encsr_control"]
merged_final_chip_df = merged_id_df.loc[:, select_cols]
merged_final_chip_df.to_excel(join(output_dir, "GSM_ENCFF_ENCSR_ID_w_TF_name_for_chip_merged.xlsx"), index=False)

### Extract GSM control, ENCFF control and ENCSR control:
select_cols = ["TF_name", "gsm_control", "encff_control", "encsr_control"]
merged_final_control_df = merged_id_df.loc[:, select_cols]
merged_final_control_df.to_excel(join(output_dir, "GSM_ENCFF_ENCSR_ID_w_TF_name_for_control_merged.xlsx"), index=False)

### Since, we care about controls, eliminate duplicated rows:
control_df = merged_final_control_df.set_index("TF_name").loc[:, ["gsm_control", "encff_control", "encsr_control"]]
unique_control_df = control_df[~control_df.duplicated()]
unique_control_sorted_df = unique_control_df.sort_values(["gsm_control"])
unique_control_sorted_df.to_excel(join(output_dir, "Unique_GSM_ENCFF_ENCSR_ID_w_TF_name_for_control_merged_and_TFs_info.xlsx"))

