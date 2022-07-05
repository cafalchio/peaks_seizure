from utils.utils import *
# from utils.config import config
import pandas as pd
from tqdm import tqdm
import logging
import warnings
warnings.filterwarnings("ignore")



# log config
logging.basicConfig(filename="batch_log.txt",
                    filemode='a',
                    format='%(message)s',
                    level=logging.INFO)


# Configuration
PATH = "C:\\Users\\cafa\\Documents\\mark_code\\laura_analysis"
DF_FILE = "rec_times.csv"
minutes = 10

df_path = os.path.join(PATH, DF_FILE)
df = read_dataframe(df_path)
df.Channel = df.Channel.apply(lambda x: int(x.split("Ch")[-1]))

results = {
    "id": [],
    "filename": [],
    "Channel": [],
    "sex": [],
    "pk_min_base": [],
    "pk_min_fluo": [],
    "pk_min_RKCN": [],
    "pk_min_levi": [],
    "pk_min_pera": [],
    "pk_min_RKCN_25_30": [],
    "pk_min_levi_25_30": [],
    "pk_min_pera_25_30": [],
    "pk_min_RKCN_55_60": [],
    "pk_min_levi_55_60": [],
    "pk_min_pera_55_60": [],
    ###
    "amp_min_base": [],
    "std_min_base": [],
    "amp_min_fluo": [],
    "std_min_fluo": [],
    "amp_min_RKCN": [],
    "std_min_RKCN": [],
    "amp_min_levi": [],
    "std_min_levi": [],
    "amp_min_pera": [],
    "std_min_pera": [],
    ##
    "amp_min_RKCN_25_30": [],
    "std_min_RKCN_25_30": [],
    "amp_min_levi_25_30": [],
    "std_min_levi_25_30": [],
    "amp_min_pera_25_30": [],
    "std_min_pera_25_30": [],
    "amp_min_RKCN_55_60": [],
    "std_min_RKCN_55_60": [],
    "amp_min_levi_55_60": [],
    "std_min_levi_55_60": [],
    "amp_min_pera_55_60": [],
    "std_min_pera_55_60": [],
    
}

print('Running batch analysis ..')
for i, value in enumerate(tqdm(df.values)):
    # try:
    meta = extract_df_value(value, PATH)
    # Load data:
    logging.info(f"\n{'#'*100}")
    logging.info(f"Working on: {meta['filename']}")
    dt = Dataset(meta)

## Logging results
    treats = meta["treatments"]
    results["id"].append(i)
    results["filename"].append(meta["filename"].split("\\")[-1])
    results["Channel"].append(meta["channel"])
    results["sex"].append(meta["sex"])
### peaks/min
    results["pk_min_base"].append(dt.get_treatment_data("Baseline", minutes))
    results["pk_min_fluo"].append(dt.get_treatment_data("Fluorocitrate", minutes))
    results["pk_min_RKCN"].append(dt.get_treatment_data("R_KCN", minutes))
    results["pk_min_levi"].append(dt.get_treatment_data("Levetiracetam", minutes))
    results["pk_min_pera"].append(dt.get_treatment_data("Perampanel", minutes))
    results["pk_min_RKCN_25_30"].append(dt.get_treatment_data("R_KCN", minutes, time="25_30"))
    results["pk_min_levi_25_30"].append(dt.get_treatment_data("Levetiracetam", minutes, time="25_30"))
    results["pk_min_pera_25_30"].append(dt.get_treatment_data("Perampanel", minutes, time="25_30"))
    results["pk_min_RKCN_55_60"].append(dt.get_treatment_data("R_KCN", minutes, time="55_60"))
    results["pk_min_levi_55_60"].append(dt.get_treatment_data("Levetiracetam", minutes, time="55_60"))
    results["pk_min_pera_55_60"].append(dt.get_treatment_data("Perampanel", minutes, time="55_60"))
### mean amplitude
    results["amp_min_base"].append(dt.get_treatment_data("Baseline", minutes, stats="mean"))
    results["amp_min_fluo"].append(dt.get_treatment_data("Fluorocitrate", minutes, stats="mean"))
    results["amp_min_RKCN"].append(dt.get_treatment_data("R_KCN", minutes, stats="mean"))
    results["amp_min_levi"].append(dt.get_treatment_data("Levetiracetam", minutes, stats="mean"))
    results["amp_min_pera"].append(dt.get_treatment_data("Perampanel", minutes, stats="mean"))
    results["amp_min_RKCN_25_30"].append(dt.get_treatment_data("R_KCN", minutes, time="25_30", stats="mean"))
    results["amp_min_levi_25_30"].append(dt.get_treatment_data("Levetiracetam", minutes, time="25_30", stats="mean"))
    results["amp_min_pera_25_30"].append(dt.get_treatment_data("Perampanel", minutes, time="25_30", stats="mean"))
    results["amp_min_RKCN_55_60"].append(dt.get_treatment_data("R_KCN", minutes, time="55_60", stats="mean"))
    results["amp_min_levi_55_60"].append(dt.get_treatment_data("Levetiracetam", minutes, time="55_60", stats="mean"))
    results["amp_min_pera_55_60"].append(dt.get_treatment_data("Perampanel", minutes, time="55_60", stats="mean"))
### std amplitude
    results["std_min_base"].append(dt.get_treatment_data("Baseline", minutes, stats="std"))
    results["std_min_fluo"].append(dt.get_treatment_data("Fluorocitrate", minutes, stats="std"))
    results["std_min_pera"].append(dt.get_treatment_data("Perampanel", minutes, stats="std"))
    results["std_min_RKCN"].append(dt.get_treatment_data("R_KCN", minutes, stats="std"))
    results["std_min_levi"].append(dt.get_treatment_data("Levetiracetam", minutes, stats="std"))
    results["std_min_RKCN_25_30"].append(dt.get_treatment_data("R_KCN", minutes, time="25_30", stats="std"))
    results["std_min_levi_25_30"].append(dt.get_treatment_data("Levetiracetam", minutes, time="25_30", stats="std"))
    results["std_min_pera_25_30"].append(dt.get_treatment_data("Perampanel", minutes, time="25_30", stats="std"))
    results["std_min_RKCN_55_60"].append(dt.get_treatment_data("R_KCN", minutes, time="55_60", stats="std"))
    results["std_min_levi_55_60"].append(dt.get_treatment_data("Levetiracetam", minutes, time="55_60", stats="std"))
    results["std_min_pera_55_60"].append(dt.get_treatment_data("Perampanel", minutes, time="55_60", stats="std"))
    savename = meta["filename"].split("\\")[-1][:-4] + "_Ch" + str(meta["channel"])
    dt.plot_peaks(f"plot_results/{savename}")
    # dt.plot_peaks(f"plot_results/test")
    # break
    # except:
    #     print(f"Failed on: {meta['filename']}")


df_res = pd.DataFrame.from_dict(results)
f = 'files'
if len(df_res) <=1:
    f = 'file'
print(f"Done! Analysed {len(df_res)} {f}.")
df_res.to_csv("test.csv")
# print(df_res.head())
    

# if __name__ == "__main__":
#     run()
