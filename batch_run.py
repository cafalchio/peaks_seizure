from utils.utils import *
# from utils.config import config
import pandas as pd

# Configuration
PATH = "C:\\Users\\cafa\\Documents\\mark_code\\laura_analysis"
DF_FILE = "rec_times.csv"
minutes = 5

df_path = os.path.join(PATH, DF_FILE)
df = read_dataframe(df_path)
df.Channel = df.Channel.apply(lambda x: int(x.split("Ch")[-1]))

results = {
    "id": [],
    "filename": [],
    "Channel": [],
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
    "sex": [],
}

   
for i, value in enumerate(df.values):
    try:
        meta = extract_df_value(value, PATH)
        # Load data:
        print(f"\n{'#'*100}")
        print(f"Working on: {meta['filename']}")
        dt = Dataset(meta)
        treats = meta["treatments"]
        results["id"].append(i)
        results["filename"].append(meta["filename"].split("\\")[-1])
        results["Channel"].append(meta["channel"])

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
    # ####
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
        results["sex"].append(meta["sex"])

#     # break
    except:
        print(f"Failed on: {meta['filename']}")
print(f"RESULTS >>> {len(results)}")
df_res = pd.DataFrame.from_dict(results)
df_res.to_csv("peak_rate_5min_amp_std.csv")
print(df_res.head())
    

# if __name__ == "__main__":
#     run()
