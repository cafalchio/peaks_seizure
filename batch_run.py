from utils.utils import *
from plots.plots import *
import pandas as pd

# Configuration
PATH = "C:\\Users\\cafa\\Documents\\mark_code\\laura_analysis"
DF_FILE = "rec_times.csv"
MINUTES = 5

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
    "sex": [],
}

for i, value in enumerate(df.values):
    # try:
    meta = extract_df_value(value, PATH)
    # Load data:
    print(f"\n{'#'*100}")
    print(f"Working on: {meta['filename']}")
    dt = Dataset(meta)
    results["id"].append(i)
    results["filename"].append(meta["filename"].split("\\")[-1])
    results["Channel"].append(meta["channel"])
    results["pk_min_base"].append(dt.peaks_minutes("Baseline", MINUTES))
    results["pk_min_fluo"].append(dt.peaks_minutes("Fluorocitrate", MINUTES))
    results["pk_min_RKCN"].append(dt.peaks_minutes("R_KCN", MINUTES))
    results["pk_min_levi"].append(dt.peaks_minutes("Levetiracetam", MINUTES))
    results["pk_min_pera"].append(dt.peaks_minutes("Perampanel", MINUTES))
    results["pk_min_RKCN_25_30"].append(dt.peaks_minutes("R_KCN", MINUTES, m25_30=True))
    results["pk_min_levi_25_30"].append(dt.peaks_minutes("Levetiracetam", MINUTES, m25_30=True))
    results["pk_min_pera_25_30"].append(dt.peaks_minutes("Perampanel", MINUTES, m25_30=True))
    results["pk_min_RKCN_55_60"].append(dt.peaks_minutes("R_KCN", MINUTES, m55_60=True))
    results["pk_min_levi_55_60"].append(dt.peaks_minutes("Levetiracetam", MINUTES, m55_60=True))
    results["pk_min_pera_55_60"].append(dt.peaks_minutes("Perampanel", MINUTES, m55_60=True))
    results["sex"].append(meta["sex"])
    # pplot(dt.get_data())
    # break
#     except:
#         print(f"Failed on: {meta['filename']}")
print(f"RESULTS >>> {len(results)}")
df_res = pd.DataFrame.from_dict(results)
df_res.to_csv("peak_rate_5min_tresh.csv")
print(df_res.head())
