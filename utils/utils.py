import neo
import numpy as np
import pandas as pd
from scipy.signal import decimate, find_peaks, savgol_filter
import os

### Collection of class and function to load and preprocess data


def read_dataframe(df_path):
    """Read the information in the dataframe that will be used in the script"""
    df = pd.read_csv(df_path)
    df = df[
        [
            "filename",
            "Channel",
            "Sex",
            "Baseline",
            "Fluorocitrate",
            "R_KCN",
            "Levetiracetam",
            "Perampanel",
            "Big_Noise_Threshold",
            "Small_Peak_Threshold",
        ]
    ]
    return df


def extract_df_value(value, PATH):
    """Extract data from a dataframe line"""
    meta = {}
    COLORS = {
        "Baseline": "green",
        "Fluorocitrate": "red",
        "R_KCN": "blue",
        "Levetiracetam": "grey",
        "Perampanel": "brown",
    }

    if value[2] == "M":
        filename = os.path.join(PATH, "data\\Male_Rats", value[0])
    if value[2] == "F":
        filename = os.path.join(PATH, "data\\Female_Rats", value[0])

    treatments = pd.DataFrame(value[3:-2]).T
    treatments.columns = list(COLORS.keys())
    treatments.dropna(inplace=True, axis=1)

    meta["treatments"] = treatments
    meta["filename"] = filename
    meta["channel"] = value[1]
    meta["big_tresh"] = float(value[-2])
    meta["small_tresh"] = float(value[-1])
    meta["colors"] = COLORS
    meta["sex"] = value[2]

    return meta


class Dataset:
    """Class to manipulate data and execute data related methods"""

    def __init__(self, meta):
        self.meta = meta
        self.data, self.sf = self.read_data()
        self.data = self.data - self.data.mean()
        self.downsample()
        if self.meta["big_tresh"]:
            print("Cleanning big tresh..")
            self.clean_noise(self.meta["big_tresh"])
        # self.peak_neg = self.get_peaks()

    def get_data(self):
        return self.data

    def get_peaks(self):
        return self.peak_neg

    def get_sample_rate(self):
        return self.sf

    def get_treatments(self):
        """Return a dataframe of treatments intervals"""
        colors = list(self.meta["colors"].keys())
        treatments = list(self.meta["treatments"].values[0])
        treatments.append(int(self.data.size / self.sf))
        t_times = []
        columns = []
        for i in range(0, len(treatments) - 1):  # get all intervals
            t_times.append((int(treatments[i]), int(treatments[i + 1])))
            columns.append(colors[i])
        res = pd.DataFrame(t_times).T
        res.columns = columns
        return res

    def read_data(self):
        """Read Spike2 file"""
        print("Reading file:")
        reader = neo.io.Spike2IO(
            filename=self.meta["filename"], try_signal_grouping=False
        )
        bl = reader.read(lazy=False, load_waveforms=True)[0]
        for seg in bl.segments:
            for i, anasig in enumerate(seg.analogsignals):
                if i == self.meta["channel"] - 1:
                    data = anasig.flatten()
                    sf = int(anasig.sampling_rate)
                    return data, sf

    def downsample(self):
        """Downsample the data to a final sample of 250 samples/s"""
        factor = int(self.sf / 250)
        self.data = decimate(self.data, factor, ftype="fir")
        self.downsampled = True
        self.sf = int(self.sf / factor)
        return None

    def clean_noise(self, noise_tresh):
        """Attempt to clean big noise from the recording by finding
        peaks above or below treshold and remove 2s of data around peak.
        """
        pos_idxs = find_peaks(self.data, threshold=noise_tresh)[0]

        for idx in pos_idxs:
            try:
                self.data[idx - int(self.sf // 2) : idx + int(self.sf // 2)] = 0
            except:
                continue

        neg_idxs = find_peaks(-self.data, threshold=noise_tresh)[0]

        for idx in neg_idxs:
            try:
                self.data[idx - int(self.sf // 2) : idx + int(self.sf // 2)] = 0
            except:
                continue

        n_removed = len(pos_idxs) + len(neg_idxs)

        # In case the find peaks miss some
        pos_idxs = list(np.where(self.data >= noise_tresh))[0]
        neg_idxs = list(np.where(self.data <= -noise_tresh))[0]
        for idx in neg_idxs:
            self.data[idx - int(self.sf // 2) : idx + int(self.sf // 2)] = 0

        for idx in pos_idxs:
            self.data[idx - int(self.sf // 2) : idx + int(self.sf // 2)] = 0

        print(f"Removed {n_removed} peaks.")
        pass

    def get_peaks(self):
        """Find sample peaks"""
        w = 33
        dist = 10  # minimal dist in samples between peaks
        width = [3]  # minimal width of a peak in samples
        # smooth the data to avoid high frequency peaks (spikes)
        f_data = savgol_filter(self.data.copy(), window_length=w, polyorder=1, deriv=0)

        peaks_neg = find_peaks(
            -f_data, height=self.meta["small_tresh"], width=width, distance=dist
        )[0]

        return np.asarray(peaks_neg)

    def peaks_minutes(self, treatment_name, minutes, m25_30=False, m55_60=False):
        """return the peaks for specified time"""
        try:  # some files dont have all treatments
            t_times = self.get_treatments()[treatment_name].values
        except:
            return None
        peaks = self.get_peaks()
        times_s = np.arange(t_times[0], t_times[1])  # get the interval between t0 - t1
        times_m = np.array_split(
            times_s, int((t_times[1] - t_times[0]) // 60)
        )  # split times in 1 min

        if treatment_name == 'Baseline': # smaller then 5m baselines
            if len(times_m) < 5:
                minutes = len(t_times)
        try:
            if m25_30:
                t = [times_m[24][0], times_m[29][-1]]
            elif m55_60:
                t = [times_m[54][0], times_m[59][-1]]
            else:
                 t = [times_m[-minutes][0], times_m[-1][-1]]
        except:
            # print(f"Maximum size: {len(times_m)}")
            return None
        peaks_min = (
            len(
                peaks[(peaks > t[0] * self.sf) & (peaks <= t[-1] * self.sf)]
            )
            / 5
        )
        # print(peaks_min)
        return peaks_min