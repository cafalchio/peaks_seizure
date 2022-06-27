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
        dist = 10
        width = [10]
        f_data = savgol_filter(self.data.copy(), window_length=w, polyorder=1, deriv=0)

        peaks_neg = find_peaks(
            -f_data, height=self.meta["small_tresh"], width=width, distance=dist
        )[0]

        return np.asarray(peaks_neg)

    def peaks_minutes(self, treatment_name, minutes, extra=False):
        """return the peaks for specified time"""
        t_times = self.treatment_start_end() # should be in seconds
        peaks = self.get_peaks()
        results = []
        for i, treat_name in enumerate(self.meta["treatments"].columns):
            if treat_name == treatment_name:
                print(treat_name)
                treat_peaks = peaks[
                    (peaks > t_times[i][0] * self.sf)
                    & (peaks < t_times[i][1] * self.sf)
                ]
                peaks_min = np.array_split(
                    treat_peaks, int((t_times[i][1] - t_times[i][0]) // 60)
                )  # spliting peaks into all minutes in treatment interval
                results = [len(pks) for pks in peaks_min]  # get peaks/min
                if extra:  # Check if there is 2 intervals
                    return sum(results[25:30]) / 5  # from 25-30 min
                else:
                    return sum(results[-minutes:]) / minutes  # last minutes

    def treatment_start_end(self):

        treatments = list(self.meta["treatments"].values[0])
        treatments.append(int(self.data.size / self.sf))  # add the last point
        t_times = []
        for i in range(0, len(treatments) - 1):  # get all intervals
            t_times.append((int(treatments[i]), int(treatments[i + 1])))

        return t_times

