import neo
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import decimate, find_peaks, savgol_filter
from plots.plots import *
import logging
import os

plt.subplots_adjust(bottom=0.3)
plt.rcParams.update({'font.size': 16})
# sns.set(font_scale=1.6, style='white')
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
    elif value[2] == "F":
        filename = os.path.join(PATH, "data\\Female_Rats", value[0])
    else: 
        filename = os.path.join(PATH, "data\\Control_ACSF", value[0])

    treatments = pd.DataFrame(value[3:-2]).T
    treatments.columns = list(COLORS.keys())

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
            logging.info("Cleanning big tresh..")
            self.clean_noise(self.meta["big_tresh"])
        logging.info(f"Sample rate: {self.sf}")
        self.base_stats = None
        self.fluo_stats = None
        self.r_KCN_stats = None
        self.leve_stats = None
        self.pera_stats = None

    # @property
    def get_data(self):
        return self.data

    # @property
    def get_peaks(self):
        return self.peak_neg

    # @property
    def get_sample_rate(self):
        return self.sf

    def get_treatments(self):
        """Return a dataframe of treatments intervals"""
        
        treatments = self.meta["treatments"]
        treatments.dropna(inplace=True, axis=1)
        columns = treatments.columns
        treatments = list(self.meta["treatments"].values[0])
        treatments.append(int(self.data.size / self.sf))
        t_times = []
        for i in range(0, len(treatments) - 1):  # get all intervals
            t_times.append((int(treatments[i]), int(treatments[i + 1])))
            # columns.append(columns[i])
        res = pd.DataFrame(t_times).T
        res.columns = columns
        logging.info(f"Treatments (start-end): \n{res}")
        return res

    def read_data(self):
        """Read Spike2 file"""
        logging.info("Reading file:")
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

        logging.info(f"Removed {n_removed} peaks.")
        pass


    def get_peaks(self):
        """Find sample peaks"""
        w = 33
        dist = 20
        width = [20]  # minimal width of a peak in samples
        # smooth the data to avoid high frequency peaks (spikes)
        f_data = savgol_filter(self.data.copy(), window_length=w, polyorder=1, deriv=0)

        peaks_neg = find_peaks(
            -f_data, height=self.meta["small_tresh"], width=width, distance=dist
        )[0]
        logging.info(f"Found total of {len(peaks_neg)} peaks")
        return np.asarray(peaks_neg)


    def get_peaks_by_interval(self, data_time):
        '''return the peaks around the interval t[0] to t[1]'''
        peaks = self.get_peaks()
        peaks_interval =  peaks[(peaks > data_time[0] * self.sf) & (peaks <= data_time[1] * self.sf)]
        logging.info(f"Found {len(peaks_interval)} peaks in the interval {data_time[0]} to {data_time[1]} => {(data_time[1]-data_time[0])/60} min")
        return peaks_interval

    def get_peaks_min(self, data_time, minutes):
        ''' Get the number of peaks per min'''
        pks = self.get_peaks_by_interval(data_time)
        pks_min =  len(pks) / minutes
        logging.info(f"Peaks per min: {pks_min} in {minutes} minutes")
        return pks_min


    def calculate_amplitudes_mean_std(self, data_time,  stats="mean"):
        """Calculate mean or Std amplitudes of peaks inside interval t"""
        pk = self.get_peaks_by_interval(data_time)
        if len(pk) < 3:
            logging.info(f"Number of peaks < 3, stats will not be calculated.")
            return 0
        amp = self.data[pk]
        if stats == "std":
            return amp.std()
        logging.info(f"Calculating amplitude {stats} for {len(pk)} peaks: {amp.mean()} in {amp.std()} minutes")
        return amp.mean()


    def get_treatment_interval(self, treatment_name):
        '''Get the entire treatment time split in seconds'''
        try:  # some files dont have all treatments
            t_times = self.get_treatments()[treatment_name].values
            logging.info(f"FOUND treatment data:  {treatment_name}.")    
        except:
            logging.info(f"NO TREATMENT DATA:  {treatment_name}.")
            return None
        
        return t_times


    # def get_baseline_data(self, minutes, stats="mean"):
    #     '''Get Baseline treatment peaks/min or mean amplitude or std on the minutes interval'''
        
    #     t_times = self.get_treatment_interval("Baseline")
    #     interval = t_times[1] - t_times[0]
    #     seconds = minutes * 60

    #     if interval < seconds:
    #         data_time = t_times
    #     else:
    #         data_time = [(t_times[1] - seconds), t_times[-1]]
    #         logging.info(f'Time for calculation {(data_time[0])/60:.1f} to {(data_time[1])/60:.1f} => {(data_time[1]-data_time[0])/60} min')

    #     if stats:
    #         return self.calculate_amplitudes_mean_std(data_time, stats)
        
    #     return self.get_peaks_min(data_time, minutes)        


    def get_treatment_data(self, treatment, minutes, stats=None, time=""):

        logging.info(f"\n\n{treatment}, stats: {stats}, time: {time}")
        t_times = self.get_treatment_interval(treatment)
        
        if t_times is None:
            return None
        interval = t_times[1] - t_times[0]
        seconds = minutes * 60
        logging.info(f"Time total: t_times {t_times[0]/60:.1f} to {t_times[1]/60:.1f} => {interval/60:.2f} min")
        
        if treatment == "Baseline":
            # return self.get_baseline_data(minutes, stats=stats)
            if interval < seconds:
                data_time = t_times
            else:
                data_time = [(t_times[1] - seconds), t_times[-1]]
                logging.info(f'Time for calculation {(data_time[0])/60:.1f} to {(data_time[1])/60:.1f} => {(data_time[1]-data_time[0])/60} min')
            if stats:
                return self.calculate_amplitudes_mean_std(data_time, stats)
            return self.get_peaks_min(data_time, minutes)        

        if time == "25_30":
            if interval >= 30:
                data_time = [(t_times[0] + (30*60 - seconds)), (t_times[0] + 30*60)]
                logging.info(f'Time for calculation = 25_30 {(data_time[0])/60:.1f} to {(data_time[1])/60:.1f} => {(data_time[1]-data_time[0])/60} min')
                if stats:
                    return self.calculate_amplitudes_mean_std(data_time, stats)
                return self.get_peaks_min(data_time, minutes)
            else:
                data_time = [(t_times[1] - seconds) , t_times[1]]
                logging.info(f'Time for calculation < 25_30 {(data_time[0])/60:.1f} to {(data_time[1])/60:.1f} => {(data_time[1]-data_time[0])/60} min')
                if stats:
                    return self.calculate_amplitudes_mean_std( data_time, stats)
                return self.get_peaks_min(data_time, minutes)

        elif time == "55_60":  # get the 55-60 min or the last 5 smaller 60
            if interval >= 60:
                data_time = [(t_times[0] + (60*60 - seconds)), (t_times[0] + 60*60)]
                logging.info(f'Time for calculation = 55_60 {(data_time[0])/60:.1f} to {(data_time[1])/60:.1f} => {(data_time[1]-data_time[0])/60} min')
                if stats:
                    return self.calculate_amplitudes_mean_std(data_time, stats)
                return self.get_peaks_min(data_time, minutes)
            elif interval > 45:
                data_time = [(t_times[1] - seconds) , t_times[1]]
                logging.info(f'Time for calculation > 55_60 {(data_time[0])/60:.1f} to {(data_time[1])/60:.1f} => {(data_time[1]-data_time[0])/60} min')
                if stats:
                    return self.calculate_amplitudes_mean_std(data_time, stats)
                return self.get_peaks_min(data_time, minutes)
            
        else:
            data_time = [(t_times[1] - seconds), t_times[1]]
            logging.info(f'Time for calculation last 5 min {(data_time[0])/60:.1f} to {(data_time[1])/60:.1f} => {(data_time[1]-data_time[0])/60} min')
            if stats:
                return self.calculate_amplitudes_mean_std(data_time, stats)
            return self.get_peaks_min(data_time, minutes)

    def plot_peaks(self, save_name=None):


        tvalues = self.meta['treatments'].values[0]
        tnames = self.meta['treatments'].columns
        tvalues = [t*self.sf for t in tvalues]
        neg = self.get_peaks()
        time = np.arange(0, self.data.size)/self.sf
        sns.set_style("whitegrid", {'axes.grid' : False})
        dy = savgol_filter(self.data, window_length=33, polyorder=1, deriv=0)
        
        _ , ax = plt.subplots(nrows=2, sharex=True, figsize=(20,12), 
                            gridspec_kw={'height_ratios': [5,1]})

        ax[0].plot(dy)
        ax[0].plot(neg, dy[neg], "x")
        ax[0].set_ylabel('Amplitude (uV)')
        ax[1].eventplot(positions=neg, data=dy, linewidths=.4, 
                        linelengths=.4)
        ax[1].set_xlim(0, dy.size)       
        ax[1].set_xlabel('Time (s)')
        ax[1].set_ylabel('Events')
        ax[1].set_xticks(ticks=tvalues, labels=tnames, rotation=-30)  

        if save_name:
            plt.savefig(f'{save_name}.jpg')
        else:
            plt.show()
        plt.close()