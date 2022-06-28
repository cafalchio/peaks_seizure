# Batch Peaks Detection in Mitocondrial model of Seizures

The script detects the peaks per minutes, peaks amplitude and plots the data for file referenced on the input csv.


## Using

A CSV with the treatment times and noise treshold should be in the same folder of **batch_run.py**
The data should be in a data folder, divided into male and female folders. ( in the future, a config file?)

## Results

A result CSV will be generated with the peaks/min and peaks amplitude for each file and treatment.


The peaks were calculated at the 5 last minutes of each treatment and also for 25-30 and 55-60 for selected treatments


![image](https://user-images.githubusercontent.com/47299428/176182606-7cbf5396-9292-44f8-9535-d0842123ea2c.png)

*todo:*
<ul>
  <li>Plot pictures</li>
  <li>Config file</li>
  <li>Calculate amplitudes</li>
</ul>


**Code style**
Black

