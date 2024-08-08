import sys
import joblib
import pandas as pd

# How to consolidate the gt distribution for a single joblib file
data = joblib.load(sys.argv[1])
programs = ['kanpig', 'sniffles', 'svjedi', 'cutesv']
rows = [data['ts'][_]['gt_dist'] for _ in programs]
rows.append(data['base_gt_dist'])
frame = pd.concat(rows)
frame['program'] = programs + ['base']
frame.reset_index(drop=True, inplace=True)
print(frame)
