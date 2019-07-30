import pandas as pd
import os
import numpy as np

os.chdir("/home/anirbandutta")
df = pd.read_csv('Data.dat', sep ='\s+', index_col = False)
df['D'] = df['D'].replace('INDEF', np.nan)

#df['E'] = df['A']-df['D']
print(df)
