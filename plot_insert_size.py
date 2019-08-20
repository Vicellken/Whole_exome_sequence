## the python version & matplotlib version is depending on the system
## in this case we use python2, since the matplotlib version on our server
## only support python2. But i recommand stick with python3

## pre-selecting insert data from .sam file
## hereby use 'awk' in linux environment

#awk '{print $9}' NexteraB.mapped.sam | awk '$1 > 0' > pyplot.insert.size.NexteraB


## below is python script in linux environment

#module load python2
#module load matplotlib

#!/usr/bin/env python2

import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors

with open('pyplot.insert.size.NexteraA', 'r') as hist:
	data = hist.read().splitlines()

value = [int(x) for x in data if int(x) <= 1000]

plt.figure(figsize = [20, 20])
plt.hist(value, bins = 1000)
plt.title('Distribution of Insert Sizes', fontsize = 18)
plt.xlabel('Insert Size', fontsize = 14)
plt.ylabel('Frequency', fontsize = 14)

#plt.axes(facecolor = '#E6E6E6')
plt.grid(color = 'r', linestyle = '-')
plt.savefig('insert_size_NexteraA.png')
