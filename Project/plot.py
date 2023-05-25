import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Code that plots a histogram of the data in the file specified

filename = input('Enter the name of the file to plot: ')
p= Path(__file__).with_name(filename)
data = np.genfromtxt(p, delimiter=',', skip_header=1, names=['lb', 'ub', 'count'])

# Plot the data in a bar chart
plt.bar(data['lb'], data['count'], width=data['ub']-data['lb'], align='edge', edgecolor='black')
plt.xlabel('Nr. of suceptible people')
plt.ylabel('Count')
plt.title('Histogram of the number of susceptible people')
plt.show()


