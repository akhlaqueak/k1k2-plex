import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams.update({'font.size': 14})
def draw(data, ds):
   x=[1,10, 100, 1000, 10000, 20000, 50000, 100000]
   plt.plot(x, data)
   plt.legend(fontsize=15)
   plt.xlabel(r'$\tau_{time}$')
   plt.ylabel('Time (sec)')
   plt.xticks(x) 
   plt.yticks()
   with PdfPages(ds+'.pdf') as pdf:
      pdf.savefig()
   plt.savefig(ds+".jpg", dpi=300)
   plt.close()

data=[
30.992, 51.82, 6.185, 2.419, 51.375, 21.279, 140.914, 76.306, 11.458, 6.304, 31.103, 51.915, 5.969, 2.498, 52.992, 23.132, 140.184, 73.459, 11.378, 6.294, 30.943, 51.608, 6.075, 2.463, 50.824, 22.662, 140.776, 71.893, 11.389, 6.284, 30.723, 51.115, 5.405, 2.503, 52.919, 22.916, 140.935, 71.826, 11.392, 6.283, 31.113, 51.93, 5.084, 2.328, 51.647, 21.853, 140.802, 72.624, 11.364, 6.287, 30.781, 51.293, 4.858, 2.336, 52.569, 22.551, 140.807, 72.577, 11.395, 6.277, 30.896, 51.978, 4.859, 2.412, 51.487, 23.771, 140.561, 71.817, 11.396, 6.281, 30.765, 51.233, 4.923, 2.314, 53.755, 22.312, 141.438, 72.196, 11.424, 6.302
]
ds = ["webbase-2001", "clue-web", "it-2004", "arabic-2005", 'uk-2005']

for i in range(len(ds)):
   data1=[data[i] for i in range(i, len(data), 10)]
   data2=[data[i] for i in range(i+1, len(data), 10)]
   draw(data1, ds[i]+"_tau_q1")
   draw(data2, ds[i]+"_tau_q2")
