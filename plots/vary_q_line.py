import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

includev0 = True

def drawfig(ds, data, x):
   versions = ["Ours", "v3", "v2", "v1", "v0"]
   dstyle = ['d', '*', 'x', 'o', '^']
   if (includev0):
      r = 5
   else:
      r = 4
   for i in range(r):
      plt.plot(x, data[i], marker=dstyle[i], linewidth=1.2, label=versions[i])

   plt.legend()
   plt.xlabel(r'$q$', fontsize=15)
   plt.ylabel('Time (sec)', fontsize=15)
   plt.xticks(x, fontsize=14) 
   plt.yticks(fontsize=14)


   with PdfPages(ds+'.pdf') as pdf:
      pdf.savefig()
   plt.savefig(ds+".jpg", dpi=300)

   plt.close()
def draw(ds, data, x1, x2):
   data1 = []
   data2=[]
   sz = len(data)
   for i in range(5):
      data1.append([data[x] for x in range(i, sz, 10)])
   for i in range(5):
      data2.append([data[x] for x in range(i+5, sz, 10)])
   drawfig(ds+"_q1", data1, x1)
   drawfig(ds+"_q2", data2, x2)

bitcoin=[
0.201, 0.54, 1.752, 0.554, 40.255, 62.329, 73.04, 78.592, 86.214, 1490.38, 0.719, 1.006, 0.744, 0.784, 28.16, 33.953, 42.962, 45.202, 50.426, 830.295, 0.509, 0.467, 0.518, 0.048, 16.921, 13.494, 19.862, 9.045, 24.042, 572.765, 0.093, 0.024, 0.329, 0.069, 5.114, 2.489, 5.859, 3.891, 6.497, 320.879, 0.007, 0.007, 0.051, 0.132, 3.522, 0.064, 0.199, 0.625, 0.992, 224.601, 0.005, 0.005, 0.005, 0.027, 1.893, 0.019, 0.039, 0.01, 0.021, 131.723, 
]
wikivote=[
0.057, 1.605, 1.001, 1.425, 60.011, 16.133, 20.516, 50.408, 60.399, 3821.49, 1.266, 1.441, 1.026, 0.916, 49.409, 0.218, 4.906, 0.59, 9.627, 2947.23, 0.098, 1.004, 0.753, 0.666, 39.798, 1.569, 1.45, 0.081, 1.235, 2298.08, 0.738, 0.711, 0.517, 0.041, 22.399, 0.936, 0.895, 0.031, 1.025, 451.929, 0.014, 0.02, 0.603, 0.557, 17.522, 0.629, 0.029, 0.446, 0.754, 776.253, 0.316, 0.257, 0.56, 0.18, 13.548, 0.722, 0.014, 0.033, 0.21, 533.504
]
ascaida=[
0.43, 2.346, 2.29, 2.374, 4.437, 0.863, 3.178, 3.126, 3.39, 16.286, 0.244, 1.208, 1.163, 1.221, 2.497, 0.406, 1.421, 1.38, 1.512, 10.452, 0.138, 0.542, 0.527, 0.547, 1.112, 0.174, 0.569, 0.554, 0.598, 6.764, 0.078, 0.215, 0.21, 0.219, 0.583, 0.066, 0.203, 0.197, 0.212, 4.168, 0.05, 0.094, 0.092, 0.097, 0.308, 0.026, 0.064, 0.063, 0.066, 3.044, 0.031, 0.044, 0.044, 0.043, 0.2, 0.016, 0.03, 0.029, 0.03, 1.88
]
mathoverflow=[
658.347, 2218.82, 2241.78, 3557.38, 11341, np.nan, np.nan, np.nan, np.nan, np.nan, 501.231, 1990.61, 1383.8, 1396.81, 8280.34, np.nan, np.nan, np.nan, np.nan, np.nan, 356.159, 1081.15, 814.304, 1105.18, 3567.22, np.nan, np.nan, np.nan, np.nan, np.nan, 288.637, 524.034, 529.533, 530.723, 3020.51, 19664, np.nan, np.nan, np.nan, np.nan, 197.458, 413.121, 351.036, 418.146, 2563.54, np.nan, np.nan, np.nan, np.nan, np.nan, 136.577, 250.031, 245.364, 251.886, 3429.29, 16630.9, np.nan, np.nan, np.nan, np.nan, 
]
x1=[10, 11, 12, 13, 14, 15]
x2=[15, 16, 17, 18, 19, 20]

draw("bitcoin", bitcoin, x1, x1)
draw("wiki-vote", wikivote, x1, x1)
draw("as-caida", ascaida, x1, x2)
draw("mathoverflow", mathoverflow, x1, x1)