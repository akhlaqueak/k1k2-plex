import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams.update({'font.size': 14})

includev0 = True

def drawfig(ds, data, x):
   plt.figure(figsize=(6, 5))
   versions = ["Ours", "Base-3", "Base-2", "Base-1", "Base-0"]
   dstyle = ['d', '*', 'x', 'o', '^']
   if (includev0):
      r = 5
   else:
      r = 4
   for i in range(r):
      plt.plot(x, data[i], marker=dstyle[i], linewidth=1.2, label=versions[i])

   plt.legend(fontsize=15)
   plt.xlabel(r'$q$')
   plt.ylabel('Time (sec)')
   plt.xticks(x) 
   plt.yticks()


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
1061.73, 3175.34, 542.362, 3229.32, 5923.51, 324.502, 5899.54, 5904.78, 1648.16, 2194.54, 791.665, 2124.93, 464.726, 2157.16, 1377.78, 324.902, 4775.04, 4765.76, 4838.88, 6132.48, 557.652, 1314.07, 1321.49, 1337.2, 672.25, 245.753, 558.164, 735.897, 3365.12, 5896.5, 72.801, 829.927, 831.231, 152.006, 4074.99, 843.888, 2168.22, 2199.88, 2221.01, 5268.76, 59.631, 102.469, 549.439, 553.464, 3567.7, 585.016, 1474.43, 1488.48, 1515.16, 4800.47, 42.36, 385.116, 379.826, 84.494, 3140.82, 390.117, 187.237, 1055.51, 1089.58, 4277.53, 
]
webgoogle=[
8.48, np.nan, np.nan, np.nan, np.nan, 9.536, np.nan, np.nan, np.nan, np.nan, 8.791, np.nan, np.nan, np.nan, np.nan, 2.823, np.nan, np.nan, np.nan, np.nan, 0.969, np.nan, np.nan, np.nan, np.nan, 4.744, np.nan, np.nan, np.nan, np.nan, 6.453, np.nan, np.nan, np.nan, np.nan, 3.087, np.nan, 18205.2, np.nan, np.nan, 5.926, np.nan, np.nan, 18083.2, np.nan, 2.758, np.nan, np.nan, np.nan, np.nan, 1.855, 15661.4, np.nan, np.nan, np.nan, 2.044, 15675.1, np.nan, np.nan, np.nan, 
]
x1=[10, 11, 12, 13, 14, 15]
x2=[15, 16, 17, 18, 19, 20]

draw("bitcoin", bitcoin, x1, x1)
draw("wiki-vote", wikivote, x1, x1)
draw("as-caida", ascaida, x1, x2)
draw("mathoverflow", mathoverflow, x1, x1)
draw("web-google", webgoogle, x2, x2)