import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


data1=[
[2.059030753,	3.840227659,	7.407488448,	14.52111242,	27.78988476],
[1.850587288,	3.5798704,	6.946726039,	13.7992117,	26.45505667],
[1.881702082,	3.643692622,	7.02600717,	13.70127756,	22.4844921],
[2.014761344,	3.930607579,	7.180600298,	12.06938273,	21.02357954],
[1.841026366,	3.401907001,	6.816278461,	14.86586758,	15.37485242],
]

data2=[
[1.853962822,	3.776112333,	7.390437719,	12.83678399,	24.04718811],
[1.846515287,	3.646647532,	7.057089163,	13.95277021,	26.64070629],
[1.883754111,	3.586748319,	7.627239792,	12.16570358,	19.75076065],
[1.938251366,	3.542435708,	6.147723485,	10.16444395,	21.07931429],
[1.666488223,	3.215518289,	6.180915674,	11.98120897,	22.29191238],
]


fig, ax = plt.subplots(figsize=(9,5))
size = 5
x = np.arange(size)
total_width, n = 0.5, 4
width = total_width / n
labels = ['arabic-2005\n($q = 1000$)', 'uk-2005\n($q = 200$)', 'it-2004\n($q = 1000$)', 'webbase-2001\n($q = 300$)', 'clue-web\n($q = 10$)']
plt.bar(x - width*1.5,[i[0] for i in data1],width,label='2 threads', color='dodgerblue', alpha = 0.5)
plt.bar(x-width*0.5,[i[1] for i in data1],width,label='4 threads', color='mediumslateblue', alpha = 0.5)
plt.bar(x+width*0.5,[i[2] for i in data1],width,label='8 threads', color='sandybrown', alpha = 0.5)
plt.bar(x+width*1.5,[i[3] for i in data1],width,label='16 threads', color='r', alpha = 0.5)
plt.bar(x+width*1.5,[i[4] for i in data1],width,label='32 threads', color='crimson', alpha = 0.5)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=15)
# plt.legend()
# plt.xlabel('network')
plt.ylabel('Speedup ratio', fontsize=15)
#plt.xticks(x)
y_ticks = [0, 1, 2, 4, 8, 16, 32]
plt.yticks(y_ticks)
plt.legend(borderpad=1, shadow=True, loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=4, fontsize=14)

fig.tight_layout()

with PdfPages('speedup_k=2.pdf') as pdf:
   pdf.savefig()
plt.savefig("speedup1.jpg", dpi=300)


plt.close()
# fig, ax = plt.subplots(figsize=(9,5))
# size = 5
# x = np.arange(size)
# total_width, n = 0.5, 4
# width = total_width / n
# labels = ['arabic-2005\n($q = 3000$)', 'uk-2005\n($q = 400$)', 'it-2004\n($q = 3000$)', 'webbase-2001\n($q = 500$)', 'clue-web\n($q = 20$)']
# plt.bar(x - width*1.5,threads_2_3,width,label='2 threads', color='dodgerblue', alpha = 0.5)
# plt.bar(x-width*0.5,threads_4_3,width,label='4 threads', color='mediumslateblue', alpha = 0.5)
# plt.bar(x+width*0.5,threads_8_3,width,label='8 threads', color='sandybrown', alpha = 0.5)
# plt.bar(x+width*1.5,threads_16_3,width,label='16 threads', color='r', alpha = 0.5)
# #plt.rcParams['text.usetex'] = True
# ax.set_xticks(x)
# ax.set_xticklabels(labels, fontsize=15)
# # plt.legend()
# # plt.xlabel('network')
# plt.ylabel('Speedup ratio', fontsize=15)
# #plt.xticks(x)
# y_ticks = [0, 1, 2, 4, 8, 16]
# plt.yticks(y_ticks)
# plt.legend(borderpad=1, shadow=True, loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=4, fontsize=14)

# fig.tight_layout()
# plt.savefig("speedup2.jpg", dpi=300)
# with PdfPages('speedup_k=3.pdf') as pdf:
#    pdf.savefig()
