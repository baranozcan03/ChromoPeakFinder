import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import scipy
import math

###SETTINGS
RESOLUTION = 1e6
JUMP = int(10e6 / RESOLUTION)
RESOLUTION_NAME = '1mb'
NORMALIZATION_NAME = 'KR'
RAW_OBSERVED_FILENAME = 'GSE63525_K562_intrachromosomal_contact_matrices_CHR21\chr21_1mb.RAWobserved'
NORMALIZATION_FILENAME = 'GSE63525_K562_intrachromosomal_contact_matrices_CHR21\chr21_1mb.KRnorm'
PEAKS_FILENAME = 'GSE63525_K562_intrachromosomal_contact_matrices_CHR21\GSE63525_K562_HiCCUPS_looplist_with_motifs.txt'
CHROMOSOME_NAME = 'chr21'
STEP_SIZE = 10 # fixed distance between different diagonals in the matrix
PLOT_AMOUNTS = 15 #chose how many plots do you want starting from the main idagonal


raw_observed_data = pd.read_csv(RAW_OBSERVED_FILENAME, delimiter='\t', names=['row', 'column', 'intensity'])
norm_coefficients = pd.read_csv(NORMALIZATION_FILENAME, delimiter='\t', names=['coefficient'])

# (row, column) -> intensity
intensity_map = {}

# build the intensity_map
for index, row in raw_observed_data.iterrows():
  r, c, intensity = row
  r = int(r / RESOLUTION)
  c = int(c / RESOLUTION)
  intensity_map[(r, c)] = intensity

number_of_rows = max(int(max(intensity_map.keys())[0]),int(max(intensity_map.keys())[1]))+3


# initialize the matrix
intensity_matrix_divide_coefficients = np.zeros([number_of_rows + 1, number_of_rows + 1])
intensity_matrix_log = np.zeros([number_of_rows + 1, number_of_rows + 1])


# build the matrix with dividing intensities by normalization coefficients
for indices, intensity in intensity_map.items():
  i, j = indices
  i = int(i)
  j = int(j)
  # handle NaN values
  if norm_coefficients.iloc[i]['coefficient'] == 'NaN' or norm_coefficients.iloc[j]['coefficient'] == 'NaN':
    intensity = 0
  else:
    # normalize by dividing
    if norm_coefficients.iloc[i]['coefficient'] != 0.0:
      intensity /= norm_coefficients.iloc[i]['coefficient']
    if norm_coefficients.iloc[j]['coefficient'] != 0.0:
      intensity /= norm_coefficients.iloc[j]['coefficient']
  # symmetric matrix, M[i][j] = M[j][i]
  intensity_matrix_divide_coefficients[i][j] = intensity
  intensity_matrix_divide_coefficients[j][i] = intensity
  # take the log if > 0
  if intensity > 0.0:
    intensity = math.log(intensity)
  # symmetric matrix, M[i][j] = M[j][i]
  intensity_matrix_log[i][j] = intensity
  intensity_matrix_log[j][i] = intensity

plt.imshow(intensity_matrix_divide_coefficients, cmap = 'YlOrBr')
plt.clim(0, 1e5) # give min and max value to colorbar range
plt.title('Intensity Matrix Divide Norm Coefficients')
plt.tight_layout()
plt.colorbar()
# plt.show()
plt.savefig('intensity_matrix_{0}_{1}.png'.format(RESOLUTION_NAME, NORMALIZATION_NAME))
plt.clf()

plt.imshow(intensity_matrix_log, cmap = 'YlOrBr')

# plt.clim(0, 1e5) # give min and max value to colorbar range
plt.title('Intensity Matrix Log')
plt.tight_layout()
plt.colorbar()
#plt.show()
plt.savefig('intensity_matrix_log_{0}_{1}.png'.format(RESOLUTION_NAME, NORMALIZATION_NAME))



def histdrawer(y):
  data = {'intensities':[]}
  # matrix[i][j]
  i = 0
  for j in range(y, number_of_rows + 1):
    intensity = intensity_matrix_divide_coefficients[i][j]
    if intensity !=0:
      data['intensities'].append(intensity)    
    i += 1
  
  #play with a to see how much of the anomalies you want to get rid of
  #you can check the list in the outpust in order to see weather if the desired 
  #values mistaken as anomalies
  a = float(len(data['intensities'])) * 0.025
  x = 0
  anom = []
  while x < int(a):
    h = max(data['intensities'])
    data['intensities'].remove(h)
    anom.append(h)
    x = x+1 


  df = pd.DataFrame(data)
  sns.histplot(data=df, x='intensities', kde=True).set(title=f'Resolution: {int(RESOLUTION)}, Starting Column Index: {int(y*RESOLUTION)}')
  print(anom)

plt.clf()
with  PdfPages('{0}_{1}.pdf'.format(RESOLUTION_NAME, NORMALIZATION_NAME)) as hists:
  t = 1
  for i in range(0,PLOT_AMOUNTS*JUMP,JUMP):
    figu = plt.figure(t)
    histdrawer(i)
    hists.savefig(figu)
    t += 1

