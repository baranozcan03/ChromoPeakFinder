import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import scipy
import math

###SETTINGS
RESOLUTION = 10000
JUMP = int(10e6 / RESOLUTION)
RESOLUTION_NAME = '10kb'
NORMALIZATION_NAME = 'VCsqrt'
RAW_OBSERVED_FILENAME = 'GSE63525_K562_intrachromosomal_contact_matrices_CHR21\chr21_10kb.RAWobserved'
NORMALIZATION_FILENAME = 'GSE63525_K562_intrachromosomal_contact_matrices_CHR21\chr21_10kb.SQRTVCnorm'
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



plt.imshow(intensity_matrix_log, cmap = 'YlOrBr')
# plt.clim(0, 1e5) # give min and max value to colorbar range
plt.title('Intensity Matrix Log')
plt.tight_layout()
plt.colorbar()
#plt.show()
plt.savefig('intensity_matrix_log_{0}_{1}.png'.format(RESOLUTION_NAME, NORMALIZATION_NAME))


plt.imshow(intensity_matrix_divide_coefficients, cmap = 'YlOrBr')
plt.clim(0, 1e3) # give min and max value to colorbar range
plt.title('Intensity Matrix Divide Norm Coefficients')
plt.tight_layout()
plt.colorbar()
# plt.show()
plt.savefig('intensity_matrix_{0}_{1}.png'.format(RESOLUTION_NAME, NORMALIZATION_NAME))
