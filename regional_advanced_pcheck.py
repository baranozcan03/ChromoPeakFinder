import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import scipy
import math
import collections

PEAKS_FILENAME = 'GSE63525_K562_intrachromosomal_contact_matrices_CHR21\GSE63525_K562_HiCCUPS_looplist_with_motifs.txt'
CHROMOSOME_NAME = '21'
RESOLUTION_NAME = 'VCsqrt 29-31 peak matrix'
NORMALIZATION_NAME = 'VCsqrt'
RAW_OBSERVED_FILENAME = 'GSE63525_K562_intrachromosomal_contact_matrices_CHR21\chr21_10kb.RAWobserved'
NORMALIZATION_FILENAME = 'GSE63525_K562_intrachromosomal_contact_matrices_CHR21\chr21_10kb.SQRTVCnorm'
RESOLUTION = 10000
MATRIX_EXTENDER = 6
STRIPE_LENGHT = 15
RANGE_START = 29e6
RANGE_END = 31e6


def peak_finder(filename,chrname,resolution): 
    with open(filename) as f:
        Peaks_raw = f.read()
    Peaks_raw_list = Peaks_raw.split("\n")
    Peaks_wanted = []
    for i in Peaks_raw_list[1::]:
        if i[0:3] == chrname + "\t":
            Peaks_wanted.append(i)
    Peaks_wanted_split = []
    for i in Peaks_wanted:
        split_but_big= i.split("\t")
        split_and_small = [[int(split_but_big[1]),int(split_but_big[2])],[int(split_but_big[4]), int(split_but_big[5])]]
        Peaks_wanted_split.append(split_and_small)
    Peaks_wanted_pure = []
    for i in Peaks_wanted_split:
        an_item = [i[0][0],i[1][0]]
        Peaks_wanted_pure.append(an_item)
    Peaks_wanted_resolutioned = []
    for i in Peaks_wanted_pure:
        if RANGE_START<=i[0]<=RANGE_END and RANGE_START<i[1]<RANGE_END:
            Peaks_wanted_resolutioned.append([math.floor(i[0]/resolution),math.floor(i[1]/resolution)])
    return Peaks_wanted_resolutioned


Peaks_wanted_resolutioned = peak_finder(PEAKS_FILENAME,CHROMOSOME_NAME,RESOLUTION)

The_Ultimate_Hash = collections.defaultdict(list)

for i in Peaks_wanted_resolutioned:
    The_Ultimate_Hash[i[1]-i[0]].append([i[0],i[1]])

List_of_start_j_and_num_of_peaks = []

for i in The_Ultimate_Hash:
    List_of_start_j_and_num_of_peaks.append([i,len(The_Ultimate_Hash[i])])

List_of_start_j_s = []

for i in The_Ultimate_Hash:
    List_of_start_j_s.append(i)

print(List_of_start_j_and_num_of_peaks)















#
#
#
#
#
#BUILD THE MATRIX
#
#
#
#
#

raw_observed_data = pd.read_csv(RAW_OBSERVED_FILENAME, delimiter='\t', names=['row', 'column', 'intensity'])
norm_coefficients = pd.read_csv(NORMALIZATION_FILENAME, delimiter='\t', names=['coefficient'])

# (row, column) -> intensity
intensity_map = {}

# build the intensity_map
for index, row in raw_observed_data.iterrows():
  r, c, intensity = row
  if (RANGE_START <= r <= RANGE_END)and(RANGE_START <= c <= RANGE_END) :
    r = int(r / RESOLUTION)
    c = int(c / RESOLUTION)
    intensity_map[(r, c)] = intensity

number_of_rows = max(int(max(intensity_map.keys())[0]),int(max(intensity_map.keys())[1]))+30


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






def histdrawer(y,list_of_peaks,len_stirpe):
    data = {'intensities':[]}
    # matrix[i][j]
    i = 0
    for j in range(y, number_of_rows + 1):
        intensity = intensity_matrix_divide_coefficients[i][j]
        if intensity !=0:
            data['intensities'].append(intensity)    
        i += 1
    print(len(data['intensities']))

    #:D
    stripes = []
    for k in list_of_peaks:
        i_value = k[0]
        j_value = k[1]
        stripes.append(intensity_matrix_divide_coefficients[i_value][j_value])
  
    
    X = []
    for i in stripes:
        if i != 0:
            X.append([i,i])

    Y = [0,len_stirpe]
  

    for i in X:
        plt.plot(i,Y)
  

    df = pd.DataFrame(data)
    sns.histplot(data=df, x='intensities', kde=True).set(title=f'Resolution: {int(RESOLUTION)}, Starting Column Index: {int(y*RESOLUTION)}')
    
    
  
  #############################print(anom) 



plt.clf()
with  PdfPages('{0}_{1}.pdf'.format(RESOLUTION_NAME, NORMALIZATION_NAME)) as hists:
  t = 1
  for i in List_of_start_j_s:
    figu = plt.figure(t)
    histdrawer(i,The_Ultimate_Hash[i],STRIPE_LENGHT)
    hists.savefig(figu)
    t += 1

