import pandas as pd
import numpy as np

###SETTINGS
RESOLUTION = 10000
RAW_OBSERVED_FILENAME = 'GSE63525_HMEC_intrachromosomal_contact_matrices_CHR2\chr2_10kb.RAWobserved'
NORMALIZATION_FILENAME = 'GSE63525_HMEC_intrachromosomal_contact_matrices_CHR2\chr2_10kb.KRnorm'
PEAKS_FILENAME = 'GSE63525_HMEC_intrachromosomal_contact_matrices_CHR2\GSE63525_mouse_lymphoblasts_HiCCUPS_looplist.txt'
CHROMOSOME_NAME = 'chr2'


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


number_of_rows = int(max(intensity_map.keys())[0])+12

# initialize the matrix
intensity_matrix_divide_coefficients = np.zeros([number_of_rows + 1, number_of_rows + 1])

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



#gives a list of peak points. Every poin is represented by another list. 
#Namely the x and y coordinates but since we dont have 1 base resolution every point is actually an interval
#peak point1 = [[x1,x2],[y1,y2]]
#NAME OF THE PEAKS LIST = Peaks_wanted_split
with open(PEAKS_FILENAME) as f:
  Peaks_raw = f.read()
Peaks_raw_list = Peaks_raw.split("\n")
Peaks_wanted = []
for i in Peaks_raw_list[1::]:
  if i[0:5] == CHROMOSOME_NAME + "\t":
    Peaks_wanted.append(i)
Peaks_wanted_split = []
for i in Peaks_wanted:
  split_but_big= i.split("\t")
  split_and_small = [[int(split_but_big[1]),int(split_but_big[2])],[int(split_but_big[4]), int(split_but_big[5])]]
  Peaks_wanted_split.append(split_and_small)


print(Peaks_wanted_split)


data = {'intensities':[]}
# matrix[i][j]
i = 0
for j in range(0, number_of_rows + 1):
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


#check if the peak point is fit with the intended anomoly points  
peaks = []
i = 0
for j in range(0, number_of_rows + 1):
    intensity = intensity_matrix_divide_coefficients[i][j]
    if intensity in anom:
        for l in Peaks_wanted_split:
            if (i*RESOLUTION == l[0][0] or i*RESOLUTION == l[0][1] ) and (j*RESOLUTION == l[1][0] or j*RESOLUTION == l[1][1]):
                data['intensities'].append(intensity)
                anom.remove(intensity)
                peaks.append(intensity)

if len(anom) == 0:
    print("YES ALL OF THEM ACTUALLY PEAK POINTS!!!!")
else:
    print("Peak points:", peaks )
    print("Real anomalies:", anom)
  
