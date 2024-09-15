import sys
import os
os.getcwd()
os.chdir("C:\\Users\\Karan\\Documents\\Project\\work")
baseline=open('HGVS_2020_big_VEP_baseline.tsv',"r")


hgvs_ids = []
ref_aas = []
mut_aas = []
for line in baseline:

    if line.startswith('#'):
        continue
    else:
        hgvs_ids.append(line.split('\t')[0])
        ref_aas.append(line.split('\t')[1][0])
        mut_aas.append(line.split('\t')[1][2])
        

Blosum_matrix=open('BLOSUM62.txt',"r")

aas=[]
aa_score=[]
for line in Blosum_matrix:
    if line.startswith('#') or line.startswith('x'):
        continue
    else:
        aas.append(line.strip('\n').strip()[0])
        aa_score.append(line.strip('\n').split()[1:len(line)])


blosum_dict={}
for aa in aas:
    blosum_dict[aa] = {}
    
for i in range(len(aas)):
    for j in range (len(aa_score)):
        score=int(aa_score[i][j])
        blosum_dict[aas[i]][aas[j]] = score


scores=[]
for i in range (len(hgvs_ids)):
    ref_aa=ref_aas[i]
    mut_aa=mut_aas[i]
    
    score=blosum_dict[ref_aa][mut_aa]
    scores.append(str(score))
    
scores
        
f=open('out_baseline.txt','w')
f.write('# ID\tScore\n')

for id, score in zip(hgvs_ids, scores):
    f.write(id + '\t' + score + '\n')
f.close()


#####################
#####################

import numbers
import os
import sys
import matplotlib.collections
import matplotlib.pyplot
import argparse
import numpy

matplotlib.use('AGG')

global type_predictor

filename = "out_baseline.txt"

if 'sift' in filename:
    type_predictor = 'sift'
elif 'polyphen' in filename:
    type_predictor = 'polyphen'
elif 'baseline' in filename:
    type_predictor = 'BLOSUM'
else:
    type_predictor = ''

predictor_dict = {}

f = open(filename, 'r')
    # Total bytes in the file (end of file)
eof = f.seek(0, 2)
print(eof)
    # Go to the beginning of the file again
f.seek(0)
    # Read the first line (should be the header)
f.readline()
    # Get the current position of the file pointer
cur = f.tell()
print(cur)
    # If the file doesn't contain predictor results (the header excluded), exit
if cur == eof:
    sys.exit('ERROR: input predictor file does not contain predictor results!')

for line in f:
    line = line.rstrip()
    print(line)
    arr = line.split("\t")
    print(arr)

    if len(arr) != 2:
        print("Warning: the following line does not have two elements separated by a tab:\n", line)

    key = (arr[0])
    value = float(arr[1])
    predictor_dict[key] = value
    print(predictor_dict)

f.close()

benchmark_dict = {}
f =  open("HGVS_2020_big_benchmark.tsv",'r')
    # Total bytes in the file (end of file)
eof = f.seek(0, 2)
    # Go to the beginning of the file again
f.seek(0)
    # Read the first line (should be the header)
f.readline()
    # Get the current position of the file pointer
cur = f.tell()
    # If the file doesn't contain predictor results (the header excluded), exit
if cur == eof:
    sys.exit('ERROR: input benchmark file does not contain benchmark results!')

for line in f:
    line = line.rstrip()
    arr = line.split("\t")

    if len(arr) < 2:
        print("Warning: the following line does not have three elements separated by a tab:\n", line)

    key = arr[0]
    value = arr[1]
    benchmark_dict[key] = value
    
    
pathogenic = 0
benign = 0
for key, value in predictor_dict.items():
    result = benchmark_dict[key]
    if result == 'Pathogenic':
        pathogenic += 1
    elif result == 'Benign':
        benign += 1
        
score_hgvs_pairs = [(v, k) for k, v in predictor_dict.items()]

sorted_score_hgvs_pairs = score_hgvs_pairs

if type_predictor == "sift":
    sorted_score_hgvs_pairs = sorted(score_hgvs_pairs)
elif type_predictor == "polyphen":
    sorted_score_hgvs_pairs = sorted(score_hgvs_pairs, reverse=True)
elif type_predictor == "BLOSUM":
    sorted_score_hgvs_pairs = sorted(score_hgvs_pairs, reverse=True)
else:
    # Handle unknown predictor types
    raise ValueError(f"Unknown predictor type: {type_predictor}")   

coordinate_score = [sorted_score_hgvs_pairs[0][0]]

tpr = [0.0]
fpr = [0.0]

total_p = pathogenic
total_n = benign





index_prebreakpoint_score = []
previous_score = sorted_score_hgvs_pairs[0][0]
for i in range(len(sorted_score_hgvs_pairs)):
    score = sorted_score_hgvs_pairs[i][0]
    if previous_score != score:
            # Add index of the score before the breakpoint
        index_prebreakpoint_score.append(i - 1)
    previous_score = score

    # Add index of the last score (for the last coordinate)
index_prebreakpoint_score.append(len(sorted_score_hgvs_pairs) - 1)

    # Iterate over HGVS IDs of SNPs and corresponding sorted predictor scores
for i in range(len(sorted_score_hgvs_pairs)):
    score = sorted_score_hgvs_pairs[i][0]
    hgvs = sorted_score_hgvs_pairs[i][1]
    
#########################
### START CODING HERE ###
#########################

num_tp = 0  # Initial count of true positives
num_fp = 0  # Initial count of false positives

# Loop through each score and its associated SNP (HGVS)
for i, (score, hgvs) in enumerate(sorted_score_hgvs_pairs):
    if score < 0:
        Predicted_label = "Pathogenic"
    else:
        Predicted_label = "Benign"
    # Check the benchmark label for the current SNP
    ##for k in 
    actual_label = benchmark_dict[hgvs]
    
    if Predicted_label == "Pathogenic":  # Pathogenic (actual positive)
        if actual_label == "Pathogenic":
            num_tp += 1  # This is a true positive
        else:
            num_fp += 1  # This is a false positive


    # At each unique score breakpoint, calculate and append TPR and FPR
    if i in index_prebreakpoint_score:  
        # Calculate TPR and FPR
        current_tpr = num_tp / total_p if total_p > 0 else 0  # Avoid division by zero
        current_fpr = num_fp / total_n if total_n > 0 else 0  # Avoid division by zero
        print(f"Threshold: {score}, TPR: {current_tpr}, FPR: {current_fpr}")
        # Append the values to the lists
        tpr.append(current_tpr)
        fpr.append(current_fpr)
        coordinate_score.append(score)
        
tpr.append(1.0)
fpr.append(1.0)
        
        
        
   
#########################
###  END CODING HERE  ###
#########################
    
f=open('output_Blosum.tsv','w')
for a,b in zip(fpr, tpr):
    f.write(str(a) + '\t' + str(b) + '\n')
    print(f"Written: {a}\t{b}")

def integrate(fpr, tpr):
    auc = 0
    last_fpr = fpr[0]
    last_tpr = tpr[0]
    
    for cur_fpr, cur_tpr in list(zip(fpr, tpr))[1:]:
            #########################
            ### START CODING HERE ###
            #########################
        delta_fpr = cur_fpr - last_fpr
    
        # Calculate the area of the trapezoid and add to auc
        auc += delta_fpr * (cur_tpr + last_tpr) / 2
    
            #########################
            ###  END CODING HERE  ###
            #########################
        last_fpr = cur_fpr
        last_tpr = cur_tpr
        
    return auc
    

auc = integrate(fpr, tpr)
print(auc)

# Draw ROC plot and write it to a file
lw = 1
figure, axes = matplotlib.pyplot.subplots(1, 1)

color=False

if color:
    lc = colorline(fpr, tpr, coordinator_score, axes=axes)
    color_bar = figure.colorbar(lc)
    colorbar_legend = type_predictor + ' score'
    color_bar.ax.set_ylabel(colorbar_legend)
else:
    axes.plot(fpr, tpr)

axes.plot((0, 1), (0, 1), '--', color='navy', lw=lw, linestyle='--', label='Random')
axes.set_xlim([-0.008, 1.008])
axes.set_ylim([-0.008, 1.008])
axes.set_xlabel('False Positive Rate')
axes.set_ylabel('True Positive Rate')
axes.set_title('AUC = %.3f' % auc)
matplotlib.pyplot.savefig("roc_plot_blosum.png")

