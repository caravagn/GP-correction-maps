import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib 
import numpy as np
import sys
import math

DOPLOTS = True
LOGTRANSFORM = True

ARGV = sys.argv[1:]

NAME = ARGV[0]

from numpy import genfromtxt
MData = genfromtxt('results/'+NAME, delimiter=',')
mData = genfromtxt('results/m_'+NAME, delimiter=',')


print "M", MData
print "m", mData

NumOfFreeVars = MData.shape[1]

RP_input = genfromtxt('results/input.csv', delimiter=',')

# Matlab fonts
csfont = {'fontname':'Arial'}
hfont = {'fontname':'Verdana'}

font = {'family' : 'Arial',
        #'weight' : 'bold',
        'size'   : 9}

matplotlib.rc('font', **font)


if(DOPLOTS):	
	# Two subplots, the axes array is 1-d
	f, axarr = plt.subplots(4, sharex=True)

	for a in range(NumOfFreeVars):
		axarr[0].scatter(RP_input, MData[:,a], c='g', alpha=0.3)
		
	axarr[0].scatter(RP_input, np.mean(MData, axis=1), c='k', alpha=0.9, label = 'mean')
	# axarr[0].set_yscale('log')
	# axarr[0].ticklabel_format(style = 'sci')


	# axarr[0].set_xticks([])         
	axarr[0].set_ylabel('Data generate\nfrom model M', **csfont)
	axarr[0].set_xlabel('')

	axarr[0].set_title('Average number of proteins', **csfont)
	
	axarr[0].legend(loc='upper left', prop={'size':9}, shadow=True, fancybox=True, scatterpoints = 1)

	axarr[1].scatter(RP_input, mData, c='r', alpha=0.5)	
	axarr[1].set_ylabel('Data generated\nfrom model m', **csfont)
	axarr[1].set_xlabel('')
	axarr[1].set_yticks(range(1,int(max(mData)),int(max(mData))/4))

## The input set for the training set (translation)
TRAINING_POINTS_INPUT = RP_input

## Correct each point
for a in range(NumOfFreeVars):
	MData[:,a] = MData[:,a] - mData

## Log-transformation of data
if(LOGTRANSFORM): 
  	MData = np.log(MData)

## As training point we pick the expected value
TRAINING_POINTS = np.mean(MData, axis=1) 

## The variance in these observations is hence
TRAINING_POINTS_VARIANCE = np.var(MData, axis=1) 

if(DOPLOTS):	
	axarr[2].scatter(RP_input, TRAINING_POINTS, alpha=0.5)
	axarr[2].set_yticks(range(1,int(max(TRAINING_POINTS)),int(max(TRAINING_POINTS))/4))
	
	if(LOGTRANSFORM): 
		axarr[2].set_ylabel('Correction\n(logscale)', **hfont)
	if(LOGTRANSFORM == False): 
		axarr[2].set_ylabel('Correction', **hfont)

# fit with np.polyfit
m, b = np.polyfit(RP_input, TRAINING_POINTS, 1)
axarr[2].plot(RP_input, m*RP_input + b, '-')

# Variance plot
if(DOPLOTS):	
	axarr[3].scatter(RP_input, TRAINING_POINTS_VARIANCE, c="#a020f0", alpha=0.5)
	axarr[3].set_xticks(range(0,110,10))
	axarr[3].set_xlabel('Protein translation rate', **hfont)
	axarr[3].set_ylabel('Correction\nVariance', **hfont)

if(LOGTRANSFORM == False): 
	axarr[3].set_yscale('log')

if(DOPLOTS):	
	f.show()

np.savetxt("results/ML_TRAINING_POINTS_INPUT_"+NAME, TRAINING_POINTS_INPUT, delimiter=",")
np.savetxt("results/ML_TRAINING_POINTS_"+NAME, TRAINING_POINTS, delimiter=",")
np.savetxt("results/ML_TRAINING_POINTS_VARIANCE_"+NAME, TRAINING_POINTS_VARIANCE, delimiter=",")
