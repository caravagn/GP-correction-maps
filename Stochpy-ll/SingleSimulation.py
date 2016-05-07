import stochpy
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import sys
import getopt

smodfull = stochpy.SSA()

ARGV = sys.argv[1:]

DOPLOTS = True

print 'ARGV      :', sys.argv[1:]

TRANSLATION = None
TRANSCRIPTION = None 

opts = None

NAME = ARGV[0]
TRANSLATION = float(ARGV[1])
if(len(ARGV) == 3):
	TRANSCRIPTION = float(ARGV[2])

## Model
smodfull.Model(NAME+'.psc')
SIM_MODE = 'time'

## Time for simulations
LENGTH = 10**2
SAMPLES = 100 

THRESHOLD_PROTEIN = 200

if(TRANSLATION != None):
	smodfull.ChangeParameter("ktranslation",TRANSLATION)
if(len(ARGV) == 3):
	smodfull.ChangeParameter("ktranscription",TRANSCRIPTION)

PROBAB = 0.0

for s in range(SAMPLES):
	smodfull.DoStochSim(trajectories=1,end=LENGTH,mode=SIM_MODE) 

	times = smodfull.data_stochsim.getTime()

	for r in range(len(times)):
		if(smodfull.data_stochsim.getDataAtTime(times[r])[0][1] > THRESHOLD_PROTEIN):
			PROBAB = PROBAB + 1
			break

print "Outcomes: ", PROBAB

PROBAB = PROBAB / SAMPLES

print "P: ", PROBAB


 
################################################# Traces Plot 
if( (NAME == 'TR-Full') & DOPLOTS ):
	print "Plotting"
	gs = gridspec.GridSpec(3, 1, width_ratios=[1],height_ratios=[0.3,0.7,1])

	ax1 = stochpy.plt.subplot(gs[0])
	smodfull.PlotSpeciesTimeSeries(species2plot = ["Gene_Active"], colors = ["red"])
	stochpy.plt.legend('',frameon=False) # remove legend
	stochpy.plt.xticks([])               # remove x ticks
	stochpy.plt.ylim([0,1.5])            # set y lim
	stochpy.plt.ylabel('Gene')
	stochpy.plt.xlabel('')
	stochpy.plt.yticks([0,1])

	ax2 = stochpy.plt.subplot(gs[1])
	smodfull.plot.ResetPlotnum()
	smodfull.PlotSpeciesTimeSeries(species2plot = ["mRNA" ], colors = ["blue"])
	stochpy.plt.legend('',frameon=False) # remove legend
	stochpy.plt.xticks([])               # remove x ticks
	stochpy.plt.ylabel('mRNA')
	stochpy.plt.xlabel('')

	ax3 = stochpy.plt.subplot(gs[2])
	smodfull.plot.ResetPlotnum()
	smodfull.PlotSpeciesTimeSeries(species2plot = ["Protein" ], colors = ["green"])

	stochpy.plt.legend('',frameon=False) # remove legend
	stochpy.plt.ylabel('Protein')

if( (NAME == 'TR-NoTranscr') & DOPLOTS):
	print "Plotting"
	gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[0.2,0.8])

	ax1 = stochpy.plt.subplot(gs[0])
	smodfull.PlotSpeciesTimeSeries(species2plot = ["Gene_Active"], colors = ["red"])
	stochpy.plt.legend('',frameon=False) # remove legend
	stochpy.plt.xticks([])               # remove x ticks
	stochpy.plt.ylim([0,1.5])            # set y lim
	stochpy.plt.ylabel('Gene')
	stochpy.plt.xlabel('')
	stochpy.plt.yticks([0,1])

	ax2 = stochpy.plt.subplot(gs[1])
	smodfull.plot.ResetPlotnum()
	smodfull.PlotSpeciesTimeSeries(species2plot = ["Protein" ], colors = ["green"])
	stochpy.plt.legend('',frameon=False) # remove legend
	stochpy.plt.ylabel('Protein')


################################################## Pdfs Plot 
if( (NAME == 'TR-Full') & DOPLOTS):
	print "Plotting"
	smodfull.PlotSpeciesDistributions(species2plot = ["mRNA"],  colors = ["blue"])
	smodfull.PlotSpeciesDistributions(species2plot = ["Protein"],  colors = ["green"])

if( (NAME == 'TR-NoTranscr') & DOPLOTS):
	print "Plotting"
	smodfull.PlotSpeciesDistributions(species2plot = ["Protein"],  colors = ["green"])


with open("results/m_PROBAB.csv", "a") as myfile:
 		print "Average logged to : ", str(PROBAB) 
 		myfile.write(str(PROBAB)  + "\n")
