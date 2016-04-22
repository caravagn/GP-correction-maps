import stochpy
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

smodfull = stochpy.SSA()

DOPLOTS = True

## Model
NAME = 'TR-Full'
smodfull.Model(NAME+'.psc')
SIM_MODE = 'time'

## Time for simulations
LENGTH = 10**3
SAMPLES = 1 # Ergodic

## Values for regression - [0,0.1,1,3,...,100]
TRANSLATION = 100 
TRANSCRIPTION = 10

smodfull.ChangeParameter("ktranslation",TRANSLATION)
smodfull.ChangeParameter("ktranslation",TRANSLATION)


smodfull.DoStochSim(trajectories=SAMPLES,end=LENGTH,mode=SIM_MODE) 

protein_avg_smodfull = smodfull.data_stochsim.species_means["Protein"]
protein_stdev_smodfull = smodfull.data_stochsim.species_standard_deviations["Protein"]
protein_var_smodfull = protein_stdev_smodfull * protein_stdev_smodfull

smodfull.ShowOverview()

smodfull.PrintAverageSpeciesTimeSeries()

################################################# Traces Plot 
gs = gridspec.GridSpec(3, 1, width_ratios=[1],height_ratios=[0.3,0.7,1])

ax1 = stochpy.plt.subplot(gs[0])
smodfull.PlotSpeciesTimeSeries(species2plot = ["Gene_Active"], colors = ["red"])
stochpy.plt.legend('',frameon=False) # remove legend
stochpy.plt.xticks([])               # remove x ticks
stochpy.plt.ylim([0,1.5])            # set y lim
#stochpy.plt.xlim([0,TIME])            # set y lim
# stochpy.plt.text(+80,0.9,'ON')      
# stochpy.plt.text(TIME+80,0,'OFF')
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

################################################## Pdfs Plot 
smodfull.PlotSpeciesDistributions(species2plot = ["Protein"],  colors = ["green"])
smodfull.PlotSpeciesDistributions(species2plot = ["mRNA"],  colors = ["blue"])


## Ergodicity test
## Simulate up to estiamtion of convergence in the first 4 moments of the empirical distributions
# smod.DoCompleteStochSim()

# protein_avg = smod.data_stochsim.species_means["Protein"]
# protein_stdev = smod.data_stochsim.species_standard_deviations["Protein"]
# protein_var = protein_stdev * protein_stdev

# print "Protein: 1st moment %s." % protein_avg
# print "Protein: 2nd moment %s." % protein_var

# Simulation for TIME = 10^6 and this procedure give consistent results
# Protein: 1st moment 34.2858754344.
# Protein: 2nd moment 186.146132686.
# Protein: 1st moment 34.2439045426.
# Protein: 2nd moment 185.431079815.




#smod.PlotSpeciesDistributions(species2plot = ["mRNA"])

#smod.PlotSpeciesDistributions(species2plot = ["mRNA"])

#smod.PlotSpeciesDistributions(species2plot = ["Gene_Active"])

#smod.ChangeInitialSpeciesCopyNumber("Protein",1000)
#smod.ChangeParameter("kunbind",0.01)


#smod.GetWaitingtimes()
#smod.PlotWaitingtimesDistributions()

#smod.GetRegularGrid(n_samples=50)

#smod.PlotAverageSpeciesTimeSeries()

# np.savetxt("averages.csv", buckets_avg, delimiter=",")
# np.savetxt("variances.csv", buckets_avg, delimiter=",")
# np.savetxt("timings.csv", buckets_times, delimiter=",")

