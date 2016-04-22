import stochpy
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

smodfull = stochpy.SSA()

DOPLOTS = True

## Model
NAME = 'TR-Full'

## transcription values -- when alpha \in F
import random
SAMPLED_VALUES_FPAR = 10
ALPHA =  random.sample(range(1, 100), SAMPLED_VALUES_FPAR)

print "Transcription values (F):"
print ALPHA


## Time for simulations
LENGTH = 10**6
SAMPLES = 1 # Ergodic

## Values for regression - [0,0.1,1,3,...,100]
SAMPLED_VALUES = 50 

Valtranslation = [0] * SAMPLED_VALUES
Valtranslation[1] = 0.1
Valtranslation[2] = 1

for x in range(3, SAMPLED_VALUES):
	Valtranslation[x] = Valtranslation[x - 1] + 2

print "Translation values (Theta):"
print Valtranslation


# This works, but let's do it parallel
# buckets_avg = [0] * SAMPLED_VALUES
# buckets_var = [0] * SAMPLED_VALUES
# for x in range(1, SAMPLED_VALUES):
# 	smodfull.ChangeParameter("ktranslation",Valtranslation[x])

# 	smodfull.DoStochSim(trajectories=SAMPLES,end=TIME,mode='time') 

# 	################################################## First and Second moments
# 	protein_avg_smodfull = smodfull.data_stochsim.species_means["Protein"]
# 	protein_stdev_smodfull = smodfull.data_stochsim.species_standard_deviations["Protein"]
# 	protein_var_smodfull = protein_stdev_smodfull * protein_stdev_smodfull

# 	print "Full Model, Protein: translation %s." % Valtranslation[x]
# 	print "Full Model, Protein: 1st moment %s." % protein_avg_smodfull
# 	print "Full Model, Protein: 2nd moment %s." % protein_var_smodfull

# 	buckets_avg[x] = protein_avg_smodfull	
# 	buckets_var[x] = protein_var_smodfull


# Creates a list containing SAMPLED_VALUES lists, each of 8 items, all set to 0
buckets_avg = np.zeros( (SAMPLED_VALUES,SAMPLED_VALUES_FPAR) )
buckets_var = np.zeros( (SAMPLED_VALUES,SAMPLED_VALUES_FPAR) )
buckets_times = np.zeros( (SAMPLED_VALUES,SAMPLED_VALUES_FPAR) )

###
print "Starting multiprocessing"
import multiprocessing as mp
import timeit

# One parameter set, one computation
def foo_pool(x):
	start_time = timeit.default_timer()

	## Set beta to x
	smodfull.ChangeParameter("ktranslation",x)
	smodfull.DoStochSim(trajectories=SAMPLES,end=LENGTH,mode='steps') 
	#smodfull.ShowOverview()

	## First and Second moments
	protein_avg_smodfull = smodfull.data_stochsim.species_means["Protein"]
	protein_stdev_smodfull = smodfull.data_stochsim.species_standard_deviations["Protein"]
	protein_var_smodfull = protein_stdev_smodfull * protein_stdev_smodfull

	elapsed = timeit.default_timer() - start_time
	print "Input ", x, " in ", (elapsed/1000), " sec(s): ", protein_avg_smodfull, "/", protein_var_smodfull

	return {'x':x, 'time':elapsed, 'avg':protein_avg_smodfull, 'var':protein_var_smodfull}


# buckets_avg = []
# buckets_var = []

# def log_result(result):
#     # This is called whenever foo_pool(i) returns a result.
#     # result_list is modified only by the main process, not the pool workers.
#     buckets_avg.append(result['avg'])
#     buckets_var.append(result['var'])



def apply_async_with_callback(ncpu, idx):
	
	alpha = ALPHA[idx]

	# Create a model
	smodfull.Model(NAME+'.psc')

	# Set alpha
	smodfull.ChangeParameter("ktranscription", alpha)

	# Pool of parallel processes
	pool = mp.Pool(processes = ncpu)
	
	print "Pooling: ", ncpu
	results = [pool.apply_async(foo_pool, args = (Valtranslation[x], )) for x in range(SAMPLED_VALUES)]
	print "Pooling: closing/join"

	pool.close()
	pool.join()
	results = [result.get() for result in results]

	avg = [result['avg'] for result in results]
	var = [result['var'] for result in results]
	times = [result['time'] for result in results]

	for z in range(SAMPLED_VALUES):
		buckets_avg[z,idx] = avg[z]
		buckets_var[z,idx] = var[z]
		buckets_times[z,idx] = times[z]

	return results
		

for a in range(SAMPLED_VALUES_FPAR):
	print "********\n******** Iteration ",(a + 1),"Alpha = ", ALPHA[a]
	results = apply_async_with_callback(25, a)
	# print(buckets_avg)

smodfull.ShowOverview()

if(DOPLOTS):
	for a in range(SAMPLED_VALUES_FPAR):
		plt.scatter(Valtranslation, buckets_avg[:, a])
	plt.figure()
	for a in range(SAMPLED_VALUES_FPAR):
		plt.scatter(Valtranslation, buckets_var[:, a])


np.savetxt("averages.csv", buckets_avg, delimiter=",")
np.savetxt("variances.csv", buckets_avg, delimiter=",")
np.savetxt("timings.csv", buckets_times, delimiter=",")

# 	plt.scatter(Valtranslation, avg)
# 	plt.figure()
# 	plt.scatter(Valtranslation, var)
# 	plt.show()



#output = mp.Queue()
#
#print "Pooling: %s" % SAMPLED_VALUES
#processes = [mp.Process(target = foo_pool, args = (Valtranslation[x],)) for x in range(SAMPLED_VALUES)]
#
#for p in processes:
#	p.start()
#
#for p in processes:
#	p.join()
#
#results = [output.get() for p in processes]
#
#buckets_avg = [result['avg'] for result in results]
#buckets_var = [result['var'] for result in results]



	# import csv

	# with open('data-FM-' + alpha +'.csv', 'w') as csvfile:
 #    	fieldnames = ['input', 'time', 'avg', 'var']
 #    	writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

	#     writer.writeheader()
 #   		writer.writerow({'input': alpha, 'last_name': 'Beans'})
 #   		writer.writerow({'first_name': 'Lovely', 'last_name': 'Spam'})
 #    	writer.writerow({'first_name': 'Wonderful', 'last_name': 'Spam'})


 
# import csv

# myfile = open("averages_"+NAME+'.csv', 'wb')
# wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
# wr.writerow(buckets_avg)


################################################## Traces Plot 
# gs = gridspec.GridSpec(3, 1, width_ratios=[1],height_ratios=[0.3,0.7,1])

# ax1 = stochpy.plt.subplot(gs[0])
# smodfull.PlotSpeciesTimeSeries(species2plot = ["Gene_Active"], colors = ["red"])
# stochpy.plt.legend('',frameon=False) # remove legend
# stochpy.plt.xticks([])               # remove x ticks
# stochpy.plt.ylim([0,1.5])            # set y lim
# #stochpy.plt.xlim([0,TIME])            # set y lim
# stochpy.plt.text(TIME+80,0.9,'ON')      
# stochpy.plt.text(TIME+80,0,'OFF')
# stochpy.plt.ylabel('Gene')
# stochpy.plt.xlabel('')
# stochpy.plt.yticks([0,1])

# ax2 = stochpy.plt.subplot(gs[1])
# smodfull.plot.ResetPlotnum()
# smodfull.PlotSpeciesTimeSeries(species2plot = ["mRNA" ], colors = ["blue"])
# stochpy.plt.legend('',frameon=False) # remove legend
# stochpy.plt.xticks([])               # remove x ticks
# stochpy.plt.ylabel('mRNA')
# stochpy.plt.xlabel('')

# ax3 = stochpy.plt.subplot(gs[2])
# smodfull.plot.ResetPlotnum()
# smodfull.PlotSpeciesTimeSeries(species2plot = ["Protein" ], colors = ["green"])
# stochpy.plt.legend('',frameon=False) # remove legend
# stochpy.plt.ylabel('Protein')

# ################################################## Pdfs Plot 
# smodfull.PlotSpeciesDistributions(species2plot = ["Protein"],  colors = ["green"])
# smodfull.PlotSpeciesDistributions(species2plot = ["mRNA"],  colors = ["blue"])


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
#smod.PrintAverageSpeciesTimeSeries()
