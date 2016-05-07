import stochpy as stochpy
import stochpy as stochpy_m

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

smodfull = stochpy.SSA()
smodred = stochpy_m.SSA()

DOPLOTS = False

## Model
NAME = 'TR-Full'
NAME_RED = 'TR-NoTranscr'

## transcription values -- when alpha \in F
import random
SAMPLED_VALUES_FPAR = 50
ALPHA =  random.sample(range(1, 100), SAMPLED_VALUES_FPAR)

print "Transcription values (F):"
print ALPHA

## Time for simulations
LENGTH = 10**6
SAMPLES = 1 # Ergodic

## Values for regression - [0,0.1,1,3,...,100]
SAMPLED_VALUES = 50 + 1

Valtranslation = np.zeros(SAMPLED_VALUES)
Valtranslation[0] = 0.1
# Valtranslation[1] = 0.1
# Valtranslation[2] = 1

for x in range(1,SAMPLED_VALUES):
	Valtranslation[x] = x * 2

print "Translation values (Theta):"
print Valtranslation

# Creates a list containing SAMPLED_VALUES lists, each of 8 items, all set to 0
buckets_avg = np.zeros( (SAMPLED_VALUES,SAMPLED_VALUES_FPAR) )
buckets_var = np.zeros( (SAMPLED_VALUES,SAMPLED_VALUES_FPAR) )
buckets_times = np.zeros( (SAMPLED_VALUES,SAMPLED_VALUES_FPAR) )

###
print "Starting multiprocessing"
import multiprocessing as mp
import timeit
import multiprocessing as mp_m

# One parameter set, one computation
def foo_pool_M(x):
	start_time = timeit.default_timer()

	## Set beta to x, simulate M
	smodfull.ChangeParameter("ktranslation", x)
	smodfull.SetQuiet()
	smodfull.DoStochSim(trajectories=SAMPLES,end=LENGTH,mode='steps') 

	## First and Second moments of M
	protein_avg_smodfull = smodfull.data_stochsim.species_means["Protein"]
	protein_stdev_smodfull = smodfull.data_stochsim.species_standard_deviations["Protein"]
	protein_var_smodfull = protein_stdev_smodfull * protein_stdev_smodfull

	elapsed = timeit.default_timer() - start_time
	print "@ M: input ", x, " in ", (elapsed/1000), " sec(s): ", protein_avg_smodfull, "/", protein_var_smodfull

	return {'x':x, 'time':elapsed, 'avg':protein_avg_smodfull, 'var':protein_var_smodfull}

def foo_pool_m(x):
	start_time = timeit.default_timer()

	## Set beta to x, simulate m
	smodred.ChangeParameter("ktranslation",x)
	smodred.DoStochSim(trajectories=SAMPLES,end=LENGTH,mode='steps') 
	
	smodred.ShowOverview()
	smodred.PlotSpeciesTimeSeries()
	smodred.PrintSpeciesMeans()

	# print smodred.data_stochsim.species_means
	print smodred.data_stochsim.species_means

	protein_avg_smodred = smodred.data_stochsim.species_means["Protein"]
	protein_var_smodred = 1

	elapsed = timeit.default_timer() - start_time
	print "@ m: input ", x, " in ", (elapsed/1000), " sec(s): ", protein_avg_smodred, "/", protein_var_smodred

	return {'x':x, 'time':elapsed, 'avg':protein_avg_smodred, 'var':protein_var_smodred}


def apply_async_with_callback_M(ncpu, idx):
	
	alpha = ALPHA[idx]

	# Create a model
	smodfull.Model(NAME+'.psc')

	# Set alpha -- M only
	smodfull.ChangeParameter("ktranscription", alpha)

	# # Pool of parallel processes
	pool = mp.Pool(processes = ncpu)
	
	print "Pooling: ", ncpu
	results = [pool.apply_async(foo_pool_M, args = (Valtranslation[x], )) for x in range(SAMPLED_VALUES)]
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
		
def apply_async_with_callback_m(ncpu):
	
	# Create a model
	smodred.Model(NAME_RED+'.psc')

	# Pool of parallel processes
	pool_m = mp_m.Pool(processes = ncpu)
	
	print "Pooling: ", ncpu
	results_m = [pool_m.apply_async(foo_pool_m, args = (Valtranslation[x], )) for x in range(SAMPLED_VALUES)]
	print "Pooling: closing/join"

	pool_m.close()
	pool_m.join()

	print results_m
	results = [x.get() for x in results_m]
	print results

	avg = [result['avg'] for result in results]
	var = [result['var'] for result in results]
	times = [result['time'] for result in results]

	for z in range(SAMPLED_VALUES):
		for w in range(SAMPLED_VALUES_FPAR):
			buckets_avg[z,w] = buckets_avg[z,w] - avg[z]
			buckets_var[z,w] = buckets_var[z,w] - var[z]
			buckets_times[z,w] = buckets_times[z,w] + times[z]

	return results		

print "******* M"
# Pool of parallel processes
for a in range(SAMPLED_VALUES_FPAR):
	print "********\n******** Iteration ",(a + 1),"Alpha = ", ALPHA[a]
	results = apply_async_with_callback_M(25, a)


print "******* m"

with open("script.sh", "a") as myfile:
	for z in range(SAMPLED_VALUES):
		command = "python SingleSimulation.py TR-NoTranscr " + str(Valtranslation[z])
		print "Command: ", command
		myfile.write(command + "\n")


if(DOPLOTS):
	for a in range(SAMPLED_VALUES_FPAR):
		plt.scatter(Valtranslation, buckets_avg[:, a])
	plt.figure()
	for a in range(SAMPLED_VALUES_FPAR):
		plt.scatter(Valtranslation, buckets_var[:, a])


# Valtranslation
np.savetxt("results/input.csv", Valtranslation, delimiter=",")
np.savetxt("results/averages.csv", buckets_avg, delimiter=",")
np.savetxt("results/variances.csv", buckets_avg, delimiter=",")
np.savetxt("results/timings.csv", buckets_times, delimiter=",")


# ##################################################  Ergodicity test
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
