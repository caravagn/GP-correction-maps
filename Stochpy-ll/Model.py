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

THRESHOLD_PROTEIN = 200

## Time for simulations
LENGTH = 10**2
SAMPLES = 100 # Ergodic

## Values for regression (50)
SAMPLED_VALUES = 50 + 1

Valtranslation = np.zeros(SAMPLED_VALUES)
Valtranslation[0] = 0.1 

for x in range(1,SAMPLED_VALUES):
	Valtranslation[x] = x * 2

print "Translation values (Theta):"
print Valtranslation

# Creates a list containing SAMPLED_VALUES lists, each of 8 items, all set to 0
buckets_PROBAB = np.zeros( (SAMPLED_VALUES,SAMPLED_VALUES_FPAR) )
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

	PROBAB = 0.0

	for s in range(SAMPLES):
		smodfull.DoStochSim(trajectories=1,end=LENGTH,mode='time') 

		times = smodfull.data_stochsim.getTime()

		for r in range(len(times)):
			if(smodfull.data_stochsim.getDataAtTime(times[r])[0][1] > THRESHOLD_PROTEIN):
				PROBAB = PROBAB + 1
				break

	# print "Outcomes: ", PROBAB

	PROBAB = PROBAB / SAMPLES

	# print "P: ", PROBAB

	elapsed = timeit.default_timer() - start_time
	print "@ M: input ", x, " in ", (elapsed/1000), " sec(s); PROBAB ", PROBAB

	return {'x':x, 'time':elapsed, 'PROBAB':PROBAB}


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

	PROBAB = [result['PROBAB'] for result in results]
	times = [result['time'] for result in results]

	for z in range(SAMPLED_VALUES):
		buckets_PROBAB[z,idx] = PROBAB[z]
		buckets_times[z,idx] = times[z]

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
np.savetxt("results/PROBAB.csv", buckets_PROBAB, delimiter=",")
np.savetxt("results/timings.csv", buckets_times, delimiter=",")

