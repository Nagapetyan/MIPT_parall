from time import time
import subprocess
import numpy as np

subprocess.call(["make"])

times = list()
for j in xrange(1,5):
	time_j = list()
	for _ in xrange(20):
		start = time()
		subprocess.check_output([ "mpirun", "-n", "%f"  % j,"./task1" , "1000000"] )  #1000000
		time_j.append(time() - start)
	times.append(np.mean(time_j))
print times