from time import time
import subprocess
from numpy import mean
from numpy.random import uniform
from tqdm import tqdm 


check = subprocess.check_output(["bash", "exist.sh", 'dataForTask2.txt'])
if check == 1:
	check = subprocess.check_output(["rm", "dataForTask2.txt"])

file = open("dataForTask2.txt", 'w+')

for _ in xrange(2):
	number = uniform(0,2,2000)
	number = number.astype(int)

	for j in xrange(2000):
		file.write(str(number[j]))
	file.write('\n')


subprocess.call(["make"])
check = subprocess.check_output(["bash", "exist.sh", "answer.txt"])
if check == 1:
	check = subprocess.check_output(["rm", "answer.txt"])


times = list()
for j in xrange(1,5):
	time_j = list()
	pbar = tqdm(total=100)
	for _ in xrange(100):
		pbar.update(1)
		start = time()
		subprocess.call([ "mpirun", "-n", "%d"  % j,"./task2" ] )  
		subprocess.call(['rm', 'answer.txt'])
		time_j.append(time() - start)
	times.append(mean(time_j))
	pbar.close()
print times