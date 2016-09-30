import shlex, subprocess
import numpy as np
import matplotlib.pyplot as plt

nprocs = 6
ngrain = 5
procs = np.array([2**n for n in range(nprocs)])
grain = np.array([4**n for n in range(ngrain)])
walltime = np.zeros((nprocs,ngrain))
for ip in range(nprocs):
	for ivec in range(ngrain):
		args = shlex.split('mpirun -np {} predprey.out {}'.format(procs[ip],grain[ivec]))
		process = subprocess.Popen(args, stdout=subprocess.PIPE)
		l = process.stdout.read().split()
		walltime[ip,ivec] = float(l[3])
speedup = walltime[0,0]/walltime

plt.imshow(speedup,interpolation='nearest')
plt.xticks(np.arange(ngrain),grain)
plt.yticks(np.arange(nprocs),procs)
plt.xlabel('Grain')
plt.ylabel('Processes')
plt.colorbar()
plt.savefig('speedup.png')