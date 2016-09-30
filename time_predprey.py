import shlex, subprocess
import numpy as np
import matplotlib.pyplot as plt

nprocs = 6
ngrain = 5
procs = 2**np.arange(nprocs)
grain = 4**np.arange(ngrain)
walltime = np.zeros((nprocs,ngrain))
for ip in range(nprocs):
	for ivec in range(ngrain):
		args = shlex.split('mpirun -np %d predprey.out %d'%(procs[ip],grain[ivec]))
		process = subprocess.Popen(args, stdout=subprocess.PIPE)
		l = process.stdout.read().split()
		walltime[ip,ivec] = float(l[3])
		print '%3d, %3d, %f'%(procs[ip], grain[ivec], walltime[ip,ivec])
speedup = walltime[0,0]/walltime

plt.imshow(speedup,interpolation='nearest')
plt.xticks(np.arange(ngrain),grain)
plt.yticks(np.arange(nprocs),procs)
plt.xlabel('Granularity')
plt.ylabel('Processes')
plt.colorbar()
plt.savefig('speedup.png')
