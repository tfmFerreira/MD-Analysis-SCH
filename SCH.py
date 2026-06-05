#!/usr/bin/env python
# Ferreira (2017)
from __future__ import print_function
import time
start = time.time()

import mdtraj as md
import math
import sys
import numpy as np

print('     ')
print('loading trajectory')

traj = md.load('traj.gro')

time1 = time.time()
totaltime = time1 - start
print('(time spent loading: %d s)' % totaltime)
print('     ')

tstep=traj.time[1]-traj.time[0]
print('simulation trajectory of %s ps with steps of %s ps' % (traj.time[-1],tstep))
print('     ')
corr=[0]*traj.n_frames
lag=traj.n_frames/2

topology = traj.topology

with open('inputSCH') as f:
	CHlabel = [line.rstrip('\n') for line in open('inputSCH')]

print('labels found in file input1:')
print(CHlabel)
print('     ')


totalH=len(CHlabel)
#print(totalH)

SCHFILE2=open('SCH.dat','w+')
nH=0
#while nH < totalH:
for line in CHlabel:
	CHpair = line.split()
	#print(CHpair)
	Ccalc=topology.select("name "+str(CHpair[0]))
	Hsurround=topology.select("name "+str(CHpair[1]))
	#print(Hcalc)
	#print(Hsurround)
	totalsurround=len(Hsurround)
	#print(totalsurround)
	nsurround=0
	OPCH=0

	while nsurround < totalsurround:

		atom2=topology.atom(Hsurround[nsurround])
		atom1=topology.atom(Ccalc[nsurround])
		start = time.time()
		d3=np.array(traj.xyz[:,Ccalc[nsurround]])-np.array(traj.xyz[:,Hsurround[nsurround]])
		d3x=d3[:,0]
		d3y=d3[:,1]
		d3z=d3[:,2]
		d=np.sqrt(d3x*d3x+d3y*d3y+d3z*d3z)
		internucvector=np.array([d3x/d, d3y/d, d3z/d])
		internucvector=internucvector.conj().transpose()
		#print(d)
		#SCHFILE=open('internuclearvector'+str(nsurround)+'.dat','w+')
		j=0
		while j<traj.n_frames:
			dx=d3x[j]
			dy=d3y[j]
			dz=d3z[j]
			#print(traj.time[j],dx,dy,dz,file=SCHFILE)
			OPCH=OPCH+0.5*(3*dz**2/d[j]**2-1)/totalsurround/traj.n_frames
			j=j+1
		#print('loop j: %s' % j)
		end = time.time()
		totaltime = end - start
		#print('total time of calculation: %d seconds' % totaltime)
		#print('     ')
		nsurround=nsurround+1
	print(CHpair[0],CHpair[1],OPCH,file=SCHFILE2)
	print(CHpair[0],CHpair[1],OPCH)
print('     ')
print('C-H bond order parameters saved!')
print('     ')
