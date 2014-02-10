#!/usr/bin/env python

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from subprocess import call
from multiprocessing import Process

beta = 1./(300*8.314e-3)

def plot_graph(*args):
	plt.figure(facecolor='white')
	for data in args[0]:
		plt.semilogx(data[0],data[1])
	plt.title(r"$%s$"%args[1],fontsize=18)
	plt.xlabel(r"$Time$",fontsize=18)
	plt.ylabel(r"$P(t)$",fontsize=18)
	plt.show()

class FourState:
	""" numerical example of a four state model"""

	def __init__(self):
		self.initialize()
		self.calc_trans()
		self.calc_rate()
		self.calc_eigs()

	def initialize(self):
		""" initialize count matrix """
		# initialize randomly
		randint = np.random.randint(1000,size=(4,4)) 
		# disconnect end states
		randint[0,3] = 0 ; randint[3,0] = 0
		# enforce detailed balance
		randint = np.triu(randint)
		for i in range(1,4):
			for j in range(i):
				randint[i,j] = randint[j,i]
		# enforce metastability of states
		for i in range(4): 
			randint[i,i] = randint[i,i]**2
		diagsort = np.sort(np.diag(randint)) # enforce FF and UU more stable than I1 and I2
		randint[0,0] = diagsort[2]
		randint[1,1] = diagsort[0]
		randint[2,2] = diagsort[1]
		randint[3,3] = diagsort[3]

		self.count = randint
		print "\n count matrix:"
		print self.count

	def calc_trans(self):
		""" calculate transition matrix (assume lag time is 1)"""
		count = self.count
		trans = np.zeros((4,4),float)
		for i in range(4):
			nsum = np.sum(count[:,i])
			for j in range(4):
				trans[j,i] = float(self.count[j,i])/nsum
		self.trans = trans
		print "\n transition matrix:"
		print self.trans
	
	def calc_rate(self):
		""" calculate rates """
		rate = self.trans
		for i in range(4):
			rate[i,i] = 0.
			rate[i,i] = -np.sum(rate[:,i])
		self.rate = rate
		print "\n rate:"
		print self.rate

	def calc_eigs(self):
		""" calculate eigenvectors """
		evals,rvecs = np.linalg.eig(self.rate)
		lvecs = np.linalg.inv(rvecs)
		# sort eigenvectors
		order = np.argsort(-evals)
		evals_sort = evals[order]
		rvecs_sort = np.zeros((4,4),float)
		for i in range(4):
			rvecs_sort[:,i] = rvecs[:,order[i]]
		lvecs_sort = np.linalg.inv(rvecs_sort)
		self.peq = rvecs_sort[:,0]/np.sum(rvecs_sort[:,0])
		self.evals = evals_sort
		self.rvecs = rvecs_sort
		self.lvecs = lvecs_sort
		print "\n eigenvalues:"
		print self.evals

	def run_simul(self,label=None):
		""" simulate a relaxation from unfolded to folded"""
		p0 = [1.,0.,0.,0.]
		pt = []
		logt = np.arange(0,6,0.1)
		time = 10**logt
		for t in time:
			expdiagkt = np.diag(np.exp(self.evals*t))
			expKt = np.dot(self.rvecs,np.dot(expdiagkt,self.lvecs))
			pt.append(np.dot(expKt,p0))
		data = []
		for i in range(4):
			data.append([time,map(lambda x: x[i],pt)])
		p = Process(target=plot_graph, args=[data,label])
		p.start()

	def run_commit(self):
		""" calculate committors and reactive flux """
		K = self.rate
		peq = self.peq
		# define end-states
		UU = [0]
		FF = [3]
		UUFF = UU+FF
		I = filter(lambda x: x not in UU+FF, range(4))
		NI = len(I)

		# calculate committors
		b = np.zeros([NI], float)
		A = np.zeros([NI,NI], float)
		for j_ind in range(NI):
			j = I[j_ind]
			sum = 0.
			for i in FF:
				sum+= K[i][j]
			b[j_ind] = -sum
			for i_ind in range(NI):
				i = I[i_ind]
				A[j_ind][i_ind] = K[i][j]		
		# solve Ax=b
		Ainv = np.linalg.inv(A)
		x = np.dot(Ainv,b)
		XX = np.dot(Ainv,A)

		pfold = np.zeros(4,float)
		for i in range(4):
			if i in UU:
				pfold[i] = 0.0
			elif i in FF:
				pfold[i] = 1.0
			else:
				ii = I.index(i)
				pfold[i] = x[ii]
		self.pfold = pfold
		print "\n pfold values:"
		print pfold
							
		# stationary distribution
		pss = np.zeros([4],float)
		for i in range(4):
			pss[i] = (1-pfold[i])*peq[i]
		self.pss = pss

		# flux matrix and reactive flux through pfold=0.5
		psepar = 0.5
		J = np.zeros([4,4],float)
		sum_flux = 0
		for i in range(4):
		    for j in range(4):
				J[j][i] = K[j][i]*peq[i]*(pfold[j]-pfold[i])
				if pfold[j] > pfold[i]:
					sum_flux += J[j][i]

		print "\n flux :"
		print J
		print "\n reactive flux: %g"%sum_flux
		self.J = J
				    
	def do_dot(self,out=None):
		""" generate network diagram using dot """
		if out==None:
			out="Flux.dot"
		fout = open(out,"w")
		JD = nx.DiGraph(self.J)
		fluxmin = np.min(np.abs(self.J[self.J.nonzero()]))
		pfold = self.pfold
		peq = self.peq
		fout.write ("strict digraph G {\n")
		peq_min = np.min(peq)
		for n in JD.nodes():
			fout.write("%i [shape=circle,width=%f];\n"%(n,peq[n]))
		fout.write("{rank=same; %i; %i;}\n"%(1,2))
		d_sum = 0
		fluxes = []
		for (u,v,d) in JD.edges(data=True):
			if d['weight'] > 0:
				fluxes.append(d['weight'])
		flux_min = np.min(fluxes)
		flux_max = np.max(fluxes)
		dflux = flux_max - flux_min
		for (u,v,d) in JD.edges(data=True):	
			if d['weight'] > 0:
				fout.write( "%i -> %i  [penwidth=%f];\n"%(v,u,9*(d['weight']-flux_min)/dflux+1))
		fout.write("}")
		fout.close()

	def mutate(self,elem,dg=None):
		""" mutate model by given amount """
		if dg is None:
			dg = 1.
		i = elem
		rate = self.rate
		for j in [x for x in range(4) if x!= i]:
			rate[j,i] = rate[j,i]*np.exp(beta*dg/2)
			rate[i,j] = rate[i,j]*np.exp(-beta*dg/2)
		for j in range(4):
			rate[j,j] = 0.
			rate[j,j] = -np.sum(rate[:,j])
		self.rate = rate

		print "\n mutated residue: %g\n change in free energy : %g kT\n"%(elem,dg*beta)
		print " mutated rate matrix"
		print self.rate

	def partial_rate(self,elem,dg):
		""" calculate derivative of rate matrix """
		rate = self.rate
		d_rate = np.zeros((4,4),float)
		for i in range(4):
			if i != elem:
				d_rate[i,elem] = beta/2.*rate[i,elem]*np.exp(beta*dg/2.);
				d_rate[elem,i] = -beta/2.*rate[elem,i]*np.exp(beta*dg/2.);
		for i in range(4):
			d_rate[i,i] = -np.sum(d_rate[:,i]))
		self.d_rate = d_rate

	def partial_peq(self,elem,dg):
		""" calculate derivative of equilibrium distribution """
		for i in range(4):
			if i != elem:
				d_peq.append(beta*self.peq[i]*self.peq[elem])
			else:
				d_peq.append(-beta*self.peq[i]*(1.-self.peq[i])
		self.d_peq = d_peq

	def partial_pfold(self,elem,dg):
		""" calculate derivative of pfold """
		K = self.rate
		peq = self.peq
		# define end-states
		UU = [0]
		FF = [3]
		UUFF = UU+FF
		I = filter(lambda x: x not in UU+FF, range(4))
		NI = len(I)
		# calculate committors
		b = np.zeros([NI], float)
		A = np.zeros([NI,NI], float)
		db = np.zeros([NI], float)
		dA = np.zeros([NI,NI], float)
		for j_ind in range(NI):
			j = I[j_ind]
			sum = 0.
			sumd = 0.
			for i in FF:
				sum+= K[i][j]
				sumd+= self.d_rate[i][j]
			b[j_ind] = -sum
			db[j_ind] = -sumd
			for i_ind in range(NI):
				i = I[i_ind]
				A[j_ind][i_ind] = K[i][j]
				dA[j_ind][i_ind] = self.d_rate[i][j]

		# solve Ax + Bd(x) = c
		Ainv = np.linalg.inv(A)
		pfold = np.dot(Ainv,b)
		dpfold = np.dot(Ainv,db - dot(dA,pfold))
		self.pfold = dpfold
