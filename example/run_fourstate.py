#!/usr/bin/env python

import sys,os,copy
import fourstate

print "\n WT MARKOV MODEL"
#bhs = fourstate.FourState("random")
bhs = fourstate.FourState()
print bhs.count
#bhs.run_simul("WT")
bhs.run_commit()
bhs.fold_rate()
#bhs.do_dot()

#dg = 0.001
#print "\n\n MUTANT MARKOV MODEL"
#mut = (copy.deepcopy(bhs)) #.mutate(1)
#mut.mutate(2,dg)
#mut.calc_eigs()
##mut.run_simul("Mutant")
#mut.run_commit()
#mut.fold_rate()
##mut.do_dot("Mut.dot")
#
#print "\n WT MARKOV MODEL: PERTURBATION"
#sum_d_flux,d_peq,d_pfold,d_J,d_K = bhs.partial_flux(2,dg)
#
#
#print "\n State      Peq     Peq_mut    Peq_mut-Peq    delta_Peq"
#for i in range(4):
#	print "%4i %12.4f %12.4e %12.4e %12.4e"%(i+1,bhs.peq[i],mut.peq[i],mut.peq[i]-bhs.peq[i],d_peq[i]*dg)
#
#print "\n State      phi     phi_mut   phi_mut-phi    delta_phi"
#for i in range(4):
#	print "%4i %12.4f %12.4e %12.4e %12.4e"%(i+1,bhs.pfold[i],mut.pfold[i],mut.pfold[i]-bhs.pfold[i],d_pfold[i]*dg)
#
#print "\n State    State      phi     phi_mut   phi_mut-phi    delta_phi"
#for i in range(4):
#	for j in range(4):
#		print "%4i %4i %12.4f %12.4e %12.4e %12.4e"%(i,j,bhs.K[i,j],mut.K[i,j],mut.K[i,j]-bhs.K[i,j],d_K[i,j]*dg)
#
#print " Difference in fluxes"
#print "    WT      Mut      Mut-WT   delta_flux"
#print "%12.4e %12.4e %12.4e %12.4e"%(bhs.sum_flux, mut.sum_flux, mut.sum_flux-bhs.sum_flux, sum_d_flux*dg)
#
#
