#!/usr/bin/env python

import sys,os,copy
import fourstate

print "\n WT MARKOV MODEL"
bhs = fourstate.FourState()
#bhs.run_simul("WT")
bhs.run_commit()
bhs.fold_rate()
#bhs.do_dot()

print "\n\n MUTANT MARKOV MODEL"
mut = (copy.deepcopy(bhs)) #.mutate(1)
mut.mutate(2,1)
mut.calc_eigs()
#mut.run_simul("Mutant")
mut.run_commit()
mut.fold_rate()
#mut.do_dot("Mut.dot")


print "\n WT MARKOV MODEL: PERTURBATION"
bhs.partial_flux(2,1)

