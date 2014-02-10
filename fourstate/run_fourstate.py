#!/usr/bin/env python

import sys,os,copy
import fourstate

print "\n WT MARKOV MODEL"
bhs = fourstate.FourState()
bhs.run_simul("WT")
bhs.run_commit()
bhs.do_dot()

print "\n\n MUTANT MARKOV MODEL"
mut = (copy.deepcopy(bhs)) #.mutate(1)
mut.mutate(2,10)
mut.calc_eigs()
mut.run_simul("Mutant")
mut.run_commit()
mut.do_dot("Mut.dot")
