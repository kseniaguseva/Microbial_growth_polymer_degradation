#!/usr/bin/env python
from math import *
from numpy import *
from pylab import *


import sys
#sys.path.append('../module_microbes')
from init_single  import *
from compounds  import *

dt = 0.005 # h
time = arange(0, 1000., dt)

s = 5
N_range = [500] #[10, 1000] #linspace(10, 200, s)
totmon_range = logspace(0, 5, 20) #[10, 1000] #linspace(1, 100, s)

efr = 0.5
fr_enz1 = 1.    # exo-chitinase  ---> 1 / endo-chitinase ----> 0
fr_enz2 = 0.    # second type of bacteria (not used here for a single grid)
 

for j, N in enumerate(N_range):
    for i, totmon in enumerate(totmon_range):
        init_PS = totmon/N
        chain_size = int(N)
        # output files (needs a data folder)
        if(fr_enz1 == 1):
            outfile = open("data/exo_PS%.3f_C%.1f.txt" %(totmon, chain_size), 'w')
        else:
            outfile = open("data/endo_PS%.3f_C%.1f.txt" %(totmon, chain_size), 'w')
        bac_time = []

        print(i, j, init_PS, chain_size)
        
        sub_PS, sub_CW, sub_PR, din, enz_exo, enz_endo, enz_CW, enz_PR, bac, MS, u_fr = init_system(init_PS, int(chain_size), fr_enz1, fr_enz2, efr, False)

        for t in time:
            bac.metabolism(sub_PS, din, sub_CW, sub_PR, enz_PR, enz_CW, enz_exo, enz_endo, dt, MS, u_fr)
        
            enz_exo.exo_activity(sub_PS, dt)
            enz_endo.endo_activity(sub_PS, dt)
            for enz in [enz_exo, enz_endo]:
                enz.decay(dt, sub_PR)

            for comp in [sub_PS, sub_CW, sub_PR, enz_PR, enz_CW, enz_endo, enz_exo]:
                comp.update()


            print(t, sum(bac.get_c()), sum(bac.get_resp()), sum(sub_PS.get_c()[:, :, 0]), sum(sub_PS.get_c()[:, :, -1]*chain_size), sum(enz_endo.get_c()), sum(enz_exo.get_c()),sum(sub_PR.get_c()[:, :, 0]), sum(sub_CW.get_c()[:, :, 0]), file = outfile)


        outfile.close()





