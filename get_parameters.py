#!/usr/bin/env python
from math import *
from numpy import *
from pylab import *
import pandas as pd
from compounds import *


def get_microorganisms(N_size, c_init, n_types, fr_enz1, fr_enz2, efr = 0):
    file_param = "./parameters.csv"
    data = pd.read_csv(file_param, index_col=0)
    rm = data["rm"].iloc[0]
    re = data["re"].iloc[0]
    rg = data["rg"].iloc[0]
    if(efr == 0):
        efr = data["efr"].iloc[0]
    u = data["u"].iloc[0]*8  # convert considering 8 C molecules in a monomer
    km_up = data["km_up"].iloc[0]
    d = data["d"].iloc[0]
    size_max = data["size_max"].iloc[0]
    size_min = data["size_min"].iloc[0]


    # ################################################################################
    #  initialization of microorganisms: Microorganisms(Size of the Grid, initial
    #  concentration of biomass, number of types, fraction of exo in type 1,
    #  fraction of exo in type 2 (not used if there is only one type))
    # ################################################################################

    bac = Microorganisms(N_size, c_init, n_types, 
                         rm, rg, re, efr,
                         double(u), double(km_up), d,
                         double(size_max), double(size_min),
                         fr_enz1, fr_enz2)
    return bac
