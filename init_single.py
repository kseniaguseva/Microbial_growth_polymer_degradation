#!/usr/bin/env python
from math import *
from numpy import *
from pylab import *
from compounds import *
from get_parameters import *

def get_total(n):
    sum_n = 0
    for i in range(0, len(n)):
        sum_n += n[i]*(i+1)
    return sum_n

 ################################################################################
 #   Units for the variables in the INPUT:
 #   Concentration -> nmol C/mm^3
 #   time -> hours
 #
 #   The program converts all units to fmol C/MS  (fmol C/microsite)
 ################################################################################

#  initial concentration of the primary substrate in (given in nmol C/mm^3) 
#  init_PS -- nmol C/mm^3
#  chain size of the chitin polymer 
def init_system(init_PS, chain_size, fr_enz1, fr_enz2, efr, output = True):

    ################################################################################
    #            GRID      
    ################################################################################
    N_size = 1        # grid size
    L = 10.         # mm --> micrometer (grid size in micrometers)
    MS = L/N_size       # N_size - number of microsites in the grid
                        #  size of the microsite in micrometers

    ################################################################################

    #chain_size = 100     # chain size of the chitin polymer 

    CN_PS  = 8           # C:N ratio for primary substrate
    CN_CW  = 150         # C:N ratio for C rich substrate (from cell wall)
    CN_PR  = 5           # C:N ratio for N rich substrate  

    u_fr = 1.            # to differentiate uptake of different substances/NOT used

    ################################################################################
    #  Initial concentration
    ################################################################################

    
    #init_PS = 0.8333           # nmol C/mm^3
    enz_u = 0.00001           # initial concentration of enzymes
    add_enz = 0.75

    # diffusion coefficient in  micrometers^2 p hour
    D_nag = 30.                # mum^2/h  
    D_din = 30.                # mum^2/h 
   # same units as MS

    ################################################################################
    #         h^{-1}
    #         nmolC/mm^3         
    #  this version has several sizes of substrate chains
    #  sub_PS.get_c().[j, i, 0]  -- monomers
    #  sub_PS.get_c().[j, i, 1]  -- dimers
    #  sub_PS.get_c().[j, i, -1]  -- chitin chain
    #  Compound(grid-size, chain-size, initial concentration C or N, CN ratio,
    #  diffusion rate, shape of the initial distribution of substrate)
    ################################################################################

    sub_PS = Compound(N_size, chain_size, init_PS, CN_PS, D_nag, 0)
    #sub_PS.add_grid(0.5, 1., 0)
    
    # other pools single size of molecules 
    sub_CW = Compound(N_size, 1,  0.1, CN_CW, 0,  0)
    sub_PR = Compound(N_size, 1, 0.1, CN_PR, 0, 0)

    din = Compound(N_size, 1, 1., 0., D_din, 0)

    ################################################################################
    #  1 mm^3 -> (1000 * micrometer)^3 = 10^9 micrometers^3 (mum^3)
    #  1 nmol -> 10^6 fmol
    #  nmol/mm^3 -> 10^6 fmol/ 10^9 mum^3 = 10^{-3} fmol/mum^3
    #  (fmol/mum^3)/(10^3) = (fmol/mum^3)/convert for:
    ################################################################################

    convert = 10**3
    for sub in [sub_PS, sub_CW, sub_PR]:
        sub.rescale_cons(MS, convert)


    for liab in [din]:
        liab.rescale_cons(MS, convert)

    # ################################################################################
    #   Enzymes(grid-size, initial concentration, km, vmax)
    # ################################################################################

    enz_exo = Enzymes(N_size, enz_u*fr_enz1, 1, 0.82)
    enz_endo = Enzymes(N_size, enz_u*(1-fr_enz1), 1, 0.82)

    enz_CW = Enzymes(N_size, 0.01, 0.50, 0.82)
    enz_PR = Enzymes(N_size, 0.01, 0.50, 0.82)

    for enz in [enz_exo, enz_endo, enz_CW, enz_PR]:
        enz.rescale_cons(MS, convert)


    # ################################################################################
    #  initialization of microorganisms: Microorganisms
    #  Parameters are taken from a file: ./parameters/*.csv
    # ################################################################################

    bac_init = 4. 
    bac_ntypes = 1

    # get_microorganisms( grid_size, initial concentration, number of types)
    bac = get_microorganisms(N_size, bac_init, bac_ntypes, fr_enz1, fr_enz2, efr)   # reads file with parameters 
    b = bac.get_c()
    for i in range(0, b.shape[0]):
        for j in range(0, b.shape[1]):
            if(b[i,j] != 0):
                enz_exo.add_enz(i, j, add_enz*fr_enz1)
                enz_endo.add_enz(i, j, add_enz*(1-fr_enz1))
    bac.rescale_cons(MS, convert, 1)          # (MS, ? , number of species)

    if(output == True):
        print("######################################################")
        print("#      PARAMETERS                                    #")
        print("######################################################")
        print("MS = %.1f mum" %MS)
        print("\n")

        print("######################################################")
        print("#      INITIAL CONDITIONS                            #")
        print("######################################################")
        print("Substrate = %.2f fmol/MS" %sub_PS.get_c()[0,0,-1])
        print("Biomass = %.1f  fmol/MS" %(mean(bac.get_c().flatten())))
        print("Enzymes:")
        print("exo-chitinase = %.3f fmol/MS" %(mean(enz_exo.get_c().flatten())))
        print("endo-chitinase = %.3f fmol/MS" %(mean(enz_endo.get_c().flatten())))

        print("######################################################")
    return sub_PS, sub_CW, sub_PR, din, enz_exo, enz_endo, enz_CW, enz_PR, bac, MS, u_fr
