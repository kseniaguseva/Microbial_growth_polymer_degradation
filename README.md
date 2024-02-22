This code produces the first part of results of the paper Guseva et al (2024),
"Bacteria face trade-offs in the decomposition of complex biopolymers". It is
written in Python, with additional functions imported from C++ using
Boost.Python. To run the code make sure to have installed Python 3.11 and
Boost.Python.

# Before running the program:
----------------------------

For installation and use of Python and Boost.Python please check the
documentation:

https://www.python.org/downloads/
https://www.boost.org/doc/libs/1_84_0/more/getting_started/index.html

Please configure the Makefile file with for the appropriate location of your
libraries (first line of the file) and appropriate version of Python installed.
Run ==make== to compile all functions written in C++. 

After that you should be ready to tun the python code.

# File list:
-----------

## Structural files
--------

1. **compunds.cc**: diffusion function and boost_python_module

2. **compunds.hh**: functions used for the substrates used in the simulation

3. **microorganisms.hh**: functions for the phyisological functions
   of microorganisms, such as uptake, maintenance, enzyme production,
   growth and mortality.

4. **enzymes.hh**: fucntions related to enzyme production, production and decay.

## Initialization files
--------

1. **get_parameters.py**: uploads parameters used to set up the
   microbial physiology and initializes microorganisms. The parameters
   should be supplied in a separate csv file.

2. **init_single.py**: Initializes a single microsite. This setting will
   correspond to a well mixed system.   

3. **init.py**: Initializes a spatial grid. The complex substrate is
   spread homogeneously on the whole surface.Microorganisms are
   initialized with an initial amount of enzymes.
   
4. **parameters.csv**: list of parameters used

## Execution files
--------

1. **main_mixed.py**: evaluates the dynamics in the well mixed
   system. The output is a txt file with the time series of concentration
   of different components of the model.


## Output
-----------
1. Please create a ./data folder for the output from main_mixed.py
