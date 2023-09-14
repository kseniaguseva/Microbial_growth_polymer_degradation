// Copyright (C) 2020 Ksenia Guseva <ksenia@skewed.de>

#ifndef COMPOUNDS_HH
#define COMPOUNDS_HH

#include <cmath>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <ext/numeric>
using __gnu_cxx::power;
#include <vector>
#include <utility>

#include <boost/multi_array.hpp>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "numpy_bind.hh"


using namespace std;
using namespace boost;

extern std::mt19937 _rng;

void diffusion(multi_array_ref<double, 3>&  M,
               multi_array_ref<double, 3>&  M_new,
               multi_array_ref<double, 3>&  D,
               double dt, double dx, double dy);



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                    Initialize
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class Compound {

public:
    // shape -- grid size
    // N -- chain size
    // c_init -- initial concentration of carbon
    // CN -- carbon to nitrogen ratio of the substance
    // D_init -- initial diffusion
    // if we start as a circle
    Compound(int shape, int N, double c_init, double CN, double D_init, int circle)
        :_l(0.0516),
         _CN(0.),
         // space grid (shape x shape) and size of the chain (N)
         _c(extents[shape][shape][N]),
         _D(extents[shape][shape][N]),
         _c_next(extents[shape][shape][N])         
    {

        _CN = CN;
        
        std::uniform_int_distribution<> sample(0, 99);
        if(circle ==1)
        {
            for (int i = 0; i < shape; i++)
            {
                for (int j = 0; j < shape; j++)
                {
                    for (int k = 0; k < N; k++)
                    {
                        _c[i][j][k] =  0;
                        _D[i][j][k] = D_init/pow((k+1), 1./3);
                    }
                    
                    _D[i][j][N-1] = 1.;///(pow(N, 2));
               
                    if(sqrt(pow(i-50,2) + pow(j-50,2)) < 5)
                    {
                        _c[i][j][0] =  c_init;   
                        
                    }
                }
            }
        }
        else
        {
            for (int i = 0; i < shape; i++)
            {
                for (int j = 0; j < shape; j++)
                {
                    for (int k = 0; k < N; k++)
                    {
                        _c[i][j][k] =  0;
                        _D[i][j][k] = D_init/pow((k+1), 1);
                   
                    }
                    
                    _c[i][j][N-1] =  c_init;
                    
                }
            }
            
        }
        
        
        _c_next = _c;
        
            
    }

    // we use the same set up: a gird of 1 mm2 area and 10 micro height
    // the grid can be divided into different number of sites
    
    void rescale_concentration(double MS, double C)
    {
        
        for (int i = 0; i < int(_c.shape()[0]); i++)
        {
            for (int j = 0; j < int(_c.shape()[1]); j++)
            {
                for (int k = 0; k <  int(_c.shape()[2]); k++)
                { 
                    _c_next[i][j][k] = _c_next[i][j][k]*pow(MS, 2)*10/C;
                    _c[i][j][k] = _c[i][j][k]*pow(MS, 2)*10/C;
                    _D[i][j][k] =_D[i][j][k]/pow(MS, 2);
                }
            }
        }
        
    
    }
    
 

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                    WRAP: PYTHON
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    python::object get_c()
    {
        return wrap_multi_array_not_owned(_c);
    }

   

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                   GRID MANIPULATION
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    
    void sub_cn(double c, int i, int j, int k)
    {
        
        _c_next[i][j][k] -= c;     
    }

    
    void add_cn(double c, int i, int j, int k)
    {
        _c_next[i][j][k] += c;
    }

    


    void add_grid(double c, double dt, int k)
    {
        for (int i = 0; i < int(_c.shape()[0]); i++)
        {
            for (int j = 0; j < int(_c.shape()[1]); j++)
            {
                _c_next[i][j][k] += c*dt;
            }
        }    
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                   DIFFUSION
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void diffuse_fixedcn(double dt, double dx, double dy)
    {
        int Nx, Ny, Nz;
        Nx =  int(_c.shape()[0]);
        Ny =  int(_c.shape()[1]);
        Nz =  int(_c.shape()[2]);
       
        multi_array<double, 3> c_new(extents[Nx][Ny][Nz]);
        multi_array<double, 3> n_new(extents[Nx][Ny][Nz]);
        c_new = _c;
        
        diffusion(_c, c_new, _D, dt, dx, dy);
        
        _c = c_new;
        _c_next = _c;
        
    }
    
   

     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                   UPTAKE
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::tuple<double, double> uptake_dom(int i, int j, double dt, double u_pot, double km_up, double u_fr,  double MS)
    {
        double c_uptake = 0;
        double n_uptake = 0;
        u_pot = u_pot /u_fr;
        
        if(_c[i][j][0] > 0)
        {
            c_uptake = u_pot*_c[i][j][0]*dt/(km_up + _c[i][j][0]);
            //cout << "upot   "<< u_pot <<" "<< km_up + _c[i][j][0] << "\n";
            //cout << c_uptake <<" "<< _c[i][j][0] << "\n";
            
            if(c_uptake < _c[i][j][0])
            {
                n_uptake = c_uptake/_CN;
                sub_cn(c_uptake, i, j, 0);
            }
            else
            {
                
                c_uptake = _c[i][j][0];    
                n_uptake = _c[i][j][0]/_CN;
                _c_next[i][j][0] = 0.00001;
                
                if(n_uptake < 0) // security check
                {
                    cout << " c0 n0 " << _c[i][j][0] <<"\n";
                    cout << " c_uptake" << c_uptake << "\n";
                    cout << " n uptake " << n_uptake << "\n";
                }
                    
                
            }
            //c_uptake = c_uptake*8;    // converts into mol C from mols of monomers of chitin
            // not necessary for nitrogen uptake (1-to-1 ratio)
        }
      
        
        return {c_uptake, n_uptake}; 
        
    }
    
    double uptake_din(int i, int j, int k, double dt, double MS)
    {
        double n_uptake = 0;

        if(_c[i][j][k] > 0)
        {
            n_uptake = 0.95*_c[i][j][k];
            _c_next[i][j][k] -= n_uptake;
          
        }

        
        return n_uptake;
        
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                   UPDATES
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void update()
    {
        _c =  _c_next;
        
    }

    void update_next()
    {
        _c_next =  _c;
    }

    
    
    
    void update_diff(Compound  sub, double D, int k)
    {
        for (int i = 0; i < int(_D.shape()[0]); i++)
        {
            for (int j = 0; j < int(_D.shape()[1]); j++)
            {
                if(sub._c[i][j][k] > 0.1)
                {
                    _D[i][j][k] = D/1;
                }
                else
                {
                    _D[i][j][k] = D;
                }
                
            }
        }
         

    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                   Variables
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
          
public:
    double _l;
    double _CN;
    multi_array<double, 3> _c;
    multi_array<double, 3> _D;
    multi_array<double, 3> _c_next;
};

#endif
