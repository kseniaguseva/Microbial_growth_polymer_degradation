// Copyright (C) 2020 Ksenia Guseva <ksenia@skewed.de>

#ifndef ENZYMES_HH
#define ENZYMES_HH

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

#include "numpy_bind.hh"

#include "compounds.hh"

using namespace std;
using namespace boost;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                  Initialize
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
class Enzymes {

public:
    Enzymes(int shape, double c_init, double km, double vmax)
        : _c(extents[shape][shape]),
          _c_next(extents[shape][shape]),
          _CN(5),
          _vmax_exo(0.82),            //   h^{-1}
          _vmax_endo(0.82),            //   h^{-1}
          _ke(0.0015),             //   h^{-1}
          _km_exo(1.),               //   nmolC
          _km_endo(1.)               //   nmolC
          
    {
        _km_exo = km;
        _km_endo = km;
        _vmax_exo = vmax;
        _vmax_endo = vmax; //vmax/2;    // we consider that endo-chitinase is slower
       
        
            
            
        for (int i = 0; i < shape; i++)
        {
            for (int j = 0; j < shape; j++)
            {
                
                _c[i][j] = c_init;
                _c_next[i][j] = c_init;
                
            }
        }
    }
    
     void rescale_concentration(double MS, double C)
    {
        _km_exo = _km_exo;//*pow(MS,3)/C;
        _km_endo = _km_exo;//*pow(MS,3)/C;
        for (int i = 0; i < int(_c.shape()[0]); i++)
        {
            for (int j = 0; j < int(_c.shape()[1]); j++)
            {
              
                _c[i][j] = _c[i][j]*pow(MS,2)*10/C;                
                _c_next[i][j] = _c_next[i][j]*pow(MS,2)*10/C;
            }
        }
    }

    void add_enz(int i, int j, double v)  // use to add enzymes at the position of microorganisms
    {
        _c[i][j] += v;
        _c_next[i][j] += v;
    }
    

     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                   WRAP: PYTHON
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    

    python::object get_c()
    {
        return wrap_multi_array_not_owned(_c);
    }

    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                   DYNAMICS
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    double produce_enzymes(double &c_uptake, double &n_uptake, int i, int j, double dt,
                         double efr, double re, double e)
    {
        // efr --> fraction of remaining C used for enzyme production
        // c_cost = c_enz + c_resp
        // the total cost of enzyme production:
        // carbon for the enzymes + carob respired in the synthesis
        double resp = 0;
        double c_ecost_pot = efr*c_uptake/(1 - re*dt);
        double n_ecost_pot = c_ecost_pot*(1 - re*dt)/_CN;

        if(c_uptake  > 0)
        {
            if(c_ecost_pot <= c_uptake)
            {
                if(n_ecost_pot <= n_uptake)
                {
                    resp += c_ecost_pot*re*dt;

                    // large uptake --> follow with enzyme production
                    add_enzymes(i, j, dt, c_ecost_pot, n_ecost_pot, c_uptake, n_uptake, re, e);
                }
                else
                {
                    resp += (n_uptake*_CN/(1-re*dt))*re*dt;
                    // n-limited 
                    add_enzymes(i, j, dt, n_uptake*_CN/(1 - re*dt), n_uptake,
                                c_uptake, n_uptake, re, e);
                }   
            }
            else
            {
                if(n_ecost_pot <= n_uptake)
                {
                    // c-limited  (this only occurs if 1 <= efr + re*dt)
                    add_enzymes(i, j, dt, c_uptake, c_uptake, c_uptake, n_uptake, re, e);
                }
                else
                {
                    //cout << "Limited by all";
                    //exit(0);
                }

            }

        }
        else
        {
            if(n_ecost_pot <= n_uptake)
            {
                add_enzymes(i, j, dt, c_ecost_pot, n_ecost_pot, c_uptake, n_uptake, re, e);
            }
            
        }        
        return resp;
    }
    
    
    void add_enzymes(int i, int j, double dt, double c_ecost, double n_ecost, double &c_uptake,
                     double &n_uptake, double re, double e)
    {
        c_uptake -= e*c_ecost;
        n_uptake -= e*n_ecost;

        double c_prod = e*c_ecost*(1 - re*dt);
        _c_next[i][j] += c_prod;
        if(c_prod < -0.0000000000000001)
        {
            exit(0);
            
        }        
        
    }
    
    void decay_enzymes(double dt, Compound &sub)
    {
         for (int i = 0; i <  int(_c.shape()[0]); i++)
        {
            for (int j = 0; j < int(_c.shape()[1]); j++)
            {
                if(_c_next[i][j] - _ke*_c[i][j]*dt > 0)
                {
                    sub.add_cn(_ke*_c[i][j]*dt, i, j, 0);
                    _c_next[i][j] -=  _ke*_c[i][j]*dt;
                }
                else
                {
                    _c_next[i][j] = 0.00;
                }
                
                
            }
        } 
    }



    void exo_activity(Compound &sub, double dt)
    {
        //double dn[shape];

        int shape = sub._c.shape()[2];
        
        for (int i = 0; i <  int(_c.shape()[0]); i++)
        {
            for (int j = 0; j < int(_c.shape()[1]); j++)
            {
                double dn[shape];
                for (int k = 0; k < shape; k++)
                {
                    dn[k] = 0;
                }
                double sum_n = 0;
                // sums all, except monomers
                for (int k = 1; k < shape; k++)
                {
                    sum_n += sub._c[i][j][k];
                }
        
                double c_u = _vmax_exo * _c[i][j]/ (sum_n + _km_exo);  // same for all sizes
                // depends on the whole pool of particles in the system (except monomers)
               
                ///////////////////////////////////////////////////////////////
                // 1. degradation into monomers
                double sum_n2 = 0;
                for (int k = 2; k < shape; k++)
                {
                    sum_n2 += sub._c[i][j][k];
                }
                
                // first term: from dimes, second term from all the rest
                dn[0] = (2 * sub._c[i][j][1] + sum_n2) * c_u;

                // 2. degradation into all other sizes, except the n-mer
                for (int k = 1; k < shape - 1; k++)
                {
                    dn[k] = c_u*(- sub._c[i][j][k]  + sub._c[i][j][k+1]);
                }

                // 3. degradation of n-mers (has only a negative term)
                dn[shape-1] = (-sub._c[i][j][shape-1] * c_u);

                ///////////////////////////////////////////////////////////////
                for (int k = 0; k < int(shape); k++)
                {
                    dn[k] = dn[k]*dt;

                    if((sub._c_next[i][j][k] + dn[k]) < 0)
                     {
                         sub._c_next[i][j][k] = 0.;
                         
                     }
                     else
                     {
                         sub.add_cn(dn[k], i, j, k);
                     }
                }
                
            }
            
        }
    }
    void endo_activity(Compound &sub, double dt)
    {
        int grid_size = int(_c.shape()[0]);
        int shape = sub._c.shape()[2];
        for (int i = 0; i <  grid_size; i ++) // c = c*e   takes into account enzyme concentration
        {
            for (int j = 0; j < grid_size; j ++)
            {
                double dn[shape];
                for (int k = 0; k < shape; k++)
                {
                    dn[k] = 0;
                }
                
                // sums all, except monomers
                double sum_ni = 0;
                for (int k = 1; k < shape; k++)
                {
                    sum_ni += sub._c[i][j][k]*k;
                }
                
                 double c_u = _vmax_endo * _c[i][j]/ (sum_ni + _km_endo);  // same for all sizes
                 for (int k = 0; k < shape; k++)
                 {
                    if(k  > 0)
                    {
                        dn[k] -= c_u*sub._c[i][j][k]*k;
                    }
                    for (int l = k + 1; l < shape; l++)                    //the smallest size produced are dimers
                    {
                        dn[k] += c_u*sub._c[i][j][l]*2;
                    }
                 }
                 for (int k = 0; k < int(shape); k++)
                 {                     
                     dn[k] = dn[k]*dt;

                     
                     if((sub._c_next[i][j][k] + dn[k]) < 0)
                     {
                         sub._c_next[i][j][k] = 0.;
                         
                     }
                     else
                     {
                         sub.add_cn(dn[k], i, j, k);
                     }
                     
                     
                 }
                 
            }
        }
    }
    


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                   UPDATE
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
  
    void update()
    {
        _c =  _c_next;
    }
    
    
    
public:
    multi_array<double, 2> _c;
    multi_array<double, 2> _c_next;
private:
    double  _CN;
    double _vmax_exo;
    double _vmax_endo;
    double _ke;
    double _km_exo;
    double _km_endo;
    
};

#endif
