// Copyright (C) 2020 Ksenia Guseva <ksenia@skewed.de>

#ifndef MICROORGANISM_HH
#define MICROORGANISM_HH


#include <cmath>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <ext/numeric>
using __gnu_cxx::power;
#include <vector>
#include <utility>
#include <random>

#include <boost/multi_array.hpp>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "numpy_bind.hh"

#include "compounds.hh"
#include "enzymes.hh"


using namespace std;
using namespace boost;
namespace bpn = boost::python::numpy;

extern std::mt19937 _rng;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                  Initialize
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class Microorganisms {

public:
    Microorganisms(int shape, double c_init, int ntypes,
                   double rm, double rg, double re,
                   double efr, double u, double km_up,
                   double d, double size_max, double size_min,
                   double fr_enz1, double fr_enz2)
        : _c(extents[shape][shape]),
          _type(extents[shape][shape]),
          _resp(extents[shape][shape]),
          _rm(ntypes, 0),        // resiration for maintenance
          _re(ntypes, 0),           // respiration enzyme production
          _rg(ntypes, 0),           // respiration growth
          _efr(ntypes, 0),         // fraction for enzyme production
          _u(ntypes, 0),         // uptake rate
          _km_up(ntypes, 0),         // half saturation uptake
          _d(ntypes, 0),          // mortality rate
          _size_max(ntypes, 0),
          _size_min(ntypes, 0),
          _cell_size(ntypes, 0)
    {
        // components of a microbial cell:
        // 0. fraction rich in carbon (CN0)
        // 1. fraction rich in proteins (CN1)
        // 2. fraction of DOM (CN2)

        int single = 1;
        double CN0 = 150;
        double CN1 = 5;
        double CN2 = 8;     // chitin monomer
        
        _cell_size[0] = _cell_size[1] = 4.;
        
        //cout << "!!!!!!!!!!!!!!!!!!!!!!!  " <<  ntypes << "\n";
        
        //boost::multi_array_ref<double, 1> get_array(rm);
        _rm[0] = _rm[1] = rm;           // resp0rat0on for ma0ntenance       h^{-1}
        _re[0] = _re[1] = re;             // resp0rat0on enzyme product0on     h^{-1}
        _rg[0] = _rg[1] = rg;             // resp0rat0on growth                h^{-1}
        _efr[0] = _efr[1] = efr;           // fract0on for enzyme product0on    
        _u[0] = _u[1] = u;         //{0.01854, 0.01854};   uptake rate      nmol C/( h^{-1} mm^2)
        _km_up[0] = _km_up[1] = km_up;         //{0.01854, 0.01854};   uptake rate      nmol C/( h^{-1} mm^2)
        
        _d[0] = _d[1] = d;            // mortal0ty rate                    h^{-1}
        _size_max[0] = size_max*_cell_size[0];        // nmolC
        _size_max[1] = size_max*_cell_size[1];
        _size_min[0] = size_min*_cell_size[0];   //         nmolC
        _size_min[1] = size_min*_cell_size[1];
        cout << size_min << "\n";
        
        // fractions correspond to {Prot, CW, DOM}
        _fr = {{0.39, 0.56, 0.05}, {0.39, 0.56, 0.05}};

        
        _CN = {{CN0, CN1, CN2, 1./(_fr[0][0]/CN0 + _fr[0][1]/CN1 + _fr[0][2]/CN2)},
               {CN0, CN1, CN2, 1./(_fr[1][0]/CN0 + _fr[1][1]/CN1 + _fr[1][2]/CN2)}};
        
        
        //_e = {{1 - fr_enz1 - 0.1, 0.05, 0.05, fr_enz1}, {1 - fr_enz2 - 0.1, 0.05, 0.05, fr_enz2}};
        _e = {{fr_enz1, 1-fr_enz1}, {fr_enz2, 1-fr_enz2}};
        

        std::uniform_int_distribution<> sample(0, 99);
        std::uniform_int_distribution<> sample_type(1, ntypes);
        
        if(shape  > 1)
        {
            if (shape == 2)
            {
                for (int i = 0; i < shape; i++)
                {
                    for (int j = 0; j < shape; j++)
                    {

                        _resp[i][j] = 0;
                        _c[i][j] = 0;
                        _type[i][j] = 0;
                        
                    }
                }
                _c[0][0] = c_init;
                _type[0][0] = 1;    
            }
                
            else
            {
                if(single == 0)
                {
                    for (int i = 0; i < shape; i++)
                    {
                        for (int j = 0; j < shape; j++)
                        {
                            _resp[i][j] = 0;
                            
                            int n_rand = sample(_rng);
                            _c[i][j] = 0;
                            _type[i][j] = 0;
                            
                            
                            if(n_rand < 25)
                                //if( i == 0 and j == 0)
                            {
                                if(ntypes == 1)
                                {
                                    _type[i][j] = 1; //rand() % (2) + 1;
                                }
                                else
                                {
                                    _type[i][j] =  sample_type(_rng);
                                    //cout<< _type[i][j] << "\n";
                                }
                            }
                            if (_type[i][j] != 0)
                            {
                                _c[i][j] = c_init;//*sample(_rng)/100.;
                            }
                        }
                    }
                }
                else
                {
                    cout << "entrei \n";
                    
                    for (int i = 0; i < shape; i++)
                    {
                        for (int j = 0; j < shape; j++)
                        {
                            _resp[i][j] = 0;
                            _c[i][j] = 0;
                            _type[i][j] = 0;
                        }
                    }

					// for (int i = int(shape/2)-4; i < int(shape/2)+4; i++)
                    // {
                    //     for (int j = int(shape/2)-4; j < int(shape/2)+4; j++)
                    //     {
					// 		//int aux = int(shape/2);
					// 		cout << "entrei \n";
					// 		_c[i][j] = c_init;//*sample(_rng)/1
					// 		_type[i][j] = rand() % (2) + 1;
					// 	}
					// }
                    
                    int aux = int(shape/2);
                    _c[aux][aux] = c_init;//*sample(_rng)/1
                    _type[aux][aux] = 1;

                    // not single anymore
                    // for (int k = -1; k < 1; k++)
                    // {
                    //     for (int l = -1; l < 1; l++)
                    //     {
                    //         _c[aux+k][aux+l] = c_init;//*sample(_rng)/1
                    //         _type[aux+k][aux+l] = 1;
                    //     }
                    // }
                    
                }    
            }
        }
        else
        {
             _resp[0][0] = 0;
            _c[0][0] = c_init;
            _type[0][0] = 1;
        }
        
    }

    void rescale_concentration(double MS, double C, int N)
    {
        
        for (int i = 0; i < N; i++)
        {
            _u[i] = _u[i]*pow(MS,2);
            // converting Volume/C , the constant comes from  conversion of [nmolC /mm3] to 
            _size_max[i] = _size_max[i]*pow(MS,2)*10/C;   // [fmol C]
            _size_min[i] = _size_min[i]*pow(MS,2)*10/C;   // [fmol C]
            cout << _size_min[i] << "\n";
            cout << _size_max[i] << "\n";
        }
        for (int i = 0; i < int(_c.shape()[0]); i++)
        {
            for (int j = 0; j < int(_c.shape()[0]); j++)
            {
                _c[i][j] = _c[i][j]*pow(MS,2)*10/C;     // [fmol C]
            }
        }
        
                          
        
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                 WRAP: PYTHON
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    
    python::object get_c()
    {
        return wrap_multi_array_not_owned(_c);
    }
    

    
    python::object get_type()
    {
        return wrap_multi_array_not_owned(_type);
    }

    python::object get_resp()
    {
        return wrap_multi_array_not_owned(_resp);
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                METABOLISM
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    
    double get_upot(int i, int j, int type, double MS)
    {

        //*************************************************************************************************/
        // The maximum possible uptake (unlimited by the substrate) is calculated according to:
        // U_pot = u*fr_ch*rsv*vol = u_pot*fr_ch*rsv*(biomass/rho)
        //
        /**************************************************************************************************/
        
        double mc = 12;                                      // [fg/fmol]
        // typical cell density  converted to appropriate units
        double rho = (1000)*pow(MS, 3);                      // [fg/ \mum^3] -->   [fg/ MS]

        // fraction of channels on the membrane (Folse et al  2012)
        double fr_ch = 0.1;                                                                 

        // considering that carbon account for 10pc (Romanov et al 2010) of the total microbial biomass,
        //  the total biomass of a microorganism is:
        double biomass = (_c[i][j])*mc/0.1;           // total biomass of a microorganism

        // considering the density of the cell given above the volume of the cell is    
        double vol = biomass/rho;                            //[fmolC/ MS][fg/fmol][1 MS]/[fg/MS] = [MS]

        // surface to volume ratio of a cell
        double rsv = 3./pow(((vol*3.)/(4*M_PI)),(1./3.));    //[1/(MS)^{1/3}]

        int n_cells = std::floor(_c[i][j]/_cell_size[type-1]) + 1;    // single cell's biomass carbon is 4 fmol C
        
        double u_pot = (pow(n_cells, 1./3)) * fr_ch * _u[type-1] * rsv * vol; // [1/(MS h) ][fmol C/MS][1/(MS)^{1/3}] = [fmol C/(MS h)]]
        //cout << "uptake pot "<< biomass  << "\n";
        //exit(0);
        
        return u_pot;
    }
    
    void maintenance(int i, int j, int type, double dt, double &c_uptake, double &n_uptake)
    {
        //cout << "c_uptake in func " << c_uptake << "\n";
        //cout << "maintenance " << _c[i][j]*_rm*dt << "\n";
        
        c_uptake = c_uptake - _c[i][j]*_rm[type-1]*dt;
        n_uptake += _c[i][j]*_rm[type-1]*dt/_CN[type -1][3]; // add the corresponding nitrogen, it can still be used in other steps before eliminated
        _resp[i][j] += _c[i][j]*_rm[type-1]*dt;
        
        
        //cout << "maintenance " << c_uptake << "\n";
        
    }

    // production of enzymes at the costs of biomass
    double produce_minenz(double upot, int i, int j, double dt)
    {
        double c_emin;
        int nt =_type[i][j]-1;

        // maximum possible enzymes produced for the given size when the uptake is not limited
        // externally available carbon (int eh best possible scenario)
        c_emin = (upot - _rm[nt]*_c[i][j])*dt; // we subtract respiration since we estimate the amount left 
        c_emin = 0.0*_efr[nt]*c_emin;         // we suppose that 10% of this amount has to be always produced by microorganisms
        c_emin = c_emin/(1 - _re[nt]*dt);          // the necessary carbon for enzymes and respiration in its production
        return c_emin;                           
    }

    void growth(double c_uptake, double n_uptake,  Compound &din , int i, int j, int type, double dt)
    {
        //cout << "c_uptake " << c_uptake << "\n";
        if(c_uptake > 0 and n_uptake > 0)
        {
            double c_g = c_uptake - c_uptake*_rg[type-1];
            double n_g = c_g/_CN[type -1][3];
            
            //cout << "growth " << c_g << " " << n_g << " " << n_uptake << "\n";

            if(n_g < n_uptake) // not limited by N
            {
                _c[i][j] += c_g;
                _resp[i][j] += c_uptake*_rg[type-1];
                din._c_next[i][j][0] += n_uptake - n_g;

                //cout << "cg: " << c_g <<" " << _CN[type -1][3]<<"\n";
                //cout << "nupt - ng "<< din_next[i][j][0] += n_uptake - c_uptake/_CN[type-1][3];
            }
            else
            {
                cout << "limited by N \n";
            }
            
        }         
    }
    
    void subtract_costs(double c_uptake, Compound &subCW, Compound &subPR,
                        Compound &subPS, int i, int j, int type)
    {
        
        _c[i][j] += c_uptake; // subtracts maintenance
        if(_c[i][j] <  _size_min[type -1]) // mortality due to minimum size
        {
            //cout << "costs_in"<< _c[i][j] << "\n";
            subCW.add_cn(_fr[type -1][0]*_c[i][j], i, j, 0);
            subPR.add_cn(_fr[type -1][1]*_c[i][j], i, j, 0);           
            subPS.add_cn(_fr[type -1][2]*_c[i][j], i, j, 0);
            
            _c[i][j] = 0;
            _type[i][j] = 0;
            
        }
        
    }

     void metabolism(Compound &subPS, Compound &din, Compound &subCW,  Compound &subPR, Enzymes &enz_PR, Enzymes &enz_CW, Enzymes &enz_exo, Enzymes &enz_endo, double dt, double MS, double u_fr)
    {
       
        
        for (int i = 0; i <  int(_c.shape()[0]); i++)
        {
            for (int j = 0; j < int(_c.shape()[1]); j++)
            {
                
                if(_c[i][j] != 0)
                {
                    int type = _type[i][j];
                    double upot = get_upot(i, j, type, MS);
                   
                    auto[c_uptake, n_uptake] = subPS.uptake_dom(i, j, dt, upot, _km_up[0], 1, MS);
                    n_uptake += din.uptake_din(i, j, 0, dt, MS);

                    //cout << "c " <<c_uptake << "\n";
                    //cout << "n " << n_uptake << "\n";
                    
                    maintenance(i, j, type, dt, c_uptake, n_uptake);
                        
                    double c_emin = produce_minenz(upot, i , j, dt); // for paper set to 0
                    //cout << "after maintenance " <<c_uptake << "\n";
                    if(c_uptake > 0)
                    {
                        if(c_uptake - c_emin > 0) // check for C
                        {

                            _resp[i][j] += enz_exo.produce_enzymes(c_uptake, n_uptake, i, j, dt, _efr[type-1], _re[type-1], _e[type-1][0]);
                            _resp[i][j] += enz_endo.produce_enzymes(c_uptake, n_uptake, i, j, dt, _efr[type-1], _re[type-1], _e[type-1][1]);
                            
                            
                        }
                        else // production at the consts of biomass
                        {
                            n_uptake += (c_emin - c_uptake)/_CN[type-1][3]; // nitrogen added from the biomass
                             _resp[i][j] += enz_exo.produce_enzymes(c_emin, n_uptake, i, j, dt, _efr[type-1], _re[type-1], _e[type-1][0]);
                            _resp[i][j] += enz_endo.produce_enzymes(c_emin, n_uptake, i, j, dt, _efr[type-1], _re[type-1], _e[type-1][1]);
                            c_uptake -= c_emin;
                            
                        }
                        
                        //cout << "before growth " <<c_uptake << "\n";

                        if(c_uptake > 0)
                        {
                            growth(c_uptake, n_uptake, din, i, j, type, dt);
                        }
                        else
                        {
                            subtract_costs(c_uptake, subCW,  subPR, subPS, i, j, type);
                            din._c_next[i][j][0] += n_uptake;
                        }
                        
                        
                    }
                    else
                    {
                        _resp[i][j] += enz_exo.produce_enzymes(c_emin, n_uptake, i, j, dt, _efr[type-1], _re[type-1], _e[type-1][0]);
                        _resp[i][j] += enz_endo.produce_enzymes(c_emin, n_uptake, i, j, dt, _efr[type-1], _re[type-1], _e[type-1][1]);
                        c_uptake -= c_emin;
                        
                        subtract_costs(c_uptake, subCW,  subPR, subPS, i, j, type);
                        if(n_uptake > 0) // eliminates excess nitrogen where there is no carbon for growth
                        {
                            // add nitrogen not used from the uptake and nitrogen
                            // left from the carbon used in the maintenance
                            // the negative sign is because c_uptake here is negative
                            din._c_next[i][j][0] += n_uptake - c_uptake/_CN[type-1][3];
                        }
                    }
                    
                 
                }
            }
        }

    }

    

    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                CELL DIVISION AND MORTALITY
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    void cell_division(Compound &sub, Enzymes &enz_PR, Enzymes &enz_CW, Enzymes &enz_exo, Enzymes &enz_endo, double dt)
    {
        int Nx, Ny, nt;
        Nx =  int(_c.shape()[0]);
        Ny =  int(_c.shape()[1]);
        
        multi_array<double, 2> c_new(extents[Nx][Ny]);
        multi_array<double, 2> exoc_new(extents[Nx][Ny]);
        multi_array<double, 2> endoc_new(extents[Nx][Ny]);
        multi_array<double, 2> type_new(extents[Nx][Ny]);
        c_new = _c;
        exoc_new = enz_exo._c_next;
        endoc_new = enz_endo._c_next;
        type_new = _type;
        
        for (int i = 0; i <  int(_c.shape()[0]); i++)
        {
            for (int j = 0; j < int(_c.shape()[1]); j++)
            {
                nt = _type[i][j];

                if(_c[i][j] > 0)
                {
                    //cout << nt <<" "<<_c[i][j] << " " << _size_max[nt -1]<< "\n";
                    
                
                    if(_c[i][j] > _size_max[nt -1])
                    {
                        cout << "divide! \n";

                        int n_empty = 0;
                        int n_inv = 0;
                        //int n_empty_sub = 0;
                        int ij_empty[2][9];
                        int ij_inv[2][9];


                        // creates a list of free cells
                        for (int aux_i = i-1; aux_i <= i+1; aux_i ++)
                        {
                            for (int aux_j = j-1; aux_j <= j+1; aux_j ++)
                            {

                                int ni = (aux_i + Nx) % Nx;
                                int nj = (aux_j + Ny) % Ny;

                                // if(sub._c[ni][nj][sub._c.shape()[2]] == 0)
                                // {
                                //     n_empty_sub += 1;
                                // }
                                if(c_new[ni][nj] == 0 )
                                {
                                    ij_empty[0][n_empty] = ni;
                                    ij_empty[1][n_empty] = nj;
                                    n_empty += 1;
                                }
                                else
                                {
                                    if(_type[ni][nj] !=  _type[i][j])
                                    {
                                        ij_inv[0][n_inv] = ni;
                                        ij_inv[1][n_inv] = nj;
                                        n_inv += 1;
                                    }

                                }

                            }
                        }




                        // chose a random empty cell
                        cout << n_empty <<"\n";
                        std::uniform_int_distribution<> sample_inv(0, 99);
                        int n_rand1 = sample_inv(_rng);
                        if(n_rand1 < 0*dt &&  n_inv > 0)
                        {
                            std::uniform_int_distribution<> sample_inv2(0, n_inv-1);
                            int n_rand2 = sample_inv2(_rng);
                            int i_inv = ij_inv[0][n_rand2];
                            int j_inv = ij_inv[1][n_rand2];

                            c_new[i_inv][j_inv] = _c[i][j]/2.;
                            c_new[i][j] = _c[i][j]/2.;

                            // enzymes taken with the cell (attached to the wall)
                            exoc_new[i_inv][j_inv] += enz_exo._c[i][j]/2.;
                            exoc_new[i][j] = enz_exo._c[i][j]/2.;
                            
                            endoc_new[i_inv][j_inv] += enz_endo._c[i][j]/2.;
                            endoc_new[i][j] = enz_endo._c[i][j]/2.;
                            
                            
                            type_new[i_inv][j_inv] = _type[i][j];
                        }
                        else
                        {

                            if(n_empty > 0) //&& n_empty_sub > 0)
                            {
                                std::uniform_int_distribution<> sample_empty(0, n_empty-1);
                                int n_rand = sample_empty(_rng); 

                                int i_empty = ij_empty[0][n_rand];
                                int j_empty = ij_empty[1][n_rand];

                                c_new[i_empty][j_empty] = _c[i][j]/2.;
                                c_new[i][j] = _c[i][j]/2.;

                                // enzymes taken with the cell (attached to the wall)
                                exoc_new[i_empty][j_empty] += enz_exo._c[i][j]/2.; 
                                exoc_new[i][j] = enz_exo._c[i][j]/2.;

                                endoc_new[i_empty][j_empty] += enz_endo._c[i][j]/2.;
                                endoc_new[i][j] = enz_endo._c[i][j]/2.;
                                
                                type_new[i_empty][j_empty] = _type[i][j];
                            }
                            // else // no empty cells
                            // {
                            //     _c[i][j] = _c[i][j]/2;


                            // }
                        }
                    }

                }
            }
            
        }
        _c = c_new;
        enz_exo._c_next = exoc_new;
        enz_endo._c_next = endoc_new;
        
        _type = type_new;
        
    }
    
    
    void mortality( Compound &subCW, Compound &subPR,
                    Compound &subPS, double dt)
    {
        std::uniform_real_distribution<> dist;
        
        for (int i = 0; i <  int(_c.shape()[0]); i++)
        {
            for (int j = 0; j < int(_c.shape()[1]); j++)
            {
                int type = _type[i][j];
                
                double n_rand = dist(_rng);
                
                
                if(n_rand < _d[type -1]*dt)
                {
                    subCW.add_cn(_fr[type-1][0]*_c[i][j], i, j, 0);
                    subPR.add_cn(_fr[type-1][1]*_c[i][j], i, j, 0);
                    //cout << "cn before " << dom._CN[i][j] << "\n";
                    subPS.add_cn(_fr[type-1][2]*_c[i][j], i, j, 0);
                    //cout << "cn after " << dom._CN[i][j] << "\n";
                    _c[i][j] = 0;
                    _type[i][j] = 0;
                }
                
            }
        }
        
    }
    
    

public:         
    multi_array<double, 2> _c;      
private:
    multi_array<double, 2> _type;
    multi_array<double, 2> _resp;

    vector<double> _rm;           // respiration for maintenance
    vector<double> _re;           // respiration enzyme production
    vector<double> _rg;           // respiration growth
    vector<double> _efr;          // fraction for enzyme production
    vector<double> _u;            // uptake rate
    vector<double> _km_up;            // uptake rate
    vector<double> _d;            // random death events
    vector<double> _size_max;
    vector<double> _size_min;
    vector<double> _cell_size;

    vector<vector<double>> _CN;
    vector<vector<double>> _fr;
    vector<vector<double>> _e;
};

#endif
