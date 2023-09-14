// Copyright (C) 2020 Ksenia Guseva <ksenia@skewed.de>

#include <random>

std::random_device _rd;
std::mt19937 _rng(_rd());

#include "compounds.hh"
#include "microorganism.hh"

void* do_import_array() {
    import_array1(NULL);
    return NULL;
}

void diffusion(multi_array_ref<double, 3>&  M,
               multi_array_ref<double, 3>&  M_new,
               multi_array_ref<double, 3>&  D,
               double dt, double dx, double dy)
{
     double Mjp, Mjm, Mip, Mim;
     int Nx =  int(M.shape()[0]);
     int Ny =  int(M.shape()[1]);
     int Nz =  int(M.shape()[2]);
     int shape = 10;

     if(Nz < 10)
     {
         shape = Nz;
     }
     
     
     for (int k = 0; k < shape; ++k)
     {        
         for (int i = 0; i < Nx; ++i)
         {
             for (int j = 0; j < Ny; ++j)
             {
                 int jp = (j < int(Ny - 1)) ? (j + 1) : 0;
                 int jm = (j > 0) ? j - 1 : Ny - 1;

                 int ip = (i < int(Nx - 1)) ? (i + 1) : 0;
                 int im = (i > 0) ? i - 1 : Ny - 1;

                 Mjp = M[i][jp][k];
                 Mjm = M[i][jm][k];
                 Mip = M[ip][j][k];
                 Mim = M[im][j][k];
                 
                 //M_new[i][j][k] = M[i][j][k];
                 // M_new[i][j][k] += dt*(D[i][j][k] + D[i][jp][k])*0.5*(Mjp - M[i][j][k])/(dy*dy);
                 // M_new[i][j][k] += dt*(D[i][jm][k] + D[i][j][k])*0.5*(Mjm - M[i][j][k])/(dy*dy);
                 // M_new[i][j][k] += dt*(D[ip][j][k] + D[i][j][k])*0.5*(Mip - M[i][j][k])/(dx*dx);
                 // M_new[i][j][k] += dt*(D[i][j][k] + D[im][j][k])*0.5*(Mim - M[i][j][k])/(dx*dx);
                 M_new[i][j][k] += dt*D[i][j][k]*(Mjp - M[i][j][k])/(dy*dy);
                 M_new[i][j][k] += dt*D[i][j][k]*(Mjm - M[i][j][k])/(dy*dy);
                 M_new[i][j][k] += dt*D[i][j][k]*(Mip - M[i][j][k])/(dx*dx);
                 M_new[i][j][k] += dt*D[i][j][k]*(Mim - M[i][j][k])/(dx*dx);

                 
             }
         }
     }
     
}


BOOST_PYTHON_MODULE(compounds)
{
    boost::python::numpy::initialize();
    do_import_array();

    using namespace boost::python;
    
    class_<Compound>("Compound", init<int, int, double, double, double, int>())
        .def("rescale_cons", &Compound::rescale_concentration)
        .def("add_cn", &Compound::add_cn)
        .def("add_grid", &Compound::add_grid)
        .def("get_c", &Compound::get_c)       
        .def("update", &Compound::update)
        .def("update_next", &Compound::update_next)
        .def("update_diff", &Compound::update_diff)
        .def("diffuse_fixedcn", &Compound::diffuse_fixedcn);
    
    
    class_<Microorganisms>("Microorganisms", init<int, double, int,
                           double, double, double,
                           double, double, double,
                           double, double, double,
                           double, double>())
        .def("rescale_cons", &Microorganisms::rescale_concentration)
        .def("get_c", &Microorganisms::get_c)
        .def("get_resp", &Microorganisms::get_resp)
        .def("get_type", &Microorganisms::get_type)
        .def("metabolism", &Microorganisms::metabolism)
        .def("mortality", &Microorganisms::mortality)
        .def("cell_division", &Microorganisms::cell_division);
    
    class_<Enzymes>("Enzymes", init<int, double, double, double>())
        .def("add_enz", &Enzymes::add_enz)
        .def("rescale_cons", &Enzymes::rescale_concentration)
        .def("get_c", &Enzymes::get_c)
        .def("update", &Enzymes::update)
        .def("exo_activity", &Enzymes::exo_activity)
        .def("endo_activity", &Enzymes::endo_activity)
        .def("decay", &Enzymes::decay_enzymes);
    
    
    
};
