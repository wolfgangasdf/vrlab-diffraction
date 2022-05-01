#include <iostream>
#include <fstream>
#include <complex>
#include <math.h>

#include "doubleslit-calc.hpp"
#include "grating-calc.hpp"

int main () {
    std::cout << "double slit demo!\n";
    Doubleslit ds;
    // this is just to demonstrate things, 
    // see docs in doubleslit-calc.hpp

    // double lambda, double beam_x, double beam_w0, double slit_w, double slit_d, int lens_exists, double lens_z, double lens_f, double screen_z
    ds.calc(0.6328, 0.0, 1500.0,  50.0,   250.0, 0, 0.0, 0.0, 600000.0, "test-ds-nolens-60cm.png"); 
    ds.calc(0.6328, 0.0, 1500.0,  50.0,   250.0, 0, 0.0, 0.0, 60000.0, "test-ds-nolens-6cm.png"); 
    ds.calc(0.6328, 0.0, 1500.0,  50.0,   250.0, 1, 100000.0, 50000.0, 150000.0, "test-ds-farfield.png"); 


    // double objlens_f, double screen_z, double collens_f, double collens_z, double grating_th0, double grating_w, double grating_d, int grating_n, int lightsource
    Grating gr;
    gr.calc(50000.0, 50000.0, 50000.0, 50000.0, 0.29893, 0.5, 2.0, 2000, 2, "test-gr-na-dlines.png"); // sodium d-lines: just visible!
    gr.calc(50000.0, 100000.0, 50000.0, 50000.0, 0.29893, 0.5, 2.0, 2000, 2, "test-gr-na-dlines-misaligned.png"); // sodium d-lines: not anymore
    gr.calc(50000.0, 50000.0, 50000.0, 50000.0, 0.29893, 0.5, 2.0, 2000, 1, "test-gr-hg-greenorange.png"); 


}
