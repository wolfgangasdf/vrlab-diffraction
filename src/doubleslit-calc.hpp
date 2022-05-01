#ifndef DOUBLESLIT_H
#define DOUBLESLIT_H

#include "helpers.hpp"
#include <complex>
#include <math.h>
#include <vector>
#include <algorithm>


class Doubleslit {
    public:
        /**
            All lengths in micrometer (um).
            double slit at z=0
            gbdx is gaussian beam displacement on double slit

            @param lambda wavelength of the light
            @param beam_x displacement of gaussian beam horizontally on double slit
            @param beam_w0 beam waist
            @param slit_w width of the slits
            @param slit_d distance between the slits
            @param lens_exists either 0 or 1
            @param lens_z z-position of the lens
            @param lens_f focal length of the lens
            @param screen_z z-position of the camera/screen.
            @return png image.
        */
        std::vector<unsigned char> calc(double lambda, double beam_x, double beam_w0, double slit_w, double slit_d, int lens_exists, double lens_z, double lens_f, double screen_z, std::string filename);
};

realfield_type gbslits(realfield_type xax, double gb0w, double gbdx, double slitd, double slitw, double sscalef) {
    realfield_type field0(xax.size());
    for (int i=0; i<xax.size(); i++) {
        if ((abs(xax[i]/sscalef)>(slitd/2+slitw/2)) || (abs(xax[i]/sscalef)<(slitd/2-slitw/2)))
            field0[i]=0.0;
        else
            field0[i] = exp(-1.0*pow(xax[i]-gbdx,2)/pow(gb0w,2));
    }
    return field0;    
}

// 2048 pixels, 14um pixel size
realfield_type alphaaxis() {
    realfield_type res(2048);
    for (int i=0; i<res.size(); i++) res[i]=-1024*14+i*14;
    return res;
}

// propagate brute force fresnel https://en.wikipedia.org/wiki/Fresnel_diffraction
realfield_type propfbf_int(realfield_type field, realfield_type xax, realfield_type txax, real_type tz, real_type lam) {
    realfield_type fieldt(txax.size());
    for (int i=0; i<fieldt.size(); i++) fieldt[i]=0;
    real_type k=2*M_PI/lam;
    for (int txi=0; txi<txax.size(); txi++) {
        real_type tx=txax[txi];
        complex_type sum = 0.0;
        for (int j=0;j<xax.size();j++) sum+=field[j]*exp(imagi*k/(2*tz)*pow(tx-xax[j],2));
        fieldt[txi]=pow(abs(sum),2);
    }
    return fieldt;
}

// makes line plot figure. y is 0..1, NOT rescaled, do before! leave filename empty to just retrieve array.
std::vector<unsigned char> plotit2(realfield_type px, realfield_type py, std::string filename) {
    // make plot manually. for now: plot autoscale. 
    // can have more simu pixels than disp pixels (inefficient...)
    int pixx=320; 
    int pixy=240;
    // walk over points and draw line.
    realfield2d_type pres(pixx, std::vector<real_type>(pixy,0));
    int pxiold=-1;
    int pyiold=-1;
    for (int xi=0; xi<px.size(); xi++) {
        int pxi=1+floor(xi*pixx/py.size());
        if (pxi >= pixx) pxi = pixx-1;
        int pyi=2+floor(py[xi]*pixy);
        if (pyi >= pixy) pyi = pixy-1;
        if (pxiold > -1) {
            // draw line
            for (int lxi=std::min(pxiold, pxi); lxi<=std::max(pxiold, pxi); lxi++) {
                for (int lyi=std::min(pyiold, pyi); lyi<=std::max(pyiold, pyi); lyi++) {
                    // std::cout << "lxi=" << lxi << " lyi=" << lyi << " \n";
                    pres[lxi][lyi]=1;
                }
            }
        }
        pxiold=pxi;
        pyiold=pyi;
    }
    return export_png(pres, pres, pres, filename);
}


std::vector<unsigned char> Doubleslit::calc(double lambda, double beam_x, double beam_w0, double slit_w, double slit_d, int lens_exists, double lens_z, double lens_f, double screen_z, std::string filename) {
    std::cout << "double slit demo from calc!\n";
    const double xaxr=500.0;
    const int xaxn=512;
    // source plane
    realfield_type xax(xaxn);
    for (int i=0; i<xaxn; i++) xax[i]=-xaxr+2.0*xaxr/xaxn * i;
    if (lens_exists==0) { // doubleslit at z=0, no lens, fresnel propagation
        if (screen_z<0) {
            std::cout << "error screen_z=" << screen_z << "\n";
            return std::vector<unsigned char>();
        }
        // gaussian beam through slits
        auto field1=gbslits(xax, beam_w0, beam_x, slit_d, slit_w, 1.0);

        // propagate to plane
        auto axax=alphaaxis();
        auto inta=propfbf_int(field1,xax,axax,screen_z,lambda);
        // normalize
        inta=int_normalize(inta);
        // real_type intamax=0.0;
        // for (int i=0; i<axax.size(); i++) if (inta[i]>intamax) intamax=inta[i];
        // for (int i=0; i<axax.size(); i++) inta[i]=inta[i]/intamax;

        // plot
        // for (int i=0; i<axax.size(); i++) if (i>100) inta[i]=(pow(sin(axax[i]/1000),2)+0.5)/1.5; else inta[i]=1.0;
        return plotit2(axax, inta, filename);
    } else { // with lens
        if (!(lens_z>0 && screen_z>0 && lens_z < screen_z)) {
            std::cout << "error parameters (with lens)!\n";
            return std::vector<unsigned char>(); 
        }
        // deviations from far field (there is always a far field)
        // in units of lensf
        real_type screen_dev_from_ff = abs(screen_z-lens_z-lens_f)/lens_f;

        // deviations from imaging, set very high if no nf plane (lens too close to slits)!
        real_type nfdist=1/(1/lens_f-1/lens_z); // from lens center...
        if (nfdist<0) nfdist=1.0e20;
        real_type screen_dev_from_nf = abs(((screen_z-lens_z)-nfdist)/lens_f);
        if (screen_dev_from_nf<0.05) screen_dev_from_nf = 0;
        std::cout << "screen_dev_from_nf=" << screen_dev_from_nf << "\n";
        std::cout << "screen_dev_from_ff=" << screen_dev_from_ff << "\n";

        // only calculate stuff if <maxdev within validity
        real_type maxdev=3; // max deviation of screen from nf/ff point in units of lensf
        auto axax=alphaaxis();

        realfield_type intff(axax.size(), 0.0);
        realfield_type intnf(axax.size(), 0.0);

        if (screen_dev_from_ff<maxdev) { // calculate far-field like field
            // at lens, it's like no lens(assume ff), at lensz+lensf: scaled ff
            real_type fffact=1/lens_f+(1/lens_z-1/lens_f)*screen_dev_from_ff;
            // calculate approx. field to consider gauss<>slit misalignment
            // gaussian beam through slits
            auto field1=gbslits(xax, beam_w0, beam_x, slit_d, slit_w, 1.0);
            realfield_type ffax(axax.size());
            for (int i=0;i<axax.size();i++) ffax[i]=axax[i]*fffact*screen_z; // fake x-axis for fresnel diffr.
            intff=propfbf_int(field1,xax,ffax,screen_z,lambda); // fake
            intff=int_normalize(intff);
            real_type frad=100*screen_dev_from_ff/maxdev;
            std::cout << "frad=" << frad << "\n";
            intff=blur(axax, intff, frad);
        }
        if (screen_dev_from_nf<maxdev) { // calculate near-field
            real_type scalef=(screen_z-lens_z)/lens_z;
            std::cout << "scalef=" << scalef << "\n";
            intnf=gbslits(axax, beam_w0, beam_x, slit_d, slit_w, scalef);
            real_type frad=100*screen_dev_from_nf/maxdev;
            std::cout << "frad=" << frad << "\n";
            intnf=blur(axax, intnf, frad);
        }
        // incoherently mix both
        real_type fffract=0.5;
        if(screen_dev_from_ff<maxdev && screen_dev_from_nf<maxdev) {
            real_type aff=1/(screen_dev_from_ff+0.01);
            real_type anf=1/(screen_dev_from_nf+0.01);
            fffract=aff/(aff+anf);
        } else if (screen_dev_from_nf<maxdev) {
            fffract=0;
        } else if (screen_dev_from_ff<maxdev) {
            fffract=1;
        }
        std::cout << "fffract=" << fffract << "\n";
        realfield_type totint(axax.size());
        for (int i=0;i<axax.size();i++) totint[i]=fffract*intff[i]+(1-fffract)*intnf[i];
        return plotit2(axax, totint, filename);
    }
}



#endif
