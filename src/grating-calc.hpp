#ifndef GRATING_H
#define GRATING_H

#include "helpers.hpp"
#include <complex>
#include <math.h>
#include <vector>
#include <algorithm>


class Grating {
    public:
        /**
            All lengths in micrometer (um).

            @param objlens_f objective lens focal length
            @param screen_z z-position of the camera/screen in um BEHIND objlens (at objlens_f in focus).
            @param collens_f collimator lens focal length
            @param collens_z collimator lens z position in um BEHIND lamp
            @param grating_th0 angle of the kijker (objective lens + screen) wrt grating. if grating rotated, add this angle!
            @param grating_w grating slits width
            @param grating_d grating slits spacing
            @param grating_n grating number of illuminated slits
            @param lightsource 0-500nm test light, 1-hg, 2-na
            @return png image
        */
        std::vector<unsigned char> calc(double objlens_f, double screen_z, double collens_f, double collens_z, double grating_th0, double grating_w, double grating_d, int grating_n, int lightsource, std::string filename);
};

// https://en.wikipedia.org/wiki/Fraunhofer_diffraction_equation#Gratings
// grw=slit width, grd=slit distance ,grn=number of slits
realfield_type gratingmonoint(realfield_type th, real_type grw, real_type grd, real_type grn, real_type lam) {
    realfield_type res(th.size());
    for (int i=0;i<th.size();i++) {
        real_type sinth=sin(th[i]);
        real_type a=M_PI*grd*sinth/lam;
        // finite slits
        res[i]=1.0/pow(grn,2) * pow(sinc(grw*sinth/lam),2) * pow(sin(grn*a),2) / pow(sin(a),2);
    }
    return res;
}


std::vector<unsigned char> Grating::calc(double objlens_f, double screen_z, double collens_f, double collens_z, double grating_th0, double grating_w, double grating_d, int grating_n, int lightsource, std::string filename) {
    std::cout << "calc gradting!\n";
    const int ccdpixx=1024;
    const int ccdpixy=768;
    const int ccdpixsize=4.0;
    realfield_type wpax(ccdpixx);
    realfield_type wax(ccdpixx);
    realfield_type wthax(ccdpixx);
    realfield_type wthaxrot(ccdpixx);
    for (int i=0;i<ccdpixx;i++) {
        wpax[i]=-512+i;
        wax[i]=ccdpixsize*wpax[i];
        wthax[i]=atan(wax[i]/objlens_f);
        wthaxrot[i]=grating_th0+wthax[i];
    }

    realfield_type lams(0);
    if (lightsource==2) { // Na
        // https://physics.nist.gov/PhysRefData/Handbook/Tables/sodiumtable2.htm
        lams.insert(lams.end(),{4392.81e-4,4405.12e-4,4455.23e-4,4490.87e-4,5889.950e-4,5895.924e-4});
        // doublet at ged=2, th=0, 0.29893
    } else if (lightsource==1) { // Hg}
    // https://physics.nist.gov/PhysRefData/Handbook/Tables/mercurytable2.htm
        lams.insert(lams.end(),{3983.931e-4,4046.563e-4,4358.328e-4,5460.735e-4,6149.475e-4,7944.555e-4});
    } else { // test
        lams.insert(lams.end(),{0.5}); // ged=2.0 => th=0, 0.25268, 0.5236
    }

    real_type blurrad=(abs(screen_z-objlens_f)/objlens_f + abs(collens_z-collens_f)/objlens_f) * 100.0;

    // ccd with auto scale
    realfield2d_type vrgb(ccdpixx, realfield_type(3));
    real_type vrgbmax = 0.0;
    for (int li=0;li<lams.size();li++) {
        realfield_type res = gratingmonoint(wthaxrot, grating_w, grating_d, grating_n, lams[li]);
        realfield_type resb = blur(wax, res, blurrad);
        auto rgb=wavelength2color(lams[li]*1000);
        for (int ci=0;ci<3;ci++) {
            for (int i=0;i<ccdpixx;i++) {
                vrgb[i][ci]=vrgb[i][ci]+resb[i]*rgb[ci];
                if (vrgbmax<vrgb[i][ci]) vrgbmax = vrgb[i][ci];
            }
        }
    }
    // std::cout << "xx: "; for (int i=0; i<50; i++) std::cout << vrgb[i][2] << ","; std::cout << "\n";
    for (int ci=0;ci<3;ci++) { // rescale
        for (int i=0;i<ccdpixx;i++) {
            vrgb[i][ci]=vrgb[i][ci]/vrgbmax;
        }
    }

    // export_pnm_line2pic(vrgb,ccdpixy,"testxxx.pnm");
    return export_png_fromline(vrgb, ccdpixy, filename);
}



#endif
