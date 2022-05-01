
#ifndef HELPERS_H
#define HELPERS_H

#include <iostream>
#include <fstream>
#include <complex>
#include <math.h>
#include <assert.h>
#include <chrono>
#include <vector>
#include <algorithm>
#include "lodepng.h"

typedef double real_type;
typedef std::complex<real_type> complex_type;
complex_type imagi(0, 1);

typedef std::vector<real_type> realfield_type;
typedef std::vector<complex_type> complexfield_type;
typedef std::vector<std::vector<real_type>> realfield2d_type;

// values of 0..1 for the 3 color channels!
void export_pnm(realfield2d_type r, realfield2d_type g, realfield2d_type b, std::string filename) {
    std::cout << "writing image " << filename << "\n";
    std::ofstream ofs;
    ofs.open(filename, std::ios::binary);
    ofs << "P6\n"
        << r.size() << " " << r[0].size() << "\n255\n";
    for (int yi = r[0].size()-1; yi >= 0; yi--) {
        for (int xi = 0; xi < r.size(); xi++) {
            ofs << static_cast<char>(std::min(255.0,r[xi][yi]*255)) << static_cast<char>(std::min(255.0,g[xi][yi]*255)) << static_cast<char>(std::min(255.0,b[xi][yi]*255));
        }
    }
    ofs.close();
}

// now for (ccdpixx,3) rgb vectors, duplicate for ccdpixy pixels...
void export_pnm_line2pic(realfield2d_type rgb, int ccdpixy, std::string filename) {
    std::cout << "writing image " << filename << "\n";
    std::ofstream ofs;
    ofs.open(filename, std::ios::binary);
    ofs << "P6\n"
        << rgb.size() << " " << ccdpixy << "\n255\n";
    for (int yi = 0; yi < ccdpixy; yi++) {
        for (int xi = 0; xi < rgb.size(); xi++) {
            ofs << static_cast<char>(std::min(255.0,rgb[xi][0]*255)) << static_cast<char>(std::min(255.0,rgb[xi][1]*255)) << static_cast<char>(std::min(255.0,rgb[xi][2]*255));
        }
    }
    ofs.close();
}

// https://github.com/lvandeve/lodepng
std::vector<unsigned char> export_png_lode(std::vector<unsigned char> image, int pixx, int pixy, std::string filename) {
    if (filename.size()>0) {
        unsigned result = lodepng::encode(filename, image, pixx, pixy);
        if(result) std::cout << "encoder error " << result << ": "<< lodepng_error_text(result) << std::endl;
    }
    std::vector<unsigned char> png;
    unsigned result = lodepng::encode(png, image, pixx, pixy);
    if(result) std::cout << "encoder error " << result << ": "<< lodepng_error_text(result) << std::endl;

    std::cout << "encoded! fn=" << filename << std::endl;
    
    return png;
}

std::vector<unsigned char> export_png(realfield2d_type r, realfield2d_type g, realfield2d_type b, std::string filename) {
    int pixx=r.size();
    int pixy=r[0].size();
    std::vector<unsigned char> image;
    image.resize(pixx * pixy * 4);
    for (int yi = 0; yi < pixy; yi++) {
        for (int xi = 0; xi < pixx; xi++) {
            image[4 * pixx * yi + 4 * xi + 0] = static_cast<char>(std::min(255.0,r[xi][pixy-yi-1]*255));
            image[4 * pixx * yi + 4 * xi + 1] = static_cast<char>(std::min(255.0,g[xi][pixy-yi-1]*255));
            image[4 * pixx * yi + 4 * xi + 2] = static_cast<char>(std::min(255.0,b[xi][pixy-yi-1]*255));
            image[4 * pixx * yi + 4 * xi + 3] = 255;
        }
    }
    return export_png_lode(image, pixx, pixy, filename);
}

std::vector<unsigned char> export_png_fromline(realfield2d_type rgb, int pixy, std::string filename) {
    int pixx=rgb.size();
    std::vector<unsigned char> image;
    image.resize(pixx * pixy * 4);
    for (int yi = 0; yi < pixy; yi++) {
        for (int xi = 0; xi < pixx; xi++) {
            for (int ci=0; ci<4; ci++) {
                image[4 * pixx * yi + 4 * xi + ci] = (ci==3) ? 255 : static_cast<char>(std::min(255.0,rgb[xi][ci]*255));
            }
        }
    }
    return export_png_lode(image, pixx, pixy, filename);
}


// based on https://github.com/razanskylab/wavelength2rgb/blob/master/wavelength2color.m
real_type adjust(real_type inputVal, real_type factor, real_type gammaVal) {
    return pow(inputVal * factor,gammaVal);
}
realfield_type wavelength2color(real_type wavelength) {
    real_type maxIntensity = 1;
    real_type gammaVal = 1;//0.8;

    real_type r,g,b;
    // real_type r = 0.0;
    if ((wavelength >= 380) && (wavelength < 440)) {
        r = -(wavelength - 440) / (440 - 380);
        g = 0;
        b = 1;
    } else if ((wavelength >= 440) && (wavelength < 490)) {
        r = 0;
        g = (wavelength - 440) / (490 - 440);
        b = 1;
    } else if ((wavelength >= 490) && (wavelength < 510)) {
        r = 0;
        g = 1;
        b = -(wavelength - 510) / (510 - 490);
    } else if ((wavelength >= 510) && (wavelength < 580)) {
        r = (wavelength - 510) / (580 - 510);
        g = 1;
        b = 0;
    } else if ((wavelength >= 580) && (wavelength < 645)) {
        r = 1;
        g = -(wavelength - 645) / (645 - 580);
        b = 0;
    } else if ((wavelength >= 645) && (wavelength < 780)) {
        r = 1;
        g = 0;
        b = 0;
    } else {
        r = 0;
        g = 0;
        b = 0;
    }

    real_type factor=0.0;
    if ((wavelength >= 380) && (wavelength < 420))
        factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380);
    else if ((wavelength >=  420) && (wavelength < 700))
        factor = 1;
    else if ((wavelength >= 700) && (wavelength < 780))
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700);

    r = adjust(r, factor, gammaVal);
    g = adjust(g, factor, gammaVal);
    b = adjust(b, factor, gammaVal);

    auto res=realfield_type({r,g,b});
    return res;
    // colorCode = colorCode * maxIntensity;
    // % colorCode = uint8(colorCode * maxIntensity * 255);

}


realfield_type blur(realfield_type ax, realfield_type field, real_type radius) {
    int rpix=floor(radius/(ax[2]-ax[1]));
    if (rpix < 2 || field.size() < 5*rpix) {
        return field;
    } else {
        realfield_type res(field.size(),0);
        for (int i=rpix+1;i<ax.size()-rpix-1;i++) {
            real_type sum=0.0;
            for (int j=i-rpix;j<i+rpix;j++) sum += field[j];
            res[i]=sum/(2*rpix);
        }
        return res;
    }
}


realfield_type int_normalize(realfield_type intf) {
    real_type intamax=0.0;
    for (int i=0; i<intf.size(); i++) if (intf[i]>intamax) intamax=intf[i];
    for (int i=0; i<intf.size(); i++) intf[i]=intf[i]/intamax;
    return intf;
}

real_type sinc(real_type x) {
    if (x==0)
        return 1.0;
    else
        return sin(x)/x;
}

inline auto timer_start() {
    auto t = std::chrono::high_resolution_clock::now();
    return t;
}

inline auto timer_printrestart(std::chrono::time_point<std::chrono::system_clock> t, std::string title) {
    auto t2 = timer_start();
    std::cout << title << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t).count() << "ms.\n";
    t = timer_start();
}





#endif