#pragma once
#include "liquiddsp.hpp"

class RealResampler
{
    float       mRate;
    resamp_rrrf mResamp;

public:

    RealResampler (float r, int d, float fc, float sbsp, int nf) {
        mRate = r;
        mResamp = resamp_rrrf_create (r,d,fc,sbsp,nf);
    }

    ~RealResampler (void) {
        resamp_rrrf_destroy (mResamp);
    }

    void reset (void) {
        resamp_rrrf_reset (mResamp);
    }

    float get_rate (void) {
        return mRate;
    }

    void set_rate (float r) {
        mRate = r;
        resamp_rrrf_set_rate (mResamp,mRate);
    }

    void print (void) {
        resamp_rrrf_print(mResamp);
    }

    py::array_t<float> execute (py::array_t<float> inp) {
        float *x = array_to_ptr<float>(inp);
        float y[py::len(inp)];
        unsigned int nd(0), nw(0); 
        for (auto n = 0; n < py::len(inp); n++) {
            resamp_rrrf_execute (mResamp, x[n], &y[nw], &nd);
            nw += nd;
        }
        py::array_t<float> ret(nw);
        float *z = array_to_ptr<float>(ret);
        std::copy (y, y+nw, z);
        return ret;    
    }

    static const char *doc;
    static const char *init_doc;
    static const char *reset_doc;
    static const char *rate_doc;
    static const char *print_doc;
    static const char *execute_doc;
};

class ComplexResampler
{
    float       mRate;
    resamp_cccf mResamp;

public:

    ComplexResampler (float r, int d, float fc, float sbsp, int nf) {
        mRate = r;
        mResamp = resamp_cccf_create (r,d,fc,sbsp,nf);
    }

    ~ComplexResampler (void) {
        resamp_cccf_destroy (mResamp);
    }

    void reset (void) {
        resamp_cccf_reset (mResamp);
    }

    float get_rate (void) {
        return mRate;
    }

    void set_rate (float r) {
        mRate = r;
        resamp_cccf_set_rate (mResamp,mRate);
    }

    void print (void) {
        resamp_cccf_print(mResamp);
    }

    py::array_t<std::complex<float>> execute (py::array_t<std::complex<float>> inp) {
        std::complex<float> *x = array_to_ptr<std::complex<float>>(inp);
        std::complex<float> y[py::len(inp)];
        unsigned int nd(0), nw(0); 
        for (auto n = 0; n < py::len(inp); n++) {
            resamp_cccf_execute (mResamp, x[n], &y[nw], &nd);
            nw += nd;
        }
        py::array_t<std::complex<float>> ret(nw);
        std::complex<float> *z = array_to_ptr<std::complex<float>>(ret);
        std::copy (y, y+nw, z);
        return ret;    
    }
};
