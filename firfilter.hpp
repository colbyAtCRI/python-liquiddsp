#pragma once
#include "liquiddsp.hpp"

// Bare Bones roll your own 
class RealFIRFilter
{
public:
    firfilt_rrrf mFilter;

    RealFIRFilter (void) : mFilter(nullptr) {
    }

    RealFIRFilter (py::array_t<float> h) {
        float *hc = array_to_ptr<float>(h);
        mFilter = firfilt_rrrf_create (hc, py::len(h));
    }

   ~RealFIRFilter (void) {
        if (mFilter)
            firfilt_rrrf_destroy (mFilter);
    }

    std::complex<float> freqresponse (float f) {
        std::complex<float> res;
        firfilt_rrrf_freqresponse (mFilter,f,&res);
        return res;
    }

    py::array_t<float> execute (py::array_t<float> inp) {
        py::array_t<float> ret(py::len(inp));
        float *x = array_to_ptr<float>(inp);
        float *y = array_to_ptr<float>(ret);
        firfilt_rrrf_execute_block (mFilter, x, py::len(inp), y);
        return ret;
    }
};

// use inheritance to save typing
class RealDCBlocker : public RealFIRFilter
{
public:
    RealDCBlocker (int m, float as) : RealFIRFilter() {
        mFilter = firfilt_rrrf_create_dc_blocker (m,as);
    }

   ~RealDCBlocker (void) {
        firfilt_rrrf_destroy (mFilter);
        mFilter = nullptr;
    }
};

class RealKaiserBessel : public RealFIRFilter 
{
public:

    RealKaiserBessel (int flen, float fc, float as, float offset) : RealFIRFilter() {
        std::complex<float> res0;
        mFilter = firfilt_rrrf_create_kaiser (flen, fc, as, offset);
        firfilt_rrrf_freqresponse(mFilter,0.0f,&res0);
        firfilt_rrrf_set_scale (mFilter,1.0/abs(res0));
    }

   ~RealKaiserBessel (void) {
        firfilt_rrrf_destroy (mFilter);
        mFilter = nullptr;
    }
};
