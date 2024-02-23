
#pragma once
#include "liquiddsp.hpp"

class Delay 
{
    int      mDelay;
    wdelayf  mRealDelay;
    wdelaycf mComplexDelay;
public:

    Delay (int nd) {
        mDelay = nd;
        mRealDelay = wdelayf_create (nd);
        mComplexDelay = wdelaycf_create (nd);
    }

   ~Delay (void) {
        wdelayf_destroy (mRealDelay);
        wdelaycf_destroy (mComplexDelay);
    }

    int get_delay (void) {
        return mDelay;
    }

    void set_delay (int nd) {
        wdelayf_destroy (mRealDelay);
        wdelaycf_destroy (mComplexDelay);
        mDelay = nd;
        mRealDelay = wdelayf_create (nd);
        mComplexDelay = wdelaycf_create (nd);
    }
 
    py::object call (py::object inp) {
        auto tp = py::getattr (inp,"dtype");
        if ( tp.is(py::dtype ("complex64")) ) {
            py::array_t<std::complex<float>> ret(py::len(inp));
            std::complex<float> *z = array_to_ptr<std::complex<float>>(inp);
            std::complex<float> *y = array_to_ptr<std::complex<float>>(ret);
            for (auto n = 0; n < (int)py::len(inp); n++) {
                wdelaycf_read (mComplexDelay, &y[n]);
                wdelaycf_push (mComplexDelay, z[n]);
            }
            return ret;
        }
        if ( tp.is(py::dtype("float32")) ) {
            py::array_t<float> ret(py::len(inp));
            float *z = array_to_ptr<float>(inp);
            float *y = array_to_ptr<float>(ret);
            for (auto n = 0; n < (int)py::len(inp); n++) {
                wdelayf_read (mRealDelay, &y[n]);
                wdelayf_push (mRealDelay,  z[n]);
            }
            return ret;
        }
        return py::none();
    }
};

py::array_t<liquid_float_complex> bytes_to_iq ( py::bytes byts )
{
    std::string data = py::cast<std::string>(byts);
    std::complex<short> *x = (std::complex<short>*) data.c_str();
    std::vector<std::complex<float>> y(data.size()/4);
    for (unsigned int n = 0; n < data.size()/4; n++) 
        y[n] = std::complex<float> ((float)x[n].real()/32767.0f,(float)x[n].imag()/32767.0f);
    return py::cast(y);
}

class HilbertTransform 
{
    firhilbf mHt_c2r;
    firhilbf mHt_r2c;

public:

    HilbertTransform (int m, float As) {
        mHt_c2r = firhilbf_create (m, As);
        mHt_r2c = firhilbf_create (m, As);
    }

   ~HilbertTransform (void) {
        firhilbf_destroy (mHt_c2r);
        firhilbf_destroy (mHt_r2c);
    }

    py::object call (py::object inp) {
        auto tp = py::getattr (inp,"dtype");
        if ( tp.is(py::dtype ("complex64")) ) {
            py::array_t<float> ret(py::len(inp));
            std::complex<float> *z = array_to_ptr<std::complex<float>>(inp);
            float *y = array_to_ptr<float>(ret);
            for (auto n = 0; n < (int)py::len(inp); n++)
                firhilbf_interp_execute (mHt_c2r, z[n], &y[n]);
            return ret;
        }
        if ( tp.is(py::dtype("float32")) ) {
            py::array_t<std::complex<float>> ret(py::len(inp));
            float *z = array_to_ptr<float>(inp);
            std::complex<float> *y = array_to_ptr<std::complex<float>>(ret);
            for (auto n = 0; n < (int)py::len(inp); n++)
                firhilbf_decim_execute (mHt_r2c, &z[n], &y[n]);
            return ret;
        }
        return py::none();
    }
};

