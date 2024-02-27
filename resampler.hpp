#pragma once
#include "liquiddsp.hpp"

class RResampler
{
    float       mRate;
    resamp_rrrf mSampler;
public:

    RResampler (float rate) : mRate(rate) 
    {
        mSampler = resamp_rrrf_create_default (mRate);
    }

   ~RResampler (void) 
    {
        resamp_rrrf_destroy (mSampler);
    }

    void reset (void) 
    {
        resamp_rrrf_reset (mSampler);
    }

    py::array_t<float> execute (py::array_t<float> inp) 
    {
        int no = (int)(py::len(inp)*mRate) + 1;
        float y[no];
        float *x = array_to_ptr<float> (inp);
        unsigned int nw;
        resamp_rrrf_execute_block (mSampler,x,py::len(inp),y,&nw);
        py::array_t<float> ret(nw);
        float *z = array_to_ptr<float>(ret);
        std::copy (y,y+nw,z);
        return ret;
    }
};

class CResampler 
{
    float       mRate;
    resamp_crcf mSampler;
public:

    CResampler (float rate) : mRate(rate)
    {
        mSampler = resamp_crcf_create_default (mRate);
    }

   ~CResampler (void) {
        resamp_crcf_destroy (mSampler);
    }

    void reset (void) {
        resamp_crcf_reset (mSampler);
    }

    py::array_t<std::complex<float>> execute (py::array_t<std::complex<float>> inp) {
        int no = (int)(py::len(inp)*mRate)+1;
        std::complex<float> y[no];
        std::complex<float> *x = array_to_ptr<std::complex<float>> (inp);
        unsigned int nw;
        resamp_crcf_execute_block (mSampler, x, py::len(inp), y, &nw);
        py::array_t<std::complex<float>> ret(nw);
        std::complex<float> *z = array_to_ptr<std::complex<float>> (ret);
        std::copy (y,y+nw,z);
        return ret;
    }

};

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
