
#pragma once
#include "liquiddsp.hpp"

std::map<std::string,liquid_iirdes_filtertype> filter_type_map =
{
    {"butter", LIQUID_IIRDES_BUTTER},
    {"cheby1", LIQUID_IIRDES_CHEBY1},
    {"cheby2", LIQUID_IIRDES_CHEBY2},
    {"ellip",  LIQUID_IIRDES_ELLIP},
    {"bessel", LIQUID_IIRDES_BESSEL}
};

std::map<std::string,liquid_iirdes_bandtype> band_type_map = 
{
    {"lowpass", LIQUID_IIRDES_LOWPASS},
    {"highpass", LIQUID_IIRDES_HIGHPASS},
    {"bandpass", LIQUID_IIRDES_BANDPASS},
    {"bandstop", LIQUID_IIRDES_BANDSTOP}
};

// I'm going to do just LIQUID_IIRDES_SOS, or the one that works
class ComplexIIRFilter
{
public:

    std::string  mFt;
    std::string  mBt;
    int          mOrder;
    float        mFc;
    float        mF0;
    float        mAp;
    float        mAs;

    iirfilt_crcf mFilter;

    ComplexIIRFilter (std::string filter_type, std::string band_type, int order, float fc, float f0, float Ap, float As) {
        mOrder = order;
        mFc = fc;
        mF0 = f0;
        mAp = Ap;
        mAs = As;

        liquid_iirdes_filtertype ft(LIQUID_IIRDES_BUTTER);
        if ( filter_type_map.find(filter_type) != filter_type_map.end() ) {
            mFt = filter_type;
            ft = filter_type_map[filter_type];
        }
        liquid_iirdes_bandtype bt(LIQUID_IIRDES_LOWPASS);
        if ( band_type_map.find(band_type) != band_type_map.end() ) {
            mBt = band_type;
            bt = band_type_map[band_type];
        }
        mFilter = iirfilt_crcf_create_prototype (ft,bt,LIQUID_IIRDES_SOS,order,fc,f0,Ap,As);
    }

    ~ComplexIIRFilter (void) {
        iirfilt_crcf_destroy (mFilter);
    }

    void print (void) {
        iirfilt_crcf_print (mFilter);
    }

    std::complex<float> response (float freq) {
        std::complex<float> ret;
        iirfilt_crcf_freqresponse (mFilter, freq, &ret);
        return ret;
    }

    py::array_t<std::complex<float>> execute (py::array_t<std::complex<float>> inp) {
        py::array_t<std::complex<float>> ret(py::len(inp));
        std::complex<float> *x = array_to_ptr<std::complex<float>>(inp);
        std::complex<float> *y = array_to_ptr<std::complex<float>>(ret);
        iirfilt_crcf_execute_block (mFilter,x,py::len(inp),y);
        return ret;
    }
};

class RealIIRFilter
{
public:

    std::string  mFt;
    std::string  mBt;
    int          mOrder;
    float        mFc;
    float        mF0;
    float        mAp;
    float        mAs;

    iirfilt_rrrf mFilter;

    RealIIRFilter (std::string filter_type, std::string band_type, int order, float fc, float f0, float Ap, float As) {
        mOrder = order;
        mFc = fc;
        mF0 = f0;
        mAp = Ap;
        mAs = As;

        liquid_iirdes_filtertype ft(LIQUID_IIRDES_BUTTER);
        if ( filter_type_map.find(filter_type) != filter_type_map.end() ) {
            mFt = filter_type;
            ft = filter_type_map[filter_type];
        }
        liquid_iirdes_bandtype bt(LIQUID_IIRDES_LOWPASS);
        if ( band_type_map.find(band_type) != band_type_map.end() ) {
            mBt = band_type;
            bt = band_type_map[band_type];
        }
        mFilter = iirfilt_rrrf_create_prototype (ft,bt,LIQUID_IIRDES_SOS,order,fc,f0,Ap,As);
    }

    ~RealIIRFilter (void) {
        iirfilt_rrrf_destroy (mFilter);
    }

    void print (void) {
        iirfilt_rrrf_print (mFilter);
    }

    std::complex<float> response (float freq) {
        std::complex<float> ret;
        iirfilt_rrrf_freqresponse (mFilter, freq, &ret);
        return ret;
    }

    py::array_t<float> execute (py::array_t<float> inp) {
        py::array_t<float> ret(py::len(inp));
        float *x = array_to_ptr<float>(inp);
        float *y = array_to_ptr<float>(ret);
        iirfilt_rrrf_execute_block (mFilter,x,py::len(inp),y);
        return ret;
    }
};

class DeemphasisFilter
{
    float        mB[1];
    float        mA[2];
    iirfilt_rrrf mFilt;

public:

    DeemphasisFilter (float sr) {
        float x = exp(-1.0/(75.0E-6 * sr));
        mA[0] = 1.0;
        mA[1] = -x;
        mB[0] = 1.0 - x;
        mFilt = iirfilt_rrrf_create (mB,1,mA,2);
    }

   ~DeemphasisFilter (void) {
        iirfilt_rrrf_destroy (mFilt);
    }

    std::complex<float> response (float freq) {
        std::complex<float> ret;
        iirfilt_rrrf_freqresponse (mFilt, freq, &ret);
        return ret;
    }

    py::array_t<float> execute (py::array_t<float> data) {
        py::array_t<float> ret(py::len(data));
        float *x = array_to_ptr<float>(data);
        float *y = array_to_ptr<float>(ret);
        for (auto n = 0; n < py::len(data); n++)
            iirfilt_rrrf_execute (mFilt, x[n], &y[n]);
        return ret;
    }
};

