#pragma once
#include "liquiddsp.hpp"

// Okay, liquid's ampmodem isn't viable as is for AM broadcast audio 
// when carriers as present. The reason stems from the DC blocker
// applied to the output of the demodulator badly distorts the frequency
// spectrum and, the results vary depending on input sample rate.
// Other than this, ampmodem_demod_dsb_pll_carrier is a fine algorithm.

class BroadcastAM
{
public:
    nco_crcf     mMixer;
    firfilt_crcf mLowpass;
    iirfilt_rrrf mDCBlock;
    wdelaycf     mDelay;

    BroadcastAM (int m) {
        float fc = 20.0f/48000.0f;
        mMixer = nco_crcf_create (LIQUID_NCO);
        nco_crcf_pll_set_bandwidth (mMixer, 0.001f);
        mLowpass = firfilt_crcf_create_kaiser (2*m+1, 0.01f, 40.0f, 0.0f);
        mDCBlock = iirfilt_rrrf_create_prototype (LIQUID_IIRDES_CHEBY2,LIQUID_IIRDES_HIGHPASS,LIQUID_IIRDES_SOS,3,fc,0.0f,0.5f,20.0f);
        mDelay = wdelaycf_create (m);
    }

   ~BroadcastAM (void) {
        nco_crcf_destroy (mMixer);
        firfilt_crcf_destroy (mLowpass);
        wdelaycf_destroy (mDelay);
    }

    void reset (void) {
        nco_crcf_reset (mMixer);
        firfilt_crcf_reset (mLowpass);
        wdelaycf_reset (mDelay);
    }

    py::array_t<float> execute (py::array_t<std::complex<float>> inp) {
        py::array_t<float> ret(py::len(inp));
        std::complex<float> *x = array_to_ptr<std::complex<float>>(inp);
        float *y = array_to_ptr<float>(ret);
        for (auto n = 0; n < py::len(inp); n++)
            y[n] = demod_one (x[n]);
        return ret;
    }

    float demod_one (std::complex<float> x) {
        std::complex<float> x0, x1;
        firfilt_crcf_push (mLowpass, x);
        firfilt_crcf_execute (mLowpass, &x0);
        wdelaycf_push (mDelay, x);
        wdelaycf_read (mDelay, &x1);

        std::complex<float> v0, v1;
        nco_crcf_mix_down (mMixer, x0, &v0);
        nco_crcf_mix_down (mMixer, x1, &v1);

        float phase_error = arg(v0);
        nco_crcf_pll_step (mMixer, phase_error);

        nco_crcf_step (mMixer);

        float outp;
        iirfilt_rrrf_execute (mDCBlock, real(v1), &outp);
        return outp;
    }
};

class SSBDemod
{
public:
    bool     mUSB;
    firhilbf mHilbert;

    SSBDemod (std::string band) {
        mUSB = band == "usb";
        mHilbert = firhilbf_create (25,60.0f);
    }

   ~SSBDemod (void) {
        firhilbf_destroy (mHilbert);
    }

    void reset (void) {
        firhilbf_reset (mHilbert);
    }

    py::array_t<float> execute (py::array_t<std::complex<float>> inp) {
        py::array_t<float> ret(py::len(inp));
        float *y = array_to_ptr<float> (ret);
        std::complex<float> *x = array_to_ptr<std::complex<float>>(inp);
        float junk;
        if (mUSB)
            for (auto n = 0; n < py::len(inp); n++)
                firhilbf_c2r_execute (mHilbert,x[n],&junk,&y[n]);
        else
            for (auto n = 0; n < py::len(inp); n++)
                firhilbf_c2r_execute (mHilbert,x[n],&y[n],&junk);
        return ret;
    }
};

class FreqDem
{
    float   mKd;
    freqdem mModem;
public:

    FreqDem (float kd) {
        mKd = kd;
        mModem = freqdem_create (mKd);
    }

    ~FreqDem (void) {
        freqdem_destroy (mModem);
    }

    void reset (void) {
        freqdem_reset (mModem);
    }

    void print (void) {
        freqdem_print (mModem);
    }

    py::array_t<float> demod (py::array_t<std::complex<float>> inp) {
        py::array_t<float> ret(py::len(inp));
        std::complex<float> *x = array_to_ptr<std::complex<float>>(inp);
        float *y = array_to_ptr<float>(ret);
        freqdem_demodulate_block (mModem, x, py::len(inp), y);
        return ret;
    }
};

std::map<std::string,liquid_ampmodem_type> ampmodem_type_map = 
{
    {"dsb",LIQUID_AMPMODEM_DSB},
    {"usb",LIQUID_AMPMODEM_USB},
    {"lsb",LIQUID_AMPMODEM_LSB}
};

class AmpModem
{
public:

    float       mModulation;
    std::string mType;
    bool        mCarrier;

    ampmodem mModem;


    AmpModem (float mod, std::string type, bool car) : mModulation(mod), mCarrier(car) {
        mModem = makeFromArgs (mod, type, car);
    }

    AmpModem (void) {
    }

   ~AmpModem (void) {
        ampmodem_destroy (mModem);
    }

    void set_type (std::string type) {
        if ( type == "dsb" or type == "usb" or type == "lsb" ) {
            mType = type;
            ampmodem_destroy (mModem);
            mModem = makeFromArgs (mModulation, mType, mCarrier);
        }
    }

    std::string get_type (void) {
        return mType;
    }

    void set_modulation (float mod) {
        mModulation = mod;
        ampmodem_destroy (mModem);
        mModem = makeFromArgs (mModulation, mType, mCarrier);
    }

    float get_modulation (void) {
        return mModulation;
    }

    void set_carrier (bool val) {
        mCarrier = val;
        ampmodem_destroy (mModem);
        mModem = makeFromArgs (mModulation, mType, mCarrier);
    }

    bool get_carrier (void) {
        return mCarrier;
    }

    void print (void) {
        ampmodem_print (mModem);
    }

    void reset (void) {
        ampmodem_reset (mModem);
    }

    py::array_t<float> demod (py::array_t<liquid_float_complex> inp) {
        py::array_t<float> ret(py::len(inp));
        liquid_float_complex *x = array_to_ptr<liquid_float_complex>(inp);
        float *y = array_to_ptr<float>(ret);
        ampmodem_demodulate_block (mModem, x, py::len(inp), y);
        return ret;
    }   
    
    ampmodem makeFromArgs (float mod, std::string type, bool car) {
        liquid_ampmodem_type mt(LIQUID_AMPMODEM_DSB);
        int suppress = (car) ? 0 : 1;
        if ( ampmodem_type_map.find(type) != ampmodem_type_map.end() ) {
            mType = type;
            mt = ampmodem_type_map[type];
        }   
        return ampmodem_create (mod,mt,suppress);
    }
};

