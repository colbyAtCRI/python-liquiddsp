#pragma once
#include "liquiddsp.hpp"

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

