#pragma once
#include "liquiddsp.hpp"

struct AGC 
{
    AGC (void) {
        mOnRise  = py::none();
        mSquelch = false;
        mLock    = false;
        mAGC     = agc_crcf_create ();
    }

   ~AGC (void) {
        agc_crcf_destroy (mAGC);
    }

    void set_bandwidth (float bw) {
        agc_crcf_set_bandwidth (mAGC,bw);
    }

    float get_bandwidth (void) {
        return agc_crcf_get_bandwidth (mAGC);
    }

    bool get_squelch (void) {
        return mSquelch;
    }

    void set_squelch (bool val) {
        mSquelch = val;
        if (mSquelch) 
            agc_crcf_squelch_enable (mAGC);
        else
            agc_crcf_squelch_disable (mAGC);
    }

    void squelch_set_threshold (double t) {
        agc_crcf_squelch_set_threshold (mAGC,t);
    }

    double squelch_get_threshold (void) {
        return agc_crcf_squelch_get_threshold (mAGC);
    }

    float get_level (void) {
        return agc_crcf_get_signal_level (mAGC);
    }

    void set_level (float level) {
        agc_crcf_set_signal_level (mAGC, level);
    }

    float get_rssi (void) {
        return agc_crcf_get_rssi (mAGC);
    }

    void set_rssi (float rssi) {
        agc_crcf_set_rssi (mAGC, rssi);
    }

    bool get_lock (void) {
        return mLock;
    }

    void set_lock (bool val) {
        mLock = val;
        if (mLock) 
            agc_crcf_lock (mAGC);
        else
            agc_crcf_unlock (mAGC);
    }

    float get_gain (void) {
        return agc_crcf_get_gain (mAGC);
    }

    void set_gain (float gain) {
        agc_crcf_set_gain (mAGC, gain);
    }

    float get_scale (void) {
        return agc_crcf_get_scale (mAGC);
    }

    void set_scale (float scale) {
        agc_crcf_set_scale (mAGC, scale);
    }

    int status (void) {
        return agc_crcf_squelch_get_status (mAGC);
    }

    void print (void) {
        agc_crcf_print (mAGC);
    }

    void reset (void) {
        agc_crcf_reset (mAGC);
    }

    py::object get_onRise (void) {
        return mOnRise;
    }

    void set_onRise (py::object clb) {
        mOnRise = clb;
    }

    py::array_t<liquid_float_complex> execute (py::array_t<liquid_float_complex> inp) {
        static int state_last = LIQUID_AGC_SQUELCH_UNKNOWN;
        py::array_t<liquid_float_complex> ret(py::len(inp));
        liquid_float_complex *x = array_to_ptr<liquid_float_complex>(inp);
        liquid_float_complex *y = array_to_ptr<liquid_float_complex>(ret);
        for (unsigned int n = 0; n < py::len(inp); n++) { 
            agc_crcf_execute (mAGC,x[n],&y[n]);
            int state = agc_crcf_squelch_get_status (mAGC);
            if ( state != state_last ) {
                state_last = state;
                if (state == LIQUID_AGC_SQUELCH_RISE) {
                    if ( not mOnRise.is(py::none()) )
                        mOnRise();
                }
            }
            if ( state == LIQUID_AGC_SQUELCH_SIGNALLO or state == LIQUID_AGC_SQUELCH_ENABLED )
                y[n] *= 0.0;
        }
        return ret;
    }

    bool       mSquelch;
    bool       mLock;
    agc_crcf   mAGC;
    py::object mOnRise;

    static const char *doc;
    static const char *threshold_doc;
    static const char *bandwidth_doc;
    static const char *squelch_doc;
    static const char *gain_doc;
    static const char *level_doc;
    static const char *level_dB_doc;
    static const char *status_doc;
    static const char *print_doc;
    static const char *lock_doc;
    static const char *scale_doc;  
    static const char *reset_doc;
    static const char *onRise_doc;
    static const char *execute_doc;
}; 
