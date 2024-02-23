#pragma once
#include "liquiddsp.hpp"

std::map<std::string,liquid_ncotype> nco_type_map =
{
    {"nco", LIQUID_NCO},    
    {"vco", LIQUID_VCO}
};

struct NCO 
{
    nco_crcf    mNCO;
    std::string mType;

    NCO (std::string type) {
        if ( type == "nco" ) { 
            mType = "nco";
            mNCO = nco_crcf_create (LIQUID_NCO);
        }
        else {
            mType = "vco";
            mNCO = nco_crcf_create (LIQUID_VCO);
        }
    }

   ~NCO (void) {
        nco_crcf_destroy (mNCO);
    }

    float frequency (void) {
        return nco_crcf_get_frequency (mNCO);
    }

    void set_frequency (float fr) {
        nco_crcf_set_frequency (mNCO,fr);
    }

    void adjust_frequency (float df) {
        nco_crcf_adjust_frequency (mNCO,df);
    }

    float phase (void) {
        return nco_crcf_get_phase (mNCO);
    }

    void set_phase (float phs) {
        nco_crcf_set_phase (mNCO,phs);
    }

    void adjust_phase (float dphs) {
        nco_crcf_adjust_phase (mNCO,dphs);
    }

    void set_pll_bandwidth (float bw) {
        nco_crcf_pll_set_bandwidth (mNCO, bw);
    }

    void pll_step (float dph) {
        nco_crcf_pll_step (mNCO, dph);
    }

    void print (void) {
        nco_crcf_print (mNCO);
    }

    py::array_t<liquid_float_complex> mix_up (py::array_t<liquid_float_complex> inp) {
        py::array_t<liquid_float_complex> ret(py::len(inp));
        liquid_float_complex *x = array_to_ptr<liquid_float_complex>(inp);
        liquid_float_complex *y = array_to_ptr<liquid_float_complex>(ret);
        nco_crcf_mix_block_up (mNCO,x,y,py::len(inp));
        return ret;
    } 

    py::array_t<liquid_float_complex> mix_down (py::array_t<liquid_float_complex> inp) {
        py::array_t<liquid_float_complex> ret(py::len(inp));
        liquid_float_complex *x = array_to_ptr<liquid_float_complex>(inp);
        liquid_float_complex *y = array_to_ptr<liquid_float_complex>(ret);
        nco_crcf_mix_block_down (mNCO,x,y,py::len(inp));
        return ret;
    } 
};

