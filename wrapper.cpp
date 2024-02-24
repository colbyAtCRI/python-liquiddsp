#include "liquiddsp.hpp"
#include "agc.hpp"
#include "resampler.hpp"
#include "iirfilter.hpp"
#include "firfilter.hpp"
#include "demod.hpp"
#include "nco.hpp"
#include "utility.hpp"

PYBIND11_MODULE (liquiddsp, m)
{

    m.def ("bytes_to_iq", bytes_to_iq);

    py::class_<Delay>(m, "Delay")
        .def (py::init<int>(), py::arg("nd")=1)
        .def_property ("delay", &Delay::get_delay, &Delay::set_delay) 
        .def ("__call__", &Delay::call);

    py::class_<CIIRFilter>(m,"CIIRFilter")
        .def (py::init<py::array_t<float>,py::array_t<float>>(),py::arg("Bc"),py::arg("Ac"))
        .def ("reset", &CIIRFilter::reset)
        .def ("freqresponse", &CIIRFilter::freqresponse)
        .def ("__call__", &CIIRFilter::execute);

    py::class_<CLowpassIIR>(m,"CLowpassIIR")
        .def (py::init<std::string,int,float,float,float>(), 
            py::arg("filter_type")="butter",
            py::arg("order"),
            py::arg("Fc"),
            py::arg("Ap")=0.5f,
            py::arg("As")=20.0f)
        .def ("reset", &CLowpassIIR::reset)
        .def ("freqresponse", &CLowpassIIR::freqresponse)
        .def ("__call__", &CLowpassIIR::execute);

    py::class_<CHighpassIIR>(m,"CHighpassIIR")
        .def (py::init<std::string,int,float,float,float>(), 
            py::arg("filter_type")="butter",
            py::arg("order"),
            py::arg("Fc"),
            py::arg("Ap")=0.5f,
            py::arg("As")=20.0f)
        .def ("reset", &CHighpassIIR::reset)
        .def ("freqresponse", &CHighpassIIR::freqresponse)
        .def ("__call__", &CHighpassIIR::execute);

    py::class_<CBandpassIIR>(m,"CBandpassIIR")
        .def (py::init<std::string,int,float,float,float,float>(), 
            py::arg("filter_type")="butter",
            py::arg("order"),
            py::arg("Fc"),
            py::arg("F0"),
            py::arg("Ap")=0.5f,
            py::arg("As")=20.0f)
        .def ("reset", &CBandpassIIR::reset)
        .def ("freqresponse", &CBandpassIIR::freqresponse)
        .def ("__call__", &CBandpassIIR::execute);

    py::class_<CBandstopIIR>(m,"CBandstopIIR")
        .def (py::init<std::string,int,float,float,float,float>(), 
            py::arg("filter_type")="butter",
            py::arg("order"),
            py::arg("Fc"),
            py::arg("F0"),
            py::arg("Ap")=0.5f,
            py::arg("As")=20.0f)
        .def ("reset", &CBandstopIIR::reset)
        .def ("freqresponse", &CBandstopIIR::freqresponse)
        .def ("__call__", &CBandstopIIR::execute);

    py::class_<ComplexIIRFilter>(m, "ComplexIIRFilter")
        .def (py::init<std::string,std::string,int,float,float,float,float>(), 
            py::arg("filter_type") = "butter",
            py::arg("band_type") = "lowpass",
            py::arg("order") = 2,
            py::arg("Fc") = 0.2f,
            py::arg("F0") = 0.3f,
            py::arg("Ap") = 0.7f,
            py::arg("As") = 60.0)
        .def_readonly ("filter_type", &ComplexIIRFilter::mFt)
        .def_readonly ("band_type", &ComplexIIRFilter::mBt)
        .def_readonly ("order", &ComplexIIRFilter::mOrder)
        .def_readonly ("Fc", &ComplexIIRFilter::mFc)
        .def_readonly ("F0", &ComplexIIRFilter::mF0)
        .def_readonly ("Ap", &ComplexIIRFilter::mAp)
        .def_readonly ("As", &ComplexIIRFilter::mAs)
        .def ("freqresponse", &ComplexIIRFilter::response)
        .def ("__call__", &ComplexIIRFilter::execute)
        .def ("print", &ComplexIIRFilter::print);

    py::class_<RealIIRFilter>(m, "RealIIRFilter")
        .def (py::init<std::string,std::string,int,float,float,float,float>(), 
            py::arg("filter_type") = "butter",
            py::arg("band_type") = "lowpass",
            py::arg("order") = 2,
            py::arg("Fc") = 0.2f,
            py::arg("F0") = 0.3f,
            py::arg("Ap") = 0.7f,
            py::arg("As") = 60.0)
        .def_readonly ("filter_type", &RealIIRFilter::mFt)
        .def_readonly ("band_type", &RealIIRFilter::mBt)
        .def_readonly ("order", &RealIIRFilter::mOrder)
        .def_readonly ("Fc", &RealIIRFilter::mFc)
        .def_readonly ("F0", &RealIIRFilter::mF0)
        .def_readonly ("Ap", &RealIIRFilter::mAp)
        .def_readonly ("As", &RealIIRFilter::mAs)
        .def ("freqresponse", &RealIIRFilter::response)
        .def ("__call__", &RealIIRFilter::execute)
        .def ("print", &RealIIRFilter::print);

    py::class_<HilbertTransform>(m, "HilbertTransform") 
        .def (py::init<int,float>(),py::arg("m")=5,py::arg("As")=60.0f)
        .def ("__call__", &HilbertTransform::call);

    py::class_<DeemphasisFilter>(m,"DeemphasisFilter")
        .def (py::init<float>(), py::arg("sample_rate")=48000) 
        .def ("freqresponse", &DeemphasisFilter::response)
        .def ("__call__", &DeemphasisFilter::execute);

    py::class_<FreqDem>(m,"FreqDem")
        .def (py::init<float>())
        .def ("reset", &FreqDem::reset)
        .def ("print", &FreqDem::print)
        .def ("__call__", &FreqDem::demod);

    py::class_<AmpModem>(m,"AmpModem")
        .def (py::init<float,std::string,bool>(),
            py::arg("modulation")=0.75, 
            py::arg("type")="dsb", 
            py::arg("carrier")=false)
        .def_property ("modulation", &AmpModem::get_modulation, &AmpModem::set_modulation)
        .def_property ("type", &AmpModem::get_type,&AmpModem::set_type)
        .def_property ("carrier", &AmpModem::get_carrier,&AmpModem::set_carrier)
        .def ("print", &AmpModem::print)
        .def ("reset", &AmpModem::reset)
        .def ("__call__", &AmpModem::demod);

    py::class_<NCO>(m,"NCO")
        .def (py::init<std::string>(), py::arg("type") = "nco")
        .def ("print", &NCO::print)
        .def_property ("freq", &NCO::frequency, &NCO::set_frequency)
        .def ("adjust_frequency", &NCO::adjust_frequency)
        .def ("adjust_phase", &NCO::adjust_phase)
        .def_property ("phase", &NCO::phase, &NCO::set_phase)
        .def ("set_pll_bandwidth", &NCO::set_pll_bandwidth)
        .def ("pll_step", &NCO::pll_step)
        .def ("__call__", &NCO::mix_up)
        .def ("mix_up", &NCO::mix_up)
        .def ("mix_down", &NCO::mix_down);

    py::class_<RealResampler>(m, "RealResampler",RealResampler::doc)
        .def (py::init<float,int,float,float,int>(),py::arg("rate"),py::arg("len")=20,py::arg("Fc"),py::arg("As")=60.0f,py::arg("nfilter")=13,RealResampler::init_doc)
        .def ("print", &RealResampler::print,RealResampler::print_doc)
        .def ("reset", &RealResampler::reset,RealResampler::reset_doc)
        .def ("__call__", &RealResampler::execute,RealResampler::execute_doc)
        .def_property ("rate", &RealResampler::get_rate, &RealResampler::set_rate,RealResampler::rate_doc);

    py::class_<ComplexResampler>(m,"ComplexResampler")
        .def (py::init<float,int,float,float,int>(),py::arg("rate"),py::arg("len")=20,py::arg("Fc"),py::arg("As")=60.0f,py::arg("nfilter")=13)
        .def ("print",&ComplexResampler::print) 
        .def ("reset", &ComplexResampler::reset)
        .def ("__call__", &ComplexResampler::execute)
        .def_property ("rate", &ComplexResampler::get_rate, &ComplexResampler::set_rate);

    py::class_<AGC>(m,"AGC",AGC::doc)
        .def (py::init())
        .def_property ("squelch", &AGC::get_squelch, &AGC::set_squelch,AGC::squelch_doc)
        .def_property ("threshold", &AGC::squelch_get_threshold, &AGC::squelch_set_threshold,AGC::threshold_doc)
        .def_property ("bandwidth", &AGC::get_bandwidth, &AGC::set_bandwidth,AGC::bandwidth_doc)
        .def_property ("level", &AGC::get_level, &AGC::set_level,AGC::level_doc)
        .def_property ("level_dB", &AGC::get_rssi, &AGC::set_rssi,AGC::level_dB_doc)
        .def_property ("lock", &AGC::get_lock, &AGC::set_lock,AGC::lock_doc)
        .def_property ("gain", &AGC::get_gain, &AGC::set_gain,AGC::gain_doc)
        .def_property ("scale", &AGC::get_scale, &AGC::set_scale,AGC::scale_doc)
        .def_property_readonly ("status", &AGC::status,AGC::status_doc)
        .def_property ("onRise", &AGC::get_onRise, &AGC::set_onRise,AGC::onRise_doc)
        .def ("print", &AGC::print,AGC::print_doc)
        .def ("reset", &AGC::reset,AGC::reset_doc)
        .def ("__call__", &AGC::execute,AGC::execute_doc);

    py::class_<RealFIRFilter>(m, "RealFIRFilter")
        .def (py::init<py::array_t<float>>())
        .def ("freqresponse", &RealFIRFilter::freqresponse)
        .def ("__call__", &RealFIRFilter::execute);

    py::class_<RealDCBlocker>(m,"RealDCBlocker")
        .def (py::init<int,float>(),py::arg("slen")=25,py::arg("As")=20.0f)
        .def ("freqresponse", &RealDCBlocker::freqresponse)
        .def ("__call__", &RealDCBlocker::execute);
        
    py::class_<RealKaiserBessel>(m,"RealKaiserBessel")
        .def (py::init<int,float,float,float>(),py::arg("flen")=25,py::arg("Fc"),py::arg("As")=20.0f,py::arg("offset")=0.0f)
        .def ("freqresponse", &RealKaiserBessel::freqresponse)
        .def ("__call__", &RealKaiserBessel::execute);

    py::class_<BroadcastAM>(m,"BroadcastAM")
        .def (py::init<int>(),py::arg("slen")=25)
        .def ("reset", &BroadcastAM::reset)
        .def ("__call__", &BroadcastAM::execute);
}

