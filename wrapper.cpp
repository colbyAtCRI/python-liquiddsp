#include "liquiddsp.hpp"
#include "agc.hpp"

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
class IIRfilter_complex
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

    IIRfilter_complex (std::string filter_type, std::string band_type, int order, float fc, float f0, float Ap, float As) {
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

    ~IIRfilter_complex (void) {
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

class IIRfilter_real
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

    IIRfilter_real (std::string filter_type, std::string band_type, int order, float fc, float f0, float Ap, float As) {
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

    ~IIRfilter_real (void) {
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
/*
struct AGC 
{
    AGC (void) {
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

    py::array_t<liquid_float_complex> execute (py::array_t<liquid_float_complex> inp) {
        py::array_t<liquid_float_complex> ret(py::len(inp));
        liquid_float_complex *x = array_to_ptr<liquid_float_complex>(inp);
        liquid_float_complex *y = array_to_ptr<liquid_float_complex>(ret);
        //agc_crcf_execute_block (mAGC, x, py::len(inp), y);
        for (unsigned int n = 0; n < py::len(inp); n++) { 
            agc_crcf_execute (mAGC,x[n],&y[n]);
            if (mSquelch and agc_crcf_squelch_get_status (mAGC) != 3)
                y[n] *= 0.0;
        }
        return ret;
    }

    bool     mSquelch;
    bool     mLock;
    agc_crcf mAGC;
}; 
*/
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

PYBIND11_MODULE (liquiddsp, m)
{

    m.def ("bytes_to_iq", bytes_to_iq);

    py::class_<Delay>(m, "Delay")
        .def (py::init<int>(), py::arg("nd")=1)
        .def_property ("delay", &Delay::get_delay, &Delay::set_delay) 
        .def ("__call__", &Delay::call);

    py::class_<IIRfilter_complex>(m, "ComplexIIRFilter")
        .def (py::init<std::string,std::string,int,float,float,float,float>(), 
            py::arg("filter_type") = "butter",
            py::arg("band_type") = "lowpass",
            py::arg("order") = 2,
            py::arg("Fc") = 0.2f,
            py::arg("F0") = 0.3f,
            py::arg("Ap") = 0.7f,
            py::arg("As") = 60.0)
        .def_readonly ("filter_type", &IIRfilter_complex::mFt)
        .def_readonly ("band_type", &IIRfilter_complex::mBt)
        .def_readonly ("order", &IIRfilter_complex::mOrder)
        .def_readonly ("Fc", &IIRfilter_complex::mFc)
        .def_readonly ("F0", &IIRfilter_complex::mF0)
        .def_readonly ("Ap", &IIRfilter_complex::mAp)
        .def_readonly ("As", &IIRfilter_complex::mAs)
        .def ("freqresponse", &IIRfilter_complex::response)
        .def ("__call__", &IIRfilter_complex::execute)
        .def ("print", &IIRfilter_complex::print);

    py::class_<IIRfilter_real>(m, "RealIIRFilter")
        .def (py::init<std::string,std::string,int,float,float,float,float>(), 
            py::arg("filter_type") = "butter",
            py::arg("band_type") = "lowpass",
            py::arg("order") = 2,
            py::arg("Fc") = 0.2f,
            py::arg("F0") = 0.3f,
            py::arg("Ap") = 0.7f,
            py::arg("As") = 60.0)
        .def_readonly ("filter_type", &IIRfilter_real::mFt)
        .def_readonly ("band_type", &IIRfilter_real::mBt)
        .def_readonly ("order", &IIRfilter_real::mOrder)
        .def_readonly ("Fc", &IIRfilter_real::mFc)
        .def_readonly ("F0", &IIRfilter_real::mF0)
        .def_readonly ("Ap", &IIRfilter_real::mAp)
        .def_readonly ("As", &IIRfilter_real::mAs)
        .def ("freqresponse", &IIRfilter_real::response)
        .def ("__call__", &IIRfilter_real::execute)
        .def ("print", &IIRfilter_real::print);

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

    py::class_<RealResampler>(m, "RealResampler")
        .def (py::init<float,int,float,float,int>(),py::arg("rate"),py::arg("len")=20,py::arg("Fc"),py::arg("As")=60.0f,py::arg("nfilter")=13)
        .def ("print", &RealResampler::print)
        .def ("reset", &RealResampler::reset)
        .def ("__call__", &RealResampler::execute)
        .def_property ("rate", &RealResampler::get_rate, &RealResampler::set_rate);

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
}
