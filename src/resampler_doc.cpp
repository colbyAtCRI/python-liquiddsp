#include "resampler.hpp"

const char *RealResampler::doc =
R"""(
Arbitrary rate resampler. While mostly used to reduce rates, may 
also be used to increase or interpolate.
)""";

const char *RealResampler::init_doc = 
R"""(
    rate    : fractional rate change.
    len     : mysterious filter length (20)
    Fc      : anti aliasing frequency cutoff
    As      : Stop band level (set to 60.0f dB)
    nfilter : number of filters to cascade (13) 
)""";

const char *RealResampler::reset_doc = 
R"""(
Reset resampler
)""";

const char *RealResampler::print_doc =
R"""(
Print resampler state to stdout
)""";

const char *RealResampler::rate_doc = 
R"""(
Access to rate. Setting replaces underlying liquid object
)""";

const char *RealResampler::execute_doc = 
R"""(
Takes an IQ (for ComplexResampler) or PCM (for RealResampler) applies 
the band pass filter an returns an array. The length of the returned 
array will in general vary from call to call such that the desired 
data rate is achived. 
)""";

