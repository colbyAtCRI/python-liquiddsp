#include "agc.hpp"

const char *AGC::doc = 
R"""(
AGC provides both scaling and squelch for complex IQ data stream.
)""";

const char *AGC::threshold_doc = 
R"""(
Squelch trigger level in dB.
)""";

const char *AGC::bandwidth_doc = 
R"""(
Adjusts AGC settling time.
)""";

const char *AGC::squelch_doc = 
R"""(
Enable and disable squelch. When enabled the IQ samples are zeroed
on status == 3.
)""";

const char *AGC::gain_doc = 
R"""(
)""";

const char *AGC::level_doc = 
R"""(
)""";

const char *AGC::level_dB_doc = 
R"""(
)""";

const char *AGC::status_doc = 
R"""(
)""";

const char *AGC::print_doc = 
R"""(
)""";

const char *AGC::lock_doc = 
R"""(
)""";

const char *AGC::scale_doc = 
R"""(
Linear output scale. This property determins the output level
of AGC. The default output scale is 1.0 or 0dBm. Setting to
0.01 the target output level is -20dBm.
)"""; 

const char *AGC::reset_doc = 
R"""(
Reset AGC to default settings. Note this will cancel lock and
squelch in the process.
)""";
