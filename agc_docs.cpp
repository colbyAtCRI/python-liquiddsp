#include "agc.hpp"

const char *AGC::doc = 
R"""(
AGC provides both scaling and squelch for complex IQ data stream.
)""";

const char *AGC::execute_doc = 
R"""(
Takes an numpy.array(...,numpy.complex64), iq, of IQ samples and returns

    iq -> self.scale * self.gain * iq

when the squelch is disabled. When squelch is enabled, entering states
LIQUID_AGC_SQUELCH_ENABLED or LIQUID_AGC_SQUELCH_SIGNALLO, cause the 
output to be zeroed,

    iq -> 0

For all other states the gain and scaled iqs are returned.  
)""";

const char *AGC::threshold_doc = 
R"""(
Squelch trigger level in dB.
)""";

const char *AGC::bandwidth_doc = 
R"""(
Adjusts AGC settling time of the gain update loop filter.
)""";

const char *AGC::squelch_doc = 
R"""(
Enable and disable squelch. 
)""";

const char *AGC::gain_doc = 
R"""(
The current linear gain. This may be set by the user to reduce gain loop settling time.
)""";

const char *AGC::level_doc = 
R"""(
The current input linear level. This may be set by the user to reduce gain loop settling time.
)""";

const char *AGC::level_dB_doc = 
R"""(
The current input level in dB.
)""";

const char *AGC::status_doc = 
R"""(
The AGC squelch logic has the following states which are returned by this function,

    LIQUID_AGC_SQUELCH_UNKNOWN   (0)
    LIQUID_AGC_SQUELCH_ENABLED,  (1)
    LIQUID_AGC_SQUELCH_RISE,     (2)
    LIQUID_AGC_SQUELCH_SIGNALHI, (3)
    LIQUID_AGC_SQUELCH_FALL,     (4)
    LIQUID_AGC_SQUELCH_SIGNALLO, (5)
    LIQUID_AGC_SQUELCH_TIMEOUT,  (6)
    LIQUID_AGC_SQUELCH_DISABLED  (7)

Since the call to AGC is vectorized, the status is polled for each IQ sample.
For a status of LIQUID_AGC_SQUELCH_ENABLED or LIQUID_AGC_SQUELCH_SIGNALLO
the IQ sample is zeroed. Otherwise, the gain and scale factors are applied.

On transition to LIQUID_AGC_SQUELCH_RISE, the callback onRise is called if
provided. The user may use this call to reset or initialize the demodulator. 
)""";

const char *AGC::print_doc = 
R"""(
Print a summary of the AGC state to stdout.
)""";

const char *AGC::lock_doc = 
R"""(
When True, The AGC is prevented from updating its gain. The AGC does continues
to update its input power levels.  
)""";

const char *AGC::scale_doc = 
R"""(
Linear output scale. This property sets the output level scaling of AGC. The 
default output scale is 1.0 or 0dBm. Setting to 0.01 the target output level 
is -20dBm. Setting scale permits the user to choose the power output of the 
AGC.
)"""; 

const char *AGC::onRise_doc = 
R"""(
Function of no arguments called when the squelch state transitions to 
LIQUID_AGC_SQUELCH_RISE. This call may be used to reset demodulators. 
)""";

const char *AGC::reset_doc = 
R"""(
Reset AGC to default settings. Note this will cancel lock and
squelch in the process.
)""";
