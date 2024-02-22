# python-liquiddsp
A vectorized pybind11 interface for liquid-dsp. 

# Install
Provides a python module liquiddsp for python 3. To build one needs 

- liquid-dsp (brew install)
- pybind11   (pip install)
- numpy      (pip install)
- setuptools (pip install)
- wheel      (pip install)

Follow the typical steps,

```
git clone https://github.com/colbyAtCRI/python-liquiddsp.git

cd python-liquiddsp

```
either
```
make
```
or do
```
pip install -e .
```
# Background

I've written several python interfaces for SDRs

- Rfspace's SDRIQ
- SDRPlay RSP2 and RSPduo (all are supported but I've only tested what I have)
- RTLSDRs 

Each of these use python callbacks of the form of a class with 
the `__call__` member defined. A typical example of an AM radio reciever 
is

```
import numpy as np
from liquiddsp import ComplexResampler, ComplexIIRFilter, AmpModem, AGC, DeemphasisFilter
from sdrplay import Radio

class AMRadio:

    def __init__(self, bandwidth = 15000, iq_rate = 2000000, pcm_rate = 48000):
        self.bandpass = ComplexIIRFilter (filter_type='cheby2', order = 8, Fc = bandwidth/iq_rate) 
        self.resample = ComplexRampler (rate = pcm_rate/iq_rate, Fc = pcm_rate/iq_rate)
        self.am = AmpModem (modulation = 0.5, type = 'dsb', carrier = True)
        self.audio_filter = DeemphasisFilter (pcm_rate)
        self.agc = AGC()
        self.agc.lock = False
        self.agc.scale = 0.01
        self.pcm = b''

    def __call__(self, iq):
        pcm = self.audio_filter ( self.am ( self.agc ( self.resample ( self.bandpass ((iq))))))
        self.pcm = self.pcm + pcm.tobytes()
        while len(self.pcm) > 4096:
            # write self.pcm[0:4096] to pyaudio stream
            self.pcm = self.pcm[4096:]

radio = Radio ()
radio.freq = 1170000
radio.onIQData = AMRadio()
radio.running = True

```
