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

Each of these use python callbacks of the form

`
def process_iqs (iq):
    # iq is an numpy array of type numpy.complex64
    pass
`
or as a class with the __call__ member defined a notional example

`
from liquiddsp import NCO, Sampler, IIRfilter
class RadioProcess:

    def __init__(self,...):
        # do useful init stuff like
        self.resample = Sampler (rate = 48000/2048000, Fc = 48000/2048000)
        self.am = AmpModem (type='dsb',carrier=True)
        self.tuner = NCO()
        
    def __call__(self, iq):
        pcm = self.resample ( self.am ( self.tuner (iq) ))
        # write pcm.tobytes() to pyaudio stream
`
