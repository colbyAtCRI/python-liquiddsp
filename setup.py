from distutils.core import setup, Extension
from sys import platform
import pybind11 as pb

if platform == 'darwin':
    lib_path = ['/opt/homebrew/lib']
    inc_path = [pb.get_include(), '/opt/homebrew/include']
    link_lib = ['liquid']
elif platform == 'linux':
    lib_path = ['/usr/local/lib']
    inc_path = [pb.get_include(), '/usr/local/include']
    link_lib = ['liquid', 'm']

liquiddsp  = Extension ( 'liquiddsp',
                    define_macros = [('MAJOR_VERSION',1), ('MINOR_VERSION',0)],
                    include_dirs = inc_path,
                    libraries = link_lib,
                    library_dirs = lib_path,
                    extra_compile_args = ['-std=c++11','-Wno-deprecated'],
                    sources = ['wrapper.cpp','agc_docs.cpp','resampler_doc.cpp'])

setup ( name = 'liquiddsp',
        version = '1.0',
        description = 'liquid-dsp wrapper',
        author = 'Paul Colby',
        author_email = 'paulccolby@earthlink.net',
        long_description = '''
Support for Liquid-dsp and liquid-sdr
''',
        ext_modules = [liquiddsp])
                    
