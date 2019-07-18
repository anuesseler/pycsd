from setuptools     import find_packages
from numpy.distutils.core import Extension, setup

ext1 = Extension(name='zuncsd',sources=['src/zuncsd.pyf', 'src/zuncsd.f'])
setup(
    name='pycsd',
    version='0.1',
    author='Alexander Nuesseler',
    py_modules=['pycsd'],
    description="Wrapper for LAPACK zuncsd function to compute cosine-sine (CS) decomposition of a partitioned unitary matrix.",
    ext_modules = [ext1]    
)