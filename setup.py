import os
import glob
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.sdist import sdist
from setuptools.command.install import install


class MyBuildExt(build_ext):

    def run(self):
        #os.system('')   # add some commands for Cython functions
        build_ext.run(self)


class MySDist(sdist):

    def run(self):
        #os.system('')   # add some commands for Cython functions
        sdist.run(self)


class MyInstall(install):

    def run(self):
        #os.system('')   # add some commands for Cython functions
        install.run(self)


setup (
    name='pyqha',
    version='0.1',
    packages=find_packages(),

    # Declare your packages' dependencies here, for eg:
    install_requires=[
        'numpy', 'scipy', 'matplotlib', 'multiprocessing', # 'wx'
    ],

    # Fill in these to make your Egg ready for upload to
    # PyPI
    author='Mauro Palumbo',
    author_email='mauropalumbo75@gmail.com',

    #summary = 'Just another Python package for the cheese shop',
    url='',
    license='MIT',
    long_description='Post processing tools for the quasi-harmonic approximation',

    cmdclass={
        'build_ext': MyBuildExt,
        'sdist': MySDist,
        'install': MyInstall
    },

    # could also include long_description, download_url, classifiers, etc.
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Physics'
    ]
)
