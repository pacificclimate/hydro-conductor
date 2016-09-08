from setuptools import setup

__version__ = '0.1-dev'

setup(
    name="hydro-conductor",
    description="",
    keywords="science climate hydrology glacier modelling",
    packages=['conductor'],
    version=__version__,
    url="http://www.pacificclimate.org/",
    author="Michael Fischer",
    author_email="mfischer@uvic.ca",
    install_requires = ['numpy', 'netCDF4'],
    tests_require = ['pytest', 'mock'],
    scripts = ['scripts/vic_rgm_conductor.py'],
    package_data = {'conductor': ['tests/input/global.txt',
                                  'tests/input/snow_band.txt',
                                  'tests/input/veg.txt',
				  'tests/input/snb_toy_64px.txt',
				  'tests/input/vpf_toy_64px.txt',
				  'tests/input/rgm_vic_map_toy_64px_auto.txt']},
    zip_safe=True,
        classifiers=[
            'Environment :: Console',
            'Natural Language :: English',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPLv3)',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Atmospheric Science',
            'Topic :: Software Development :: Libraries :: Python Modules'
        ]
)
