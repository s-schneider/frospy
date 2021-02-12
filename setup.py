import pathlib
from setuptools import setup
import os

"""
:copyright:
    Simon Schneider, Sujania Talavera-Soza, Lisanne Jagt
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
INSTALL_REQUIRES = [
    'pytorch torchvision -c pytorch',
    'obspy -c conda-forge',
    'basemap cartopy -c conda-forge',
    'mechanize',
    'natsort',
    'numpy',
    'pandas',
    'pandoc',
    'pip',
    'scipy',
    'sqlite',
    'tk=8.6.8']
EXTRAS_REQUIRE = {}

MIN_PYTHON_VERSION = (3, 6)

def find_packages():
    """
    Simple function to find all modules under the current folder.
    """
    modules = []
    for dirpath, _, filenames in os.walk(os.path.join(SETUP_DIRECTORY,
                                                      "frospy")):
        if "__init__.py" in filenames:
            modules.append(os.path.relpath(dirpath, SETUP_DIRECTORY))
    return [_i.replace(os.sep, ".") for _i in modules]


setup(
    name="frospy",
    version="0.8.0",
    description="Python toolbox for normal mode seismologists",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/s-schneider/frospy",
    author="Simon Schneider",
    author_email="s.a.schneider@mailbox.org",
    license='GNU Lesser General Public License, Version 3 (LGPLv3)',
    platforms='OS Independent',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU Library or ' +
            'Lesser General Public License (LGPL)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics'],
    packages=find_packages(),
    python_requires=f'>={MIN_PYTHON_VERSION[0]}.{MIN_PYTHON_VERSION[1]}',
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    # entry_points={
    #     "console_scripts": [
    #         "realpython=reader.__main__:main",
    #     ]
    },
)
