# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
# Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland and Adam Goldstein
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Daniel Kocevski
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
import shutil
from pathlib import Path

from setuptools import setup, find_namespace_packages
import gdt.core as core

if __name__ == '__main__':
    setup(
        name="astro-gdt-fermi",
        version=core.__version__,
        description="Gamma-ray Data Tools: Fermi Mission",
        author='Cleveland, Goldstein, Kocevski',
        url='https://github.com/USRA-STI/gdt-fermi',
        packages=find_namespace_packages(where='src', include=["*"]),
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: Apache Software License",
            "Operating System :: MacOS",
            "Operating System :: POSIX :: Linux",
            "Topic :: Scientific/Engineering :: Astronomy",
            "Topic :: Software Development :: Libraries",
        ],
        license_files=['license.txt'],
        keywords=['astronomy', 'gammaray', 'gamma-ray', 'usra'],
        package_dir={"": "src"},
        package_data={
            'gdt.missions.fermi': ['data/McIlwainL_Coeffs.npy']
        },
        include_package_data=True,
        python_requires='>=3.8',
        install_requires=[
            'astro-gdt>=2.0.0',
            'pyproj>=1.9.6',
            'numpy>=1.17.3',
            'scipy>=1.1.0',
            'matplotlib>=3.7.1',
            'astropy>=3.1',
            'healpy>=1.12.4',
            'cartopy>=0.21.1',
        ],
        data_files=[('config', ['config/test_files.urls'])]
    )

    # Copy test_files.url
    src = Path(__file__).parent / 'config' / 'test_files.urls'
    dest = core.data_path / 'fermi-gbm' / 'test_files.urls'
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dest)
