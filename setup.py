import shutil
from setuptools import setup, find_namespace_packages
from gdt.core import data_path
from src.gdt.missions.fermi import __version__
from tests import tests_path

if __name__ == '__main__':
    setup(
        name='gdt-fermi',
        version=__version__,
        packages=find_namespace_packages(where='src/', include=['gdt.missions.*']),
        package_dir={'': 'src'},
        url='github.com/USRA-STI/gdt',
        license='Apache 2.0',
        author='Cleveland, Goldstein, Kocevski',
        description='The Gamma-ray Data Tools',
        package_data={
            'gdt.missions.fermi': ['data/McIlwainL_Coeffs.npy']
        },
        include_package_data=True,
        python_requires='>=3.8',
        install_requires=[
            'gdt-core',
            'pyproj>=1.9.6',
            'numpy>=1.17.3',
            'scipy>=1.1.0',
            'matplotlib',
            'astropy>=3.1',
            'healpy>=1.12.4',
        ],
        extras_require={
            'docs': [
                'Sphinx',
                'astropy_sphinx_theme',
                'nbsphinx',
                'ipython',
                'sphinx_automodapi',
                'notebook'
            ],
            'test': [
                'pytest'
            ]
        }
    )
    # Copy test_files.url
    src = tests_path / 'test_files.urls'
    dest = data_path / 'fermi-gbm' / 'test_files.urls'
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dest)
