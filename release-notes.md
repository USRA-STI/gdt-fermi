# Release Notes for Gamma-ray Data Tools: Fermi
## Version 2.1.1 (Released Aug 16, 2024)  

This release included the following updates from pull requests:  

- Updated to use the finder API in gdt-core 2.1.0
- Jupiter Notebooks updated and organized [#41](https://github.com/USRA-STI/gdt-fermi/pull/41)
- Handle recent matplotlib and numpy API changes [#42](https://github.com/USRA-STI/gdt-fermi/pull/42)
- Updating to new BaseFinder from gdt-core [#38](https://github.com/USRA-STI/gdt-fermi/pull/38)
- Updated the GbmSaa class to include the new SAA polygon [#36](https://github.com/USRA-STI/gdt-fermi/pull/36)

## Version 2.1.0 (Released Jun 3, 2024)

This release included the following updates from pull requests:  

- Added Fermi GBM DoL (Daughter of LocBurst) [#33](https://github.com/USRA-STI/gdt-fermi/pull/33)
- Fix GbmHealPix multiply [#30](https://github.com/USRA-STI/gdt-fermi/pull/30)
- Update DATE-OBS and DATE-END to use "ISOT" [#29](https://github.com/USRA-STI/gdt-fermi/pull/29)
- More Trigdat fixes [#28](https://github.com/USRA-STI/gdt-fermi/pull/28)
- Fix for GbmHealpix.write() [#23](https://github.com/USRA-STI/gdt-fermi/pull/23)
- Trigdat fix for missing data columns in OB_CALC [#22](https://github.com/USRA-STI/gdt-fermi/pull/22)
- Documentation fixes [#20](https://github.com/USRA-STI/gdt-fermi/pull/20)
- Bug fix for McIlwain L reported by standard_title() [#17](https://github.com/USRA-STI/gdt-fermi/pull/17)
- Added link to documentation [#8](https://github.com/USRA-STI/gdt-fermi/pull/8)
- Copy TTE header info when converting to PHAII [#6](https://github.com/USRA-STI/gdt-fermi/pull/6)

## Version 2.0.0 (Released Apr 12, 2023)

Initial release of this library.

We started the version at 2.0.0 to indicate that this is our major API update from GBM Data Tools which is version 1.1.1

Changes from GBM Data Tools include:

- Making the library more generalized so that it can be used for other gamma-ray observatories besides GBM.
- Functions that is common to all gamma-ray observatories is released as Gamma-ray Data Tools: Core
- Functions that is specific to Fermi-GBM is released as Gamma-ray Data Tools: Fermi
- Mission Elapsed Time is stored within AstroPy Time which handles the conversion to other time systems.
- Spacecraft coordinates is now using Astropy Coordinate Frames.
- Quaternions use Numpy to vectorize operations and SciPy Rotation to translate between Quaternions and Direction Cosine Matrixes.
- Fits files now support context management and perform header verification.

