# _WARNING: UNDER CONSTRUCTION, NOT A WORKING PACKAGE_
# Astrophysical Rubble piles Simulation Software

This is a software package for modeling and simulating small planetary bodies, as collections of discrete, rigid-body grains of irregular shape. The DEM heavy-lifting is provided by NVIDIA's [PhysX](https://developer.nvidia.com/technologies/physx) engine, and the self-gravity module requires a [CUDA](https://developer.nvidia.com/category/zone/cuda-zone) enabled GPU to run. ARSS was developed and used on a Windows PC, exclusively, although some care has been taken to allow porting to other platforms. This package ('ARSS_linux') is a port attempt to unix-like environments. An `ARSS_OSX` will also exist some day. The PhysX SDK version >3.0 is required for linux and MacOSX, and the crucial binaries are bundled with ARSS, but you may want to obtain the full PhysX SDK to look at the samples etc.

Contents of this package
+ `ARSS_linux/`  -  The top-level directory
+ `ARSS_linux/README.md` - This file
+ `ARSS_linux/install.md` - Instructions for building and installing the package using the provided scripts. _Try this first._ If the scripts fail, they will hopefully provide helpful error messages. More detailed build instructions are provided in `install.linux.txt`.
+ `ARSS_linux/configure*`, `ARSS_linux/Makefile.*`, `ARSS_linux/install-sh`, `ARSS_linux/missing` - files used by the build system.
+ `ARSS_linux/install.linux.md` - A guide to the various prerequisites and dependencies required to build an ARSS project. Use this to configure your own build environment, if the provided scripts fail.
+ `ARSS_linux/doc/` - Will one day contain documentation and user's guide.
+ `ARSS_linux/common/` - Contains source and header files used in all ARSS projects.
+ `ARSS_linux/extern/` - Contains the headers and binaries of external libraries required to run ARSS projects. Some of these are open source, so I could have potentially saved some disk space by bundling the source and build files. But the amount of space/bandwidth saved is negligible compared with the size of the PhysX binaries anyway.
+ `ARSS_linux/_projectname_/` - Each project in the ARSS family resides in a separate folder. This folder usually contains just the project's source and header files, and the `Makefile`s for the build system. Some projects also have data files and/or MATLAB script files. Most projects have a `_projectname_.ini` file that serves as a template for a run time configuration file.

