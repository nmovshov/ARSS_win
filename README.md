# _WARNING: UNDER CONSTRUCTION, NOT A WORKING PACKAGE_
# Astrophysical Rubble piles Simulation Software

This is a software package for modeling and simulating small planetary bodies as collections of discrete, rigid-body grains of irregular shape. The DEM heavy-lifting is provided by NVIDIA's [PhysX](https://developer.nvidia.com/technologies/physx) engine, and the self-gravity module requires a [CUDA](https://developer.nvidia.com/category/zone/cuda-zone) enabled GPU to run. ARSS was developed and used on a Windows PC, exclusively, although some care has been taken to allow porting to other platforms. The PhysX SDK version >3.0 is required for linux and MacOSX, and the crucial binaries are bundled with ARSS, but you may want to obtain the full PhysX SDK to look at the samples etc.

Contents of this package
+ `ARSS_win/`  -  The top-level directory
+ `ARSS_win/README.md` - This file
+ `ARSS_win/install.md` - A brief guide to what needs to happen before you can attempt to build an ARSS project.
+ `ARSS_win/ARSS.sln` - This is a Visual Studio 2010 solution file. If you have VS2010 Pro or better (no Express, sorry, get a free copy at Microsoft's Dreamspark if you're a college student) and have installed NSIGHT VS Edition, you should be able to double click on this file to load all the projects and build with one click, after installing external dependencies (see below). A "solution" is a Visual Studio concept for organizing multiple "projects" that share code. A "project" is a set of source files that will compile and link with dependencies to create an executable. ARSS contains several projects under the ARSS solution to maximize code reuse. But the projects all create standalone executables that never call each other. For non VS users the solution concept and ARSS.sln file are irrelevant. Build each project independently.
+ `ARSS_win/common/` - Contains source and header files used in all ARSS projects.
+ `ARSS_win/extern/` - Contains the headers and binaries of external libraries required to run ARSS projects. Some of these are open source, so I could have potentially saved some disk space by bundling the source and build files. But the amount of space/bandwidth saved is negligible compared with the size of the PhysX binaries anyway.
+ `ARSS_win/_projectname_/` - Each project in the ARSS family resides in a separate folder. This folder usually contains just the project's source and header files, and the Visual Studio `.proj` file. Some projects also have data files and/or MATLAB script files. Most projects have a `_projectname_.ini` file that serves as a template for a run time configuration file.

