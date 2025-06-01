idock
=====

idock is a standalone tool for structure-based [virtual screening] powered by fast and flexible ligand docking. It was inspired by [AutoDock Vina], and is hosted on GitHub at https://GitHub.com/HongjianLi/idock under [Apache License 2.0]. idock is also available as a web server at [istar].


Features
--------

* idock invents a lightweight and efficient io service pool to reuse threads and maintain a high CPU utilization throughout the entire docking procedure.
* idock recognizes as many as 28 chemical elements, i.e. H, C, N, O, S, Se, P, F, Cl, Br, I, Zn, Fe, Mg, Ca, Mn, Cu, Na, K, Hg, Ni, Co, Cd, As, Sr, U, Cs, Mo.
* idock outputs atomwise free energy values for subsequent identification of intermolecular interaction hotspots.
* idock writes the scores of each predicted conformation into a CSV file for subsequent sorting and analysis.
* idock provides precompiled 64-bit executables for Linux, Mac and Windows.


Supported operating systems and compilers
-----------------------------------------

* Linux x86_64 and clang 3.7.0
* Mac OS X x86_64 and clang 3.7.0
* Windows x64 and Visual Studio 2015 Update 1

Pre-built 64-bit executables can be found in the `bin` directory. The executable was statically compiled for Linux and Windows, and dynamically compiled for Mac OS X.


Compilation from source code
----------------------------

idock depends on the [Boost C++ Libraries]. The Boost libraries required by idock are `System`, `Filesystem`, `Program Options` and `Thread`. Boost 1.60.0 was tested.

### Compilation on Linux

The Makefile uses clang as the default compiler. To compile, simply run

    make

One may modify the Makefile to use a different compiler or different compilation options.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `bin` folder.

### Compilation on Windows

Visual Studio 2015 solution and project files are provided. To compile, simply run

    msbuild /t:Build /p:Configuration=Release

Or one may open `idock.sln` in Visual Studio 2015 and do a full rebuild.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `bin` folder.


Usage
-----

First add idock to the PATH environment variable.

To display a full list of available options, simply run the program without arguments

    idock

The `examples` folder contains several use cases. For example, to dock the ligand TMC278 against HIV-1 RT of PDB ID 2ZD1,

    cd examples/2ZD1/T27

One can supply the options from command line arguments

    idock --receptor ../../../receptors/2ZD1.pdbqt --ligand ../../../ligands/T27 --center_x 49.712 --center_y -28.923 --center_z 36.824 --size_x 18 --size_y 18 --size_z 20

Or one can instruct idock to load the options from a configuration file

    idock --config idock.conf


Documentation
-------------

Documentations in both HTML and LaTeX formats can be esaily created by running [doxygen]

    doxygen idock.dox

The created documents will be placed in `doc` folder. To compile LaTeX files into PDF, one must have `pdflatex` installed.

    cd doc/latex
    make

The generated PDF will be `refman.pdf`.


Change Log
----------

### 2.2.1 (2016-01-07)

* Used the working directory as the default output directory.
* Forced docking if the two program options `ligand` and `out` are equivalent.
* Upgraded boost from 1.59.0 to 1.60.0.
* Upgraded RF-Score trained on PDBbind v2015 refined set.

### 2.2.0 (2015-12-05)

* Used program options `ligand` and `out` to specify input ligand(s) and output folder.
* Added the `score_only` option to allow scoring without docking.
* Filtered input ligand filenames with either .pdbqt or .PDBQT extension.
* Sorted input ligands alphabetically by their filename before docking.
* Added support for two new atom types Mo and Si.
* Ignored empty branches silently instead of throwing an exception.

### 2.1.4 (2015-02-01)

* Treated non-supported atom types as carbons.
* Wrote unsorted predicted free energy of both successfully and unsuccessfully docked ligands to log output.
* Changed default value of granularity from 0.1 to 0.125.
* Provided precompiled executable for 64-bit Mac OS X.
* Cleared all compilation warnings.
* Upgraded RF-Score trained on PDBbind v2014 refined set.

### 2.1.3 (2014-06-17)

* Bypassed already docked ligands by detecting file existence in output folder.
* Numbered REMARK records of output PDBQT files.
* Shortened the decimal digits from 3 to 2 in the REMARK records of output PDBQT files.
* Flushed ligand file stem before parsing.

### 2.1.2 (2014-04-06)

* Fixed a data race bug.

### 2.1.1 (2014-04-04)

* Fixed a segmentation fault in the precompiled executable for Linux.
* Enlarged the default number of Monte Carlo tasks from 32 to 64.
* Refined the default grid map granularity from 0.15625 to 0.1.
* Added 3 new examples: 1AQ1, 1PKD and 4MBS.

### 2.1.0 (2014-03-18)

* Conformed to [semantic versioning] using major.minor.patch.
* Added support for new chemical elements U and Cs.
* Added RF-Score trained on PDBbind v2013 refined set using 42 features for prospective rescoring.
* Updated the grid map creation algorithm.
* Updated the option `ligand_folder` to `input_folder`.
* Updated the output format.
* Updated the VC project files to Visual Studio 2013.
* Updated doxygen to idock.dox with version 1.8.6.
* Updated the extension name of configuration files in all examples from .cfg to .conf.
* Fixed a data race bug in the original thread pool by substituting the brandnew io service pool.
* Fixed an assertion bug caused by ligands with zero rotatable bond.
* Removed support for Mac, FreeBSD and Solaris.
* Removed precompiled 32bit executables for Linux and Windows.
* Removed support for gzip and bzip2.
* Removed the output of ligand efficiency and putative hydrogen bonds.
* Removed vina.cfg and vina.sh in all examples.

### 2.0 (2012-11-02)

* Output ligand efficiency.
* Updated Visual Studio 2012 project settings.
* Fixed a bug of repeatly executing the same Monte Carlo tasks.
* Updated the logo.

### 1.6 (2012-08-21)

* Added a new example 2VQZ.
* Output putative inter-molecular hydrogen bonds for each predicted conformation.
* Precompiled idock for Windows using Visual Studio 2012.
* Upgraded Visual Studio project from 2010 to 2012.
* Fixed a bug of aligning columns in log.csv.
* Fixed a bug of writing repeated energies to log.csv when compiling idock with clang 3.1 on Mac OS X and FreeBSD.
* Supported CentOS 6.3.
* Upgraded boost from 1.50.0 to 1.51.0.

### 1.5 (2012-06-13)

* Added a new example 2ZNL.
* Supported a new chemical element strontium (Sr).
* Updated clang from 3.0 to 3.1.
* Supported file error detection in output folder.
* Supported reading and writing ligands in gzip and/or bzip2 format.
* Output the number of hydrogen bonds for each conformation.

### 1.4 (2012-04-16)

* Fixed a segmentation fault bug when the number of heavy atoms exceeds 100.
* Added two new examples 2IQH and 1HCL.
* Reverted to file stem only in the ligand column in log.csv to shrink size.
* Added sufficient commas in log.csv to align rows.
* Skipped already docked ligands.
* Prevented dead loop by limiting the number of initial conformation trials.

### 1.3 (2012-03-05)

* Used a more compact and constant data structure for ligand representation.
* Refactored program option `conformations` to `max_conformations`.
* Output full path to docked ligands to csv.
* Removed boost::math::quaternion and implemented a lightweight quaternion class.
* Added BibTeX citation to the idock paper accepted and to be published in CIBCB 2012.
* Added bash scripts for running AutoDock Vina for docking ZINC clean drug-like ligands.
* Output predicted total free energy, predicted inter-ligand free energy and predicted intra-ligand free energy to docked PDBQT files.
* Output predicted free energy for each heavy atom to docked PDBQT files.
* Updated Boost from 1.48.0 to 1.49.0.
* Supported compilation on Windows 8 Consumer Preview x64 with Visual Studio 11 Ultimate Beta.
* Added a new example with PDB code 1V9U.
* Supported compilation on Solaris 11 11/11 with GCC 4.5.2.

### 1.2 (2012-02-06)

* Added program option `csv` for dumping docking summary sorted in the ascending of predicted free energy.
* Profiled by the Valgrind tool suite to ensure zero memory leak.
* Replaced a switch statement by table lookup to decrease indirect branch misprediction rate.
* Added move constructors for several classes to boost performance.
* Revised the precision of coordinates and free energy to be 3 digits.
* Parallelized the precalculation of scoring function.
* Fixed a numerical bug when docking a ligand of only one single heavy atom.
* Added support for Mac OS X 10.7.2 and FreeBSD 9.0.
* Added support for docking ligands created by igrow.

### 1.1 (2011-12-20)

* Changed the version control system from TFS to Git.
* Project migrated from CodePlex to GitHub.
* Tested Solaris 11, clang 3.0, and Intel C++ Compiler v11.
* Provided Visual C++ solution, project and bat files to ease recompilation on Windows.
* Added precompiled executables for both 32-bit and 64-bit Linux and Windows.
* Added program option `config` to allow users to specify a configuration file.
* Added thread-safe progress bar.
* Output predicted free energy of the top 5 conformations.
* Reverted the evaluation of intra-molecular free energy to Vina's implementation to obtain better RMSD for certain cases.

### 1.0 (2011-07-20)

* Initial release at [CodePlex].


Reference
---------

Hongjian Li, Kwong-Sak Leung, and Man-Hon Wong. idock: A Multithreaded Virtual Screening Tool for Flexible Ligand Docking. 2012 IEEE Symposium on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB), pp.77-84, San Diego, United States, 9-12 May 2012. [DOI: 10.1109/CIBCB.2012.6217214]


Author
--------------

[Jacky Lee]


Logo
----

![idock logo](https://github.com/HongjianLi/idock/raw/master/logo.png)


[virtual screening]: http://en.wikipedia.org/wiki/Virtual_screening
[docking]: http://en.wikipedia.org/wiki/Docking_(molecular)
[AutoDock Vina]: http://vina.scripps.edu
[Apache License 2.0]: http://www.apache.org/licenses/LICENSE-2.0
[istar]: http://istar.cse.cuhk.edu.hk/idock
[Boost C++ Libraries]: http://www.boost.org
[doxygen]: http://www.doxygen.org
[semantic versioning]: http://semver.org
[CodePlex]: http://idock.codeplex.com
[DOI: 10.1109/CIBCB.2012.6217214]: http://dx.doi.org/10.1109/CIBCB.2012.6217214
[Jacky Lee]: http://www.cse.cuhk.edu.hk/~hjli
