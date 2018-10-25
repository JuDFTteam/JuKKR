project: The Jülich KKR code
summary: Source code documentation of the Jülich KKR code
author: The Jülich KKR team
email: p.ruessmann@fz-juelich.de
license: by-nc
version: 2.4
project_website: https://iffgit.fz-juelich.de/kkr/kkrjm
project_download: https://iffgit.fz-juelich.de/kkr/kkrjm
src_dir: ./ALL_SOURCE_FILES
page_dir: ./
output_dir: ./doc
media_dir: ./media
docmark: !
predocmark: >
docmark_alt: <
predocmark_alt: #
display: public
         protected
         private
extra_filetypes: 1 #
		 inc !
		 mk #
source: true
search: true 
graph: true 
coloured_edges: true
warn: false 
dbg: true


#### Welcome to the Jülich KKR code for bulk and interfaces!

Our Jülich KKR package allows to perform all electron density functional theory calculations to analyze

* bulk, film and semi infinite systems
* concentrated disordered and dilute alloys (in the virtual crystal approzimation or the coherent potential approximation)
* single atoms and finite-sized clusters in a host-system

employing

* non-relativistic, scalar relativistic or relativistic calculations including spin-orbit coupling effects
* embedding-technique in real space (e.g. for single impurities without the need for supercells) with the [KKRimp code](https://iffgit.fz-juelich.de/kkr/kkrimp)
* treatment of the full potential (using the [Voronoi code](https://iffgit.fz-juelich.de/kkr/voronoi)) as well as using the atomic sphere approximation (ASA)
* Boltzmann and Landauer-Büttiker transport formalisms (using the [Pkkprime code](https://iffgit.fz-juelich.de/kkr/pkkr)) 
* time-dependent density functional theory (using the [KKRsusc extension](https://iffgit.fz-juelich.de/kkr/kkrsusc))

to predict

* scattering effects (with the [Pkkprime code](https://iffgit.fz-juelich.de/kkr/pkkr))
* electronic transport properties (e.g. conductivities, spin relaxation, family of Hall and Nernst effects with the [Pkkprime code](https://iffgit.fz-juelich.de/kkr/pkkr))
* magnetic response functions (with the [KKRsusc code](https://iffgit.fz-juelich.de/kkr/kkrsusc))
* magnetic (exchange) parameters for micromagnetic and atomistic spin models
* quasiparticle interference (QPI) spectra

#### Code characteristics:

* linear scaling for layered systems
* highly parallelized, hybrid parallelization using 2-levels of MPI and an OpenMP level
* The KKR package is part of the [JuDFT code family](http://www.judft.de/pm/index.php)
* Our KKR code is interfaced with the [AiiDA framework](http://www.aiida.net)


#### Code location and further reading
  * Down load the code via the [gitlab link](https://iffgit.fz-juelich.de/kkr/kkrjm)
  * See examples of the code usage on the [code's wiki page](https://iffwiki.fz-juelich.de/kkr/doku.php)
  * Find the latest [online version of the documentation](https://kkr.iffgit.fz-juelich.de/kkrjm)
  * Download the [`aiida-kkr` plugin](https://github.com/broeder-j/aiida-kkr)
  * Check the [online documentation of `aiida-kkr`](https://aiida-kkr.readthedocs.io)

@note
Please check the ***Readme*** page of this documentation for further information on:

  * how to install the code
  * how to contribute
  * the code's license
  * the changelog
  
@endnote

@bug
If you find any bugs, please file a new issue on the [gitlab page](https://iffgit.fz-juelich.de/kkr/kkrjm/issues)
@endbug

