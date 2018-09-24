project: Jülich KKR code for bulk and interfaces
src_dir: ./
exclude_dir: doc
             libs/FEAST/3.0/src/banded
             libs/FEAST/3.0/src/kernel
             libs/FEAST/3.0/src/sparse
             libs/FEAST/3.0/src/dense
             libs/FEAST/3.0/src/banded
             libs/FEAST/3.0/src/banded/spike-smp
             libs/FEAST/3.0/src
output_dir: ./doc
summary: Source code documentation of the Jülich KKR code for bulk and interfaces
author: The Jülich KKR team
email: p.ruessmann@fz-juelich.de
docmark: !
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: true
graph: true
extra_filetypes: 1 #
		 inc !
		 mk #
coloured_edges: true
search: true 
warn: true 
license: by-nc
version: 1.1.1
project_website: https://iffgit.fz-juelich.de/kkr/kkrjm
dbg: true

This is the introductory text ...
code locations [gitlab link](https://iffgit.fz-juelich.de/kkr/kkrjm)
[code wiki](https://iffwiki.fz-juelich.de/kkr/doku.php)

Some latex example:
\begin{equation}
\frac{\partial H}{\partial t} - \kappa\frac{\partial^{2} H}{\partial x^{2}} = f(x)
\end{equation}

@Note
This is an example note


@Bug
If you find any bugs, please email [me](mailto:p.ruessmann@fz-juelich.de)
or file a new issue on the [gitlab page](https://iffgit.fz-juelich.de/kkr/kkrjm/issues)
