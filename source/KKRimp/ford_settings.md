project: KKRimp
src_dir: ./SOURCE/
exclude_dir: doc
             SOURCE/deprecated
             SOURCE/deprecated/test_Rllsll
output_dir: ./doc
summary: Source code documentation of the KKRimp code for impurity calculation in a real-space cluster embedded into infinite host systems
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
version: 0.0.1
project_website: https://iffgit.fz-juelich.de/kkr/kkrimp


This is the introductory text ...
code locations [gitlab link](https://iffgit.fz-juelich.de/kkr/kkrimp)
[code wiki](https://iffwiki.fz-juelich.de/kkr/doku.php)

Some latex example:
\begin{equation}
\frac{\partial H}{\partial t} - \kappa\frac{\partial^{2} H}{\partial x^{2}} = f(x)
\end{equation}

@Note
This is an example note


@Bug
If you find any bugs, please email [me](mailto:p.ruessmann@fz-juelich.de)
or file a new issue on the [gitlab page](https://iffgit.fz-juelich.de/kkr/kkrimp/issues)
