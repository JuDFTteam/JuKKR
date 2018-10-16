project: Jülich KKR code for bulk and interfaces
summary: Source code documentation of the Jülich KKR code for bulk and interfaces
author: The Jülich KKR team
email: p.ruessmann@fz-juelich.de
license: by-nc
version: 1.1.1
project_website: https://iffgit.fz-juelich.de/kkr/kkrjm
project_download: https://iffgit.fz-juelich.de/kkr/kkrjm
src_dir: ./ALL_SOURCE_FILES
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
warn: true 
dbg: true

This is the introductory text ...

Code location:
 - [gitlab link](https://iffgit.fz-juelich.de/kkr/kkrjm)
 - [code wiki](https://iffwiki.fz-juelich.de/kkr/doku.php)
 - [online version of the documentation](https://kkr.iffgit.fz-juelich.de/kkrjm)

Some latex example:
\begin{equation}
\frac{\partial H}{\partial t} - \kappa\frac{\partial^{2} H}{\partial x^{2}} = f(x)
\end{equation}

@note
This is an example note
@endnote


@bug
If you find any bugs, please email [me](mailto:p.ruessmann@fz-juelich.de)
or file a new issue on the [gitlab page](https://iffgit.fz-juelich.de/kkr/kkrjm/issues)
@endbug

@warning
Some warning
@endwarning