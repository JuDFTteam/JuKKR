verdi code setup --non-interactive \
                 --label="voronoi" \
                 --on-computer \
                 --computer="slurmcontrol" \
                 --remote-abs-path="/builds/kkr/voronoi.exe" \
                 --input-plugin="kkr.voro" \
                 --prepend-text="ln -s /builds/kkr/ElementDataBase ."

