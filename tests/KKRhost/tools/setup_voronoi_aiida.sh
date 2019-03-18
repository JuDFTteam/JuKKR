verdi code setup \
    --non-interactive \
    --label="voronoi" \
    --input-plugin="kkr.voro" "slurmcontrol" \
    --code-folder="/builds/kkr" \
    --code-rel-path="voronoi.exe" \
    --prepend-text="ln -s /builds/kkr/ElementDataBase ."
