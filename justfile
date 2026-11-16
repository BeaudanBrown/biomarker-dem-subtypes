sync-m3:
    rsync -avrz --progress --delete m3:bc41_scratch/biomarker-dem-subtypes/.targets_beau/ ~/.config/targets/addf_targets/

sync-collab:
    rsync -avrz --progress --delete ~/.config/targets/addf_targets/ ~/collab/addf_targets/
