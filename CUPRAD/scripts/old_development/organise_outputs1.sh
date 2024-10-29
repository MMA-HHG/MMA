#!/bin/bash


mkdir dPhi_maps
mv abs_dPhi_dz_* dPhi_maps/
mv dPhi_dz_* dPhi_maps/

mkdir Curvature
mv Curvature_* Curvature/
mv ddr_Curvature_* Curvature/

mkdir Cutoffs_slices
mv Cutoff_t* Cutoffs_slices/

mkdir Fluences_and_max_cutoffs
mv Fluence_* Fluences_and_max_cutoffs/
mv Cutoff_max_* Fluences_and_max_cutoffs/

mkdir Intensity
mv Intens_* Intensity/

mkdir Efield
mv Efield_* Efield/

mkdir Lcoh
mv Lcoh_* Lcoh/

mkdir Plasma
mv Plasma_* Plasma/

mkdir phases_onax
mv phase_onax_* phases_onax/

mkdir misc
mv phase_in_* misc/

echo "done" 
