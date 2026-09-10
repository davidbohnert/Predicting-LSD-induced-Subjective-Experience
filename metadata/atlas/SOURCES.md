# Atlas metadata sources

| Distributed metadata | Derivation and use |
|---|---|
| `Schaefer416_8networks.txt` | The ordered 400 Schaefer cortical parcels and 16 AAL3 subcortical parcels, assigned to seven cortical networks plus a subcortical group |
| `atlas_parcel_metadata.csv` | Parcel labels, MNI centroids, and approximate AAL overlap descriptors for descriptive edge/node rankings |
| `Schaefer439_structure_membership.csv` | Cortical, subcortical, or cerebellar membership in the ordered Schaefer400 + Tian Scale II + Buckner7 atlas |

Parcel order must match the connectomes. AAL overlap labels are approximate
anatomical descriptions, not replacements for the parcel definitions.
The seven cerebellar parcels are identified from structure membership rather
than a hardcoded index boundary. Atlas images are not included.

## References and source distributions

- Schaefer A et al. (2018). *Local-Global Parcellation of the Human Cerebral
  Cortex from Intrinsic Functional Connectivity MRI*. Cerebral Cortex 28,
  3095–3114. [Paper](https://academic.oup.com/cercor/article/28/9/3095/3978804);
  [source distribution](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal).
- AAL/AAL3: [GIN atlas distribution and references](https://www.gin.cnrs.fr/fr/outils/aal/)
  and [AAL3 user guide](https://www.gin.cnrs.fr/wp-content/uploads/AAL3_UserGuide_April2024.pdf).
  AAL3 is described by Rolls and colleagues (2020), *Automated anatomical
  labelling atlas 3*.
- Tian Y et al. (2020). *Topographic organization of the human subcortex
  unveiled with functional connectivity gradients*.
  [Author distribution](https://github.com/yetianmed/subcortex).
- Buckner RL et al. (2011). *The organization of the human cerebellum
  estimated by intrinsic functional connectivity*.
  [Author distribution and references](https://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011).

The repository's MIT license covers its analysis code; upstream atlas resources
retain their own terms. Consult the linked distributions for source licenses.
