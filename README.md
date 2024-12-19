# PET DRS
PET DRS algorithm to assess targeting of lateralized hypometabolism for epilepsy surgery planning.

# Prerequisites
SPM12 is a prerequesite to run this script and can be downloaded here for free: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
*Make sure that SPM12 is in your path when running this script
Note that this script was developed and tested using MATLAB 2021a.

# Inputs
Imaging

This script takes T1-weighted MRI and FDG-PET images in NIFTI format as inputs. They should be placed in the "parent directory" folder and named "T1.nii" and "PET.nii", respectively. All outputs will be placed in this parent directory folder.

Resection

The regions in the planned resection should be defined using the AAL atlas (see AAL_regions.txt). Users can also use one of the predefined resections within the code ('Right ATL', 'Left ATL', 'Right SAH', or 'Left SAH').

Registration

Once this script is run the first time on a patient's data with reg = 1, the registration information has been saved and reg can be set to 0 for subsequent tests of alternative resections to decrease processing time.

# Outputs
A swarmchart will be outputted with the PET DRS score, along with a PET_LI.nii image. The PET_LI.nii image can be overlaid on the wT1.nii image (patient's T1 in MNI space) or the MNI template for visualization (see examples in Sainburg et al., 2024). Note that as per Sainburg et al., 2024, a PET DRS < 0.22 is suggestive of an Engel I outcome after surgery.

# Example Usage (put T1.nii and PET.nii in "test_folder"):
Resection of right amygdala and hippocampus (right selective amygdalohippocampectomy):

pet_drs('test_folder',[4102,4202],'rsah',1)

pet_drs('test_folder','Right SAH','rsah',1)

# Notes
We would like to emphasize that this tool should be used to supplement other clinical data and not to guide surgery on its own.

We encourage others to edit and adjust this algorithm to their own needs. For example, other *anatomically symmetric* atlases could be substituted for the AAL atlas
