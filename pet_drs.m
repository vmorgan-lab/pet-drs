%% Calculation of PET DRS given PET data and planned resection

% L Sainburg October 2024

% SPM12 is a prerequisite for this script

% This script runs the PET DRS algorithm, as outlined in Sainburg et al.,
% 2025. The algorithm takes presurgical FDG-PET and T1-weighted images along
% with regions that are planned to be resected and returns a region-wise 
% PET laterality map along with a PET DRS measure. This algorithm can be
% used to obtain PET DRS measures during epilepsy surgery planning to
% assess the likelihood of seizure freedom given a planned surgery and the
% patient's PET scan.

% Once this script is run initially with reg = 1, reg can be set to 0
% when rerunning to test new sets of regions for planned resections.

function pet_drs(parent_dir,res_labels,sim_res_fname,reg)

% INPUTS:

% parent_dir: Parent directory folder where inputs and outputs will be
% placed. This folder should contain a presurgical T1-weighted MRI named 
% 'T1.nii' and a presurgical PET scan named 'PET.nii'

% res_labels: numbers of regions in planned resection based on provided AAL
% atlas region-number mapping (third column of AAL_regions.txt). If a
% standard selective amygdalohippocampectomy (SAH) or anterior temporal
% lobectomy (ATL) is to be tested, the user can alternatively input one of
% the following: 'Right ATL','Left ATL','Right SAH','Left SAH' (case
% sensitive).

% sim_res_fname: Name of folder to put results from simulated resection.
% This folder will be placed within parent_dir, which allows for multiple
% resections to be tested, given different folder names.

% reg: 1 or 0 (default = 1); whether to register the PET to MNI space. This
% can be set to 0 if PET was already warped to MNI space by running this 
% code and the user only wants to recalculate PET DRS.

% OUTPUTS:

% PET_LI.nii: AAL parcellation of PET Laterality map (overlay on wT1.nii or
% MNI template to visualize)
% PET DRS Value report in sim_res_fname (pet DRS and resected region names)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE

% Example use for planned resection of right hippocampus and amygdala:
% pet_drs('test_folder',[4102,4202],'rsah',1)

% Example use for planned resection of right hippocampus and amygdala using 
% preset regions in standard SAH:
% pet_drs('test_folder','Right SAH','rsah',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(reg)
    reg = 1;
end

% STEP 1: REGISTER PET IMAGE TO MNI TEMPLATE (SPM12)
if reg == 1

% find spm
spm_path = which('spm');
spm_path = spm_path(1:end-5);

% reset origin of PET image to center (sometimes needed for registration)
V = spm_vol(fullfile(parent_dir,'PET.nii'));
pet_orig = spm_read_vols(V);
% find center of mass of image to set as origin
sumTotal = sum(pet_orig(:));
coivox = ones(4,1);
coivox(1) = sum(sum(sum(pet_orig,3),2)'.*(1:size(pet_orig,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(pet_orig,3),1).*(1:size(pet_orig,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(pet_orig,2),1))'.*(1:size(pet_orig,3)))/sumTotal; %dimension 3
XYZ_mm = V.mat * coivox; %convert from voxels to millimeters
V.mat(1,4) =  V.mat(1,4) - XYZ_mm(1);
V.mat(2,4) =  V.mat(2,4) - XYZ_mm(2);
V.mat(3,4) =  V.mat(3,4) - XYZ_mm(3);
spm_write_vol(V,pet_orig);

clear matlabbatch
    
% Segment MRI to get Grey Matter Map and MNI warping
matlabbatch{1}.spm.spatial.preproc.channel.vols = {[fullfile(parent_dir,'T1.nii'),',1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[fullfile(spm_path,'tpm','TPM.nii'),',1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[fullfile(spm_path,'tpm','TPM.nii'),',2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[fullfile(spm_path,'tpm','TPM.nii'),',3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[fullfile(spm_path,'tpm','TPM.nii'),',4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[fullfile(spm_path,'tpm','TPM.nii'),',5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[fullfile(spm_path,'tpm','TPM.nii'),',6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];

% Register PET to MRI Grey Matter
matlabbatch{2}.spm.spatial.coreg.estwrite.ref = {[fullfile(parent_dir,'c1T1.nii'),',1']};
matlabbatch{2}.spm.spatial.coreg.estwrite.source = {[fullfile(parent_dir,'PET.nii'),',1']};
matlabbatch{2}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Register T1 and PET to MNI
matlabbatch{3}.spm.spatial.normalise.write.subj.def = {fullfile(parent_dir,'y_T1.nii')};
matlabbatch{3}.spm.spatial.normalise.write.subj.resample = {[fullfile(parent_dir,'T1.nii'),',1']
                                                            [fullfile(parent_dir,'rPET.nii'),',1']
                                                            [fullfile(parent_dir,'c1T1.nii'),',1']};
matlabbatch{3}.spm.spatial.normalise.write.woptions.bb = [Inf Inf Inf
                                                          Inf Inf Inf];
matlabbatch{3}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
matlabbatch{3}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{3}.spm.spatial.normalise.write.woptions.prefix = 'w';

spm_jobman('run',matlabbatch)

end

% STEP 2: CALCULATE PET DRS

% get grey matter mask in MNI space
gm_info = spm_vol(fullfile(parent_dir,'wc1T1.nii'));
gm = spm_read_vols(gm_info);
gm_msk = gm>0.5;

% load pet image
pet_info = spm_vol(fullfile(parent_dir,'wrPET.nii'));
pet = spm_read_vols(pet_info);

% normalize PET intensity based on grey matter
pet_m = mean(pet(gm_msk));
pet_std = std(pet(gm_msk));
pet_z = (pet-pet_m)/pet_std;

% load in AAL labels
roi_labels = readtable('AAL_regions.txt');
rois = roi_labels.Var2;
roi_vals = roi_labels.Var3;
roi_labels = addvars(roi_labels,zeros(length(rois),1),[ones(90,1);zeros(length(rois)-90,1)],'NewVariableNames',{'ROI_side','use_pet'});

% get ROI sides
for ii = 1:length(rois)
    if contains(rois{ii}(end),'L')==1
        roi_labels.ROI_side(ii) = 1;
    elseif contains(rois{ii}(end),'R')==1
        roi_labels.ROI_side(ii) = 2;
    end
end

% load in AAL
aal_info = spm_vol('AAL_3mm.nii');
aal = spm_read_vols(aal_info);

% get ROI-wise PET
pet_z_roi = nan(size(roi_labels,1),1);
for ii = 1:size(roi_labels,1)
    pet_z_roi(ii,1) = mean(pet_z(aal==roi_vals(ii) & gm_msk==1));
end

% Calculate PET Laterality Map
pet_z_roi_side_diff = pet_z_roi(roi_labels.ROI_side==1)-pet_z_roi(roi_labels.ROI_side==2);
pet_z_roi_li = nan(size(pet_z_roi));
pet_z_roi_li(roi_labels.ROI_side==1) = pet_z_roi_side_diff;
pet_z_roi_li(roi_labels.ROI_side==2) = -pet_z_roi_side_diff;

% Save PET Laterality Map
pet_li_map = nan(size(pet));
for ii = 1:numel(pet_li_map)
    if aal(ii)>0 && sum(roi_labels.Var3==aal(ii),'all')>0
        pet_li_map(ii) = pet_z_roi_li(roi_labels.Var3==aal(ii));
    end
end
pet_li_map_info = pet_info;
pet_li_map_info.fname = fullfile(parent_dir,'PET_LI.nii');
spm_write_vol(pet_li_map_info,pet_li_map);

% Calculate PET DRS with PET Laterality Map and Resection Labels
rois_use = [1:70,79:90]; % exclude subcortical and cerebellar regions that wouldn't be resected from analysis

% Pre-determined standard resections vs user-defined resection
if ischar(res_labels)
    if contains(res_labels,'Right ATL')
        labels = ismember(roi_vals,[4102,4112,4202,8122,8202,8212,8302]);
    elseif contains(res_labels,'Left ATL')
        labels = ismember(roi_vals,[4101,4111,4201,8121,8201,8211,8301]);
    elseif contains(res_labels,'Right SAH')
        labels = ismember(roi_vals,[4102,4202]);
    elseif contains(res_labels,'Left SAH')
        labels = ismember(roi_vals,[4101,4201]);
    else
        error('Resection not defined')
    end
else
    labels = ismember(roi_vals,res_labels);
end
rois_res = rois(labels); % names of resected regions
mkdir(fullfile(parent_dir,sim_res_fname))
writecell(rois_res,fullfile(parent_dir,sim_res_fname,'resected_regions.txt'))

% Compute PET DRS
[~,~,~,auc] = perfcurve(labels(rois_use),-pet_z_roi_li(rois_use),1);
pet_drs = 1-auc;

% Plot PET Laterality values of resected and spared regions
rng('default')
pltsz = [350 350];
fig = figure;fig.Position(3:4) = pltsz;
swarmchart(double(labels(rois_use)),pet_z_roi_li(rois_use),75,'k','MarkerFaceColor',[.5 .5 .5]);hold on
xticks([0 1]);xticklabels({'Spared','Resected'});axis square
ylabel('PET Laterality');axis square;xlabel('Regions')
ax = gca;ax.FontSize = 16;ax.FontWeight = 'bold';
annotation('textbox', [0.45, 0.8, 0.4, 0.1], 'String', ['PET D_R_S = ' num2str(round(pet_drs,3))],...
    'LineStyle','none','FontSize',12,'FontWeight','bold')
saveas(fig,fullfile(parent_dir,sim_res_fname,'PET_DRS_scatter.jpg'))

% Display PET DRS
disp(['PET DRS = ' num2str(pet_drs)])

end