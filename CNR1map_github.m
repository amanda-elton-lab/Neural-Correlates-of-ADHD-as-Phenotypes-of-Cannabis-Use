% This code is used in MATLab to create a CNR1 gene expression map from Allen Human Brain Atlas microarray data and covary expression with functional activation from stop signal task fMRI data.
% Authors: Jacqueline Aloumanis, Amanda Elton

% Download normalized microarray datasets from AHBA 

% Load and merge data from whole-brain donors
donor_ids = ["H0351.2001","H0351.2002"];

all_expr_data = [];
all_gene_symbols = [];
all_sample_cords = [];

for i = 1:length(donor_ids)
    donor = donor_ids(i);

    expr_file = fullfile('downloadpath');
    probe_file = fullfile('downloadpath');
    sample_file = fullfile('downloadpath');

    expr_data = readtable(expr_file);
    probes = readtable(probe_file);
    samples = readtable(sample_file);

    if i == 1
        all_gene_symbols = probes.gene_symbol;
    end

    all_expr_data = [all_expr_data, expr_data{:, 2:end}];

    all_sample_cords = [all_sample_cords; samples.mni_x, samples.mni_y, samples.mni_z];
end 

% Extracting gene of interest (CNR1)
gene_symbol = "CNR1";

gene_idx = strcmp(all_gene_symbols, gene_symbol);
gene_expression = mean(all_expr_data(gene_idx, :), 1);

%Map MNI coordinates into whole-brain grid
[xq, yq, zq] = ndgrid(-80:2:80, -110:2:90, -60:2:80);

expression_map = griddata(all_sample_cords(:,1), all_sample_cords(:,2), all_sample_cords(:,3), gene_expression', xq, yq, zq, 'natural');

nii = make_nii(expression_map);
save_nii(nii, 'CNR1_expression_map_fullbrains.nii');

% Load Gray Matter Mask of choice (we used a binary mask created in AFNI)
gm_mask_nii = load_nii('GM_mask_afni_2x2x2.nii');
disp(gm_mask_nii.hdr.dime)
gm_mask = double(gm_mask_nii.img>0);
target_size = size(gm_mask);

%Resample map to match GM Mask size (if needed)
cnr1_nii = load_nii('CNR1_expression_map_fullbrains.nii');
cnr1_map = double(cnr1_nii.img);
resampled_cnr1_map = imresize3(cnr1_map, target_size, 'cubic');

resampled_map = make_nii(resampled_cnr1_map);
resampled_map_nii.hdr.dime.pixdim = gm_mask_nii.hdr.dime.pixdim;
resampled_map_nii.hdr.dime.dim = gm_mask_nii.hdr.dime.dim;
resampled_map_nii.hdr.hist.originator = gm_mask_nii.hdr.hist.originator;

%Apply GM Mask to Gene Expression Map
masked_expression_map = resampled_cnr1_map .* gm_mask;
masked_nii = make_nii(masked_expression_map);
masked_nii.hdr.dime.pixdim = gm_mask_nii.hdr.dime.pixdim;
masked_nii.hdr.dime.dim = gm_mask_nii.hdr.dime.dim;
masked_nii.hdr.hist.originator = gm_mask_nii.hdr.hist.originator;
save_nii(masked_nii, 'CNR1_expressionmap_fullbrains_afnimaskedbin.nii');


%Check overlap of GM Mask and Map
cnr1_nii = masked_nii;
cnr1_map = double(cnr1_nii.img);

disp(cnr1_nii.hdr.dime.pixdim);  
disp(cnr1_nii.hdr.hist.originator);

disp(gm_mask_nii.hdr.dime.pixdim);
disp(gm_mask_nii.hdr.hist.originator);

overlap = double(cnr1_nii.img > 0) .* double(gm_mask_nii.img > 0);
disp(['Number of overlapping voxels: ', num2str(sum(overlap(:)))]);


%Create CNR1 Vector
cnr1_vector = cnr1_map(gm_mask > 0);

% Set up applying GM Mask to functional images and covariance
input_dir = 'yourpath';
output_dir = 'yourpath';

valid_ids = ['yourdataids'];

subject_covariances = zeros(length(valid_ids), 1);

%Apply GMM to each participant's SST contrast, create SST vector and covary with CNR1 Vector
for i = 1:length(valid_ids)
    sub = valid_ids(i);
    participant_id = sprintf('%03d', sub);
    contrast_filename = fullfile('S:\Research\Elton\Studies\LIA Study\MRI data\sst\',['sub-' participant_id '_task-sst_reml.nii']);

    if isfile(contrast_filename)
        contrast_nii = load_nii(contrast_filename);
        contrast_data = double(contrast_nii.img);

        volume = squeeze(contrast_data(:,:,:,1,14));
        masked_volume = volume .* double(gm_mask);

        %Creating SST vector
        sst_vector = masked_volume(gm_mask > 0);

        %Winsorize outliers
        sst_vector(find(sst_vector>prctile(sst_vector,99)))=prctile(sst_vector,99);
        sst_vector(find(sst_vector<prctile(sst_vector,1)))=prctile(sst_vector,1);

        %Calculate vectors for valid (non-NaN) voxels
        valid_voxels = ~isnan(cnr1_vector) & ~isnan(sst_vector) & std([cnr1_vector sst_vector], 0, 2) > 0;

        %Covary vectors
        cov_values = cov(cnr1_vector(valid_voxels), sst_vector(valid_voxels))
        subject_covariances(i) = cov_values(1,2)

    else
        fprintf('Subject %s not found, skipping.\n', participant_id);
    end
end





