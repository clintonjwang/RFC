function img = display_masked_img(img_name, mask_names, data_dir, slice_fraction)
%DISPLAY_MASKED_IMG Displays MRI slice with binary mask overlay
%   Mask is colored

    if nargin < 4
        slice_fraction = .5;
    end
    nii_ext = {'*.nii; *.hdr; *.img; *.nii.gz'};
    
    data = load_nii(try_find_file(data_dir, ['**/', img_name, '.nii'],...
            'Select the image.', nii_ext));
    img = double(flip_image(data.img));
    [N1,N2,N3] = size(img);
    
    img_slice = img(:,:,round(N3*slice_fraction));
    img_slice = img_slice/max(max(img_slice));
    img_slice(:,:,2) = img_slice(:,:,1);
    img_slice(:,:,3) = img_slice(:,:,1);
    
    if ~isempty(mask_names)
        h = zeros(length(mask_names), 1);
        
        f = try_find_file(data_dir, ['**/', mask_names{1}, '.ids'],...
                    'Select the first mask to apply', '*.ids');
        mask = get_mask(f, N1,N2,N3);
        mask_slice = logical(mask(:,:,round(N3*slice_fraction)));
        img_slice(:,:,2) = img_slice(:,:,2) .* ~mask_slice;
        img_slice(:,:,3) = img_slice(:,:,3) .* ~mask_slice;

        if length(mask_names) > 1
            f = try_find_file(data_dir, ['**/', mask_names{2}, '.ids'],...
                        'Select the second mask to apply', '*.ids');
            mask = get_mask(f, N1,N2,N3);
            mask_slice = logical(mask(:,:,round(N3*slice_fraction)));

            img_slice(:,:,1) = img_slice(:,:,1) .* ~mask_slice;
            img_slice(:,:,2) = img_slice(:,:,2) .* ~mask_slice;
        %     img_slice(:,:,3) = img_slice(:,:,3);

            if length(mask_names) > 2
                f = try_find_file(data_dir, ['**/', mask_names{3}, '.ids'],...
                            'Select the third mask to apply', '*.ids');
                mask = get_mask(f, N1,N2,N3);
                mask_slice = logical(mask(:,:,round(N3*slice_fraction)));

                img_slice(:,:,1) = img_slice(:,:,1) .* ~mask_slice;
                img_slice(:,:,3) = img_slice(:,:,3) .* ~mask_slice;
            end
        end
    end

    img = figure('Name', [img_name, ', slice fraction ', num2str(slice_fraction)]);
    image(img_slice,'CDataMapping','scaled');
    hold on;
    title([img_name, ', slice fraction ', num2str(slice_fraction)]);
    
    if ~isempty(mask_names)
        h(1) = plot(NaN,NaN,'or');
        if length(mask_names) > 1
            h(2) = plot(NaN,NaN,'ob');
            if length(mask_names) > 2
                h(3) = plot(NaN,NaN,'og');
            end
        end
        legend(h, mask_names, 'Interpreter', 'none');
    end
end