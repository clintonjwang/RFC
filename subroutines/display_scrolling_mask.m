function display_scrolling_mask( )
%DISPLAY_SCROLLING_MASK Summary of this function goes here
%   Detailed explanation goes here

    data_dir = '../small_data/0347479';
    mask_names = {'vasculature_mask', 'necrosis_mask', 'viable_tumor_mask'};
    img_name = '20s';
    img = obtain_3d_mask(img_name, mask_names, data_dir);
    [~,~,max_slice,~] = size(img);
    slice_num = round(max_slice/2);

    img_slice = squeeze(img(:,:,slice_num,:));
    f = figure('Name', [img_name, ' masked']);
    image(img_slice,'CDataMapping','scaled');
    hold on;
    title([img_name, ' masked']);

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


    % Create pop-up menu
    % popup = uicontrol('Style', 'popup',...
    %        'String', {'vasculature_mask','necrosis_mask','viable_tumor_mask'},...
    %        'Position', [20 340 100 50],...
    %        'Callback', @setmap);

    uicontrol('Style', 'checkbox', 'String', 'Show vasculature',...
        'Position', [20 70 50 20], 'Value', 1,...
        'Callback', @toggle_vasc);
    uicontrol('Style', 'checkbox', 'String', 'Show necrosis',...
        'Position', [20 45 50 20], 'Value', 1,...
        'Callback', @toggle_nec);
    uicontrol('Style', 'checkbox', 'String', 'Show viable tumor',...
        'Position', [20 20 50 20], 'Value', 1,...
        'Callback', @toggle_tumor);

    % Create slider
    uicontrol('Style', 'slider',...
        'Min',1,'Max',max_slice,'Value',slice_num,...
        'Position', [400 20 120 20],...
        'Callback', @shift_slice);

    % Add a text uicontrol to label the slider.
    uicontrol('Style','text',...
        'Position',[400 45 120 20],...
        'String','Slice (fraction)');

    % Make figure visble after adding all components
%     f.Visible = 'on';
    % This code uses dot notation to set properties. 
    % Dot notation runs in R2014b and later.
    % For R2014a and earlier: set(f,'Visible','on');

    % function setmap(source,event)
    %     val = source.Value;
    %     maps = source.String;
    %     newmap = maps{val};
    %     colormap(newmap);
    % end

    function shift_slice(source,~)
        old_slice = slice_num;
        slice_num = round(source.Value);
        if old_slice ~= slice_num
            replot();
        end
    end

    function toggle_vasc(source,~)
        toggle_mask('vasculature_mask', source);
    end

    function toggle_nec(source,~)
        toggle_mask('necrosis_mask', source);
    end

    function toggle_tumor(source,~)
        toggle_mask('viable_tumor_mask', source);
    end

    function toggle_mask(mask_name, source)
        if source.Value == source.Max
            mask_names{length(mask_names) + 1} = mask_name;
        else
            todelete = false(size(mask_names));
            for c = 1:length(mask_names)
                todelete(c) = strcmp(mask_names{c}, mask_name) == 1;
            end
            mask_names(todelete) = [];
        end

        replot();
    end

    function replot()
        img = obtain_3d_mask(img_name, mask_names, data_dir);
        img_slice = squeeze(img(:,:,slice_num,:));
        image(img_slice,'CDataMapping','scaled');
        
        h = zeros(length(mask_names), 1);
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
end


function img = obtain_3d_mask(img_name, mask_names, data_dir)
%DISPLAY_MASKED_IMG Displays MRI slice with binary mask overlay
%   Mask is colored

    nii_ext = {'*.nii; *.hdr; *.img; *.nii.gz'};
    
    data = load_nii(try_find_file(data_dir, ['**/', img_name, '.nii'],...
            'Select the image.', nii_ext));
    img = double(flip_image(data.img));
    [N1,N2,N3] = size(img);
    
    img = img/max(max(max(img)));
    img(:,:,:,2) = img(:,:,:,1);
    img(:,:,:,3) = img(:,:,:,1);
    
    if ~isempty(mask_names)
        f = try_find_file(data_dir, ['**/', mask_names{1}, '.ids'],...
                    'Select the first mask to apply', '*.ids');
        mask = logical(get_mask(f, N1,N2,N3));
%         img(:,:,:,2) = img(:,:,:,2) .* ~mask;
%         img(:,:,:,3) = img(:,:,:,3) .* ~mask;
        img(:,:,:,1) = img(:,:,:,1) + (1 - img(:,:,:,1)) .* (mask * 0.1);

        if length(mask_names) > 1
            f = try_find_file(data_dir, ['**/', mask_names{2}, '.ids'],...
                        'Select the second mask to apply', '*.ids');
            mask = logical(get_mask(f, N1,N2,N3));

%             img(:,:,:,1) = img(:,:,:,1) .* ~mask;
%             img(:,:,:,2) = img(:,:,:,2) .* ~mask;
            img(:,:,:,3) = img(:,:,:,3) + (1 - img(:,:,:,3)) .* (mask * 0.1);

            if length(mask_names) > 2
                f = try_find_file(data_dir, ['**/', mask_names{3}, '.ids'],...
                            'Select the third mask to apply', '*.ids');
                mask = logical(get_mask(f, N1,N2,N3));

%                 img(:,:,:,2) = img(:,:,:,2) .* ~mask + mask * 255;
                img(:,:,:,2) = img(:,:,:,2) + (1 - img(:,:,:,2)) .* (mask * 0.1);
            end
        end
    end
end