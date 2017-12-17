function [ header, metadata, voxel_data ] = liver_seg_interface( input_args )
%LIVER_SEG_INTERFACE Summary of this function goes here
%   Detailed explanation goes here

voxel_data = user_main(false);

header.HeaderVersion = '20171216';
header.NumberOfVoxels = size(voxel_data);
header.VolumeExtent = size(voxel_data);
header.SliceThickness = VolumeExtent(3);
header.RescaleSlope = 1;
header.RescaleIntercept = 0;
header.TransformationMatrix = zeros(12,1);
header.Type = 0; %1 for binary
header.Description = 'Segmentation';
end

