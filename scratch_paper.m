isoname_map = containers.Map;
isoname_map('pre') = '/temp/pre_reg_isotropic.nii';
isoname_map('art') = '/temp/20s_isotropic.nii';
isoname_map('ven') = '/temp/70s_reg_isotropic.nii';
isoname_map('t2') = '/temp/t2_bfc_reg_isotropic.nii';
for a = keys(isoname_map)
    disp(a);
end
