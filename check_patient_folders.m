function check_patient_folders( train_dir )
%CHECK_PATIENT_FOLDERS Summary of this function goes here
%   Detailed explanation goes here

    filename_map = containers.Map;
    filename_map('pre') = '**/pre_reg.nii*';
    filename_map('art') = '**/20s.nii*';
    filename_map('pv') = '**/70s_reg.nii*';
    filename_map('t2') = '**/t2_bfc_reg.nii*';
    filename_map('t1-bfc') = '**/bias_field_isotropic.nii';
    filename_map('liver_seg') = '**/*liver.ids';
    filename_map('tumor_seg') = '**/*tumor*.ids';
    filename_map('vasc_seg') = '**/*vessel*.ids';
%     filename_map('necro_seg') = '**/*nec*.ids';

    patients = dir(train_dir);
    filenames = {patients.name};
    patients = filenames([patients.isdir]);
    patients = patients(3:end);
    
    if isempty(patients)
        disp('bad directory');
    end
    
    for patient = patients
        patient = char(patient);
        for k = values(filename_map)
            k = char(k);
            mult = check_for_file([train_dir,'/',patient], k);
            if mult == 0
                disp([patient, ' does not contain file matching ', k]);
            elseif mult == 2
                disp([patient, ' has too many files matching ', k]);
            end
        end
    end
end

