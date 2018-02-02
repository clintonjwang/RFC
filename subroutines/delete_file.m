function success = delete_file( filepath )
%DELETE_FILE Summary of this function goes here
%   Detailed explanation goes here

    if exist(filepath, 'file') == 2
        delete(filepath);
        success = true;
    else
        success = false;
    end
end