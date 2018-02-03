function data = load_wrapper( path )
%LOAD_WRAPPER Summary of this function goes here
%   Detailed explanation goes here
    data = load( path );
    data = data.data;
end

