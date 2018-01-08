function [ output_args ] = delete_file( input_args )
%DELETE_FILE Summary of this function goes here
%   Detailed explanation goes here

if exist('ResultFile.txt', 'file') == 2
  delete('ResultFile.txt');
end

end

