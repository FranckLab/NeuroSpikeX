function [data] = LoadData(fileName, pathName)

%% FUNCTION TO LOAD FDATA INTO THE CURRENT WORKSPACE

%%
selectedDataFile = fullfile(pathName,fileName);
struct_data = load(selectedDataFile);

if strcmp(fileName, 'fdata.mat')
    data = struct_data.fdata;
elseif strcmp(fileName, 'result.mat')
    data = struct_data.result;
end

end