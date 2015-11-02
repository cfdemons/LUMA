% Get file name
input_filename = 'tecplotout.Lev0.Reg0.5000.dat';
output_filename = 'tecplotsorted.Lev0.Reg0.5000.dat';

% Read file
data = importdata(input_filename,'\t',7);

% Sort
sorteddata = sortrows(sortrows(sortrows(data.data,1),2),3);

% Write header
tbl = cell2table(data.textdata);
writetable(tbl, output_filename);

% Write data
dlmwrite(output_filename',sorteddata,'-append', 'delimiter', '\t');