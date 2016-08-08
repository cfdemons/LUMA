function tecplotsorterdebug(t,tstep, l, r, precision)

for tt = 0 : tstep : t
    
    for lev = 0 : 1 : l
        
        for reg = 0 : 1 : r

            try
            % Get file name
            input_filename = ['tecplotdebug.Lev' num2str(lev) '.Reg' num2str(reg) '.' num2str(tt) '.dat'];
            output_filename = ['tecplotsrteddebug.Lev' num2str(lev) '.Reg' num2str(reg) '.' num2str(tt) '.dat'];

            % Read file
            data = importdata(input_filename,'\t',7);

            % Sort
            sorteddata = sortrows(sortrows(sortrows(data.data,1),2),3);

            % Write header
            data.textdata(7,:) = [];
            tbl = cell2table(data.textdata);
            writetable(tbl, output_filename);

            % Write data
            dlmwrite(output_filename',sorteddata,'-append',  ...
            'delimiter', '\t',  ...
            'precision', ['%.' num2str(precision) 'f'] );
            end
        end
    end
    
end


end