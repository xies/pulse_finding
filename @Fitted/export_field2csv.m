function export_field2csv(fits,filepath,fieldname)
%Exports the given FIELDNAME of a FIT array to CSV file
%
% USAGE: export_field2csv(fits,filepath,fieldname);

fitIDs = [fits.fitID]';
manual_flag = [fits.manually_added]';

mat2write = cat(2,fitIDs,manual_flag);

data = cat(1,fits.(fieldname));

mat2write = cat(2,mat2write,data);

header = [NaN NaN 1:numel(fits(1).(fieldname))];

mat2write = cat(1,header,mat2write);

csvwrite([filepath '/fits_' fieldname,'.csv'] , mat2write);

end %export_field2csv