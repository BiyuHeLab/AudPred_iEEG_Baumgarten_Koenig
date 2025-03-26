function data_filtered = myft_trialfilter(data, filter, field_length)
% data_filtered = myft_trialfilter(data, filter, field_length)
% 
% Filter the fields of a fieldtrip data struct.
% 
% inputs
% ------
% data         - fieldtrip data struct
% filter       - the filter to apply to the relevant fields of "data". 
%                can be a vector indices to select, or a vector of logical values
% field_length - the length of the fields in data that will undergo filtering.
%                all other fields will remain as originally specified in data.
% outputs
% -------
% data_filtered - the filtered version of the data


data_filtered = data;

fields = fieldnames(data);
for i = 1:length(fields)
    if eval(['length(data.' fields{i} ') == field_length'])      
        eval(['data_filtered.' fields{i} ' = data_filtered.' fields{i} '(filter);'])
    end
end