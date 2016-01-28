function process_data

% This function processes the computed ECG parameters for each
% ddbb contained in data_path (vfdb_4, cudb_4, ahadb_4, ohca_4, 
% vfdb_8, cudb_8, ahadb_8, ohca_8) to create a matlab Table, 
% includes an patient_id, and stores all data in a single .mat file, 
% for each windown size (4/8 seconds). Results: two .mat files saved 
% in data_res (data_4; data_8). 

% this code, by Felipe Alonso-Atienza
% felipe.alonso@urjc.es

clear all; close all; clc;

VarNames = {'tci', 'tcsc', 'exp', 'expmod', 'cm', 'cvbin', 'frqbin',...
    'abin', 'kurt','vfleak', 'M', 'A1', 'A2', 'A3', 'mav', 'psr',...
    'hilb', 'SamEn', 'x1', 'x2', 'x3', 'x4','x5','bCP', 'bWT', 'bW', 'Li',...
    'count1', 'count2', 'count3','labels','y','patients_id','record_idx'};
              
data_path = '../data/computed parameters/';
data_res = '../data/';

db_names = {'vfdb', 'cudb', 'ahadb','ohcadb'};
w_length = [4, 8];


for j = 1:length(w_length)
    
    data        = [];
    data_labels = [];
    data_y      = [];
    data_id     = [];
    data_index  = [];
    
    offset    = 1;
    offset_id = 1;
    
    samples_for_dbs = struct;

    for i = 1:length(db_names)
        r_filename  = sprintf('%s%s_%d', data_path , db_names{i} ,w_length(j));
        
        load(r_filename);
        
        % Include a unique patient id in the Table
        if strcmp(db_names{i},'ohcadb') %OHCA
            
            %load patients_id
            id_filename  = sprintf('%sohca_id', data_path);
            load(id_filename)
                        
            BeginE = cumsum(Lr) - Lr + 1;
            EndE   = cumsum(Lr);
                        
            patients_id = set_patients_id_ohca(BeginE,EndE,ohca_id,offset_id);
                        
        else % the rest of the ddbbs
            
            BeginE = cumsum(Lr) - Lr + 1;
            EndE   = cumsum(Lr);
            
            [patients_id,offset_id] =...
                set_patients_id(BeginE,EndE,offset_id);
        end
        
        [X,labels,y,patients_id,no_r_t] = clean_data(X,y,patients_id);
        
        data        = [data; X]; 
        data_labels = [data_labels; labels];
        data_y      = [data_y; y];
        data_id     = [data_id; patients_id];
        data_index  = [data_index; find(no_r_t)];
        
        samples_for_dbs = setfield(samples_for_dbs,...
            db_names{i},offset : offset + length(y) - 1) ;  
        
        offset = offset + length(y);
    end
    
    data_ddbb = cell(size(data_y));
    data_ddbb(samples_for_dbs.vfdb) = {'vfdb'};
    data_ddbb(samples_for_dbs.cudb) = {'cudb'};
    data_ddbb(samples_for_dbs.ahadb) = {'ahadb'};
    data_ddbb(samples_for_dbs.ohcadb) = {'ohcadb'};
    
    Tabla = array2table([data data_labels data_y data_id data_index],...
        'VariableNames',VarNames);
    Tabla.ddbb = data_ddbb;
    
    output_file = sprintf('%sdata_%d',data_res, w_length(j));
    save(output_file,'Tabla','samples_for_dbs')
    
end

end

function patients_id = set_patients_id_ohca(BeginE,EndE,ohca_id,offset)

number_of_patients = size(BeginE,1); 
patients_id = zeros(EndE(end),1);

for i=1:number_of_patients
    id = ohca_id(i).patID;
    patients_id(BeginE(i):EndE(i)) = id + offset;
end

end

function [patients_id,offset] = set_patients_id(BeginE,EndE,offset)

number_of_patients = size(BeginE,1); 
patients_id = zeros(EndE(end),1);

id = offset;

for i=1:number_of_patients
    patients_id(BeginE(i):EndE(i)) = id;
    id = id + 1;
end
offset = id;
end

function [X,y,blabel,patients_id,no_r_t] = clean_data(X,y,patients_id)

% Remove ASY, fine VF and slow VT
r_t = (y==4) | (y==24) | (y==22);

no_r_t = ~r_t;

y(r_t) = [];
X(r_t,:)=[];
patients_id(r_t) = [];

% Select shock and noshock
r_s  = (y==19) | (y==20) | (y==21);

blabel = ones(size(y));
blabel(~r_s) = -1;

end