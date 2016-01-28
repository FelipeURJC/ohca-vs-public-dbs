function show_types_of_rhythms_for_each_ddbb

% This function shows the different types of rhythms for 
% analyzed database
%
% this code, by Felipe Alonso-Atienza
% felipe.alonso@urjc.es

close all; clear all; clc;

db_names = {'vfdb', 'cudb', 'ahadb','ohcadb'};
w_length = [4,8];

for j = 1:length(w_length)
    
    % load data
    data_path = '../data/';
    filename  = sprintf('%sdata_%d',data_path, w_length(j));
    load(filename);
    
    for i = 1:length(db_names)
        
        msg = sprintf ('\nRESULTS FOR %s-%d seconds',...
            db_names{i},w_length(j));
        disp(msg)

        idx = strcmp(Tabla.ddbb,db_names{i});
        
        rhythms         = Tabla.labels(idx);
        type_of_rhythms = unique(rhythms);
        N               = length(rhythms);
                
        fprintf('%s\t%s\t%s\n','Name','Count','Percent')
        for k = 1:numel(type_of_rhythms)
            
            this_rhythm = type_of_rhythms(k);
           
            name    = label_names(this_rhythm,db_names{i});
            count   = sum( rhythms == this_rhythm );
            percent = (count / N) * 100;
            
            msg = sprintf('%s\t%d\t%2.2f',name,count,percent);
            disp(msg)
            
        end
        
    end
    
end


end


function label_name = label_names(id_label,db_name)

cases = {...
    '(AB';...       %1
    '(AFIB';...     %2
    '(AFL';...      %3
    '(ASYS';...     %4      
    '(B';...        %5
    '(BI';...       %6
    '(BII';...      %7
    '(HGEA';...     %8
    '(IVR';...      %9
    '(N';...        %10
    '(NOD';...      %11
    '(NOISE';...    %12     
    '(P';...        %13
    '(PREX';...     %14
    '(SBR';...      %15
    '(SVTA';...     %16
    '(T';...        %17
    '(VER';...      %18
    '(VF';...       %19
    '(VFL';...      %20
    '(VT';...       %21
    '(sTV';...      %22
    '(others';...   %23     
    '(fineVF'};     %24

if strcmp(db_name,'ohcadb') && (id_label == 10)
    label_name = '(PR'; 
elseif strcmp(db_name,'ohcadb') && (id_label == 23)
    label_name = '(PEA';
else
    label_name = cases{id_label};
end


end
