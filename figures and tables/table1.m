function table1

% This function provides the information of table 1 of the paper
%
% this code, by Felipe Alonso-Atienza
% felipe.alonso@urjc.es

close all; clear all; clc;

data_path = '../data/';

db_names = {'vfdb', 'cudb', 'ahadb','ohcadb'};
w_length = [4,8];

k=1;
for j = 1:length(w_length)

    % load data
    filename  = sprintf('%sdata_%d',data_path, w_length(j));
    load(filename);
    
    for i = 1:length(db_names)
        idx = strcmp(Tabla.ddbb,db_names{i});
        Sh = sum(Tabla.y(idx)==1);
        NSh = sum(Tabla.y(idx)==-1);
        
        table_of_results(i,k) = Sh;
        table_of_results(i,k+1) = NSh;
    end
    
    k=3;
end

msg = sprintf('\n\t\t4-s segments\t8-s segments');
disp(msg)
msg = sprintf('Database\tSh\tNSh\tSh\tNSh');
disp(msg)
msg = sprintf('Public\t\t%d\t%d\t%d\t%d',...
    sum(table_of_results(1:3,1)),...
    sum(table_of_results(1:3,2)),...
    sum(table_of_results(1:3,3)),...
    sum(table_of_results(1:3,4)));
disp(msg)
msg = sprintf('  vfdb\t\t%d\t%d\t%d\t%d',table_of_results(1,:));
disp(msg)
msg = sprintf('  cudb\t\t%d\t%d\t%d\t%d',table_of_results(2,:));
disp(msg)
msg = sprintf('  ahadb\t\t%d\t%d\t%d\t%d',table_of_results(3,:));
disp(msg)
msg = sprintf('ohcadb\t\t%d\t%d\t%d\t%d',table_of_results(4,:));
disp(msg)