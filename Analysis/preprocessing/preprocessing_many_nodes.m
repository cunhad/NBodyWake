function [ nodes_list redshift_list  ] = preprocessing_many_nodes( path,spec,aux_path)
%   This functions pre selects part of the data from the cubic particle
%   distribution made by the cubep3m code. It also insert the wake at the
%   center if there is a wake. The last entry specifies the percentage of
%   the data we what to analyse (edges are removed)
%(example ) [ nodes_list redshift_list ] = preprocessing_many_nodes( '/home/acer/Documents/storage/guillimin/','32Mpc_96c_zi63_nowakes','/' )


 path_in=strcat(path,spec,aux_path);
% files_list1 = dir(strcat(path_in,'*xv*'));
% sorted_files_list1={files_list1.name};
% filename=char(sorted_files_list1(1));
% pat=char(filename);
% pat=pat(1:2);

files_list = dir(strcat(path_in,'*xv0.dat'));
sorted_files_list={files_list.name};
redshift_list=cellfun(@(x) x(1:end-7),sorted_files_list,'UniformOutput', false);
%this is not good for nodes with more than one decimal character?



files_list2 = dir(strcat(path_in,char(redshift_list(1)),'xv*','.dat'));
%display(strcat(path_in,char(redshift_list(1)),'xv*','.dat'));
sorted_files_list2={files_list2.name};
%nodes_list=cellfun(@(x) x(end-4:end-4),sorted_files_list2,'UniformOutput', false);
nodes_list=cellfun(@(x) x(2+cell2mat(strfind(sorted_files_list2(1:1), 'xv')):-1+cell2mat(strfind(sorted_files_list2(1:1), '.dat'))),sorted_files_list2,'UniformOutput', false);
cd('../processing');

nodes_list=sort_nat(nodes_list);
redshift_list=sort_nat(redshift_list);

cd('../preprocessing');

end

