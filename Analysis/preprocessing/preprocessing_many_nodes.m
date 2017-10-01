function [ nodes_list redshift_list  ] = preprocessing_many_nodes( path,spec,aux_path)
%   This function takes the phase space output from CUBEP3M and returns
%   all the list of nodes and redshifts produced by that simulation


%(example ) [ nodes_list redshift_list ] = preprocessing_many_nodes( '/home/asus/Dropbox/extras/storage/guillimin/old/','32Mpc_96c_48p_zi63_nowakes','/','63.000xv0.dat' )
%(example ) [ nodes_list redshift_list ] = preprocessing_many_nodes( '/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','65.000xv0.dat')


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

