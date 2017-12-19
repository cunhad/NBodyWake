function [ i_node,j_node,k_node,number_node_dim,z ] = preprocessing_filename_info( path,spec,aux_path,filename)
%   This function takes the phase space output from CUBEP3M and returns
%   all relevant information but the global positions and velocity of all particles at the
%   corresponding redshift and node volume

%(example)[ i_node,j_node,k_node,number_node_dim,z ] = preprocessing_filename_info('/home/asus/Dropbox/extras/storage/graham/small_res/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','63.000xv0.dat' );
%(example)[ i_node,j_node,k_node,number_node_dim,z ] = preprocessing_filename_info('/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','65.000xv0.dat');



%   Detailed explanation:


%node structure

 node=char(filename);
 node=str2num(node(strfind(filename, 'xv')+2:strfind(filename,'.dat')-1));
 [ nodes_list redshift_list ] = preprocessing_many_nodes(path,spec,aux_path);
 
 number_node_dim=nthroot(numel(nodes_list), 3);
 k_node=floor(node/number_node_dim^2);
 res=mod(node,number_node_dim^2);
 j_node=floor(res/number_node_dim);
 i_node=mod(res,number_node_dim);
 
   %extracts the redshift of the file to be analised
 
 z=char(filename);
 z=str2num(z(1:end-7));
 


end

