function [  ] = slices_curv_2d( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,numb_rand,slice,NSIDE,lev,sigma )

%(example) [ ] = slices_curv_2d('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,1,1,4,2,5 );

%(example) [ ] = slices_curv_2d('/home/asus/Dropbox/extras/storage/graham/ht/','/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024/','4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m','/sample3001/','','3.000xv0.dat',1,1,1,32,8,2,5 );


%addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
%addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));

 addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
 addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));


cd('../../preprocessing');

[~,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info(root,spec,aux_path );

% [ size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw,z,path_file_in,~ ] = preprocessing_part(root,spec,aux_path,filename,1024,1);


z_glob=str2num(filename(1:end-7))


nb=np*resol_factor/lenght_factor;

display(nb)

% filename='_2dproj_z3_data_sl';    

F = ones(nb);
X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
%X = F * sqrt(prod(size(F)));
%C = fdct_wrapping(X,0);
C = fdct_wrapping(X,0,2);
%C = fdct_wrapping(F,0);
E = cell(size(C));
for s=1:length(C)
    E{s} = cell(size(C{s}));
    for w=1:length(C{s})
        A = C{s}{w};
        E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
    end
end

display(E{1}{1})


cd('../wake_detection/slices_curv_2d/');


end

