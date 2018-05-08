function [  ] = healpix_angles_largearrays( nside )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



% addpath('/home/asus/Programs/s2let/src/main/matlab','/home/asus/Programs/ssht/src/matlab','/home/asus/Programs/so3/src/matlab','/home/asus');

% addpath('/home/asus/Programs/s2let/src/main/matlab','/home/asus/Programs/ssht/src/matlab','/home/asus/Programs/so3/src/matlab');
tic;

% addpath('/home/asus/Programs/s2let/src/main/matlab');



% %prelocate memory not recomended, because if one uses a function we endup
% alocating twice
% 
% thetas=zeros(N_side);
% 
% phis=zeros(N_side);
% 

% [thetas, phis] = s2let_hpx_sampling_ring(N_side);

% csvwrite(strcat('/home/asus/Documents/angles',num2str(N_side),'_t.cvs'),thetas,'single');
% csvwrite(strcat('/home/asus/Documents/angles',num2str(N_side),'_p.cvs'),phis,'single');

% csvwrite(strcat('/home/asus/Documents/angles',num2str(N_side),'_t.cvs'),thetas);
% csvwrite(strcat('/home/asus/Documents/angles',num2str(N_side),'_p.cvs'),phis);
% 
% % %problem with nside=16384
% % 

npix = 12 * nside^2;


% % %thetas = zeros(npix,1);
% % % thetas = sparse(npix,1);
% % %     ipix = 1:npix;
% % 
% % %this saves to a file, but it is very slow
% % 
% % % pw=java.io.PrintWriter(java.io.FileWriter('/home/asus/Documents/angles.txt'));
% % % line=num2str(0:size(thetas,2)-1);
% % % pw.println(line);
% % % for index=1:length(thetas)
% % %     disp(index);
% % %     line=num2str(full(thetas(index,:)));
% % %     pw.println(line);
% % % end
% % % pw.flush();
% % % pw.close();
% 
% %We could also save everithing in mat format
% % 
% % fileName = 'ipix.mat';
% % matObj = matfile(fileName);
% % matObj.Properties.Writable = true;
% % 
% % size = npix;
% % chunk = 1000000;
% % 
% % nout = 0;
% % 
% % while(nout < size)
% %     fprintf('Writing %d of %d\n',nout,size);
% %     chunkSize = min(chunk,size-nout);
% %     data = (nout+1):(nout+chunkSize);
% %     matObj.data(1,(nout+1):(nout+chunkSize)) = data;
% %     nout = nout + chunkSize;
% % end

%try saving in  text format

nout = 0;
chunk=1000000;
size=npix;

while(nout < size)
    fprintf('Writing %d of %d\n',nout,size);
    chunkSize = min(chunk,size-nout);  
    dlmwrite('ipix.txt',[(nout+1):(nout+chunkSize)]','-append') ;
    nout = nout + chunkSize;
end

%now hte processing starts

% ds = datastore('ipix.mat');

% fds = fileDatastore('ipix.mat','ReadFcn',@load,'FileExtensions','.mat');
% tt = tall(fds);

% T = table([1:192]')
% writetable(T)


toc;
end

            