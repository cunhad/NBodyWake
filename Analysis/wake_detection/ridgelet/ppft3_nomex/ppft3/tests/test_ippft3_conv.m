% function testippft3_gpu
%
% Compare the performance for the functions ippft3_conv and ippft3_gpu.
%
% Yoel Shkolnisky 10/04/2014

function test_ippft3_conv

sz=[16 32 50 64 80 100 128];

% Precompute all filters, so time measurements of the inversion won't be
% affected by the preprocessing of generating the filters.
for n=sz
    
    t_single=-1; % In case we don't need to generate the filters, this 
    t_double=-1; % indicates that filter was loaded and not generated.
    
    precision='single';
    filtname=getppfiltname(n,precision);
    H=loadppftfilter(filtname);
    if ~isstruct(H)
        tic;
        fprintf('Filter for n=%d (%s) not found. Generating filter.\n',n,precision);
        H=makeppftfilter(n,1,precision);
        saveppftfilter(H,filtname);
        t_single=toc;
    end
    
    precision='double';
    filtname=getppfiltname(n,precision);
    H=loadppftfilter(filtname);
    if ~isstruct(H)
        tic;
        fprintf('Filter for n=%d (%s) not found. Generating filter.\n',n,precision);
        H=makeppftfilter(n,1,precision);
        saveppftfilter(H,filtname);
        t_double=toc;
    end
    
   if t_single~=-1 || t_double~=-1
    fprintf('Generating filters for n=%d  t_single=%6.3f  t_double=%6.3f\n',n,t_single,t_double);    
   end
end

% Run timing tests
fprintf('n \t err_conv \t \t t_conv \n');
for n=sz
    % Test the function ppft3 by comparing it to ppft3_ref.
    im = rand(n,n,n);
    pp = ppft3(im);
    runTest(n,pp,1.e-8,100,im,0);
end


%%%%%%%%%%%%%%%
% Sub functions
%%%%%%%%%%%%%%%

% Execute a single inversion test.
% pp           The 3-D Radon sectors.
% ErrTol       Residual error required from the inversion algorithm.
% MaxIts       Number of iterations of the inversion algorithm.
% ref          The original image. Used as a reference to check the absolute error.
% description  Test description for printing purposes.
% verbose      If true, print the inversion log.
function runTest(n,pp,ErrTol,MaxIts,ref,verbose)


tic;
[Y1,~,~,~] = ippft3_conv(pp,ErrTol,MaxIts,verbose);
t1=toc;

err_conv=norm(Y1(:)-ref(:))/norm(ref(:));

fprintf('%d \t %e \t %5.2f \n', n,err_conv, t1);

   