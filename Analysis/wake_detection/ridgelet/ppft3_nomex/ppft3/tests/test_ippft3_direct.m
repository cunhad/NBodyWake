% function testippftd3_direct
%
% Tests the function ippftd3_direct.
%
% Tests the correctness and performance of the inversion algorithm of the 3-D discrete ppft transform.
% The function tests also the inversion of the 3-D pseudo-polar Fourier transform since 
% inverting the 3-D discrete ppft transform requirers inverting the 3-D pseudo-polar Fourier transform.
%
% Yoel Shkolnisky 01/03/03

function test_ippft3_direct



for n = [8 32 64 128]
    a = randn(n,n,n);
    pp = ppft3(a);
    runTest(pp,1e-15,a,sprintf('Real random cube %dx%dx%d from [0,1]',n,n,n),1,0);
end


%%%%%%%%%%%%%%%
% Sub functions
%%%%%%%%%%%%%%%

% Execute a single inversion test.
% pp           The 3-D ppft sectors.
% ErrTol       Residual error required from the inversion algorithm.
% MaxIts       Number of iterations of the inversion algorithm.
% ref          The original image. Used as a reference to check the absolute error.
% description  Test description for printing purposes.
% verbose      If true, print the inversion log.
% trueimage    True if the input represents the discrete ppft transform of an integer image.
%              In this case it is possible to exactly inverting the transform
%              by truncating any floating parts of the result.
function runTest(pp,ErrTol,ref,description,verbose,trueimage)

fprintf('Test name : %s\n',description);
fprintf('Requested error = %-2.5e\n',ErrTol);
if verbose
    fprintf('Inversion log:\n');
   fprintf('--------------\n');
end
tic;
Y = ippft3_direct(pp,ErrTol);
t = toc;

fprintf('\nResults:\n');
fprintf('---------\n');

diff = abs(Y-ref);
fprintf('L_infinity error = %e\n',max(diff(:)));

if trueimage
    diff = abs(round(Y)-ref);
    fprintf('L_infinity error of reconstructed image = %-1.3f\n',max(diff(:)));
end
   
fprintf('Computation time = %-3.2f secs\n',t);
fprintf('--------------------------------------------------------\n\n\n');
