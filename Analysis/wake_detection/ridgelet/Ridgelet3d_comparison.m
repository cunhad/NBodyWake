function [] = Ridgelet3d_comparison()


% Two methods will be used


clearvars;

% Original data

nc=32;
ncx=nc;ncy=nc;ncz=nc;

% ncx=64;ncy=16;ncz=64;

X=randn([ncx,ncy,ncz]);


Y=X;

% Y(:,nc/2)=1;

% Y(:,:,ncz/2)=1;

% Parallel

for i=1:ncx
    for j=1:ncy
    z = ncz/2;
    Y(i,j,z)=1+Y(i,j,z);
    end
end

% % Diagonal
% 
% for i=1:ncx
%     for j=1:ncy
%     z = (i+j)*(ncz/(ncx+ncy));
%     Y(i,j,floor(z))=1+Y(i,j,floor(z));
%     end
% end

% Conserve the total amount

Y = Y-(1/ncz); 



volshow(X,'Renderer','MaximumIntensityProjection');
volshow(Y,'Renderer','MaximumIntensityProjection');


% Test on original data to see if it comes from a normal distribution


data = real(X(:));
dataW = real(Y(:));


P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN



figure; histogram(data);
figure; histogram(dataW);



figure;
% h1 = histogram(signal_sample_nw(:),'BinWidth',10);
h1 = histogram(data);
hold on
% h2 = histogram(signal_sample_w(:),'BinWidth',10);
h2 = histogram(dataW,'BinWidth',h1.BinWidth);
% xlabel('$S$ value', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('histogram', 'interpreter', 'latex', 'fontsize', 20);
% legend('G\mu=0','G\mu=1 \times 10^{-7}','location','northeast')
set(gca, 'YScale', 'log')



[h,p,ksstat,cv] = kstest(X(:))

% p =  0.8367

figure;
cdfplot(X(:))
hold on
x_values = linspace(min(X(:)),max(X(:)));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')




% Test on radon-transformed (cubic interpolation) data  
% to see if it comes from a normal distribution

% [ Radon_Z_ones,interp_Z_info,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(ones(size(X)));


% [ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(X);
[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(repelem(X,2,2,2));
[Radon_Z_] = RemoveLowVolRad_dev3(Radon_Z_,interp_Z_info,1/2);
% Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
% [ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,interp_Z_info,1);
% data = real(Radon_Z_norm_(1,:));
data = real(Radon_Z_(1,:));
% data = (Radon_Z_(1,:));


% [ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(Y);
[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(repelem(Y,2,2,2));
[Radon_Z_] = RemoveLowVolRad_dev3(Radon_Z_,interp_Z_info,1/2);
% Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
% [ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,interp_Z_info,1);
% dataW = real(Radon_Z_norm_(1,:));
dataW = real(Radon_Z_(1,:));
% dataW = real(Radon_Z_(1,:));



P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN
COMP2 = P_W/P
COMP3 = sum(dataW(dataW>4*STD_W))/sum(data(data>4*STD))


% COMP =   1.4589
% COMP2 =   1.7288

figure; histogram(data);
figure; histogram(dataW);



figure;
% h1 = histogram(signal_sample_nw(:),'BinWidth',10);
h1 = histogram(data);
hold on
% h2 = histogram(signal_sample_w(:),'BinWidth',10);
h2 = histogram(dataW,'BinWidth',h1.BinWidth);
% xlabel('$S$ value', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('histogram', 'interpreter', 'latex', 'fontsize', 20);
% legend('G\mu=0','G\mu=1 \times 10^{-7}','location','northeast')
set(gca, 'YScale', 'log')


figure;
plot(data(:))
hold on
plot(dataW(:))


data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    8.1654e-06     (not passed)
    
figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')


% Test on radon-transformed (cubic interpolation, normalized) data  
% to see if it comes from a normal distribution

[ Radon_Z_ones,interp_Z_info,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(ones(size(X)));


[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(X);
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
% [ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,interp_Z_info,1);
data = real(Radon_Z_norm_(1,:));
% data = real(Radon_Z_(1,:));


[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(Y);
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
% [ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,interp_Z_info,1);
dataW = real(Radon_Z_norm_(1,:));
% dataW = real(Radon_Z_(1,:));



P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN
COMP2 = P_W/P


% COMP =   1.7372
% COMP2 =   2.0445


figure; histogram(data);
figure; histogram(dataW);



figure;
% h1 = histogram(signal_sample_nw(:),'BinWidth',10);
h1 = histogram(data);
hold on
% h2 = histogram(signal_sample_w(:),'BinWidth',10);
h2 = histogram(dataW,'BinWidth',h1.BinWidth);
% xlabel('$S$ value', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('histogram', 'interpreter', 'latex', 'fontsize', 20);
% legend('G\mu=0','G\mu=1 \times 10^{-7}','location','northeast')
set(gca, 'YScale', 'log')



data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    0.8257    (passed)
    
figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')





% Test on ridgelet-transformed (cubic interpolation) data to see if it
% comes from a normal distribution 

% % % passed the test
% % 
% % % COMP =    1.6708
% % % p =    0.2924
% % 
% % 

[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3(X,1);
data = real(Radon_Z_(1,:));

% Ridgelet3d_interp_visual_dev3(Radon_Z_,interp_Z_info);


[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3(Y,1);
dataW = real(Radon_Z_(1,:));

% Ridgelet3d_interp_visual_dev3(Radon_Z_,interp_Z_info);




P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN
COMP2 = P_W/P

% COMP =    1.8984
% COMP2 =    2.0271



figure; histogram(data);
figure; histogram(dataW);


figure;
% h1 = histogram(signal_sample_nw(:),'BinWidth',10);
h1 = histogram(data);
hold on
% h2 = histogram(signal_sample_w(:),'BinWidth',10);
h2 = histogram(dataW,'BinWidth',h1.BinWidth);
% xlabel('$S$ value', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('histogram', 'interpreter', 'latex', 'fontsize', 20);
% legend('G\mu=0','G\mu=1 \times 10^{-7}','location','northeast')
set(gca, 'YScale', 'log')





data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    3.2735e-07 (nope)

figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')



% Test on ridgelet-transformed (cubic interpolation and normalized - rid/rad1) data  
% to see if it comes from a normal distribution

lev_3drig = 1;

[ Radon_Z_ones,interp_Z_info,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(ones(size(X)));


[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3(X,lev_3drig);
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
data = real(Radon_Z_norm_(1,:));


[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3(Y,lev_3drig);
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
dataW = real(Radon_Z_norm_(1,:));



P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN
COMP2 = P_W/P

% COMP =   1.4224
% COMP2 =   1.5148


figure; histogram(data);
figure; histogram(dataW);



figure;
% h1 = histogram(signal_sample_nw(:),'BinWidth',10);
h1 = histogram(data);
hold on
% h2 = histogram(signal_sample_w(:),'BinWidth',10);
h2 = histogram(dataW,'BinWidth',h1.BinWidth);
% xlabel('$S$ value', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('histogram', 'interpreter', 'latex', 'fontsize', 20);
% legend('G\mu=0','G\mu=1 \times 10^{-7}','location','northeast')
set(gca, 'YScale', 'log')





data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    1.2860e-06 (nope)   
    
figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')




% test to see if radon + radon_to_rid is equal to rid (OK)

% [ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3(X,1);
% data = real(Radon_Z_(1,:));
% 
% 
% [ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(X);
% [ Ridgelet_Z__ ] = Ridgelet3d_fromRadon_dev3( Radon_Z_,interp_Z_info,1);
% data_comp = real(Ridgelet_Z__(1,:));
% 
% max(abs(data(:)-data_comp(:)))





% Test on ridgelet-transformed (cubic interpolation and normalized - rid(normalz_rad)) data  
% to see if it comes from a normal distribution

lev_3drig = 1;


[ Radon_Z_ones,interp_Z_info,sample_points_Z_id,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(ones(size(X)));


[ Radon_Z_,interp_Z_info,sample_points_Z_id,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(X);
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
[ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,interp_Z_info,lev_3drig);
% [ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,sample_points_Z_id,interp_Z_info,lev_3drig);
% Ridgelet_Z_norm(isnan(Ridgelet_Z_norm)) = [];
data = real(Ridgelet_Z_norm(1,:));


[ Radon_Z_,interp_Z_info,sample_points_Z_id,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(Y);
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
[ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,interp_Z_info,lev_3drig);
% [ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,sample_points_Z_id,interp_Z_info,lev_3drig);
% Ridgelet_Z_norm(isnan(Ridgelet_Z_norm)) = [];
dataW = real(Ridgelet_Z_norm(1,:));



P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN
COMP2 = P_W/P

% COMP =   1.4147
% COMP2 =   1.5081

figure; histogram(data);
figure; histogram(dataW);



figure;
% h1 = histogram(signal_sample_nw(:),'BinWidth',10);
h1 = histogram(data);
hold on
% h2 = histogram(signal_sample_w(:),'BinWidth',10);
h2 = histogram(dataW,'BinWidth',h1.BinWidth);
% xlabel('$S$ value', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('histogram', 'interpreter', 'latex', 'fontsize', 20);
% legend('G\mu=0','G\mu=1 \times 10^{-7}','location','northeast')
set(gca, 'YScale', 'log')



figure;
plot(data(:))
hold on
plot(dataW(:))



data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    1.1643e-06   
    
figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')







% Test on ridgelet-transformed (cubic interpolation and normalized - rid(normalz_rad)) data (Oversampled) 
% to see if it comes from a normal distribution

lev_3drig = 2;

% repelem(X,2,2,2);
% repelem(Y,2,2,2);

[ Radon_Z_ones,interp_Z_info,sample_points_Z_id,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(ones(size(repelem(X,2,2,2))));


[ Radon_Z_,interp_Z_info,sample_points_Z_id,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(repelem(X,2,2,2));
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
[ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,interp_Z_info,lev_3drig);
% [ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,sample_points_Z_id,interp_Z_info,lev_3drig);
% Ridgelet_Z_norm(isnan(Ridgelet_Z_norm)) = [];
data = real(Ridgelet_Z_norm(1,:));


[ Radon_Z_,interp_Z_info,sample_points_Z_id,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(repelem(Y,2,2,2));
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
[ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,interp_Z_info,lev_3drig);
% [ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,sample_points_Z_id,interp_Z_info,lev_3drig);
% Ridgelet_Z_norm(isnan(Ridgelet_Z_norm)) = [];
dataW = real(Ridgelet_Z_norm(1,:));



P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN
COMP2 = P_W/P
COMP3 = sum(dataW(dataW>4*STD_W))/sum(data(data>4*STD))



% COMP =   1.4147
% COMP2 =   1.5081

figure; histogram(data);
figure; histogram(dataW);



figure;
% h1 = histogram(signal_sample_nw(:),'BinWidth',10);
h1 = histogram(data);
hold on
% h2 = histogram(signal_sample_w(:),'BinWidth',10);
h2 = histogram(dataW,'BinWidth',h1.BinWidth);
% xlabel('$S$ value', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('histogram', 'interpreter', 'latex', 'fontsize', 20);
% legend('G\mu=0','G\mu=1 \times 10^{-7}','location','northeast')
set(gca, 'YScale', 'log')



data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    1.1643e-06   
    
figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')




% Test on ridgelet-transformed (cubic interpolation and normalized - rid(normalz_rad), closest to the paper as possible) data  
% to see if it comes from a normal distribution

lev_3drig = 2;
frac_cut=1/2;

[ Radon_Z_ones,interp_Z_info,sample_points_Z_id,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(ones(size(repelem(X,2,2,2))));


[ Radon_Z_,interp_Z_info,sample_points_Z_id,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(repelem(X,2,2,2));
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
[ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,interp_Z_info,lev_3drig);
[Ridgelet_Z_norm] = RemoveLowVolRad_dev3(Ridgelet_Z_norm,interp_Z_info,frac_cut);
% [ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,sample_points_Z_id,interp_Z_info,lev_3drig);
% Ridgelet_Z_norm(isnan(Ridgelet_Z_norm)) = [];
data = real(Ridgelet_Z_norm(1,:));


[ Radon_Z_,interp_Z_info,sample_points_Z_id,ncx,ncy,ncz] = Radon3d_interp_forwards_dev3(repelem(Y,2,2,2));
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
[ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,interp_Z_info,lev_3drig);
[Ridgelet_Z_norm] = RemoveLowVolRad_dev3(Ridgelet_Z_norm,interp_Z_info,frac_cut);
% Ridgelet_Z_norm(real(Radon_Z_ones)<ncx*ncy*frac_cut)=[];
% [ Ridgelet_Z_norm ] = Ridgelet3d_fromRadon_dev3( Radon_Z_norm_,sample_points_Z_id,interp_Z_info,lev_3drig);
% Ridgelet_Z_norm(isnan(Ridgelet_Z_norm)) = [];
dataW = real(Ridgelet_Z_norm(1,:));



P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN
COMP2 = P_W/P
COMP3 = sum(dataW(dataW>4*STD_W))/sum(data(data>4*STD))

% COMP =   1.4147
% COMP2 =   1.5081

figure; histogram(data);
figure; histogram(dataW);



figure;
% h1 = histogram(signal_sample_nw(:),'BinWidth',10);
h1 = histogram(data);
hold on
% h2 = histogram(signal_sample_w(:),'BinWidth',10);
h2 = histogram(dataW,'BinWidth',h1.BinWidth);
% xlabel('$S$ value', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('histogram', 'interpreter', 'latex', 'fontsize', 20);
% legend('G\mu=0','G\mu=1 \times 10^{-7}','location','northeast')
set(gca, 'YScale', 'log')



figure;
plot(data(:))
hold on
plot(dataW(:))



data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    1.1643e-06   
    
figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')




end
