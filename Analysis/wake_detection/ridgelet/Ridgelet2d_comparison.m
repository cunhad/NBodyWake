function [ output_args ] = Ridgelet2d_comparison( )


% Comparison of the various methods for the ridgelet transformation
% Only rectopolar apears to preserve gaussianity

% the cubic ridgelet with cubic interpolation appears to the the one with
% best discriminator for a line discontinuity on top of gaussian
% distribution

clearvars;

% Original data

nc=64;
ncx=nc;ncy=nc;
X=randn([ncx,ncy]);

Y=X;
% Y(:,nc/2)=1;

for i=1:floor(ncx/2)
%     Y(i+floor(ncx/2),i)=1;
%     Y(-1+i+floor(ncx/2),i)=1;
%     Y(-2+i+floor(ncx/2),i)=1;
    Y(i+floor(ncx/2),i)=1+Y(i+floor(ncx/2),i);
end

figure; imagesc(X); colorbar;
figure; imagesc(Y); colorbar;


figure; histogram(X(:));


% Test on original data to see if it comes from a normal distribution

% passed the test (as it should)

% COMP =    1.0016
% p =  0.5155

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


% COMP =    1.0016


[h,p,ksstat,cv] = kstest(X(:))

% p =  0.5155

figure;
cdfplot(X(:))
hold on
x_values = linspace(min(X(:)),max(X(:)));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')

% Test on ridgelet-transformed (rectopolar interpolation) data (passed the
% test) to see if it comes from a normal distribution

% passed the test

% COMP =    1.5014
% p =  0.8985

[ Radon_hor_, Radon_vert_] = Ridgelet2d_RP_crude_forwards_dev1(X);
data = real([Radon_hor_(:);Radon_vert_(:)]);

figure;imagesc(real(Radon_hor_)); colorbar;
figure;imagesc(real(Radon_vert_)); colorbar;

[ Radon_hor_, Radon_vert_] = Ridgelet2d_RP_crude_forwards_dev1(Y);
dataW = real([Radon_hor_(:);Radon_vert_(:)]);

figure;imagesc(real(Radon_hor_)); colorbar;
figure;imagesc(real(Radon_vert_)); colorbar;




P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN

% COMP =    1.5014

    
figure; histogram(data);

data_norm=(data-mean(data))/std(data);

% [h,p,ksstat,cv] = kstest(data)
[h,p,ksstat,cv] = kstest(data_norm)

% p =  0.8985

figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')



% Test on ridgelet-transformed (rectopolar interpolation, with zero padding) data  
% to see if it comes from a normal distribution

% passed the test

% COMP =    1.6459
% p =    0.6872

[ Radon_hor_, Radon_vert_] = Ridgelet2d_RP_crude_forwards_dev2(X);
% data = real([Radon_hor_(:);Radon_vert_(:)]);
data = real([Radon_hor_(:,ncy+1:2*ncx);Radon_vert_(:,ncx+1:2*ncx)]);
data = data(:);

% figure;imagesc(real(Radon_vert_)); colorbar;
figure;imagesc(real(Radon_hor_(:,ncy+1:2*ncx))); colorbar;
figure;imagesc(real(Radon_vert_(:,ncx+1:2*ncx))); colorbar;



[ Radon_hor_, Radon_vert_] = Ridgelet2d_RP_crude_forwards_dev2(Y);
% dataW = real([Radon_hor_(:);Radon_vert_(:)]);
dataW = real([Radon_hor_(:,ncy+1:2*ncx);Radon_vert_(:,ncx+1:2*ncx)]);
dataW = dataW(:);


% figure;imagesc(real(Radon_vert_)); colorbar;
figure;imagesc(real(Radon_hor_(:,ncy+1:2*ncx))); colorbar;
figure;imagesc(real(Radon_vert_(:,ncx+1:2*ncx))); colorbar;


P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN

% COMP =    1.6459


figure; histogram(data);

data_norm=(data-mean(data))/std(data);

% [h,p,ksstat,cv] = kstest(data)
[h,p,ksstat,cv] = kstest(data_norm)

% p =    0.6872
    
figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')


% Test on ridgelet-transformed (rectopolar interpolation, with zero padding, and normalized) data  
% to see if it comes from a normal distribution

% passed the test

% COMP =   1.3048
% p =    0.3469


[ Radon_hor_ones_, Radon_vert_ones_] = Ridgelet2d_RP_crude_forwards_dev2(ones(size(X)));


[ Radon_hor_, Radon_vert_] = Ridgelet2d_RP_crude_forwards_dev2(X);
% data = real([Radon_hor_(:);Radon_vert_(:)]);
Radon_hor_norm_ = Radon_hor_./Radon_hor_ones_;
Radon_vert_norm_ = Radon_vert_./Radon_vert_ones_;
data = real([Radon_hor_norm_(:,ncy+1:2*ncx);Radon_vert_norm_(:,ncx+1:2*ncx)]);
data = data(:);

% figure;imagesc(real(Radon_vert_)); colorbar;
figure;imagesc(real(Radon_hor_norm_(:,ncy+1:2*ncx))); colorbar;
figure;imagesc(real(Radon_vert_norm_(:,ncx+1:2*ncx))); colorbar;



[ Radon_hor_, Radon_vert_] = Ridgelet2d_RP_crude_forwards_dev2(Y);
% dataW = real([Radon_hor_(:);Radon_vert_(:)]);
Radon_hor_norm_ = Radon_hor_./Radon_hor_ones_;
Radon_vert_norm_ = Radon_vert_./Radon_vert_ones_;
dataW = real([Radon_hor_norm_(:,ncy+1:2*ncx);Radon_vert_norm_(:,ncx+1:2*ncx)]);
dataW = dataW(:);


% figure;imagesc(real(Radon_vert_)); colorbar;
figure;imagesc(real(Radon_hor_norm_(:,ncy+1:2*ncx))); colorbar;
figure;imagesc(real(Radon_vert_norm_(:,ncx+1:2*ncx))); colorbar;


P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN

% COMP =   1.3048

figure; histogram(data);


data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    0.3469
    
    
figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')




% Test on ridgelet-transformed (cubic interpolation) data to see if it
% comes from a normal distribution 

% passed the test

% COMP =    1.6708
% p =    0.2924


[ Radon_hor_, Radon_vert_,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy] = Ridgelet2d_RP_crude_forwards_dev3(X);
data = real([Radon_hor_(1,:)';Radon_vert_(1,:)']);

Ridgelet2d_RP_crude_visual_dev3(Radon_hor_, Radon_vert_,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy);


[ Radon_hor_, Radon_vert_,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy] = Ridgelet2d_RP_crude_forwards_dev3(Y);
dataW = real([Radon_hor_(1,:)';Radon_vert_(1,:)']);

Ridgelet2d_RP_crude_visual_dev3(Radon_hor_, Radon_vert_,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy);





P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN

% COMP =    1.6708

figure; histogram(data);


data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    0.2924

figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')




% Test on ridgelet-transformed (cubic interpolation and normalized) data  
% to see if it comes from a normal distribution

% passed the test

% COMP =   1.5628
% p =    0.6721    


[ Radon_hor_ones_, Radon_vert_ones_,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy] = Ridgelet2d_RP_crude_forwards_dev3(ones(size(X)));



[ Radon_hor_, Radon_vert_,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy] = Ridgelet2d_RP_crude_forwards_dev3(X);
Radon_hor_norm_ = Radon_hor_./Radon_hor_ones_;
Radon_vert_norm_ = Radon_vert_./Radon_vert_ones_;
data = real([Radon_hor_norm_(1,:)';Radon_vert_norm_(1,:)']);

Ridgelet2d_RP_crude_visual_dev3(Radon_hor_norm_, Radon_vert_norm_,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy);



[ Radon_hor_, Radon_vert_,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy] = Ridgelet2d_RP_crude_forwards_dev3(Y);
Radon_hor_norm_ = Radon_hor_./Radon_hor_ones_;
Radon_vert_norm_ = Radon_vert_./Radon_vert_ones_;
dataW = real([Radon_hor_norm_(1,:)';Radon_vert_norm_(1,:)']);

Ridgelet2d_RP_crude_visual_dev3(Radon_hor_norm_, Radon_vert_norm_,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy);


P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN

% COMP =   1.5628

figure; histogram(data);


data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    0.6721    
    
figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')





end