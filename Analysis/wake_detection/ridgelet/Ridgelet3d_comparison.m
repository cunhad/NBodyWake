function [] = Ridgelet3d_comparison()


% Two methods will be used


clearvars;

% Original data

nc=8;
ncx=nc;ncy=nc;ncz=nc;
X=randn([ncx,ncy,ncz]);


Y=X;
% Y(:,nc/2)=1;

for i=1:ncx
    for j=1:ncy
    z = (i+j)*(ncz/(ncx+ncy));
    Y(i,j,floor(z))=1+Y(i,j,floor(z));
    end
end

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




[h,p,ksstat,cv] = kstest(X(:))

% p =  0.8367

figure;
cdfplot(X(:))
hold on
x_values = linspace(min(X(:)),max(X(:)));
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

[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3(X);
data = real(Radon_Z_(1,:));

Ridgelet3d_interp_visual_dev3(Radon_Z_,interp_Z_info);


[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3(Y);
dataW = real(Radon_Z_(1,:));

Ridgelet3d_interp_visual_dev3(Radon_Z_,interp_Z_info);




P = max(data)
P_W = max(dataW)


M = mean(data)
M_W = mean(dataW)


STD = std(data)
STD_W = std(dataW)


STN = P/STD
STN_W = P_W/STD_W

COMP = STN_W/STN

% COMP =    1.3578



figure; histogram(data);
figure; histogram(dataW);

data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    0.9189

figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')



% Test on ridgelet-transformed (cubic interpolation and normalized) data  
% to see if it comes from a normal distribution

[ Radon_Z_ones,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3(ones(size(X)));


[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3(X);
Radon_Z_norm_ = Radon_Z_./Radon_Z_ones;
data = real(Radon_Z_norm_(1,:));


[ Radon_Z_,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3(Y);
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

% COMP =   1.2903

figure; histogram(data);
figure; histogram(dataW);


data_norm=(data-mean(data))/std(data);

[h,p,ksstat,cv] = kstest(data_norm)

% p =    0.8257    
    
figure;
cdfplot(data_norm)
hold on
x_values = linspace(min(data_norm),max(data_norm));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')



end
