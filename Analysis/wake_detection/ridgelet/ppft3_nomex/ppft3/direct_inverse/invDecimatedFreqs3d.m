% function rim=invDecimatedFreqs3d(fim)
%
% This function restores the initial image from it's FFT, given on carthesian grid
%   
%                                           n/2-1   n/2-1   n/2-1
%             FI(ox,oy,oz)=  sum     sum     sum  I(u,v,w)exp(-2*pi*i*3*(u*ox+v*oy+w*oz)/m)   m=3n+1
%                                           u=-n/2  v=-n/2  w=-n/2
% where
%               ox, oy, oz = -n/2...n/2
%
% Oren 11/12/2005
function rim=invDecimatedFreqs3d(fim)
n = size(fim);

if (n(1)~=n(2))|(n(1)~=n(3))|(length(n)~=3)
    error('Input matrix must be of size (n+1)x(n+1)x(n+1)');
end

n=n(1)-1;
c=zeros(n,1);
r=zeros(1,n);
m = 3*n+1;
for k=-n/2:n/2-1
    for l=-n/2:n/2
        c(k+n/2+1)=c(k+n/2+1)+exp(6*pi*i*l*(-n/2-k)/m);
        r(k+n/2+1)=r(k+n/2+1)+exp(6*pi*i*l*(k+n/2)/m);
    end
end

[m1,m2,m3,m4]=topinv(c,r);
[D1,D2,D3,D4]=topprepinv(m1,m2,m3,m4);

aa=zeros(n+1,n+1,n);
aaa=zeros(n+1,n,n);
rim=zeros(n,n,n);

% apply inverse once
for k=1:n+1
    for l =1:n+1
        v=fim(k,l,:);
        v=adjF3d(v);
        u=topinvmul(D1,D2,D3,D4,v);
        aa(k,l,:)=reshape(u(:),1,n);
    end
end

% apply inverse twice
for k=1:n+1
    for l =1:n
        v=aa(k,:,l);
        v=adjF3d(v);
        u=topinvmul(D1,D2,D3,D4,v);
        aaa(k,:,l)=reshape(u(:),1,n);
    end
end

% apply inverse third time
for k=1:n
    for l =1:n
        v=aaa(:,k,l);
        v=adjF3d(v);
        u=topinvmul(D1,D2,D3,D4,v);
        rim(:,k,l)=reshape(u(:),1,n);
    end
end