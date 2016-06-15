%required vars is id and kmax
xFileName= strcat('x_',id, '.csv');

x = csvread(xFileName);

DEMEAN=2;
jj=1;

T=size(x,1);
N=size(x,2);
NT=N*T;
NT1=N+T;
CT=zeros(1,kmax);
ii=1:1:kmax;
if jj ==1;CT(1,:)=log(NT/NT1)*ii*NT1/NT;end;
if jj==2; CT(1,:)=(NT1/NT)*log(min([N;T]))*ii;end;
GCT=min([N;T]);
if jj==3; CT(1,:)=ii*log(GCT)/GCT; end;
if jj==4; CT(1,:)=2*ii/T; end;
if jj==5; CT(1,:)=log(T)*ii/T;end;
if jj==6; CT(1,:)=2*ii*NT1/NT; end;
if jj==7; CT(1,:)=log(NT)*ii*NT1/NT;end;

 if DEMEAN ==2;
 X=standard(x);
 end;

if DEMEAN ==1;
 X=demean(x);
 end;
if DEMEAN==0;
  X=x;
  end;

IC1=zeros(size(CT,1),kmax+1);
Sigma=zeros(1,kmax+1);
XX=X*X';
[Fhat0,eigval,Fhat1]=svd(XX');
for i=kmax:-1:1;
Fhat=Fhat0(:,1:i);
lambda=Fhat'*X;
chat=Fhat*lambda;
ehat=X-chat;
E=diag(eigval);
E=sort(E,'ascend');
E=E/N;
m=i;
p=T;
c=p/N;
sigma1=sum(E(1:(p-m)))/(p-m);
	
	alpha1=zeros(m,1);
	for j=1:m
	    alpha1(j)=0.5*(sqrt((sigma1*(c+1)-E(p-j+1))^2-4*c*sigma1^2)-sigma1*(c+1)+E(p-j+1));
	end
	
	I1=zeros(m,1);
	for k=1:m
	    I1(k)=1/alpha1(k);
	end
	
	b=sqrt(c/2)*(m+sigma1*sum(I1));
	
	Sigma(i)=sigma1+(b*sigma1*sqrt(2*c))/(p-m);

IC1(:,i)=Sigma(i)+CT(:,i)*Sigma(kmax);
end;
Sigma(kmax+1)=mean(sum(X.*X/T));
IC1(:,kmax+1)=Sigma(kmax+1);
ic1=minindc(IC1')';
ic1=ic1 .*(ic1 <= kmax);
Fhat=[];
Fhat=Fhat0(:,1:ic1);
lambda=Fhat'*X;
chat=Fhat*lambda;

%%
resFileName= strcat('res_',id, '.csv');
fileID = fopen(resFileName,'w');
k = ic1;
fprintf(fileID, 'k, sigma\n');
fprintf(fileID, '%f, %f\n', [k sigma1]);
fclose(fileID);
