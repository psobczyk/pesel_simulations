xFileName= strcat('temp_data/x_',id, '.csv');

Y = csvread(xFileName);

p = size(Y,2);
n = size(Y,1);

if p>n
  Y = Y';
  p = size(Y,2);
n = size(Y,1);
end

[U,d,V] = svd(Y,0);
d=diag(d);
%%
R_W=zeros(maxPC,1);
k=1;
integrand = @(z) exp(-(z.^2/(2*sigma^2))).*(z.^(n-p)).*prod(abs(repmat(z.^2,length(d)-1,1)-repmat(d(1:end ~= k).^2, 1, length(z))));

a = integral(integrand, d(1), Inf);
c = integral(integrand, d(2), Inf);

R_W(k)=a/c;

for k=2:(maxPC)

  integrand = @(z) exp(-(z.^2/(2*sigma^2))).*(z.^(n-p)).*prod(abs(repmat(z.^2,length(d)-1,1)-repmat(d(1:end ~= k).^2, 1, length(z))));

  a = integral(integrand, d(k), d(k-1));
  c = integral(integrand, d(k+1), d(k-1));
    R_W(k)=a/c;
end

%%
resFileName= strcat('temp_data/res_',id, '.csv');
fileID = fopen(resFileName,'w');
k = find(R_W>0.05, 1, 'first')-1;
fprintf(fileID, 'k, sigma\n');
fprintf(fileID, '%f, %f\n', [k sigma]);
fclose(fileID);
