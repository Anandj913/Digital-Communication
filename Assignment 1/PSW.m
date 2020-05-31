function [psw, lambda] = PSW(c, N ,to , factor, order)
matrix = zeros(N,N);
%defining matrix for eigen value calculation as given in paper
for i=0:N-1
    matrix(i+1,i+1)=i*(i+1)+(2*i*(i+1)-1)/((2*i-1)*(2*i+3))*c*c;
    if i+3<=N
        matrix(i+1,i+3)=(i+2)*(i+1)/((2*i+3)*(2*i+5))*c*c;
        matrix(i+3,i+1)=(i+2)*(i+1)/((2*i+3)*(2*i+1))*c*c;
    end
end
[eigen_vector, eigen_value] = eig(matrix);
[~, sorted_inx] = sort(diag(eigen_value));
eigen_vector_sorted = eigen_vector(:,sorted_inx);

%calculating normalization factor
deno = 0;
if mod(order,2)==0
    num=(-1)^(order/2)*factorial(order)/(2^order*(factorial(order/2))^2);
    for i=0:2:N-1
        deno=deno+(-1)^(i/2)*factorial(i)*eigen_vector_sorted(i+1,order+1)/(2^i*factorial(i/2)*factorial(i/2));
    end          
else
    num=(-1)^((order-1)/2)*factorial(order+1)/(2^order*(factorial((order-1)/2))*(factorial((order+1)/2)));
    for i=1:2:N-1
        deno=deno+(-1)^((i-1)/2)*factorial(i+1)*eigen_vector_sorted(i+1,order+1)/(2^i*factorial((i-1)/2)*factorial((i+1)/2));
    end           
end   

norm = num/deno;
eigen_vector_sorted(:,order+1) = norm*eigen_vector_sorted(:,order+1);

if mod(order,2)==0
    lambda=2*c/pi*(2^order*eigen_vector_sorted(1,order+1)*(factorial(order/2))^2/(factorial(order)))^2;
    k=factorial(order)/(2^order*(factorial(order/2))^2*eigen_vector_sorted(1,order+1));
else
    lambda=2*c/pi*(2^order*eigen_vector_sorted(2,order+1)*c*(factorial((order-1)/2))*(factorial((order+1)/2))/(3*factorial(order+1)))^2;
    k=3*factorial((order+1))/(2^order*(factorial((order-1)/2))*(factorial((order+1)/2))*c*eigen_vector_sorted(2,order+1));
end    
N0=0;
if mod(order,2)==0
    for i=0:2:N-1
        N0=N0+2*(1/(2*i+1))*(eigen_vector_sorted(i+1,order+1))^2;
    end    
else
    for i=1:2:N-1
        N0=N0+2*(1/(2*i+1))*(eigen_vector_sorted(i+1,order+1))^2;
    end        
end

t=-factor*to:1/100:factor*to;
SIZE=int32(size(t));

psw=zeros(1,SIZE(2));
if mod(order,2)==0
    for i=0:2:N-1
        psw=psw+(lambda/(to*N0))^(0.5)*k*(-1)^((i-order)/2)*eigen_vector_sorted(i+1,order+1)*sphbes(i,c/to*t);
    end
else
    for i=1:2:N-1
        psw=psw+(lambda/(to*N0))^(0.5)*k*(-1)^((i-order)/2)*eigen_vector_sorted(i+1,order+1)*sphbes(i,c/to*t);
    end        
end 
psw=real(psw);
for i=1:(SIZE(2)+1)/2
    if i == (SIZE(2)+1)/2
        if mod(order,2)==0
            psw(i)=psw(i-1);
        else
            psw(i)=0;
        end
    else
        psw(i)=-psw(i);
    end
end    

function result = sphbes(n, x)
[lx ln] = size(n);
xm = repmat(x, 1, ln);
result = sqrt(pi ./(2* xm)) .* besselj(n + 0.5, x);
end
end