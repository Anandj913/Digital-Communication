t0 =1;
c = [0.5, 1, 2, 4];
N = 20;
order = 6;
%---------------Finding inftinite time PSW---------------------------
scale = 20;
psw_inftime = [;;];
eigen_value_inftime = [;];
t_inftime =-scale*t0:1/100:scale*t0;
for j=1:length(c)
    for i=1:order
        [psw_inftime(j,i,:), eigen_value_inftime(j,i)] = PSW(c(j), N, t0, scale, i-1);
    end
end
%-----------|Plotting PSW in range of [-20,20] for c = 0.5|-----------------------
%------------------Use for plotting functions-------------------------
% for i=1:6
%  subplot(3,2,i);
%    k = psw_inftime(4,i,:);
%    k = k(1,:);
%    plot(t_inftime/t0,k);
%    xlim([-scale,scale]);
%    xlabel('t/t0');
%    ylabel('Ampltiude');
%    title(sprintf('PSW: Order:%d | c=%0.1f | eig-value=%0.2f/10^6',i-1, c(4), eigen_value_inftime(4,i)*10^6));
% end

eigen_value_inftime = eigen_value_inftime';
%--------------| for c=4 |--------------------------
%-------------Proving orthonormality--------------------
orthonormal_mat_inftime = zeros(order, order);
for j=1:order
    for i=1:order
        %k = dot(psw(4,j,:), psw(4,i,:));
        k = trapz(t_inftime,psw_inftime(4,j,:).*psw_inftime(4,i,:));
        if k < 10^(-1)
            orthonormal_mat_inftime(j,i) = 0;
        else
            orthonormal_mat_inftime(j,i) = k;
        end
    end
end
%---------------Finding ftinite time PSW---------------------------
scale_finite = 1;
psw_finitetime = [;;];
eigen_value_finitetime = [;];
t_finitetime =-scale_finite*t0:1/100:scale_finite*t0;
for j=1:length(c)
    for i=1:order
        [psw_finitetime(j,i,:), eigen_value_finitetime(j,i)] = PSW(c(j), N, t0, scale_finite, i-1);
    end
end
eigen_value_finitetime = eigen_value_finitetime';
%---------------|for c = 4|----------------------------
%-------------Proving orthogonality--------------------
orthogonal_mat = zeros(order, order);
for j=1:order
    for i=1:order
        k = trapz(t_finitetime,psw_finitetime(4,j,:).*psw_finitetime(4,i,:));
        if k < 10^(-7)
            orthogonal_mat(j,i) = 0;
        else
            orthogonal_mat(j,i) = k;
        end   
    end
end
%----------Printing Results-------------------------
disp('Eigen_value for infinite time support for all c');
disp(eigen_value_inftime);
disp('Orthonormal matrix for infinite time support for c = 4');
disp(orthonormal_mat_inftime);
disp('Eigen_value for finite time support for all c');
disp(eigen_value_finitetime);
disp('Orthogonal matrix for infinite time support for c = 4');
disp(orthogonal_mat);

%-----------L2 function approximation using finite PSWs with c =4--------------
coff_finite = zeros(5);
for i=1:5
    k = psw_finitetime(4,i,:);
    k = k(1,:);
    coff_finite(i)= trapz(t_finitetime, sinc(t_finitetime).*k);
end

reconstructed_sinc_finite = zeros(1,length(k));
for i=1:5
    k = psw_finitetime(4,i,:);
    k = k(1,:);
    reconstructed_sinc_finite = reconstructed_sinc_finite + k.*coff_finite(i);
end
%----------Plotting---------------
subplot(2,2,1)
plot(t_finitetime, sinc(t_finitetime));
xlim([-scale_finite,scale_finite]);
xlabel('t/t0');
ylabel('Ampltiude');
title('Original Sinc Function');

subplot(2,2,2)
plot(t_finitetime, reconstructed_sinc_finite );
xlim([-scale_finite,scale_finite]);
xlabel('t/t0');
ylabel('Ampltiude');
title('Reconstructed Sinc Function with finite c=4 PSWs');

error_finite = (rms(sinc(t_finitetime))- rms(reconstructed_sinc_finite))/rms(sinc(t_finitetime));
fprintf('\n Percentage Energy loss due to finite approximation = %.3f', error_finite*100);

%-----------L2 function approximation using infinite PSWs with c =4--------------
coff_infinite = zeros(5);
for i=1:5
    k = psw_inftime(4,i,:);
    k = k(1,:);
    coff_infinite(i)= trapz(t_inftime, sinc(t_inftime).*k);
end

reconstructed_sinc_infinite = zeros(1,length(k));
for i=1:5
    k = psw_inftime(4,i,:);
    k = k(1,:);
    reconstructed_sinc_infinite = reconstructed_sinc_infinite + k.*coff_infinite(i);
end
%----------Plotting---------------
subplot(2,2,3)
plot(t_inftime, sinc(t_inftime));
xlim([-scale,scale]);
xlabel('t/t0');
ylabel('Ampltiude');
title('Original Sinc Function');

subplot(2,2,4)
plot(t_inftime, reconstructed_sinc_infinite );
xlim([-scale,scale]);
xlabel('t/t0');
ylabel('Ampltiude');
title('Reconstructed Sinc Function with infinite c=4 PSWs');

error_infinite = (rms(sinc(t_inftime))- rms(reconstructed_sinc_infinite))/rms(sinc(t_inftime));
fprintf('\n Pecentage Energy loss due to infinite approximation = %.3f', error_infinite*100);

%---------------Complete ----------------------------





