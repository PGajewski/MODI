function [W, a,b, Error] = getDynamicModel(na, nb, degree x_learning, y_learning)
%UNTITLED Summary of this function goes here

%Construct M matrix.
M = zeros(length(y_learning)-nb,nb-tau+na+2);
for i=(nb+1):length(y_learning)
    temp = 1;
    for j=1:degree
        for k=1:nb
           temp = [temp u_learning(i-k)^j]; 
        end
    end
    for j=1:degree
        for k=1:na
           temp = [temp y_learning(i-k)^j]; 
        end
    end    
    
    M(i-nb,:) = temp;
    %M(i-nb,:) = [1 flip(u_learning(i-nb:i-tau)) flip(y_learning(i-na:i-1))];
end

W =M\y_learning(nb+1:end)';

y_vector_poly = zeros(1,length(u_val));
y_vector_poly(1:nb) = y_learning(1:nb);
for i=(nb)+1:length(u_val)
    y_vector_poly(1,i) = W'*[1 flip(u_learning(i-nb:i-tau)) flip(y_vector_poly(i-na:i-1))]';
end

%Compare model for learning data.
figure;
hold on;
plot(1:length(y_vector_poly),y_vector_poly);
plot(1:length(y_learning),y_learning);
plot(1:length(u_learning),u_learning);
hold off;
title('Symulacja modelu MNK dla danych ucz�cych');
legend('y_{learning}','y_{model}','u');
xlabel('t');
ylabel('y,u');

figure;
h = scatter(y_learning,y_vector_poly,0.1);
h.Marker='.';
title('Relacja modelu MNK dla danych ucz�cych');
xlabel('y_{learning}');
ylabel('y_{mod}');


y_vector_poly = zeros(1,length(u_val));
y_vector_poly(1:nb) = y_val(1:nb);
for i=(nb)+1:length(u_val)
    y_vector_poly(1,i) = W'*[1 flip(u_val(i-nb:i-tau)) flip(y_vector_poly(i-na:i-1))]';
end

%Separate a and b coefficients.
b = [2: 2+nb];
a = [3+nb:end];

%Count MSE.
Error = immse(y_vector_poly,y_learning);