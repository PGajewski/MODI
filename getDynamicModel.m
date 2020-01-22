function [W, a,b, Error] = getDynamicModel(na, nb, degree, u_learning, y_learning, print_mode)
%UNTITLED Summary of this function goes here

%Construct M matrix.
M = zeros(length(y_learning)-max(nb,na),nb*degree+na*degree+1);
for i=(max(nb,na)+1):length(y_learning)
    temp = 1;
    for j=1:degree
        for k=1:nb
           temp = [temp u_learning(i-k).^j]; 
        end
    end
    for j=1:degree
        for k=1:na
           temp = [temp y_learning(i-k).^j]; 
        end
    end    
    
    M(i-nb,:) = temp;
end

W =M\y_learning(nb+1:end)';

y_vector_poly = zeros(1,length(u_learning));
y_vector_poly(1:nb) = y_learning(1:nb);
for i=max(na,nb)+1:length(u_learning)
    u = 1;
    for j=1:degree
        for k=1:nb
           u = [u u_learning(i-k).^j]; 
        end
    end
    for j=1:degree
        for k=1:na
           u = [u y_learning(i-k).^j]; 
        end
    end      
    y_vector_poly(1,i) = W'*u';
end

%Compare model for learning data.
if strcmp(print_mode,'yes')
    figure;
    hold on;
    plot(1:length(y_vector_poly),y_vector_poly);
    plot(1:length(y_learning),y_learning);
    plot(1:length(u_learning),u_learning);
    hold off;
    title('Symulacja modelu MNK dla danych ucz¹cych');
    legend('y_{learning}','y_{model}','u');
    xlabel('t');
    ylabel('y,u');

    figure;
    h = scatter(y_learning,y_vector_poly,0.1);
    h.Marker='.';
    title('Relacja modelu MNK dla danych ucz¹cych');
    xlabel('y_{learning}');
    ylabel('y_{mod}');
end

%Separate a and b coefficients.
b = W(2: 1+nb*degree);
a = W(2+nb*degree:end);

%Count MSE.
Error = immse(y_vector_poly,y_learning);