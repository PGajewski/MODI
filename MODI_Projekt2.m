close all;
clear all;


% Read data.
[u_stat, y_stat] = readData('danestat4/danestat4.txt');

figure(1);
h = scatter(u_stat,y_stat,0.1);
h.Marker='.';
title('Data for static object');
xlabel('u_{stat}');
ylabel('y_{stat}');

% Divide data for learning and validating.
divide_coef = 0.3;
divide_index = floor(divide_coef*length(u_stat));
u_stat_learning = u_stat(1:divide_index);
u_stat_valid = u_stat(divide_index+1:end);

y_stat_learning = y_stat(1:divide_index);
y_stat_valid = y_stat(divide_index+1:end);


%% Prepare mean square linear model.
%Construct M matrix.
N = 3;
for j=1:N
    M = zeros(length(y_stat_learning),j+1);
    for i=1:j+1
        M(:,i) = u_stat_learning.^(i-1);
    end

    W =M\y_stat_learning';

    %Count error for learning data.
    y_vector_learning = zeros(1,length(y_stat_learning));
    for i=1:j+1
        y_vector_learning = y_vector_learning + W(i)*u_stat_learning.^(i-1);
    end
    e_learning = immse(y_stat_learning, y_vector_learning);
    
    %Count error for validating data.
    y_vector_valid = zeros(1,length(y_stat_valid));
    for i=1:j+1
        y_vector_valid = y_vector_valid + W(i)*u_stat_valid.^(i-1);
    end
    e_valid = immse(y_stat_valid, y_vector_valid);   
    
    % Plot for 
    figure(2*j);
    hold on;
    h = scatter(u_stat_learning,y_vector_learning,0.1);
    h.Marker='.';
    h = scatter(u_stat_learning,y_stat_learning,0.1);
    h.Marker='.';
    hold off;
    legend('y_{mod}','y_{data}');
    title(sprintf('Error for learning data: %3.6f', e_learning));
    xlabel('u');
    ylabel('y');
    
    figure(2*j+1);
    hold on;
    h = scatter(u_stat_valid,y_vector_valid,0.1);
    h.Marker='.';
    h = scatter(u_stat_valid,y_stat_valid,0.1);
    h.Marker='.';
    hold off;
    legend('y_{mod}','y_{data}');
    title(sprintf('Error for validating data: %3.6f', e_valid));
    xlabel('u');
    ylabel('y');
    
end

%% Dynamic models identification.
na_vector = [1 2 3 4];
nb_vector = [1 2 3 4];
model_degrees = [1 2 3 4];

% Read data.
[u_val, y_val] = readData('dane_wer.txt');
[u_learning, y_learning] = readData('dane.txt');

for degree=model_degrees
    pritnf('Models of %i', degree)
    for na=na_vector
        for nb=nb_vector
            printf('na=%i, nb=%i',na, nb);
            getDynamicModel(na, nb, degree x_learning, y_learning)
        end
    end
end