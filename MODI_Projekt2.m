close all;
clear all;


% Read data.
[u_stat, y_stat] = readData('danestat4/danestat4.txt');

figure(1);
h = scatter(u_stat,y_stat);
h.Marker='.';
title('Data of static object');
xlabel('u_{stat}');
ylabel('y_{stat}');

% Divide data for learning and validating.
divide_coef = 0.3;
divide_index = floor(divide_coef*length(u_stat));
u_stat_learning = u_stat(1:divide_index);
u_stat_valid = u_stat(divide_index+1:end);

y_stat_learning = y_stat(1:divide_index);
y_stat_valid = y_stat(divide_index+1:end);


% Wizualize data.
figure(30);
h = scatter(u_stat_valid,y_stat_valid);
hold on;
h.Marker='.';
h = scatter(u_stat_learning,y_stat_learning);
h.Marker='.';
hold off;
title('Static data after divide.');
xlabel('u_{stat}');
ylabel('y_{stat}');
legend('Validate data','Learning data')


%% Prepare mean square linear model.
%Construct M matrix.
N = 6;
learning_err_vector = zeros(N,1);
valid_err_vector = zeros(N,1);
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
    learning_err_vector(j) = e_learning;
    
    %Count error for validating data.
    y_vector_valid = zeros(1,length(y_stat_valid));
    for i=1:j+1
        y_vector_valid = y_vector_valid + W(i)*u_stat_valid.^(i-1);
    end
    e_valid = immse(y_stat_valid, y_vector_valid);   
    valid_err_vector(j) = e_valid;
    
    % Plot for
    temp = [u_stat_learning; y_vector_learning];
    temp = sort(temp');
    figure(2*j);
    hold on;
    plot(temp(:,1),temp(:,2));
    h = scatter(u_stat_learning,y_stat_learning);
    h.Marker='.';
    hold off;
    legend('y_{mod}','y_{data}');
    title(sprintf('N=%i Error for learning data: %3.6f', j, e_learning));
    xlabel('u');
    ylabel('y');
    
    temp = [u_stat_valid; y_vector_valid];
    temp = sort(temp');
    figure(2*j+1);
    hold on;
    plot(temp(:,1),temp(:,2));
    h = scatter(u_stat_valid,y_stat_valid);
    h.Marker='.';
    hold off;
    legend('y_{mod}','y_{data}');
    title(sprintf('N=%i Error for validating data: %3.6f', j, e_valid));
    xlabel('u');
    ylabel('y');
    
end

%Plot error functions
figure(100);
plot(1:N,learning_err_vector);
hold on;
plot(1:N,valid_err_vector);
hold off;
xlabel('N');
ylabel('error');
title('Static models error');
legend('Learning set error','Validating set error');
%% Dynamic models identification.
na_vector = [1 2 3 4 5 6 7];
nb_vector = [1 2 3 4 5 6 7];
model_degrees = [1 2 3 4 5 6 7];

% Read data.
[u_val, y_val] = readData('danedynwer4/danedynwer4.txt');
[u_learning, y_learning] = readData('danedynucz4/danedynucz4.txt');

arx_learning_error_vector = zeros(length(model_degrees),length(na_vector));
oe_learning_error_vector = zeros(length(model_degrees),length(na_vector));

arx_valid_error_vector = zeros(length(model_degrees),length(na_vector));
oe_valid_error_vector = zeros(length(model_degrees),length(na_vector));


% Wizualize data.
figure(40);
plot(1:length(u_val),u_val);
hold on;
plot(1:length(y_val),y_val);
hold off;
title('Dynamic model - validating data');
xlabel('k');
ylabel('u,y');
legend('u','y')

figure(41);
plot(1:length(u_learning) ,u_learning);
hold on;
plot(1:length(y_learning),y_learning);
hold off;
title('Dynamic model - learning data');
xlabel('k');
ylabel('u,y');
legend('u','y')

for degree=model_degrees
    fprintf('Models of %i\n', degree)
    for na=na_vector
        for nb=nb_vector
            fprintf('na=%i, nb=%i\n',na, nb);
            fprintf('Parameter number=%i',na*nb*degree);
            [W, a,b, Error] = getDynamicModel(na, nb, degree, u_learning, y_learning,'no') %yes, no
            
            %Validating for validatiing set.
            figure;
            hold on;
            plot(1:length(u_val),u_val)
            plot(1:length(y_val),y_val)
            title(sprintf('Symulacja modelu [%i,%i,%i] dla danych weryfikujących', na, nb, degree));
            xlabel('t');
            ylabel('y,u');
            for mode=["ARX","OE"]
                %Verify model.
                y_vector_poly = zeros(1,length(u_val));
                y_vector_poly(1:nb) = y_val(1:nb);
                for i=max(na,nb)+1:length(u_val)
                    u = 1;
                    for j=1:degree
                        for k=1:nb
                           u = [u u_val(i-k).^j]; 
                        end
                    end
                    for j=1:degree
                        for k=1:na
                            switch(mode)
                                case 'ARX'
                                    u = [u y_val(i-k).^j]; 
                                case 'OE'
                                    u = [u y_vector_poly(i-k).^j]; 
                            end
                        end
                    end      
                    y_vector_poly(1,i) = W'*u';
                end
                plot(1:length(y_vector_poly),y_vector_poly);
                if na == nb
                    switch(mode)
                        case 'ARX'
                            arx_valid_error_vector(degree,na) = immse(y_vector_poly, y_learning);
                        case 'OE'
                            oe_valid_error_vector(degree,na) = immse(y_vector_poly, y_learning);
                    end                   
                end
            end
            legend('u','y_{val}','y_{ARX}','y_{OE}');
            hold off;
            
            %Validating for learning set.
            figure;
            hold on;
            plot(1:length(u_learning),u_learning)
            plot(1:length(y_learning),y_learning)
            title(sprintf('Symulacja modelu [%i,%i,%i] dla danych uczących', na, nb, degree));
            xlabel('t');
            ylabel('y,u');
            for mode=["ARX","OE"]
                %Verify model.
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
                            switch(mode)
                                case 'ARX'
                                    u = [u y_learning(i-k).^j]; 
                                case 'OE'
                                    u = [u y_vector_poly(i-k).^j]; 
                            end
                        end
                    end      
                    y_vector_poly(1,i) = W'*u';
                end
                plot(1:length(y_vector_poly),y_vector_poly);
                if na == nb
                    switch(mode)
                        case 'ARX'
                            arx_learning_error_vector(degree,na) = immse(y_vector_poly, y_learning);
                        case 'OE'
                            oe_learning_error_vector(degree,na) = immse(y_vector_poly, y_learning);
                    end                   
                end
            end
            legend('u','y_{learning}','y_{ARX}','y_{OE}');
            hold off;
            
            %Display results for werifying.
% 
%             figure;
%             h = scatter(y_val,y_vector_poly,0.1);
%             h.Marker='.';
%             title('Relacja modelu MNK dla danych uczących');
%             xlabel('y_{val}');
%             ylabel('y_{mod}');
        end
    end
end

%Plot error surfs.
figure(200);
surf(arx_learning_error_vector);
hold on;
surf(arx_valid_error_vector);
hold off;
title('ARX error surface');
xlabel('N');
ylabel('na,nb');
zlabel('error');
legend('Learning data set', 'Validate data set');

figure(201);
surf(oe_learning_error_vector);
hold on;
surf(oe_valid_error_vector);
hold off;
title('OE error surface');
xlabel('N');
ylabel('na,nb');
zlabel('error');
legend('Learning data set', 'Validate data set');

%% Experimental static characteristic degree
Wdyn = [0.000664964992374448;0.00840022541185356;-0.00385019962425829;0.0126918890956658;0.0297501377147536;-0.00908219924044421;-0.0209870455236019;-0.0347898616553102;0.00389270250560253;-0.0147251628540515;0.0379405433374021;0.0113439870870900;0.0264050893831101;0.0129786889612299;0.0399999350273695;0.0234411739204131;0.00691202657199602;0.545394830849041;0.326464311053444;0.133234846180000;-0.0990295278644668;0.191613623406062;-0.0220781684059992;-0.114388769814720;-0.0281154786455814;-0.0267948523934677;0.0150910218935099;-0.0177829463551522;0.0356502576404348;-0.0384919731640682;-0.00508822854579525;0.0342361381018197;-0.00178940157318698];
na_dyn = 4;
nb_dyn = 4;
degree_dyn = 4;

Wstat = [];

y0 = 0;
u0= -1:0.1:1;
y_dyn_vector = zeros(length(u0),1);
y_stat_vector = zeros(length(u0),1);

for i=1:length(u0)
    i
    
    %Found static value form dynamic model.
    static_ch_fun = @(x) x - countStaticFromModel(u0(i), x, W, na_dyn, nb_dyn, degree_dyn);
    y_dyn_vector(i) = fsolve(static_ch_fun,y0);
    
    %Count static value from static model.
    for j=1:length(Wstat)
        y_stat_vector = y_stat_vector + W(j)*u0(i).^(j-1);
    end
end

figure;
plot(u0,y_dyn_vector);
hold on;
plot(u0,y_stat_vector);

title('Static model from dynamic model');
xlabel('u_{stat}');
ylabel('y_{stat}');
legend('from dynamic model', 'from static model')

