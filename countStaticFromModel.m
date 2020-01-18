function [y] = countStaticFromModel(u0, y0, W,na, nb, degree)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    %Prepare u vector.
    u = 1;
    for j=1:degree
        for k=1:nb
           u = [u u0.^j]; 
        end
    end
    for j=1:degree
        for k=1:na
           u = [u y0.^j]; 
        end
    end
    y = W'*u';
end

