function [y_q] = interp1_better(x_vec,y_vec,x_new)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x_q = x_new(:); % force soc to be col-vector

x_vec = x_vec(:); % force to be col vector... 03/24/10
[x_vec,IX] = sort(x_vec);
y_vec = y_vec(:); % force to be col vector... 03/24/10
y_vec=y_vec(IX);


y_q=zeros(size(x_q));
I1=find(x_q <= x_vec(1)); %if ~isempty(I1), disp('low x_vec'); end
I2=find(x_q >= x_vec(end)); %if ~isempty(I2), disp('high x_vec'); end
I3=find(x_q > x_vec(1) & x_q < x_vec(end));
I6=isnan(x_q);

% for voltages less than 0% x_vec... 07/26/06
% extrapolate off low end of table (for x_vec(1) < 0... 03/23/10)
if ~isempty(I1),
  diffx_vec=x_vec(2)-x_vec(1);
  dv = (y_vec(2)) - (y_vec(1));
  y_q(I1)= (x_q(I1)-x_vec(1)).*dv/diffx_vec + y_vec(1);
end

% for voltages greater than 100% x_vec... 07/26/06
% extrapolate off high end of table (for x_vec(end) > 1... 03/23/10)
if ~isempty(I2),
  diffx_vec=x_vec(end)-x_vec(end-1);
  dv = (y_vec(end)) - (y_vec(end-1));
  y_q(I2) = (x_q(I2)-x_vec(end)).*dv/diffx_vec + y_vec(end);
end

% for normal x_vec range...
% manually interpolate (10x faster than "interp1")

for ii=1:length(I3)
    if sum(x_vec==x_q(I3(ii)))>0
        y_value=y_vec(x_vec==x_q(I3(ii)));
        y_q(I3(ii))=y_value(1);
    else
    x_high=x_vec(x_vec>x_q(I3(ii)));
    x_high=x_high(1);
    y_high=y_vec(x_vec>x_q(I3(ii)));
    y_high=y_high(1);
    x_low=x_vec(x_vec<x_q(I3(ii)));
    x_low=x_low(end);
    y_low=y_vec(x_vec<x_q(I3(ii)));
    y_low=y_low(end);
    dy_dx=(y_high-y_low)/(x_high-x_low);
    y_b=y_high-dy_dx*x_high;
    y_q(I3(ii))=y_b+x_q(I3(ii))*dy_dx;
    end
end
% diffx_vec=x_vec(2)-x_vec(1);
% I4=(x_q(I3)-x_vec(1))/diffx_vec; % for x_vec(1) < 0... 03/23/10
% I5=floor(I4); I45 = I4-I5; omI45 = 1-I45;
% y_q(I3)=y_vec(I5+1).*omI45 + y_vec(I5+2).*I45;
y_q(I6)=0; % replace NaN x_vecs with zero voltage... 03/23/10
y_q = reshape(y_q,size(x_new));

end

