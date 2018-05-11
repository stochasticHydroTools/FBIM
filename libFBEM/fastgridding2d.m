function [a,ix,iy] = fastgridding2d(x_sc,alpha,P,M,h)

% a = zeros(M,M);
a = zeros(P,P);

E0 = zeros(1,P);
Ex = E0; Ey = E0;

j0 = floor([x_sc(1)/h x_sc(2)/h]);

tx_m = floor(x_sc(1)/h)*h - x_sc(1);
tx_p = -tx_m - h;
ty_m = floor(x_sc(2)/h)*h - x_sc(2);
ty_p = -ty_m - h;

e0    = exp(-alpha*h^2);

e1x_m = exp(2*alpha*h*tx_m);
e1x_p = exp(2*alpha*h*tx_p);
e1y_m = exp(2*alpha*h*ty_m);
e1y_p = exp(2*alpha*h*ty_p);

e2x_m = exp(-alpha*tx_m^2);
e2x_p = exp(-alpha*tx_p^2);
e2y_m = exp(-alpha*ty_m^2);
e2y_p = exp(-alpha*ty_p^2);

lstart = (P-1)/2+1;
for l = -(P-1)/2 : 0
   E0(l+lstart) = e0^(l^2);
   Ex(l+lstart) = (e1x_m)^(abs(l)) * e2x_m;
   Ey(l+lstart) = (e1y_m)^(abs(l)) * e2y_m;
end

for l = 1 : (P-1)/2
   E0(l+lstart) = e0^((l-1)^2);
   Ex(l+lstart) = (e1x_p)^(l-1) * e2x_p;
   Ey(l+lstart) = (e1y_p)^(l-1) * e2y_p;
end


% istart = (P-1)/2+1;
% for i2 = -(P-1)/2 : (P-1)/2
%         iy = mod(j0(2) + i2,M)+1;
%         Vy = E0(i2+istart) * Ey(i2+istart);
%         for i1 = -(P-1)/2 : (P-1)/2
%             ix = mod(j0(1)+i1,M) + 1;
%             a(ix,iy) = Vy* E0(i1+istart) * Ex(i1+istart);
%         end
% end

for i2 = 1 : P
        Vy = E0(i2) * Ey(i2);
        for i1 = 1:P
            a(i1,i2) = Vy* E0(i1) * Ex(i1);
        end
end

idx = -(P-1)/2 : (P-1)/2;
ix = mod(idx + j0(1), M) + 1;
iy = mod(idx + j0(2), M) + 1;

