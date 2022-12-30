hold off;
L=1;
x     = linspace(-3*L, 3*L, 300); % Enlarge the interval to 3 periods
Const = -2*L/pi;                  % Constant in the expression for B_n
Sn    = zeros(size(x));           % Initialize vector of series sum values
Cn    = zeros(size(x));
Cn = Cn + 1/2;


for n = 1 : 3
    Const = -Const;       % Efficient way to implement alternating sign
    Bn = Const/n*((-1)^n-1/2*(1+(-1)^n)*(-1)^((n+4)/2));         % Coefficients inversely proportional to n
    Fn = Bn * sin(n*pi*x);
    Sn = Sn + Fn;
    figure(1)
    plot(x, Sn, 'color', [1 0 0] + (3-n)/2.5*[0 1 1], 'linewidth', 4-n);
    hold on;
end

for n = 1 : 5
    Const = -Const;       % Efficient way to implement alternating sign
    An = Const/n*(1/2*(1+(-1)^(n+1))*(-1)^((n+3)/2));
    Fn2= An * cos(n*pi*x);
    if ( Fn2 == 0)
        continue;
    end
    Cn = Cn + Fn2;
    figure(2)
    plot(x, Cn, 'color', [1 0 0] + (5-n)/4.5*[0 1 1], 'linewidth', 6-n);
    hold on;
end

figure(1)
xlabel('x'); ylabel('Sum(B_nsin(n\pix/L))');
title('First three terms in sum of B_nsin(n\pix/L)');
legend('1st harmonic', 'Two terms', 'Three terms', 'Location', 'Southeast');
figure(2)
xlabel('x'); ylabel('Sum(A_ncos(n\pix/L))');
title('First three terms in sum of A_ncos(n\pix/L)');
legend('1st harmonic', 'Two terms', 'Three terms', 'Location', 'Southeast');

% Note how the sum approaches the odd periodic extension of f(x)