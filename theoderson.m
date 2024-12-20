
close all;
clc;

usydcolours();

w = 2*pi;
amp = deg2rad(5);
c = 1;
b = 0.5*c;


v = 1;

k = w*b/v;

p = -1*c;

t = linspace(0,2,1000);

cl = 2 * pi * (a(amp,w,t) + h_dot(amp,w,t) + a_dot(amp,w,t)*(0.5 - p))*c_k(k) + ...
         pi * (h_ddot(amp,w,t) + a_dot(amp,w,t) - p*a_ddot(amp,w,t));

% cl1 = pi*b^2*(v*a_dot(amp,w,t) - p*b*a_ddot(amp,w,t)) + 2*pi*b*v*c_k(k)*(v*a(amp,w,t) + b*(0.5-p)*a_dot(amp,w,t));

aoa_t = amp*sin(2*pi*t);
qs_cl = aoa_t * 2*pi;



figure;
plot(t,cl,'LineWidth',2)
hold on;
plot(deg3.e1(4:end),deg3.e00(4:end),'LineWidth',2)
% plot(t,50*a(amp,w,t),'LineWidth',2)
xlim([0,2])
ylim([-20,20])
grid on
grid minor
fontname(gcf,"Times New Roman")
fontsize(gcf,12,'points')
xlabel('Time (s)','FontSize',12)
ylabel('C_L','FontSize',12)
legend('Theodorsen', 'LAUTAT')


figure;
plot(a(amp,w,t), cl)

figure;
plot(t,pi * (h_ddot(amp,w,t) + a_dot(amp,w,t) - p*a_ddot(amp,w,t)))

function[c] = c_k(k)
    c = besselh(1,2,k)./(besselh(1,2,k)+1j.*besselh(0,2,k));
    disp(c)
end

function[x] = a(amp,w, t)
    
    x = amp*sin(w*t);

end

function[x] = a_dot(amp,w, t)
    
    x = w*amp*cos(w*t);

end

function[x] = a_ddot(amp,w, t)
    
    x = -w*w*amp*sin(w*t);

end

function[x] = h(amp,w, t)
    
    x = 0.5*sin(amp) * sin(w*t);

end

function[x] = h_dot(amp,w, t)
    
    x = w*0.5*sin(amp) * cos(w*t);

end

function[x] = h_ddot(amp,w, t)
    
    x = -w*w*0.5*sin(amp) * sin(w*t);

end