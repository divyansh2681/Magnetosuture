syms theta I_1 I_2 I_3 I_4 y v_d omega_d x y xd yd kr rdd vf lam vdc wdc;


% Inverse Kinematics
x1=0.01;        % Feedback from vision sensor
y1=0.01;      
theta11 = pi/4;

vd = 0.00001;   % Input desired velocity
wd = 0;

rk_1 = [0.08; 0];   % Vectors to the centers of the coils
rk_2 = [0; 0.08];
rk_3 = [-0.08; 0];
rk_4 = [0; -0.08];

rc_1 = [1; 0];  % Unit vectors to the centre of the coils
rc_2 = [0; 1];
rc_3 = [-1; 0];
rc_4 = [0; -1];

rc_1t = transpose(rc_1); % Done to match the dimension to perform matrix/vector multiplications
rc_2t = transpose(rc_2);
rc_3t = transpose(rc_3);
rc_4t = transpose(rc_4);

sigma_r = 2 * pi * 0.43 * 18 * 10^(-3) * log(18/0.8);                       % Constants
sigma_t = 2 *pi * 0.43 * 0.8 * 10^(-3) * (18 * 10^(-3))^2;
c = pi * 0.8*10^(-3) * 0.8*10^(-3) * 18 * 10^(-3) * 1.33 * 0.0915 * 12 * 54;

r = [x1; y1];   % Position vector from the centre of the petri dish to the centre of the needle
h = [cos(theta11); sin(theta11)];   % Heading Vector
ht = transpose(h);

dk_1 = r - rk_1;    % Vector from the centre of the coils to the centre of the needle
dk_2 = r - rk_2;
dk_3 = r - rk_3;
dk_4 = r - rk_4;

mod_1 = norm(dk_1);
mod_2 = norm(dk_2);
mod_3 = norm(dk_3);
mod_4 = norm(dk_4);

dc_1 = dk_1/mod_1;  % Unit vector from the centre of the coils to the centre of the needle
dc_2 = dk_2/mod_2;
dc_3 = dk_3/mod_3;
dc_4 = dk_4/mod_4;

dc_1t = transpose(dc_1);
dc_2t = transpose(dc_2);
dc_3t = transpose(dc_3);
dc_4t = transpose(dc_4);

s=[0,-1;1,0];  % 2x2 Skew symmetric matrix for calculating cross product

vv1 = 3*c*(2*ht*dc_1*rc_1t*h + rc_1t*dc_1 - 5*rc_1t*dc_1*dc_1t*h)/(sigma_r*mod_1^4);    % Parametric velocity component due to magnetic feilds each coil
vv2 = 3*c*(2*ht*dc_2*rc_2t*h + rc_2t*dc_2 - 5*rc_2t*dc_2*dc_2t*h)/(sigma_r*mod_2^4);
vv3 = 3*c*(2*ht*dc_3*rc_3t*h + rc_3t*dc_3 - 5*rc_3t*dc_3*dc_3t*h)/(sigma_r*mod_3^4);
vv4 = 3*c*(2*ht*dc_4*rc_4t*h + rc_4t*dc_4 - 5*rc_4t*dc_4*dc_4t*h)/(sigma_r*mod_4^4);

w1 = c*(ht*s*rc_1 - 3*ht*s*dc_1*dc_1t*rc_1)/(sigma_t*mod_1^3);  % Parametric angular velocity component due to magnetic feilds of each coil
w2 = c*(ht*s*rc_2 - 3*ht*s*dc_2*dc_2t*rc_2)/(sigma_t*mod_2^3);
w3 = c*(ht*s*rc_3 - 3*ht*s*dc_3*dc_3t*rc_3)/(sigma_t*mod_3^3);
w4 = c*(ht*s*rc_4 - 3*ht*s*dc_4*dc_4t*rc_4)/(sigma_t*mod_4^3);

gx_r = [vv1 ,vv2 , vv3 ,vv4; w1 ,w2, w3, w4 ];

% disp(size(gx));

gi = pinv(gx_r); % MoorePenrose Inverse for non square matrices

% disp(size(gi));

m = [vd ; wd];  % Desired velocities

I = gi*m;   % Control currents for each coil

% disp(size(I));

% disp(I);
disp('current in coil 1');
display(I(1,1));
disp('current in coil 2');
display(I(2,1));
disp('current in coil 3');
display(I(3,1));
disp('current in coil 4');
display(I(4,1));



% Forward Kinematics

r_k_1 = [0.08; 0; 0];
r_k_2 = [0; 0.08; 0];
r_k_3 = [-0.08; 0; 0];
r_k_4 = [0; -0.08; 0];

sigma_r = 2 * pi * 0.43 * 18 * 10^(-3) * log(18/0.8);
sigma_theta = 2 *pi * 0.43 * 0.8 * 10^(-3) * (18 * 10^(-3))^2;

I_1 = I(1,1);
I_2 = I(2,1);
I_3 = I(3,1);
I_4 = I(4,1);

r = [x; y; 0];
rd = [xd;yd];
heading = [cos(theta); sin(theta); 0]; % 3 cross 1

heading_T = transpose(heading);  % 1 cross 3
d_k_1 = r - r_k_1;
d_k_2 = r - r_k_2;
d_k_3 = r - r_k_3;
d_k_4 = r - r_k_4;

mod1 = norm(d_k_1);
mod2 = norm(d_k_2);
mod3 = norm(d_k_3);
mod4 = norm(d_k_4);

M_needle = pi * 0.8*10^(-3) * 0.8*10^(-3) * 18 * 10^(-3) * 1.33 / (4 * pi * 10^(-7))  * heading;

B_constant = 4 * pi * 10^(-7) * 0.0915 * 12 * 54;

I = [I_1; I_2; I_3; I_4];

B1 = (-B_constant / (4 * mod1^3)) * ((r_k_1/norm(r_k_1)) - 3 * (d_k_1 / norm (d_k_1)) * transpose(d_k_1 / norm (d_k_1)) * (r_k_1/norm(r_k_1)))*I_1; % Magnetic feild due to each coil
B2 = (-B_constant / (4 * mod2^3)) * ((r_k_2/norm(r_k_2)) - 3 * (d_k_2 / norm (d_k_2)) * transpose(d_k_2 / norm (d_k_2)) * (r_k_2/norm(r_k_2)))*I_2;
B3 = (-B_constant / (4 * mod3^3)) * ((r_k_3/norm(r_k_3)) - 3 * (d_k_3 / norm (d_k_3)) * transpose(d_k_3 / norm (d_k_3)) * (r_k_3/norm(r_k_3)))*I_3;
B4 = (-B_constant / (4 * mod4^3)) * ((r_k_4/norm(r_k_4)) - 3 * (d_k_4 / norm (d_k_4)) * transpose(d_k_4 / norm (d_k_4)) * (r_k_4/norm(r_k_4)))*I_4;

F1 = -gradient(dot(transpose(M_needle), B1)); % Forces due to each coil on the needle
F2 = -gradient(dot(transpose(M_needle), B2));
F3 = -gradient(dot(transpose(M_needle), B3));
F4 = -gradient(dot(transpose(M_needle), B4));

tau_1 = cross (M_needle, B1); % Torque due to each coil on the needle
tau_2 = cross (M_needle, B2);
tau_3 = cross (M_needle, B3);
tau_4 = cross (M_needle, B4);

r_dot_1 = F1 / sigma_r; % Magnitude of velocities
r_dot_2 = F2 / sigma_r;
r_dot_3 = F3 / sigma_r;
r_dot_4 = F4 / sigma_r;

v1 = transpose(heading) * r_dot_1;  % Components of velocity due to each coil
v2 = transpose(heading) * r_dot_2;
v3 = transpose(heading) * r_dot_3;
v4 = transpose(heading) * r_dot_4;

omega1 = tau_1 / sigma_theta;   % Angular velocity due to each coil
omega2 = tau_2 / sigma_theta;
omega3 = tau_3 / sigma_theta;
omega4 = tau_4 / sigma_theta;

g1 = [vpa(norm(v1)) ; vpa(norm(omega1))];
g2 = [vpa(norm(v2)) ; vpa(norm(omega2))];
g3 = [vpa(norm(v3)) ; vpa(norm(omega3))];
g4 = [vpa(norm(v4)) ; vpa(norm(omega4))];
gg = g1 + g2 + g3 + g4;                     % Resultant velocities
disp(gg);
theta = pi/4;
x = 0.05;
y = 0.05;
r=[x;y];
vi = [subs(gg);0];
disp(vi);

zd = kr*(rd - r) + rdd;                                             % PD control equation
t = [cos(theta),sin(theta); -sin(theta)/lam , cos(theta)/lam ] * zd;
vdc = t(1,1);
wdc = t(2,1);


