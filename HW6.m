%% MAE472 HW6
%Sam Masten


%% Problem 1

%Given Values
E1 = 145e9; %Pa
E2 = 10.5e9; %Pa
G12 = 7.0e9; %Pa
v12 = 0.28;
v21 = v12 * (E2/E1);

%Inputted Values
n = input('\nNumber of layers:');

t = .001; %m, thickness of one ply
h = .008; %m, height
theta = [-45 20 45 -45 20 45];

%Code to Solve

%% Calculating the C Matrix

C11 = E1/(1-(v12*v21));
C12 = (v12*E2)/(1-(v12*v21));
C16= 0;
C22 = E2/(1-(v12*v21));
C26 = 0;
C66 = G12;

C = [C11 C12 C16; C12 C22 C26; C16 C26 C66];

fprintf('\nThe C matrix in GPa is\n');
fprintf('%14.4f %14.4f %14.4f\n', C/10^9)


%% Calculating the C Bar Matrix

Cbar = {};

for i = 1:length(theta)
    c = cosd(theta(i));
    s = sind(theta(i));
    
    C_11(i) = C11*c^4+C22*s^4+(2*C12+4*C66)*c^2*s^2;
    C_12(i) = C12*(c^4+s^4)+(C11+C22-4*C66)*c^2*s^2;
    C_22(i) = C11*s^4+C22*c^4+(2*C12+4*C66)*c^2*s^2;
    C_16(i) = (C11-C12-2*C66)*c^3*s-(C22-C12-2*C66)*c*s^3;
    C_26(i) = (C11-C12-2*C66)*c*s^3-(C22-C12-2*C66)*c^3*s;
    C_66(i) = (C11+C22-2*C12-2*C66)*c^2*s^2+C66*(c^4+s^4);

    Cbar{i} = [C_11(i) C_12(i) C_16(i); C_12(i) C_22(i) C_26(i); C_16(i) C_26(i) C_66(i)];

    fprintf('\nC Bar Matrix (%.1f degrees) in GPa\n',theta(i))
    fprintf('%14.4f %14.4f %14.4f\n', Cbar{i}/(10^9))
   
end



%% Calculating the z values

%Accounting for varying thickness
for i = 1:n+1
    z(i) = -h/2 + (t*(i-1));
end

for i = 3:n+1
    z(i) = z(i) +.001;
end

for i = 6:n+1
    z(i) = z(i) +.001;
end

%% Creating each Z vector

z1 = [];
z2 = [];
z3 = [];

for i = 1:n
    z1_loop = z(i+1)-z(i);
    z2_loop = z(i+1)^2 - z(i)^2;
    z3_loop = z(i+1)^3 - z(i)^3;
    z1 = [z1 z1_loop];
    z2 = [z2 z2_loop];
    z3 = [z3 z3_loop];
end

%% Calculating the A,B,D Matrices

for i = 1:n
    coeff = 1/3;
    A11(i) = C_11(i).*z1(i);
    A12(i) = C_12(i).*z1(i);
    A16(i) = C_16(i).*z1(i);
    A22(i) = C_22(i).*z1(i);
    A26(i) = C_26(i).*z1(i);
    A66(i) = C_66(i).*z1(i);
    B11(i) = 0.5*C_11(i)*z2(i);
    B12(i) = 0.5*C_12(i)*z2(i);
    B16(i) = 0.5*C_16(i)*z2(i);
    B22(i) = 0.5*C_22(i)*z2(i);
    B26(i) = 0.5*C_26(i)*z2(i);
    B66(i) = 0.5*C_66(i)*z2(i);
    D11(i) = coeff*C_11(i)*z3(i);
    D12(i) = coeff*C_12(i)*z3(i);
    D16(i) = coeff*C_16(i)*z3(i);
    D22(i) = coeff*C_22(i)*z3(i);
    D26(i) = coeff*C_26(i)*z3(i);
    D66(i) = coeff*C_66(i)*z3(i);
end

A = [sum(A11) sum(A12) sum(A16); sum(A12) sum(A22) sum(A26); sum(A16) sum(A26) sum(A66)];
B = [sum(B11) sum(B12) sum(B16); sum(B12) sum(B22) sum(B26); sum(B16) sum(B26) sum(B66)];
D = [sum(D11) sum(D12) sum(D16); sum(D12) sum(D22) sum(D26); sum(D16) sum(D26) sum(D66)];

%% Outputs

fprintf('A Matrix in N/m\n')
A

fprintf('B Matrix in N\n')
B

fprintf('D Matrix in N-m\n')
D


%% Problem 2

%Given Values
E1 = 145e9; %Pa
E2 = 10.5e9; %Pa
G12 = 7.0e9; %Pa
v12 = 0.28;
v21 = v12 * (E2/E1);

%Inputted Values
n = input('\nNumber of layers:');

t = .00025; %m, thickness of one ply
h = t*3; %m, height
theta = [0 45 -45];

%Code to Solve

%% Calculating the C Matrix

C11 = E1/(1-(v12*v21));
C12 = (v12*E2)/(1-(v12*v21));
C16= 0;
C22 = E2/(1-(v12*v21));
C26 = 0;
C66 = G12;

C = [C11 C12 C16; C12 C22 C26; C16 C26 C66];

fprintf('\nThe C matrix in GPa is\n');
fprintf('%14.4f %14.4f %14.4f\n', C/10^9)


%% Calculating the C Bar Matrix

Cbar = {};

for i = 1:length(theta)
    c = cosd(theta(i));
    s = sind(theta(i));
    
    C_11(i) = C11*c^4+C22*s^4+(2*C12+4*C66)*c^2*s^2;
    C_12(i) = C12*(c^4+s^4)+(C11+C22-4*C66)*c^2*s^2;
    C_22(i) = C11*s^4+C22*c^4+(2*C12+4*C66)*c^2*s^2;
    C_16(i) = (C11-C12-2*C66)*c^3*s-(C22-C12-2*C66)*c*s^3;
    C_26(i) = (C11-C12-2*C66)*c*s^3-(C22-C12-2*C66)*c^3*s;
    C_66(i) = (C11+C22-2*C12-2*C66)*c^2*s^2+C66*(c^4+s^4);

    Cbar{i} = [C_11(i) C_12(i) C_16(i); C_12(i) C_22(i) C_26(i); C_16(i) C_26(i) C_66(i)];

    fprintf('\nC Bar Matrix (%.1f degrees) in GPa\n',theta(i))
    fprintf('%14.4f %14.4f %14.4f\n', Cbar{i}/(10^9))
   
end



%% Calculating the z values

%Accounting for varying thickness
for i = 1:n+1
    z(i) = -h/2 + (t*(i-1));
end

%% Creating each Z vector

z1 = [];
z2 = [];
z3 = [];

for i = 1:n
    z1_loop = z(i+1)-z(i);
    z2_loop = z(i+1)^2 - z(i)^2;
    z3_loop = z(i+1)^3 - z(i)^3;
    z1 = [z1 z1_loop];
    z2 = [z2 z2_loop];
    z3 = [z3 z3_loop];
end

%% Calculating the A,B,D Matrices

for i = 1:n
    coeff = 1/3;
    A11(i) = C_11(i).*z1(i);
    A12(i) = C_12(i).*z1(i);
    A16(i) = C_16(i).*z1(i);
    A22(i) = C_22(i).*z1(i);
    A26(i) = C_26(i).*z1(i);
    A66(i) = C_66(i).*z1(i);
    B11(i) = 0.5*C_11(i)*z2(i);
    B12(i) = 0.5*C_12(i)*z2(i);
    B16(i) = 0.5*C_16(i)*z2(i);
    B22(i) = 0.5*C_22(i)*z2(i);
    B26(i) = 0.5*C_26(i)*z2(i);
    B66(i) = 0.5*C_66(i)*z2(i);
    D11(i) = coeff*C_11(i)*z3(i);
    D12(i) = coeff*C_12(i)*z3(i);
    D16(i) = coeff*C_16(i)*z3(i);
    D22(i) = coeff*C_22(i)*z3(i);
    D26(i) = coeff*C_26(i)*z3(i);
    D66(i) = coeff*C_66(i)*z3(i);
end

A = [sum(A11) sum(A12) sum(A16); sum(A12) sum(A22) sum(A26); sum(A16) sum(A26) sum(A66)];
B = [sum(B11) sum(B12) sum(B16); sum(B12) sum(B22) sum(B26); sum(B16) sum(B26) sum(B66)];
D = [sum(D11) sum(D12) sum(D16); sum(D12) sum(D22) sum(D26); sum(D16) sum(D26) sum(D66)];

%% Outputs

fprintf('A Matrix in N/m\n')
A

fprintf('B Matrix in N\n')
B

fprintf('D Matrix in N-m\n')
D

%% Applied Loads

Nx = 5400; %kN/m
Ny = -3200; %kN/m
Ns = -1900; %kN/m
Mx = 4500; %kN
My = 2900; %kN
Ms = -2300; %kN

%Load Vector/Matrix
NM_mat = [Nx; Ny; Ns; Mx; My; Ms];
NM_mat = NM_mat * (10^3); %converting kN to N

%% Finding Strain Values

%Solving for Midplane Strain & Curvatures
F = [A B; B D]^(-1);
%changed N to kN

midplane_mat = F*NM_mat;

%Extracting each variable
epsilon_x0 = midplane_mat(1);
epsilon_y0 = midplane_mat(2);
epsilon_s0 = midplane_mat(3);
k_x = midplane_mat(4);
k_y = midplane_mat(5);
k_s = midplane_mat(6);

%Forming Vectors
strain0_mat = [epsilon_x0; epsilon_y0; epsilon_s0];
k_mat = [k_x; k_y; k_s]; %convert to /mm

%% Calculating Strains at Other Distances from Midplane

epsilon_x = strain0_mat(1) + z*k_mat(1);
epsilon_y = strain0_mat(2) + z*k_mat(2);
gamma_s = strain0_mat(3) + z*k_mat(3);


%% Calculating Stresses at Other Distances from Midplane

z_new = [z(1) z(2) z(2) z(3) z(3) z(4)];

for j = 1:n
    
    spot = 2*j-1;
    spot2 = 2*j;
    
    c11 = C_11(j)*(10^-9);
    c12 = C_12(j)*(10^-9);
    c16 = C_16(j)*(10^-9);
    c22 = C_22(j)*(10^-9);
    c26 = C_26(j)*(10^-9);
    c66 = C_66(j)*(10^-9);
    
    ex = epsilon_x(j);
    ey = epsilon_y(j);
    gs = gamma_s(j);
    
    ex2 = epsilon_x(j+1);
    ey2 = epsilon_y(j+1);
    gs2 = gamma_s(j+1);
    
    sigma_x(spot) = [c11*ex + c12*ey + c16*gs];
    sigma_y(spot) = [c12*ex + c22*ey + c26*gs];
    tau_s(spot) = [c16*ex + c26*ey + c66*gs];

    sigma_x(spot2) = [c11*ex2 + c12*ey2 + c16*gs2];
    sigma_y(spot2) = [c12*ex2 + c22*ey2 + c26*gs2];
    tau_s(spot2) = [c16*ex2 + c26*ey2 + c66*gs2];    
    
end

%% Plots

%% Strain Field Plot

figure(1)
plot(z,epsilon_x,'b','LineWidth',1.5)
hold on;
plot(z,epsilon_y,'r','LineWidth',1.5)
hold on;
plot(z,gamma_s,'g','LineWidth',1.5)
title('Laminate Strain Field')
xlabel('z')
ylabel('Strain')
legend('{\epsilon}_{x}','{\epsilon}_{y}','{\gamma}_{s}')


%% Stress Field Plot

figure(2)
plot(z_new,sigma_x,'b','LineWidth',1.5)
hold on;
plot(z_new,sigma_y,'r','LineWidth',1.5)
hold on;
plot(z_new,tau_s,'g','LineWidth',1.5)
title('Laminate Stress Field')
xlabel('z')
ylabel('Stress (GPa)')
legend('{\sigma}_{x}','{\sigma}_{y}','{\tau}_{s}')