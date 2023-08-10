% P-V DATA ANALYSIS
%Solving the 2nd & 3rd-order Birch-Murnaghan EoS 
%% CODE INFORMATION 
% Input arguments 
%P and error if available
%V and error if available
% Version 1.0 2020
% @E.K.Tripoliti ucfbetr@ucl.ac.uk
%
% INPUT VARIABLES
rng(9845, 'twister')                 % For reproductivity; controls random number generation. 
 P_exp = data(:,1); 
 V_exp = data(:,2);        
error_V = data(:,3);    
V_pred = min(V_exp)-100:1:max(V_exp)+10;     V_pred = V_pred';  
%% SOLVE MODEL FOR THE 2ND-ORDER BMEoS
% initial values for the unknown model parameters: V0, K0
% V0=x(1);                      K0=x(2);                       
x0 = [170 1600];
xl = [100  0];
xu = [1500 1800];
options = optimset( 'algorithm', 'trust-region-reflective', ...
    'MaxIter', 3000, ...
    'maxFunEvals', 3000, ...
    'FinDiffType', 'central',...
    'TolX', 1e-30, ...
    'TolFun',1e-30, ...
    'Display', 'final');
% Create function for 2nd-order BMEoS
fun2BMEoS = @(x) ((((((3/2)*x(1)).*(((x(2)./V_exp).^(7/3))-((x(2)./V_exp).^(5/3)))- P_exp)).^2)./error_V);
[x,resnorm,residual,exitflag,output,lambda, jacobian]=  lsqnonlin(fun2BMEoS,x0, xl,xu, options);
%
% Calculated Pressures for finer volume steps
% Data
P2BMEoS_data = ((3/2)*x(1)).*(((x(2)./V_exp).^(7/3))-((x(2)./V_exp).^(5/3)));                  
% Extrapolated 
P2BMEoS = ((3/2)*x(1)).*(((x(2)./V_pred).^(7/3))-((x(2)./V_pred).^(5/3)));                     % Extrapolated Pressures     
%
SR2 = (P_exp - P2BMEoS_data).^2; SSR2 = sum(SR2);
R2_2 = regress( P_exp, P2BMEoS_data);
%%  SOLVE MODEL FOR THE 3RD ORDER BMEoS
% unknown fitted parameters : K0, K'0, V0
% initial values for the unknowns
% V0=x(1);             K'0=x(2);           KT0=x(3);    
m0 = [150 4 1500];
mlb = [100 1 0];
mub = [1500 6 1600];
options = optimset( 'algorithm', 'trust-region-reflective', ...
    'MaxIter', 3000, ...
    'maxFunEvals', 3000, ...
    'FinDiffType', 'central',...
    'TolX', 1e-30, ...
    'TolFun',1e-30, ...
    'Display', 'final');
% Create function for 3nd-order BMEoS
fun3BMEoS = @(m) ((((((3/2)*m(1)).*(((m(3)./V_exp).^(7/3))-((m(3)./V_exp).^(5/3))).*(1+(3/4).*(m(2)-4).*(((m(3)./V_exp).^(2/3))-1))-P_exp)).^2)./error_V);
% Export fitted patameters  
[m,resnorm,residual,exitflag,output,lambda, jacobian] = lsqnonlin(fun3BMEoS,m0,mlb,mub, options);
%
% Calculated Pressures
% Data
P3BMEoS_data = ((3/2)*m(1)).*(((m(3)./V_exp).^(7/3))-((m(3)./V_exp).^(5/3))).*(1+(3/4).*(m(2)-4).*(((m(3)./V_exp).^(2/3))-1));      
% Extrapolated
P3BMEoS = ((3/2)*m(1)).*(((m(3)./V_pred).^(7/3))-((m(3)./V_pred).^(5/3))).*(1+(3/4).*(m(2)-4).*(((m(3)./V_pred).^(2/3))-1));      
%
SR3 = (P_exp - P3BMEoS_data).^2; SSR3 = sum(SR3);
R2_3 = regress( P_exp, P3BMEoS_data);
%% CALCULATE BULK MODULUS 
% from BMEoS models
Kt_3rd = - V_pred(2:end).*(diff(P3BMEoS))./(diff(V_pred));
Kt_2nd = - V_pred(2:end).*(diff(P2BMEoS))./(diff(V_pred));
% From Experimental Data
BulkM = - V_exp(2:end).*(diff(P_exp))./(diff(V_exp));   
% 
%% CALCULATE  STRESS & STRAIN
% P2BMEoS
stress2 = -diff(V_pred)./(V_pred(2:end));
strain2 = Kt_2nd.*stress2;
f2 = (((x(2)./V_pred).^(2/3)) -1).^(1/2);                  % Eulerian srain
F2 = P2BMEoS./3.*f2.*((2.*f2+1).^(5/2));                   % Normalized stress
% P3BMEoS
stress3 = -diff(V_exp)./(V_exp(2:end));
strain3 = BulkM.*stress3;
f3 =0.5.*(((m(3)./V_exp).^(2/3)) -1);                   % Eulerian finite strain
F3 = P_exp./((3.*f3).*(((2.*f3)+1).^(5/2)));            % Normalized pressure
% extrapolations
f33 =0.5.*(((m(3)./V_pred).^(2/3)) -1);                   % Eulerian finite strain
F33 = P3BMEoS./((3.*f33).*(((2.*f33)+1).^(5/2)));         % Normalized pressure
%K''0 from Anderson (1995)
K_double_prime3 = (-1/m(1)).*((3-m(2)).*(4-m(2))+(35/9)); %for the 3rd-order
K_double_prime2 = (-1/m(1)).*((3-4).*(4-4)+(35/9)); %for the 2nd-order
display([' K''03:'     num2str(K_double_prime3)]);   display([' K''02:'      num2str(K_double_prime2)]); 
%% CALCULATE DEFORMATION ENERGY
%data
DE = m(3).*((9/2).*(f3.^2)).*(m(1)+((2.*((3/2).*(m(1).*(m(2)-4))))./(3.*f3))+((3/2).*(m(1).^2).*K_double_prime2+((3/2).*(m(2).^2).*m(1))-((21/2).*(m(1).*m(2))) + (143/6.*m(1))));
%extrapolations
DE2 = m(3).*((9/2).*(f33.^2)).*(m(1)+((2.*((3/2).*(m(1).*(m(2)-4))))./(3.*f33))+((3/2).*(m(1).^2).*K_double_prime2+((3/2).*(m(2).^2).*m(1))-((21/2).*(m(1).*m(2))) + (143/6.*m(1))));
%% CALCULATE BULK MODULUS
%extrapolations
Kt2 = - V_pred(2:end).*(diff(P2BMEoS))./(diff(V_pred));
Kt3 = - V_pred(2:end).*(diff(P3BMEoS))./(diff(V_pred));
%data
BulkM = - V_exp(2:end).*(diff(P_exp))./(diff(V_exp));   

%% STANDARD ERRORS
% Partial Derivatives and Weights 
% for P2BMEoS
D2_1 = ((3/2)*x(1).*(1+0.001)).*(((x(2)./V_exp).^(7/3))-((x(2)./V_exp).^(5/3)))./(x(1)*0.001); 
D2_2 = ((3/2)*x(1)).*((((x(2).*(1+0.001))./V_exp).^(7/3))-(((x(2).*(1+0.001))./V_exp).^(5/3)))./(x(2)*0.001);
% for P3BMEoS
D3_1 = ((3/2)*(m(1)*(1+0.001))).*(((m(3)./V_exp).^(7/3))-((m(3)./V_exp).^(5/3))).*(1+(3/4).*(m(2)-4).*(((m(3)./V_exp).^(2/3))-1))./(m(1)*0.001);
D3_2 = ((3/2)*m(1)).*(((m(3)./V_exp).^(7/3))-((m(3)./V_exp).^(5/3))).*(1+(3/4).*((m(2).*(1+0.001))-4).*(((m(3)./V_exp).^(2/3))-1))./(m(2)*0.001);
D3_3 = ((3/2)*m(1)).*((((m(3)*(1+0.001))./V_exp).^(7/3))-(((m(3)*(1+0.001))./V_exp).^(5/3))).*(1+(3/4).*(m(2)-4).*((((m(3)*(1+0.001))./V_exp).^(2/3))-1))./(m(3)*0.001);
%
% For weighted error use P_error = experimental, for equally weighted use P_error = 1.
P_error = 1;       % Change this if P_error from experiments is available
% for P2BMEoS
D2_11 = P_error.*(D2_1.^2);          SD2_11 = sum(D2_11,'all');
D2_12 =  P_error.*D2_1.*D2_2;        SD2_12 = sum(D2_12,'all');
%
D2_22 = P_error.*(D2_2.^2);          SD2_22 = sum(D2_22,'all');
%
% aij & bij matrices
aij2 = zeros(2,2);
aij2(1,1) = SD2_11; aij2(1,2) = SD2_12; 
aij2(2,1) = SD2_12; aij2(2,2) = SD2_22;
bij2 = inv(aij2);                                        % Inversion matrix
Cij2 = corrcoef(bij2);                                   % Correlation matrix 
NoU2 = 2;                                                % Number of unknown variables in the model
% e.s.d 
Error_2V0 = sqrt(bij2(1,1).*SSR2./(length(P_exp) - NoU2));
Error_2KT0 = sqrt(bij2(2,2).*SSR2./(length(P_exp) - NoU2));
% for P3BMEoS
D3_11 = P_error.*(D3_1.^2);             SD3_11 = sum(D3_11,'all');
D3_12 =  P_error.*D3_1.*D3_2;           SD3_12 = sum(D3_12,'all');
D3_13 =  P_error.*D3_1.*D3_3;           SD3_13 = sum(D3_13,'all');
%
D3_22 = P_error.*(D3_2.^2);             SD3_22 = sum(D3_22,'all');
D3_23 = P_error.*D3_2.*D3_3;            SD3_23 = sum(D3_23,'all');
%
D3_33 = P_error.*D3_3.*D3_3;            SD3_33 = sum(D3_33,'all');
%
% aij & bij matrices
aij3 = zeros(3,3);
aij3(1,1) = SD3_11; aij3(1,2) = SD3_12; aij3(1,3) = SD3_13;
aij3(2,1) = SD3_12; aij3(2,2) = SD3_22; aij3(2,3) = SD3_23;
aij3(3,1) = SD3_13; aij3(3,2) = SD3_23; aij3(3,3) = SD3_33;
bij3 = inv(aij3);                           % Inversion matrix
Cij3 = corrcoef(bij3);                      % Correlation matrix 
NoU3 = 3;                                   % Number of unknown variables in the model
% e.s.d 
Error_3KT0 = sqrt(bij3(1,1).*SSR3./(length(P_exp) - NoU3));
Error_KT0_prime = sqrt(bij3(2,2).*SSR3./(length(P_exp) - NoU3));
Error_3VT0 = sqrt(bij3(3,3).*SSR3./(length(P_exp) - NoU3));
%%   PRINT ALL FITTED VALUES
% 2nd
 display([' 2K0:'     num2str(x(1))]);   display([' 2V0:'      num2str(x(2))]);  
 display([' Error 2V0:'     num2str(Error_2V0)]);             display([' Error 2KT0:'   num2str(Error_2KT0)]);
% 3rd
display([' 3K0:'     num2str(m(1))]); display([' 3K_prime0:'   num2str(m(2))]); display([' 3V0:'      num2str(m(3))]);
display([' Error 3KT0:'     num2str(Error_3KT0)]);             display([' Error 3KT0_prime:'   num2str(Error_KT0_prime)]);    display([' Error 3V0:'     num2str(Error_3VT0)]);
%% PLOTTING 
% Pressure & Volume
figure(1)
errorbar(P_exp, V_exp, error_V,  'ko','MarkerSize',5, 'linewidth',0.8, 'MarkerFaceColor','k')
hold on
plot(P2BMEoS, V_pred, 'k--','LineWidth',0.5)
plot(P3BMEoS, V_pred, 'k-','LineWidth',0.5)
xlim([min(P_exp) max(P_exp)])
xlabel('Pressure (GPa)', 'FontWeight','bold','fontsize', 12)
ylabel('V (Å^3)', 'FontWeight','bold','fontsize', 12)
set(gca,'fontsize', 13,  'FontWeight','bold')
hold off

% Stress & Volume
figure(2)
plot(f3, F3, 'ko','MarkerSize',5, 'linewidth',0.8, 'MarkerFaceColor','k')
hold on
plot(f33, F33, 'k--','LineWidth',0.5)
xlabel('Eulerian Finite Strain, \itf_E', 'FontWeight','bold','fontsize', 12)
ylabel('Normalized Stress, \itF_E', 'FontWeight','bold','fontsize', 12)
set(gca,'fontsize', 13,  'FontWeight','bold')
hold off

%Deformation Energy
figure(3)
plot(P_exp, DE, 'ko','MarkerSize',5, 'linewidth',0.8, 'MarkerFaceColor','k')
hold on
plot(P2BMEoS, DE2, 'k--','LineWidth',0.5)
xlabel('Pressure (GPa)', 'FontWeight','bold','fontsize', 12)
ylabel('Deformation Energy', 'FontWeight','bold','fontsize', 12)
set(gca,'fontsize', 13,  'FontWeight','bold')
hold off

% Pressure & Bulk Modulus
figure(4)
plot(P2BMEoS(2:end), Kt2,  'ko','MarkerSize',5, 'linewidth', 0.8, 'MarkerFaceColor','k')
hold on
plot(P3BMEoS(2:end), Kt3,  'ko','MarkerSize',5, 'linewidth', 0.8, 'MarkerFaceColor','b')
xlim([min(P_exp) max(P_exp)])
xlabel('Pressure (GPa)', 'FontWeight','bold','fontsize', 12)
ylabel('K_T (GPa)', 'FontWeight','bold','fontsize', 12)
set(gca,'fontsize', 13,  'FontWeight','bold')
hold off
%--------------------------------------------------------------------------
