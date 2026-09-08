%% ME/RBE 4322 Homework 1 - Six-Bar Linkage Analysis
% Correct topology: AB, BC, CDE, EF, GF. Units: m, kg, s, N, N*m.
clc; clear; close all;

%% OPERATING ASSUMPTIONS
g=9.81; partsRequired=12500; operatingHours=9; partsPerRevolution=1;
inputRPM=partsRequired/(operatingHours*60*partsPerRevolution);
omega1=2*pi*inputRPM/60; alpha1=0;
artifactWeight=200; gripperDistance=1.843;

%% SOLIDWORKS PROPERTIES: [AB, BC, CDE, EF, GF], 6061 ALUMINUM
m=[8.02698,19.40883,33.48428,16.16478,35.18223];
J=[0.255848,3.389539,17.813936,1.970070,19.996261];
% Local CoM [x,y] from origins A, B, D, E, G respectively.
comLocal=[0,0.28630; 0,0.70785; 0.00006,1.24630; 0,0.58770; 0,1.29205];

%% FIRST-POSITION JOINTS
A=[1.400,0.485]; 
B0=[1.670,0.990]; 
C0=[0.255,1.035];
D=[0.285,0.055]; 
E0=[0.195,2.540]; 
F0=[-0.980,2.570]; 
G=[0.050,0.200];
lAB=norm(B0-A); 
lBC=norm(C0-B0); 
lCD=norm(C0-D);
lDE=norm(E0-D); 
lCE=norm(E0-C0); 
lEF=norm(F0-E0); 
lGF=norm(F0-G);
fprintf('Input speed = %.4f rad/s = %.4f rpm\n',omega1,inputRPM);
fprintf('Lengths (m): AB %.4f, BC %.4f, CD %.4f, DE %.4f, CE %.4f, EF %.4f, GF %.4f\n', lAB,lBC,lCD,lDE,lCE,lEF,lGF);

%% STORAGE
thetaDeg=0:1:360; n=numel(thetaDeg);
pos=nan(n,7,2); 
vel=nan(n,7,2); 
acc=nan(n,7,2);
omega=nan(n,5); 
alpha=nan(n,5); 
comPos=nan(n,5,2); 
comAcc=nan(n,5,2);
P=nan(n,2); 
Fstatic=nan(n,15); 
Fdynamic=nan(n,15);
thetaInput0=atan2(B0(2)-A(2),B0(1)-A(1)); 
Cprev=C0; 
Fprev=F0;

%% ALL POSITIONS
for i=1:n
    th=thetaInput0+deg2rad(thetaDeg(i));
    B=A+lAB*[cos(th),sin(th)];
    if i==1
        C=C0;
    else
        candidates=circleIntersections(B,lBC,D,lCD);
        if isempty(candidates), warning('No C solution at %.1f deg.',thetaDeg(i)); break; end
        [~,q]=min(vecnorm(candidates-Cprev,2,2)); C=candidates(q,:);
    end

    % CDE is one rigid ternary link grounded at D.
    phi3=atan2(C(2)-D(2),C(1)-D(1))-atan2(C0(2)-D(2),C0(1)-D(1));
    R3=rot2(phi3); E=D+(R3*(E0-D)')';
    if i==1
        F=F0;
    else
        candidates=circleIntersections(E,lEF,G,lGF);
        if isempty(candidates), warning('No F solution at %.1f deg.',thetaDeg(i)); break; end
        [~,q]=min(vecnorm(candidates-Fprev,2,2)); F=candidates(q,:);
    end
    Pnow=F+gripperDistance*(F-G)/norm(F-G); % beyond F, away from G

    %% VELOCITIES: link order [AB, BC, CDE, EF, GF]
    vB=omega1*kcross(B-A);
    M1=[kcross(C-B)',-kcross(C-D)'];
    if rcond(M1)<1e-10, warning('Loop 1 singular at %.1f deg.',thetaDeg(i)); break; end
    w23=M1\(-vB'); 
    w2=w23(1); 
    w3=w23(2);
    vC=w3*kcross(C-D); 
    vE=w3*kcross(E-D);
    M2=[kcross(F-E)',-kcross(F-G)'];
    if rcond(M2)<1e-10, warning('Loop 2 singular at %.1f deg.',thetaDeg(i)); break; end
    w45=M2\(-vE'); 
    w4=w45(1); 
    w5=w45(2); 
    vF=w5*kcross(F-G);

    %% ACCELERATIONS
    aB=alpha1*kcross(B-A)-omega1^2*(B-A);
    rhs1=-aB+w2^2*(C-B)-w3^2*(C-D);
    al23=M1\rhs1'; 
    al2=al23(1); 
    al3=al23(2);
    aC=al3*kcross(C-D)-w3^2*(C-D);
    aE=al3*kcross(E-D)-w3^2*(E-D);
    rhs2=-aE+w4^2*(F-E)-w5^2*(F-G);
    al45=M2\rhs2'; 
    al4=al45(1); 
    al5=al45(2);
    aF=al5*kcross(F-G)-w5^2*(F-G);

    %% CENTERS OF MASS FROM SOLIDWORKS
    S1=localPoint(A,B,comLocal(1,:)); 
    S2=localPoint(B,C,comLocal(2,:));
    S3=localPoint(D,E,comLocal(3,:)); 
    S4=localPoint(E,F,comLocal(4,:));
    S5=localPoint(G,F,comLocal(5,:));
    aS1=alpha1*kcross(S1-A)-omega1^2*(S1-A);
    aS2=aB+al2*kcross(S2-B)-w2^2*(S2-B);
    aS3=al3*kcross(S3-D)-w3^2*(S3-D);
    aS4=aE+al4*kcross(S4-E)-w4^2*(S4-E);
    aS5=al5*kcross(S5-G)-w5^2*(S5-G);

    joints=[A;B;C;D;E;F;G]; 
    jointV=[0,0;vB;vC;0,0;vE;vF;0,0];
    jointA=[0,0;aB;aC;0,0;aE;aF;0,0];
    pos(i,:,:)=joints; 
    vel(i,:,:)=jointV; 
    acc(i,:,:)=jointA;
    omega(i,:)=[omega1,w2,w3,w4,w5]; 
    alpha(i,:)=[alpha1,al2,al3,al4,al5];
    comPos(i,:,:)=[S1;S2;S3;S4;S5]; 
    comAcc(i,:,:)=[aS1;aS2;aS3;aS4;aS5]; 
    P(i,:)=Pnow;

    geom={A,B,C,D,E,F,G,S1,S2,S3,S4,S5,Pnow};
    Fstatic(i,:)=forceSolve(geom,m,J,zeros(5,2),zeros(1,5),artifactWeight,g);
    Fdynamic(i,:)=forceSolve(geom,m,J,[aS1;aS2;aS3;aS4;aS5],[alpha1,al2,al3,al4,al5],artifactWeight,g);
    Cprev=C; 
    Fprev=F;
end

%% KEEP VALID ROWS
valid=find(~isnan(omega(:,1))); 
thetaDeg=thetaDeg(valid);
pos=pos(valid,:,:); 
vel=vel(valid,:,:); 
acc=acc(valid,:,:); 
omega=omega(valid,:);
alpha=alpha(valid,:); 
comPos=comPos(valid,:,:); 
comAcc=comAcc(valid,:,:);
P=P(valid,:); 
Fstatic=Fstatic(valid,:); 
Fdynamic=Fdynamic(valid,:);

%% FIRST-POSITION TABLES
jointNames={'A';'B';'C';'D';'E';'F';'G'}; 
linkNames={'AB';'BC';'CDE';'EF';'GF'};
FirstPositionJoints=table(jointNames,squeeze(pos(1,:,1))',squeeze(pos(1,:,2))',squeeze(vel(1,:,1))',squeeze(vel(1,:,2))',squeeze(acc(1,:,1))',squeeze(acc(1,:,2))','VariableNames',{'Joint','x_m','y_m','vx_mps','vy_mps','ax_mps2','ay_mps2'});
disp(FirstPositionJoints)
FirstPositionLinks=table(linkNames,omega(1,:)',alpha(1,:)',m',J','VariableNames',{'Link','omega_rad_s','alpha_rad_s2','mass_kg','J_kg_m2'});
disp(FirstPositionLinks)
FirstPositionCoM=table(linkNames,squeeze(comPos(1,:,1))',squeeze(comPos(1,:,2))',squeeze(comAcc(1,:,1))',squeeze(comAcc(1,:,2))''VariableNames',{'Link','CoM_x_m','CoM_y_m','CoM_ax_mps2','CoM_ay_mps2'});
disp(FirstPositionCoM)
forceNames={'Ax';'Ay';'Bx';'By';'Cx';'Cy';'Dx';'Dy';'Ex';'Ey';'Fx';'Fy';'Gx';'Gy';'Tin'};
FirstPositionForces=table(forceNames,Fstatic(1,:)',Fdynamic(1,:)','VariableNames',{'Quantity','Static','Dynamic'});
disp(FirstPositionForces)
SolidWorksProperties=table(linkNames,m',(m*g)',comLocal(:,1),comLocal(:,2),J','VariableNames',{'Link','Mass_kg','Weight_N','Local_CoM_x_m','Local_CoM_y_m','Jzz_kg_m2'});
disp(SolidWorksProperties)

%% PLOTS
figure('Name','Joint Positions vs Input Rotation');
tiledlayout(2,1);
movingJoints = [2,3,5,6];
movingNames = {'B','C','E','F'};

nexttile;
hold on;
grid on;
for j = 1:length(movingJoints)
    plot(thetaDeg,squeeze(pos(:,movingJoints(j),1)),'LineWidth',1.2);
end
xlabel('Input Rotation (deg)');
ylabel('X Position (m)');
legend(movingNames,'Location','best');
title('Joint X Positions');

nexttile;
hold on;
grid on;
for j = 1:length(movingJoints)
    plot(thetaDeg,squeeze(pos(:,movingJoints(j),2)),'LineWidth',1.2);
end
xlabel('Input Rotation (deg)');
ylabel('Y Position (m)');
legend(movingNames,'Location','best');
title('Joint Y Positions');

figure('Name','Joint Linear Velocity Magnitudes');
hold on;
grid on;
for j = 1:length(movingJoints)
    jointNumber = movingJoints(j);
    velocityMagnitude = hypot(squeeze(vel(:,jointNumber,1)),squeeze(vel(:,jointNumber,2)));
    plot(thetaDeg,velocityMagnitude,'LineWidth',1.2);
end
xlabel('Input Rotation (deg)');
ylabel('Absolute Velocity (m/s)');
legend(movingNames,'Location','best');
title('Joint Linear Velocity Magnitudes');

figure('Name','Joint Linear Acceleration Magnitudes');
hold on;
grid on;
for j = 1:length(movingJoints)
    jointNumber = movingJoints(j);
    accelerationMagnitude = hypot(squeeze(acc(:,jointNumber,1)),squeeze(acc(:,jointNumber,2)));
    plot(thetaDeg,accelerationMagnitude,'LineWidth',1.2);
end
xlabel('Input Rotation (deg)');
ylabel('Absolute Acceleration (m/s^2)');
legend(movingNames,'Location','best');
title('Joint Linear Acceleration Magnitudes');

figure('Name','Joint paths'); hold on; axis equal; grid on;
for j=[2,3,5,6], plot(squeeze(pos(:,j,1)),squeeze(pos(:,j,2)),'LineWidth',1.3); end
plot(P(:,1),P(:,2),'k--','LineWidth',1.5); 
xlabel('x (m)'); 
ylabel('y (m)');
legend('B','C','E','F','Gripper','Location','best'); title('Joint and Gripper Paths');
figure('Name','Angular velocities'); 
plot(thetaDeg,omega,'LineWidth',1.2); 
grid on;
xlabel('Input rotation (deg)'); 
ylabel('\omega (rad/s)'); 
legend(linkNames,'Location','best');
title('Link Angular Velocities');
figure('Name','Angular accelerations'); 
plot(thetaDeg,alpha,'LineWidth',1.2); 
grid on;
xlabel('Input rotation (deg)'); 
ylabel('\alpha (rad/s^2)'); 
legend(linkNames,'Location','best');
title('Link Angular Accelerations');
figure('Name','Input torque'); 
plot(thetaDeg,Fstatic(:,15),'LineWidth',1.3); 
hold on;
plot(thetaDeg,Fdynamic(:,15),'LineWidth',1.3); 
grid on; xlabel('Input rotation (deg)');
ylabel('Input torque (N*m)'); 
legend('Static','Dynamic','Location','best'); 
title('Motor Input Torque');
figure('Name','Dynamic joint forces'); 
hold on; 
grid on;
for j=1:7, plot(thetaDeg,hypot(Fdynamic(:,2*j-1),Fdynamic(:,2*j)),'LineWidth',1.1); end
xlabel('Input rotation (deg)'); 
ylabel('Force magnitude (N)'); 
legend(jointNames,'Location','best');
title('Dynamic Joint Reaction Forces');

figure('Name','Static joint forces');
hold on;
grid on;
for j = 1:7
    staticMagnitude = hypot(Fstatic(:,2*j-1),Fstatic(:,2*j));
    plot(thetaDeg,staticMagnitude,'LineWidth',1.1);
end
xlabel('Input Rotation (deg)');
ylabel('Force Magnitude (N)');
legend(jointNames,'Location','best');
title('Static Joint Reaction Forces');

figure('Name','Center of Mass Accelerations');
hold on;
grid on;
for j = 1:5
    comAccelerationMagnitude = hypot(squeeze(comAcc(:,j,1)),squeeze(comAcc(:,j,2)));
    plot(thetaDeg,comAccelerationMagnitude,'LineWidth',1.2);
end
xlabel('Input Rotation (deg)');
ylabel('CoM Acceleration Magnitude (m/s^2)');
legend(linkNames,'Location','best');
title('Link Center-of-Mass Accelerations');

%% MAXIMUM RESULTS
[v,ix]=max(abs(Fstatic(:,15))); 
fprintf('\nMax absolute static torque = %.4f N*m at %.1f deg\n',v,thetaDeg(ix));
[v,ix]=max(abs(Fdynamic(:,15))); 
fprintf('Max absolute dynamic torque = %.4f N*m at %.1f deg\n',v,thetaDeg(ix));
for j=1:7
    staticMag = hypot(Fstatic(:,2*j-1),Fstatic(:,2*j));
    [vStatic,ixStatic] = max(staticMag);
    fprintf('Max static force at %s = %.3f N at %.1f deg\n',jointNames{j},vStatic,thetaDeg(ixStatic));
    mag=hypot(Fdynamic(:,2*j-1),Fdynamic(:,2*j)); 
    [v,ix]=max(mag);
    fprintf('Max dynamic force at %s = %.3f N at %.1f deg\n',jointNames{j},v,thetaDeg(ix));
end

%% LOCAL FUNCTIONS
function v = kcross(r)
    v = [-r(2),r(1)];
end

function R = rot2(a)
    R = [cos(a),-sin(a);sin(a),cos(a)];
end

function p = localPoint(origin,axisPoint,xy)
    ey = (axisPoint-origin)/norm(axisPoint-origin);
    ex = [ey(2),-ey(1)];
    p = origin+xy(1)*ex+xy(2)*ey;
end

function pts = circleIntersections(c1,r1,c2,r2)
    d = norm(c2-c1);
    tol = 1e-10;
    if d > r1+r2+tol || d < abs(r1-r2)-tol || d < tol
        pts = [];
        return;
    end
    a = (r1^2-r2^2+d^2)/(2*d);
    h = sqrt(max(r1^2-a^2,0));
    p = c1+a*(c2-c1)/d;
    q = h*[-(c2(2)-c1(2)),c2(1)-c1(1)]/d;
    pts = [p+q;p-q];
end

function solution = forceSolve(geom,m,J,aS,alphas,loadWeight,g)
    x0 = zeros(15,1);
    r0 = forceResidual(x0,geom,m,J,aS,alphas,loadWeight,g);
    K = zeros(15);
    for k = 1:15
        ek = zeros(15,1);
        ek(k) = 1;
        K(:,k) = forceResidual(ek,geom,m,J,aS,alphas,loadWeight,g)-r0;
    end
    if rcond(K) < 1e-12
        solution = nan(1,15);
    else
        solution = (K\(-r0))';
    end
end

function r = forceResidual(x,Q,m,J,aS,al,loadWeight,g)
% Links: AB, BC, CDE, EF, GF
    A=Q{1}; 
    B=Q{2}; 
    C=Q{3}; 
    D=Q{4}; 
    E=Q{5}; 
    F=Q{6}; 
    G=Q{7};
    
    S1=Q{8}; 
    S2=Q{9}; 
    S3=Q{10}; 
    S4=Q{11}; 
    S5=Q{12}; 
    
    P=Q{13};
    FA=x(1:2)'; 
    FB=x(3:4)'; 
    FC=x(5:6)'; 
    FD=x(7:8)';
    FE=x(9:10)'; 
    FF=x(11:12)'; 
    FG=x(13:14)'; 
    Tin=x(15);
    
    W = @(mass)[0,-mass*g];
    L = [0,-loadWeight];
    mom = @(ra,f)ra(1)*f(2)-ra(2)*f(1);
    r = zeros(15,1);
    k = 0;

    f = FA+FB+W(m(1))-m(1)*aS(1,:);
    z = mom(A-S1,FA)+mom(B-S1,FB)+Tin-J(1)*al(1);
    r(k+(1:3)) = [f,z]; 
    k=k+3;

    f = -FB+FC+W(m(2))-m(2)*aS(2,:);
    z = mom(B-S2,-FB)+mom(C-S2,FC)-J(2)*al(2);
    r(k+(1:3)) = [f,z]; 
    k=k+3;

    f = -FC+FD+FE+W(m(3))-m(3)*aS(3,:);
    z = mom(C-S3,-FC)+mom(D-S3,FD)+mom(E-S3,FE)-J(3)*al(3);
    r(k+(1:3)) = [f,z]; 
    k=k+3;

    f = -FE+FF+W(m(4))-m(4)*aS(4,:);
    z = mom(E-S4,-FE)+mom(F-S4,FF)-J(4)*al(4);
    r(k+(1:3)) = [f,z]; 
    k=k+3;

    f = -FF+FG+W(m(5))+L-m(5)*aS(5,:);
    z = mom(F-S5,-FF)+mom(G-S5,FG)+mom(P-S5,L)-J(5)*al(5);
    r(k+(1:3)) = [f,z];
end
