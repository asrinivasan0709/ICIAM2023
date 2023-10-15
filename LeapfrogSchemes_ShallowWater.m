clc; clear; close all;
% Leapfrog schemes + Mimetic
% ASrinivasan, 29Apr23, Shallow Water eqns

addpath('...\mole-master\mole_MATLAB')

% for energy evolution
Ne = 5*[64, 128]'; 
xmin = -50; xmax = 50; 
tEnd = 25;  
CFL = 0.3;  cV = 1;  

[uOutWFL4, tOutWFL4, eOutWFL4, xGridWFL4, nWFL4, tExecWFL4] = WFL4(Ne, tEnd, CFL, cV, xmin, xmax); 
[uOutHORA, tOutHORA, eOutHORA, xGridHORA, nHORA, tExecHORA] = LiTren(Ne, tEnd, CFL, cV, xmin, xmax); 
[uOutRK4, tOutRK4, eOutRK4, xGridRK4, nRK4, tExecRK4] = RK3(Ne, tEnd, CFL, cV, xmin, xmax); 
[uOutRRK4, tOutRRK4, eOutRRK4, xGridRRK4, gamRRK, nRRK, tExecRRK] = RelaxRK4(Ne, tEnd, CFL, cV, xmin, xmax); 
[uOutRRK4r, tOutRRK4r, eOutRRK4r, xGridRRK4r, gamRRKr, nRRKr, tExecRRKr] = RelaxRK4root(Ne, tEnd, CFL, cV, xmin, xmax); 

NeOut = Ne(end); 

% Plots
figure
hold on; plot(xGridWFL4, uOutWFL4(1:NeOut+2, end-1), 'k', 'LineWidth', 1.0, 'MarkerSize', 2)
hold on; plot(xGridHORA, uOutHORA(1:NeOut+2, end-1), '-.r', 'LineWidth', 1.0, 'MarkerSize', 2)
hold on; plot(xGridRK4, uOutRK4(1:NeOut+2, end), 'ob', 'LineWidth', 1.0, 'MarkerSize', 1)
hold on; plot(xGridRRK4, uOutRRK4(1:NeOut+2, end), 'sr', 'LineWidth', 1.0, 'MarkerSize', 3)
hold on; plot(xGridRRK4r, uOutRRK4r(1:NeOut+2, end), '*k', 'LineWidth', 2.0, 'MarkerSize', 3)
legend( ' WFL4', 'HORA', 'RK3', 'RRK4', 'RRK4-root', 'Location', 'best') %
title({'Shallow Water Equations'},{'\eta_t + \nabla \cdot ((D+\eta)u) = 0, u_t + \nabla \cdot \eta + u \nabla \cdot u = 0'})
xlabel('x'); ylabel('u(x, t)')

figure
hold on; plot(xGridWFL4, uOutWFL4(NeOut+3:end, end-1), 'k', 'LineWidth', 1.0, 'MarkerSize', 2)
hold on; plot(xGridHORA, uOutHORA(NeOut+3:end, end-1), '-.r', 'LineWidth', 1.0, 'MarkerSize', 2)
hold on; plot(xGridRK4, uOutRK4(NeOut+3:end, end), 'ob', 'LineWidth', 1.0, 'MarkerSize', 1)
hold on; plot(xGridRRK4, uOutRRK4(NeOut+3:end, end), 'sr', 'LineWidth', 1.0, 'MarkerSize', 3)
hold on; plot(xGridRRK4r, uOutRRK4r(NeOut+3:end, end), '*k', 'LineWidth', 2.0, 'MarkerSize', 3)
legend( ' WFL4', 'HORA', 'RK3', 'RRK4', 'RRK4-root', 'Location', 'best') %
title({'Shallow Water Equations'},{'\eta_t + \nabla \cdot ((D+\eta)u) = 0, u_t + \nabla \cdot \eta + u \nabla \cdot u = 0'})
xlabel('x'); ylabel('v(x, t)')


figure
hold on; plot(tOutWFL4, abs( eOutWFL4./eOutWFL4(1) - 1. ), '-.k', 'LineWidth', 1.5)
hold on; plot(tOutHORA, abs( eOutHORA./eOutHORA(1) - 1. ), '-.r', 'LineWidth', 1.5)
hold on; plot(tOutRK4, abs( eOutRK4./eOutRK4(1) - 1. ), 'ob', 'LineWidth', 2, 'MarkerSize', 3) % )
hold on; plot(tOutRRK4, abs( eOutRRK4./eOutRRK4(1) - 1), '-sr', 'LineWidth', 2, 'MarkerSize', 1)
hold on; plot(tOutRRK4r, abs( eOutRRK4r./eOutRRK4r(1) - 1 ), '*k', 'LineWidth', 2, 'MarkerSize', 1)
legend( ' WFL4', 'HORA', 'RK3', 'RRK4', 'RRK4-root', 'Location', 'best') %
set(gca, 'YScale', 'log');  
title({'Shallow Water Equations'},{'\eta_t + \nabla \cdot ((D+\eta)u) = 0, u_t + \nabla \cdot \eta + u \nabla \cdot u = 0'})
xlabel('time [s]'); ylabel('Relative error, |(E_n - E_0)/E_0|'); 
grid on; %ylim([10E-8, 10E-1]); 



function [uOut, tOut, eOut, xGrid, normVal, tExec] = RK4(Ne, tEnd, CFL, cV, xmin, xmax)
% RK4 

    for j = 1:size(Ne, 1)
        fprintf('Executing RK4, iteration %i of %i ... \n', j, size(Ne, 1)); 

        tic; 
        NElem = Ne(j, 1); 
        dh = (xmax-xmin)/NElem; 
        dt =  CFL*dh/cV; 
        iter = floor(tEnd/dt); 
        
        xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
        u0 = u0Func(xGrid); 
    
        C(1,1) = 0;        C(2,1) = 1/2;
        C(3,1) = 1/2;        C(4,1) = 1;
    
        A(2,1) = 1/2;
        A(3,1) = 0;   A(3,2) = 1/2;
        A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;     
    
        B(1,1) = 1/6;        B(2,1) = 1/3;
        B(3,1) = 1/3;        B(4,1) = 1/6;

        D = div(4, NElem, dh);    
        IntD = interpC2N(NElem, 4); 
        DIntD = D*IntD; 

        uOut = zeros(2*(NElem+2), iter); tOut = zeros(iter, 1); 
        uOut(:,1) = u0;
        tOut(1,1) = 0; 
        y = u0; 
            
        eOut(1,1) = EnergyEval(u0); 
        t = dt;
        
        for i = 2:iter
            z1 = y; 
            [k1] = FuncEval(DIntD, z1); 
                       
            % Step 2        
            z2 =  y + dt*A(2,1)*k1;       
            [k2] = FuncEval(DIntD, z2); 
    
            % Step 3    
            z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
            [k3] = FuncEval(DIntD, z3); 
            
            % Step 4    
            z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
            [k4] = FuncEval(DIntD, z4); 
            
            BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
                B(4,1)*k4; 
                

            uNew = y + dt*BKsum;  

            uOut(:,i) = uNew;
            y = uNew;  
            t = t + dt;
            tOut(i,1) = t;         

            eOut(i,1) = EnergyEval(uNew); 

        end
        tExec(j, 1) = toc; 
        uEx = uExact(xGrid, xmin, xmax, cV, tOut(end-1)); 
        normVal(j, 1) = norm(uEx - uOut(:, end), 'inf'); 
         
    end



end

function [uOut, tOut, eOut, xGrid, gamRRK, normVal, tExec] = RelaxRK4(Ne, tEnd, CFL, cV, xmin, xmax)
% Relaxation RK4 

    for j = 1:size(Ne, 1)
        fprintf('Executing Relaxation RK4, iteration %i of %i ... \n', j, size(Ne, 1)); 
        tic; 
        NElem = Ne(j, 1); 
        dh = (xmax-xmin)/NElem; 
        dt =  CFL*dh/cV; 
        iter = floor(tEnd/dt); 
        
        xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
        u0 = u0Func(xGrid); 
    
        C(1,1) = 0;
        C(2,1) = 1/2;
        C(3,1) = 1/2;
        C(4,1) = 1;
    
        A(2,1) = 1/2;
        A(3,1) = 0;   A(3,2) = 1/2;
        A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;     
    
        B(1,1) = 1/6; 
        B(2,1) = 1/3;
        B(3,1) = 1/3;
        B(4,1) = 1/6;

        D = div(4, NElem, dh);    
        IntD = interpC2N(NElem, 4); 
        DIntD = D*IntD; 

        uOut = zeros(2*(NElem+2), iter); tOut = zeros(iter, 1); 
        uOut(:,1) = u0;
        tOut(1,1) = 0; 
        y = u0; 
            
        eOut(1,1) = EnergyEval(u0); 
        E0 = eOut(1,1); 
        t = dt;
        
        for i = 2:iter
            z1 = y; 
            [k1] = FuncEval(DIntD, z1); 
            
    
            % Step 2        
            z2 =  y + dt*A(2,1)*k1;       
            [k2] = FuncEval(DIntD, z2); 

            % Step 3    
            z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
            [k3] = FuncEval(DIntD, z3); 

            % Step 4    
            z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
            [k4] = FuncEval(DIntD, z4); 
            
            BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + B(4,1)*k4; 
               
            EE = y(1:NElem+2,1); UU = y(NElem+3:end,1);
            dEE = BKsum(1:NElem+2,1); dUU = BKsum(NElem+3:end,1);
            
            BB = 2*EE'*dEE + 2*UU'*dUU + 2*UU'*(EE.*dUU) + dEE'*(UU.*UU);                
            CC = dEE'*dEE + dUU'*dUU + EE'*(dUU.*dUU) + 2*dEE'*(UU.*dUU);
            DD = dEE'*(dUU.*dUU);

            gam = (-CC + (CC^2 - 4*DD*BB)^(1/2))/(2*DD*dt);                


            gamRRK(i,1) = gam; 
            uNew = y + gam*dt*BKsum;              

            uOut(:,i) = uNew;
            y = uNew;  
            t = t + gam*dt;
            tOut(i,1) = t;                  
            eOut(i,1) = EnergyEval(uNew); 

        end
        tExec(j, 1) = toc; 
        uEx = uExact(xGrid, xmin, xmax, cV, tOut(end-1)); 
        normVal(j, 1) = norm(uEx - uOut(:, end), 'inf'); 

    end



end

function [uOut, tOut, eOut, xGrid, normVal, tExec] = WFL4(Ne, tEnd, CFL, cV, xmin, xmax)
% Leapfrog Williams modified 4th order filter


    for j = 1:size(Ne, 1)
        fprintf('Executing WFL 4th order, iteration %i of %i ... \n', j, size(Ne, 1)); 

        tic; 
        NElem = Ne(j, 1); 
        dh = (xmax-xmin)/NElem; 
        dt =  CFL*dh/cV; 
        iter = floor(tEnd/dt); 
        
        xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
        u0 = u0Func(xGrid); 
        
        D = div(4, NElem, dh);    
        IntD = interpC2N(NElem, 4); DIntD = D*IntD; 
        
        uOut = zeros(2*(NElem+2), iter); tOut = zeros(iter, 1); 
        uOut(:,1) = u0;
        u1 = u0 + dt*FuncEval(DIntD, u0); 
        u2 = u1 + dt*FuncEval(DIntD, u1); 
        u3 = u2 + dt*FuncEval(DIntD, u2); 
        uOut(:,2:4) = [u1, u2, u3]; 
        tOut(1,1) = 0; tOut(2:4,1) = [dt, 2*dt, 3*dt]'; 
        eu0 = EnergyEval(u0); eu1 = EnergyEval(u1); 
        eu2 = EnergyEval(u2); eu3 = EnergyEval(u3);
        eOut = [eu0; eu1; eu2; eu3]; 
        t = 4*dt;
        v1 = u1; w1 = u1; 
        v2 = u2; w2 = u2; 
        v3 = u3; w3 = u3; 

        nu = 0.1; alpha = 0.5; 
        gamma = (5-9*nu)/(2*(4-7*nu));


        for i = 5:iter
            
            Fv = cV*FuncEval(DIntD, v3);
            Fw = cV*FuncEval(DIntD, w3);
     
            w4 = u2 + 2*dt*(gamma*Fv + (1-gamma)*Fw); 
            fDisp = nu*( w4 - 4*v3 + 6*u2 -4*u1 + u0 );
            u3 = v3 + fDisp*alpha;  
            v4 = w4 + fDisp*(alpha - 1);            

            uOut(:,i-1) = u3;
            u0 = uOut(:, i-3); u1 = uOut(:, i-2); u2 = u3; 
            v1 = v2; v2 = v3; v3 = v4; 
            w1 = w2; w2 = w3; w3 = w4; 

            t = t + dt;
            tOut(i,1) = t;      
            eu2 = u1; 
            eOut(i,1) = EnergyEval(eu2); 
    

        end
        tExec(j,1) = toc; 
        uEx = uExact(xGrid, xmin, xmax, cV, tOut(end-2)); 
        normVal(j, 1) = norm(uEx - uOut(:, end-1), 'inf'); 
    
    end

end

function [uOut, tOut, eOut, xGrid, gamRRK, normVal, tExec] = RelaxRK4root(Ne, tEnd, CFL, cV, xmin, xmax)
% Relaxation RK4 

    for j = 1:size(Ne, 1)
        fprintf('Executing Relaxation RRK4 root solver, iteration %i of %i ... \n', j, size(Ne, 1)); 
        tic; 
        NElem = Ne(j, 1); 
        dh = (xmax-xmin)/NElem; 
        dt =  CFL*dh/cV; 
        iter = floor(tEnd/dt); 
        
        xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
        u0 = u0Func(xGrid); 
    
        D = div(4, NElem, dh);    
        IntD = interpC2N(NElem, 4); 
        DIntD = D*IntD; 

        uOut = zeros(2*(NElem+2), iter); tOut = zeros(iter, 1); 
        uOut(:,1) = u0;
        tOut(1,1) = 0; 
        y = u0; 
            
        eOut(1,1) = EnergyEval(u0); 
        E0 = eOut(1,1); 
        t = dt;
        MaxIt = 30; MinErr = 1E-10;
        
        for i = 2:iter
        x0 = 0.99; x1 = 1.1; %1.000001;                 

        % Secant root finder for gamma
        for k = 1:MaxIt
            u_x0 = RKSolver(y, DIntD, dt, x0, NElem);
            u_x1 = RKSolver(y, DIntD, dt, x1, NElem);
            
            ut0 = u_x0;  ut1 = u_x1;
                        
            E_x0 = abs(EnergyEval(ut0) - E0);
            E_x1 = abs(EnergyEval(ut1) - E0); 

            x2 = x1 - E_x1*(x1 - x0)/(E_x1 - E_x0); 
            x0 = x1; x1 = x2;         
            if abs(E_x1 - E_x0) <= MinErr, break, end
    
        end

               
        MaxItCount(i, 1) = k; 
        gam = x0; 
        uNew = u_x0; 

            gamRRK(i,1) = gam; 

            uOut(:,i) = uNew;
            y = uNew;  
            t = t + gam*dt;
            tOut(i,1) = t;                  
            eOut(i,1) = EnergyEval(uNew); 

        end
        tExec(j,1) = toc; 
        uEx = uExact(xGrid, xmin, xmax, cV, tOut(end-1)); 
        normVal(j, 1) = norm(uEx - uOut(:, end), 'inf'); 

    end



end

function [uOut, tOut, eOut, xGrid, normVal, tExec] = LiTren(Ne, tEnd, CFL, cV, xmin, xmax)
% hoRA, Li & Trenchea

    for j = 1:size(Ne, 1)
        fprintf('Executing Li/Trenchea hoRA, iteration %i of %i ... \n', j, size(Ne, 1)); 
        tic; 
        NElem = Ne(j, 1); 
        dh = (xmax-xmin)/NElem; 
        dt =  CFL*dh/cV; 
        iter = floor(tEnd/dt); 
        
        xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
        u0 = u0Func(xGrid); 
        
        D = div(4, NElem, dh);    
        IntD = interpC2N(NElem, 2); DIntD = D*IntD; 
        
        uOut = zeros(2*(NElem+2), iter); tOut = zeros(iter, 1); 
        uOut(:,1) = u0;
        u1 = u0 + dt*FuncEval(DIntD, u0); 
        u2 = u1 + dt*FuncEval(DIntD, u1); 

        uOut(:,2:3) = [u1, u2]; 
        tOut(1:3,1) = [0; dt; 2*dt];  
        eu0 = EnergyEval(u0); eu1 = EnergyEval(u1); eu2 = EnergyEval(u2); 
        eOut = [eu0; eu1; eu2]; 
        t = 3*dt;
        v1 = u1; v2 = u2; 
        beta = 0.4; 

        for i = 4:iter
            
 
            v3 = u1 + 2*dt*FuncEval(DIntD, v2); 
            u2 = v2 + beta/2*( v3 - 2*v2 + u1 ) - beta/2*(v2 - 2*u1 + u0);
 
            uOut(:,i-1) = u2;
            u0 = u1; u1 = u2; 
            v1 = v2; v2 = v3; 
            t = t + dt;
            tOut(i,1) = t;      
            eu2 = u2; 
            eOut(i,1) = EnergyEval(eu2); 
    

        end
        tExec(j,1) = toc; 
        uEx = uExact(xGrid, xmin, xmax, cV, tOut(end-2)); 
        normVal(j, 1) = norm(uEx - uOut(:, end-1), 'inf'); 
    
    end

end


function [uOut, tOut, eOut, xGrid, normVal, tExec] = RK3(Ne, tEnd, CFL, cV, xmin, xmax)
% RK4 

    for j = 1:size(Ne, 1)
        fprintf('Executing RK3, iteration %i of %i ... \n', j, size(Ne, 1)); 

        tic; 
        NElem = Ne(j, 1); 
        dh = (xmax-xmin)/NElem; 
        dt =  CFL*dh/cV; 
        iter = floor(tEnd/dt); 
        
        xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
        u0 = u0Func(xGrid); 
    
        D = div(4, NElem, dh);    
        IntD = interpC2N(NElem, 4); 
        DIntD = D*IntD; 

        uOut = zeros(2*(NElem+2), iter); tOut = zeros(iter, 1); 
        uOut(:,1) = u0;
        tOut(1,1) = 0; 
        y = u0; 
            
        eOut(1,1) = EnergyEval(u0); 
        t = dt;
        
        for i = 2:iter
            [k1] = y + 1/3*dt*FuncEval(DIntD, y); %k1 = periodicBC(k1);
                       
            % Step 2        
            [k2] = y + 1/2*dt*FuncEval(DIntD, k1); %k2 = periodicBC(k2);
    
            % Step 3    
            [k3] = y + dt*FuncEval(DIntD, k2); %k3 = periodicBC(k3);
            
                
            uNew = k3;  %uNew = periodicBC(uNew);

            uOut(:,i) = uNew;
            y = uNew;  
            t = t + dt;
            tOut(i,1) = t;                  
            eOut(i,1) = EnergyEval(uNew); 

        end
        tExec(j, 1) = toc; 
        uEx = uExact(xGrid, xmin, xmax, cV, tOut(end-1)); 
        normVal(j, 1) = norm(uEx - uOut(:, end), 'inf'); 
         
    end



end





function uOut = FuncEval(DIntD, u1)
% evaluate the RHS function
    NElem = size(u1, 1)/2; 
    z4u = u1(1:NElem, 1); z4v = u1(NElem+1:end, 1);
    uOut = -[DIntD*((1 + z4u).*z4v); 
            DIntD*z4u + z4v.*DIntD*z4v]; 


end

function uNew = RKSolver(y, DIntD, dt, gam, NElem)

C(1,1) = 0;   C(2,1) = 1/2;
C(3,1) = 1/2; C(4,1) = 1;

A(2,1) = 1/2;
A(3,1) = 0;   A(3,2) = 1/2;
A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;    

B(1,1) = 1/6; B(2,1) = 1/3;
B(3,1) = 1/3; B(4,1) = 1/6;

    % Step 1
    z1 = y; 
    [k1] = FuncEval(DIntD, z1); 
           

    % Step 2        
    z2 =  y + dt*A(2,1)*k1;       
    [k2] = FuncEval(DIntD, z2); 

    % Step 3    
    z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
    [k3] = FuncEval(DIntD, z3); 
    
    % Step 4    
    z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
    [k4] = FuncEval(DIntD, z4); 
    
    
    uNew = y  + gam*dt*(B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
        B(4,1)*k4);          

end

function u = periodicBC(u)
    
    u(1) = 0.5*(u(end-1) + u(2));
    u(end) = u(1); 

end

function u0 = u0Func(xGrid)
    u = 1 + 0.5/5*exp(-1*xGrid.^2);    
    v = zeros(size(xGrid,1),1);
    
    u0 = [u; v]; 

end

function uEx = uExact(xGrid, xmin, xmax, cV, i)

    t1 = mod(i, xmax - xmin); 
    if mod(cV*i, xmax - xmin) == 0
        xt = 0; 
        uEx = u0Func(xGrid - cV*xt);
    else
        xt = floor(t1/xmax)*xmax - mod(cV*t1, xmax);       
        uEx = u0Func(xGrid + cV*xt); 
    end



end

function E = EnergyEval(u0)

    NN = size(u0,1)/2; 

    u = u0(1:NN, 1);
    v = u0(NN+1:end, 1);
    E = u'*u + v'*v + u'*v.^2; 



end
