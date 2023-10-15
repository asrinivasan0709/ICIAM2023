clc; clear; close all;
% Leapfrog 2nd order, Mimetic
% ASrinivasan, 9Jan23
addpath('C:\Anand\Acer_Data\SDSU\MOLE\mole-master\mole_MATLAB')
addpath('C:\Anand\Acer_Data\SDSU\Mimetic\Quadratures\RelaxationRK\Matlab\WaveEqn\Rev2')

Ne = 10*[8, 16]'; 
tEnd = 10;  
CFL = 0.95;  cV = 1; % Lax W conserves energy only for CFL = 1

[uEx, uOut, tOut, eOut, normVal, xGrid] = LEAP(Ne, tEnd, CFL, cV); 
[uExW, uOutW, tOutW, eOutW, normValW, xGridW] = WEND(Ne, tEnd, CFL, cV); 


% Convergence
if size(Ne, 1) > 1
    for i = 1:size(Ne, 1) - 1
        Cv(i, 1) = 1/log(2)*log(normVal(i, 1)/normVal(i+1, 1));         
    end
end


% Plots
figure
plot(xGrid, uOut(:, end), '-.k', 'LineWidth', 1.5)
hold on; plot(xGridW, uOutW(:, end), '--or', 'LineWidth', 1.0, 'MarkerSize', 2)
hold on; plot(xGrid, uEx, '-b', 'LineWidth', 1.5)
legend('LF-MIM2', 'Lax-Wendroff', 'Exact')
title({'1D Advection Equation, Solution at t = 10 s'},{'u_t + \nabla \cdot u = 0'})
xlabel('x'); ylabel('u(x, t)')

figure
plot(tOut, abs( (eOut - eOut(1))./eOut(1) ), '-.k', 'LineWidth', 1)
hold on; plot(tOutW, abs( (eOutW - eOutW(1))./eOutW(1) ), '-r', 'LineWidth', 1)
legend('LF-MIM2', 'Lax-Wendroff')
set(gca, 'YScale', 'log');  
title({'Energy Evolution, 1D Advection Equation'}, {'u_t + \nabla \cdot u = 0'})
xlabel('time [s]'); ylabel('Relative error, |(E_n - E_0)/E_0|'); 
grid on; ylim([10E-6, 10E-1]); 


function [uEx, uOut, tOut, eOut, normVal, xGrid] = LEAP(Ne, tEnd, CFL, cV)

    for j = 1:size(Ne, 1)
        fprintf('Executing Leapfrog, iteration %i of %i ... \n', j, size(Ne, 1)); 
        NElem = Ne(j, 1); 
        xmin = -1; xmax = 1; 
        dh = (xmax-xmin)/NElem; 
        dt =  CFL*dh/cV; 
        iter = floor(tEnd/dt); 
        
        xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
        u0 = u0Func(xGrid); 
        
        D = div(2, NElem, dh);    
        IntD = interpC2N(NElem, 2);
        
        % Convert IntD to periodic
        IntD(2,:) = zeros(1,NElem+2); 
        IntD(end-1,:) = zeros(1,NElem+2);     
        IntD(2,2:4) = IntD(3,3:5); IntD(2,end-1) = IntD(3,2); 
        IntD(end-1,:) = fliplr(IntD(2,:));    
        
        uOut = zeros(NElem+2, iter); tOut = zeros(iter, 1); 
        uOut(:,1) = u0;
        u1 = u0 - dt*D*IntD*u0; 
        uOut(:,2) = u1; 
        tOut(1,1) = 0; tOut(2,1) = dt; 
        eu0 = u0; eu1 = u1;
        eOut = [eu0'*eu0; eu1'*eu1]; 
        t = dt;
        
        for i = 3:iter
            u2 = u0 - 2*dt*cV*D*IntD*u1; 
    
            % Periodic BC
            u2(1) = u2(2:3,1)'*IntD(3,4:5)' + ...
                    u2(end-2:end-1)'*IntD(3,2:3)'; 
            u2(end) = u2(1);     
    
            uOut(:,i) = u2;
            u0 = u1; u1 = u2; 
            t = t + dt;
            tOut(i,1) = t;      
    %         eu2 = G*u2; 
            eu2 = u2; 
            eOut(i,1) = eu2'*eu2; 
        end
        
        uEx = uExact(xGrid, xmin, xmax, cV, tOut(end)); 
        normVal(j, 1) = norm(uEx - uOut(:, end), 'inf'); 
    
    end

end

function [uEx, uOut, tOut, eOut, normVal, xNod] = WEND(Ne, tEnd, CFL, cV)
% Lax Wendroff scheme, 

    for j = 1:size(Ne, 1)
        fprintf('Executing Lax Wendroff, iteration %i of %i ... \n', j, size(Ne, 1)); 
        NElem = Ne(j, 1); 
        xmin = -1; xmax = 1; 
        dh = (xmax-xmin)/(NElem-1); 
        dt =  CFL*dh/cV; 
        iter = floor(tEnd/dt); 
        
        xNod = [xmin xmin+dh:dh:xmax-dh xmax]';  % Nodal grid
        u0 = u0Func(xNod); 

        % A matrix, with periodic BC
        bm1 = CFL/2*(CFL + 1);
        b0 = 1 - CFL^2;
        b1 = CFL/2*(CFL - 1); 

        A = 0*eye(NElem);
        for k = 2:NElem-1
            A(k, k-1:k+1) = [bm1, b0, b1];             
        end
        A(1, 1:2) = [b0, b1]; A(1, end) = bm1;
        A(end, end-1:end) = [bm1, b0]; A(end, 1) = b1; 

        A = sparse(A); 

        uOut = zeros(NElem, iter); tOut = zeros(iter, 1); 
        uOut(:,1) = u0;
        tOut(1,1) = 0; 
            
        eOut(1,1) = u0'*u0; 
        t = dt;
        
        for i = 2:iter
            u1 = A*u0; 
            uOut(:,i) = u1;
            u0 = u1;  
            t = t + dt;
            tOut(i,1) = t;                  
            eOut(i,1) = u1'*u1; 

        end

        uEx = uExact(xNod, xmin, xmax, cV, tOut(end)); 
        normVal(j, 1) = norm(uEx - uOut(:, end), 'inf'); 

    end



end

function u0 = u0Func(xGrid)
    muu = 0.; sg = 0.01;
    u0 = 1/sqrt(sg*2*pi)*exp(-(xGrid-muu).^2/2/sg); %

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