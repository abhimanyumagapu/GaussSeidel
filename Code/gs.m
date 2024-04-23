
%%
clc 
clear
Ybus = [3-9i,      -2+6i,        -1+3i,       0     ;
        -2+6i,      3.666-11i,   -0.666+2i,   -1+3i ;
        -1+3i,      -0.666+2i,   3.666-11i    -2+6i ;
        0,          -1+3i,       -2+6i,        3-9i  ];


P = [0; 0.5; -1; 0.3];
Q = [0; 0; 0.5i; -0.1i];

PV_buses = [2];
PV_voltage_magnitudes = [1.04]; %#ok<*NBRAK2>
PV_angle = [0];

n = size(Ybus, 1);

tolerance = 1e-4;

max_iter = 100;
alpha = 1.2;
Qlim_low = -800;
Qlim_high = 1000;  


V = ones(n, 1);

% Slack bus settings
slack_bus = 1;
V(slack_bus) = 1.04;  


% Gauss-Seidel Iteration
for iter = 1:max_iter
    V_prev = V;
    Q_prev = Q;
    

    for i = 1:n
        if i ~= slack_bus

                if ismember(i, PV_buses)
                    V(i) = PV_voltage_magnitudes(PV_buses == i) * exp(1j * PV_angle(PV_buses == i));
                    Q(i) = 0; 
                    for q = 1:n
                        Q(i) = Q(i) - abs(V(i)) * abs(V(q)) * abs(Ybus(i, q)) * sin(angle(Ybus(i, q)) + angle(V(q)) - angle(V(i)));
                        
                    end
                    disp(Q(i));

                    if Q(i) < Qlim_low
                        Q(i) = Qlim_low;
                        V(i) = (1/Ybus(i, i)) * (((P(i) - 1j*Q(i)) / conj(V(i))) - sum([Ybus(i, 1:i-1) * V(1:i-1); Ybus(i, i+1:n) * V_prev(i+1:n)]));
                        PV_angle(PV_buses == i) = angle(V(i));
                    elseif Q(i) > Qlim_high
                        Q(i) = Qlim_high;
                        
                        V(i) = (1/Ybus(i, i)) * (((P(i) - 1j*Q(i)) / conj(V(i))) - sum([Ybus(i, 1:i-1) * V(1:i-1); Ybus(i, i+1:n) * V_prev(i+1:n)]));
                        PV_angle(PV_buses == i) = angle(V(i));
                    else
                        
                        V_new = (1/Ybus(i, i)) * (((P(i) - 1i*Q(i)) / conj(V(i))) - sum([Ybus(i, 1:i-1) * V(1:i-1); Ybus(i, i+1:n) * V_prev(i+1:n)]));
                        V(i) = abs(V(i)) * exp(1j * angle(V_new));
                        PV_angle(PV_buses == i) = angle(V(i));
                    end  
                    
                else 
                    V_new = (1/Ybus(i, i)) * (((P(i) - Q(i)) / conj(V(i))) - sum([Ybus(i, 1:i-1) * V(1:i-1); Ybus(i, i+1:n) * V_prev(i+1:n)]));     
                    

                    V(i) = V(i) + alpha * (V_new - V(i)); 
                end

        end

    end
    
    if max(abs(V - V_prev)) < tolerance 
        if (Q(PV_buses) - Q_prev(PV_buses)) < tolerance
            fprintf('Converged in %d iterations\n', iter);
            break;
        end
    end

end




disp('Voltage Magnitude at each bus:');
disp(abs(V));
disp('Angles of each voltage (degrees, radians):');
disp(angle(V) * 180/pi + " ; " + angle(V));
disp('Voltage at each bus:')
disp(V);

 

%%