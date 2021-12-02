% MT-Praktikum
% Versuch:  MT_Scope
% Ahmed Habib, 02.12.2021 
clc;
clear;
%% 3.2.1 Nominal characteristic curves of the 4 Thermistors

    devices = {'NTC-2K2','NTC-100K','Pt100','KTY-10-6'};
    T       = (20:1:85);                                % °C

    % NTC-2k2
    R_25    = 2200;                                     % Ohms
    B       = 3977;                                     % K

    R_1     =  arrayfun(@(t) R_NTC(R_25,B,t),T);        % Ohms 
    %rn     =  R_25 / R_1(6);                           % 1 
    R_1n    =  arrayfun(@(t) R_NTC(1,B,t),T);           % 1
    
    

    % NTC-100K
    R_25    = 100000;                                   % Ohms
    B       = 4190;                                     % K

    R_2     =  arrayfun(@(t) R_NTC(R_25,B,t),T);        % Ohms
    %rn      =  R_25 / R_1(6);                          % 1 
    R_2n    =  arrayfun(@(t) R_NTC(1,B,t),T);           % 1 


    % Pt100
    R_0     = 100;                                      % Ohms       
    A       = 3.9083*10^-3;                             % 1/k¹ 
    B       = -5.775*10^-7;                             % 1/k²
    C       = -4.183*10^-12;                            % 1/k²

    R_3     = arrayfun(@(t) R_platin(R_0,A,B,C,t),T);   % Ohms
    %rn      = R_0 / R_3(6);                            % 1
    R_3n    = arrayfun(@(t) R_platin(1,A,B,C,t),T);     % 1

    % KTY-10-6
    R_25    = 2 * 10^3;                                 % Ohms
    A       = 7.88 * 10^-3;                             % 1/K
    B       = 1.973 * 10^-5;                            % 1/K²

    R_4     = arrayfun(@(t) R_silizium(R_25,A,B,t),T);  % Ohms 
    %rn      = R_25 / R_4(6);                           % 1 
    R_4n    = arrayfun(@(t) R_silizium(1,A,B,t),T);     % Ohms 

   
%% Write tables as csv 
    headers = {'T', 'R_1','R_2','R_3','R_4'};
    tbl_1   = table(T',R_1',R_2',R_3',R_4','VariableNames', headers);
    tbl_1n  = table(T',R_1n',R_2n',R_3n',R_4n','VariableNames', headers);

    base    = "table_321";
    ext     = ".csv";
    fname   = base + ext;
    fname_n = base + "_n" + ext;
    writetable(tbl_1,fname,'WriteVariableNames', true);
    writetable(tbl_1n,fname_n,'WriteVariableNames', true);

%% Plotting 

% Kennlinien in getrennten Diagrammen

    fig1 = figure;
    set(0, 'DefaultLineLineWidth', 2);
    x_label = "Temperatur °C";
    y_label = "Widerstand Ω";
    
    subplot(2,2,1);
    plot(T,R_1)
    legend(devices(1))
    xlabel(x_label)
    ylabel(y_label)
    
    subplot(2,2,2)
    plot(T,R_2)
    legend(devices(2))
    xlabel(x_label)
    ylabel(y_label)
    
    subplot(2,2,3)
    plot(T,R_3)
    legend(devices(3))
    xlabel(x_label)
    ylabel(y_label)
    
    subplot(2,2,4)
    plot(T,R_4)
    legend(devices(4))
    xlabel(x_label)
    ylabel(y_label)
    sgtitle("Kennlinien in getrennten Diagrammen")

 % Kennlinien in gemeinsamem Diagramm
 
    fig2 = figure;
    plot(T,R_1)
    hold on 
    plot(T,R_2)
    plot(T,R_3)
    plot(T,R_4)
    hold off 
    legend(devices)
    title("Kennlinien in gemeinsamem Diagramm")
 % Normierten Kennlinien in gemeinsamem Diagramm

    fig3 = figure;
    plot(T,R_1n)
    hold on 
    plot(T,R_2n)
    plot(T,R_3n)
    plot(T,R_4n)
    hold off 
    legend(devices)
    title("Normierten Kennlinien in gemeinsamem Diagramm")

%% 3.2.3
R = [485.3 40158 131.89 2110.6];
T = [T_NTC(2200,3977,R(1)) T_NTC(100000,4190,R(2)),T_platin(100,3.9083*10^-3,-5.775*10^-7,-4.183*10^-12,R(3)),T_silizium(R(4))];
sol = table(devices',R',T','VariableNames',{'thermistor','R','T'});
display(sol)

%% 4.1 
R_vi        = [770 33*10^3 220 2.2*10^3];
Al_0        = [2.9355 3.1705 3.3835 3.5595];
Al_1        = [3.8965 2.8183 1.848 1.2343];
Al_2        = [3.9556 2.8383 1.8180 1.1820];
Al_3        = [1.6468 1.7374 1.8152 1.8758];
Al_4        = [2.3369 2.5634 2.7555 2.9026];
Al_7        = [5.003 5.003 5.004 5.003];

u           = Al_7;
ui          = [Al_1; Al_2; Al_3; Al_4];
R           = ones(size(ui));
for i=1:4
   for j=1:4
      R(i,j) = ui(i,j) * R_vi(i) / (u(i)-ui(i,j));
   end
end

R_ptat      = 10 * 10^3;             %ohm
k           = 1 * 10^-6;                  % Amp / K
T           = ( (Al_0 ./ (R_ptat * k)) - 273.15);
ans4_1      = table(round(T',2),round(R(1,:)',2),round(R(2,:)',2),round(R(3,:)',2),round(R(4,:)',2),'VariableNames',[{'T_ptat'}, devices]');
writetable(ans4_1,"table_41.csv",'WriteVariableNames', true);
display(ans4_1)
%% Functions' definitions
%% Characteristic curve  

function R = R_platin(R0,A,B,C,T)
%R_platin Finds the thermistor resistance corresponding to a temperature T (Platinum Resistance thermistor)
% inputs:
%   R0 : Ice point resistance of the device.
%   A  : Device parameter 
%   B  : Device parameter 
%   C  : Device parameter 
%   T  : temperatur in celsius
% outputs:
%   R  : Resistance of the thermistor in ohms 

    if  T <= 0 
        R =R0 *( 1 + A*T + B * T^2 + C*(T-100)*(T)^3);
    elseif T > 0 
        R = R0 * (1+A*T + B*T^2);
    end
end


function R = R_NTC(R_25,B,T)
%R_NTC Finds the thermistor resistance corresponding to a temperature T (Negative Temperature Coefficient)
% inputs:
%   R_25 : Reference point resistance of the device.
%   B  : Device parameter (material constant) 
%   T  : temperatur in celsius
% outputs:
%   R  : Resistance of the thermistor in ohms 

    
    Tref = 298.15;
    t = 273.15 + T;
    R = R_25 * exp(B*((1/t)-(1/Tref)));
 
end
      

function R = R_silizium(R_25,A,B,T)
%R_platin Finds the thermistor resistance corresponding to a temperature T (Infineon Silicon Temperature Sensors)
% inputs:
%   R_25 : Reference point resistance of the device.
%   A    : Device parameter 
%   B    : Device parameter 
%   T    : temperatur in celsius
% outputs:
%   R  : Resistance of the thermistor in ohms 

    R = R_25 *( 1 + A*(T-25) + B * (T-25)^2);
end     

%% (Umkehrfunktionen) Inverse Characteristic curve using simple brute force iterative approach 
% time complexity ≈ O(n) where n is range of accurace of the thermistor
% binary search with simple recursion is also possible (since all curves
% are monotonic -> which is basically a sorted array with our desired resolution) -> O(log(n))
function T  = T_platin(R_0,A,B,C,R)
%T_platin Finds the thermistor temperature T given its resistance R
%(Platinum Resistance thermistor) with a deviation of 1°C for each 1000°C
% inputs:
%   R0 : Ice point resistance of the device.
%   A  : Device parameter 
%   B  : Device parameter 
%   C  : Device parameter 
%   R  : Resistance in ohms
% outputs:
%   T  : temperature of the thermistor in °C 

    t           = 0;                                        % Start value 
    dev         = 1;                                        % deviation                                
    n           = 0;                                        % num iterations
    dt          = 1;
while dev > 0.001
    n           = n +1;
    r           = R_platin(R_0,A,B,C,t);
    dev         = abs(r-R)/R;
    
    
    % the function is linear, therefore simple inc and dec is sufficent
    if r > R
        t       = t -dt;
    else
        t       = t + dt;
    end
    
    % finer increments following lower deviations  
    if dev < 0.2
        dt = 0.1;
    elseif dev < 0.01
        dt = 0.01;
    else 
        dt = 1;
    end
    
    if n > 1000
        break
    end

end

    T = t;
end

function T = T_NTC(R_25,B,R)
%T_NTC Finds the thermistor temperature T given its resistance R
%(Negative Temperature Coefficient) with a deviation of 0.1% °C
% inputs:
%   R_25    : reference point resistance of the device.
%   B       : Device parameter 
%   R       : Resistance in ohms
% outputs:
%   T  : temperature of the thermistor in °C 
    t           = 0;
    dev         = 1;
    n           = 0;
    dt          = 1;
    while dev > 0.001
         n      = n +1;
         r      = R_NTC(R_25,B,t);
         dev    = abs(r-R)/R;
         if r > R
             t  = t +dt ;
         else
             t  = t -dt ;
         end
         if n > 1000
             break         
         end 
         % naive approach to handle the the non linear nature of the device. 
        if dev < 0.1
            dt = 0.1;
        end
    end
    T = t;
end

function T  = T_silizium(R)
%T_silizium Finds the thermistor temperature T given its resistance R
%(Infineon Silicon Temperature Sensors) with a deviation of 1°C for each 1000°C
% inputs:
%   R  : Resistance in ohms
% outputs:
%   T  : temperature of the thermistor in °C 
    R_25        = 2 * 10^3;                                 % Ohms
    A           = 7.88 * 10^-3;                             % 1/K
    B           = 1.973 * 10^-5;                            % 1/K²
    t           = 0;                                        % Start value 
    dev         = 1;                                        % deviation                                    
    n           = 0;                                        % num iterations
    dt          = 1;
while dev > 0.001
    n           = n +1;
    r           = R_silizium(R_25,A,B,t);
    dev         = abs(r-R)/R;
    
    
    % the function is linear, therefore simple inc and dec is sufficent
    if r > R
        t       = t -dt;
    else
        t       = t + dt;
    end
    
    % finer increments following lower deviations  
    if dev < 0.2
        dt = 0.1;
    elseif dev < 0.01
        dt = 0.01;
    else 
        dt = 1;
    end
    
    if n > 1000
        break
    end

end

    T = t;
end


    