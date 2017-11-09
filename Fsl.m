function [F,R4] = Fsl(phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Short-wave radiation distribution scheme to estimate potential solar
% radiation on mountain hillslopes. The algorithm used is that demonstrated 
% by Swift (1976), based on the equations of Lee (1963) and Frank and 
% Lee (1966). See README File for complete variable descriptions.
% Reggie Walters, February 2011, Boise State Unversity
% Output 'F' is ratio of potential solar radiation on hillslope to that of
% an identical location with a 0deg slope.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% START %
%%%%%%%%%

%%% INPUTS: requires scalar inputs, use a looping driver function to
%%% evaluate over a larger mesh

% phi is vector of [L0 I A J Lon R2]

L0 = phi(1);                  % Latitude, N is +, S is -, [degrees]
I = phi(2);                   % Slope inclination angle [degrees]
A = phi(3);                   % Slope azimuth angle, CW from N(0) [degrees]
J = phi(4);                   % Julian Day (*since Jan 1) [day]
Lon = phi(5);                 % Longitude, West of Greenwich is -                   
R2 = phi(6);                  % *Measured solar radiation [cal/cm^2/day]


R0 = 1.962;                    % Solar flux Constant [cal/cm^2/min] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isfinite(I)&&isfinite(A))   % Operate PSR on finite values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Latitude of equivalent slope (see Lee (1963))
L1 = asind( cosd(I)*sind(L0) + sind(I)*cosd(L0)*cosd(A) );

% Time offset in hour angle between actual and equivalent slopes 
L2 = gOffset(I,L0,A);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% BEGIN MAIN ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[D,E] = Func1(J);
R1 = 60*R0/(E*E);
T = Func2(L1,D);   
T7 = T - L2;
T6 = -T - L2;
T = Func2(L0,D);
T1 = T;
T0 = -T;

if(T7 < T1)
    T3 = T7;
else
    T3 = T1;
end
if(T6 > T0)
    T2 = T6;
else
    T2 = T0;
end

%%% R4  Step, alternative routine from Swift 1974 %%%%%%%%%%%%%%%%%%%%%%%%%
if(T3 < T2)
    T2 = 0;
    T3 = 0;
end
T6 = T6 + 360;
if(T6 < T1)
    T8 = T6;
    T9 = T1;
    R4 = Func3(L2,L1,T3,T2,R1,D) + Func3(L2,L1,T9,T6,R1,D);
else
    T7 = T7 - 360;
    if(T7 > T0)
        T8 = T0;
        T9 = T7;
        R4 = Func3(L2,L1,T3,T2,R1,D) + Func3(L2,L1,T9,T6,R1,D);
    else
        R4 = Func3(L2,L1,T3,T2,R1,D);
    end
end
%%% End R4 Alternative Routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T4 = T2 / 15;
T5 = T3 / 15;
R3 = Func3(0.0,L0,T1,T0,R1,D);
F = R4/R3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    D  = NaN;
    F  = NaN;
    R3 = NaN;
    R4 = NaN;
    T  = NaN;
    T2 = NaN;
    T3 = NaN;
    T4 = NaN;
    T5 = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEFINE SUBROUTINE FUNCTIONS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [D,E] = Func1(J)
        % Func1 estimates declination angle and radius vector
        % Sun Declination
        D = asin(0.39785*sin(4.868961 + 0.017203*J+0.033446*sin(6.224111 + 0.017202*J)));
        D = D*180/pi; % convert to degrees
        %Z = W - X*cosd( (J+Y)*0.986 );
        % Radius Vector
        E = 1.0 - 0.0167*cos((J-3)*0.0172);
        %NOTE: these equations were given in radians, so the trig function usage is correct
    end

    function T = Func2(Y,D)
        % Func2 begins calculation of sunrise/sunset times
        T = acosd( -tand(Y)*tand(D) );
        T = real(T);
    end

    function Z = Func3(V,W,X,Y,R1,D)
        % Func3 gives potential irradiation for a horizontal
        Z = R1*( sind(D)*sind(W)*(X-Y)/15 + cosd(D)*cosd(W)*(sind(X+V) - sind(Y+V))*12/pi);
    end

    function L2 = gOffset(I,L0,A)
        % gOffset gets the hour angle offset between actual and equivalent slopes
        D1 = cosd(I)*cosd(L0) - sind(I)*sind(L0)*cosd(A);
        if(D1 == 0)
            D1 = 1.0e-10;
        end
        L2 = atand( sind(I)*sind(A)/D1 );   % This is the time offset
        if (D1 < 0)
            L2 = L2 + 180;                  % ""  "" if D1 is negative
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end             % End of Function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EOF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



