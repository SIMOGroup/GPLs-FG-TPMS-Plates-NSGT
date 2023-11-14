function [gg,dgg] = Dist_Func(h,gpt,imix)

switch imix
    case 0000  % Hoang Nguyen
        gg = -gpt -8*gpt + 10*gpt^3/h^2 + 6/5.*gpt^5/h^4 + 8/7.*gpt^7/h^6;
        dgg= -1 -8 + 30*gpt^2/h^2 + 6*gpt^4/h^4 + 8*gpt^6/h^6;
        ddgg = 60*gpt^1/h^2 + 24*gpt^3/h^4 + 48*gpt^5/h^6;
        
   case 0  % Tuan Nguyen (Model 3)
        gg = gpt - 87/20*gpt^3/h^2 +   169/10*gpt^5/h^4 -   138/5*gpt^7/h^6;    %%% Tuan
        dgg= 1 - 3*87/20*gpt^2/h^2 + 5*169/10*gpt^4/h^4 - 7*138/5*gpt^6/h^6;        %%% Tuan 
        
    case 1
        gg = gpt - 4/3*gpt^3/h^2;      %%% Model Reddy
        dgg= 1 - 4*gpt^2/h^2;         %%% Model Reddy
        ddgg = -2*4*gpt^1/h^2;
   case 2
        gg = gpt - 17/10*gpt^3/h^2 + 22/25*gpt^5/h^4;    %%% Tuan
        dgg= 1 - 51/10*gpt^2/h^2 + 22/5*gpt^4/h^4;        %%% Tuan
    case 22222
        gg = sin(pi*gpt/h)-gpt;    %%% Model Arya
        dgg= pi/h*cos(pi*gpt/h)-1;        %%% Model Arya
    case 3
        gg = (h/pi)*sin(pi*gpt/h)-gpt;      %%% Model Touratier
        dgg= cos(pi*gpt/h)-1;             %%% Model Touratier
        
    case 4
        gg = gpt*exp(-2*(gpt/h)^2)-gpt;             %%% Model Karama
        dgg= (1-4*(gpt/h)^2)*exp(-2*(gpt/h)^2)-1; %%% Model Karama
        
    case 5
        gg = -h*sinh(gpt/h)+gpt*cosh(1/2)-gpt;   %%% Model Soldatos (chua dung)
        dgg= -cosh(gpt/h)+cosh(1/2)-1;       %%% Model Soldatos
        
    case 6
        alpha=3;
        somu =-2*(gpt/h)^2/log(alpha);
        gg   = gpt*alpha^somu-gpt;               %%% Model Aydogdu
        dgg  = alpha^somu*(1-4*(gpt/h)^2)-1;   %%% Model Aydogdu
        
    case 7
        gg  = -2*gpt+(h^2+4)/4*atan(gpt);          %%% Model Chien
        dgg = -2+(h^2+4)/4/(1+gpt^2);            %%% Model Chien
        
    case 8
        gg = 2*gpt-sqrt(h^2+4)*asinh(gpt);          %%% Model Chien
        dgg= 2-sqrt(h^2+4)*log(gpt+sqrt(gpt^2+1));  %%% Model Chien
        
    case 9
        gg = gpt*2.85^(-2*(gpt/h)^2)-gpt;         %%% Model Muntari
        dgg= 2.85^(-2*(gpt/h)^2)*(1-4*(gpt/h)^2*log(2.85))-1;
    case 10
        gg = h*atan(2*gpt/h)-2*gpt;              % % chien  3
        dgg= 2/(1+(2*gpt/h)^2)-2;
        
    case 11
        gg = atan(sin(pi*gpt/h));             % % chien 4
        dgg= pi/h*cos(pi*gpt/h)/(1+(sin(pi*gpt/h))^2);

        ddgg = - (pi^2*sin((pi*gpt)/h))/(h^2*(sin((pi*gpt)/h)^2 + 1)) - ...
            (2*pi^2*cos((pi*gpt)/h)^2*sin((pi*gpt)/h))/(h^2*(sin((pi*gpt)/h)^2 + 1)^2);
        
    case 12
        gg = asinh(sin(pi*gpt/h))-gpt;             % % chien 5
        dgg= pi/h*cos(pi*gpt/h)/sqrt(1+(sin(pi*gpt/h))^2)-1;
        
    case 13
        gg = -1/8*gpt-2/h^2*gpt^3+2/h^4*gpt^5;             % % hung
        dgg= -1/8-6/h^2*gpt^2 + 10/h^4*gpt^4;
        ddgg = - (12*gpt)/h^2 + (40*gpt^3)/h^4 ;
        
    case 14 %Hoang Nguyen
        gg = -gpt -8*gpt + 10*gpt^3/h^2 + 6/5.*gpt^5/h^4 + 8/7.*gpt^7/h^6;
        dgg= -1 -8 + 30*gpt^2/h^2 + 6*gpt^4/h^4 + 8*gpt^6/h^6;
        ddgg = 60*gpt^1/h^2 + 24*gpt^3/h^4 + 48*gpt^5/h^6;
        
    case 15 % %
        gg =  1/4*gpt - 5/3*gpt^3/h^2;
        dgg=  1/4 - 3*5/3*gpt^2/h^2;
        ddgg = - 2*3*5/3*gpt^1/h^2;
        
    case 16     %Zenkour 2013
        gg = - gpt + h*sinh(gpt/h)-4*gpt^3/(3*h^2)*cosh(1/2);
        dgg = - 1 + cosh(gpt/h) - (317398492294797*gpt^2)/(70368744177664*h^2);
        ddgg = sinh(gpt/h)/h - (317398492294797*gpt)/(35184372088832*h^2);
        %
end