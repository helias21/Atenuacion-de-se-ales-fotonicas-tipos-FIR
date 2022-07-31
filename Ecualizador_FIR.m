% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%   Calcula los valores de un filtro FIR óptimo para ecualizar un EDFA
% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

clear all
clc
L=20;                   %Se asume 20Kms de recorrido de fibra.  Este valor
                        % pueda cambiar de acuerdo a la distancia
carga=load('ASE.txt');  %Se importan lo valores de potencia sin atenuación
                        %todas las frecuencias desde 1510 a 1580nm

lamb=carga(251:324,1);  %se cargan valores en los rangos de interes
pot=carga(251:324,2);   %1545-1555nm
potdB=10*log10(pot);    %Potencia en dB
fre=1e9*299792458./lamb;%conversión de long de onda a frecuencia

%figure; plot (lamb, pot, 'm')   %grafica potencia vs. long de onda
%figure; plot (lamb, potdB, 'b') %grafica potencia en dB vs. long
                                %de onda

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Atenución de con aproximación parabolica  Y=ax2+bx+c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphaDBKm=2.15e-6.*lamb.^2-6.665e-3.*lamb+5.355375;
coefatt=alphaDBKm./4.343 ;  %ecuacion 3.7 de Le Nguyen Binh pagina 59
potatt=pot.*exp(coefatt.*(-L)); %ecuacion 3.6 de Le Nguyen Binh pagina 59

figure; plot(lamb,potatt, lamb, pot) %ver potencia inicial y potencia atenuada

ecu=(potatt./min(potatt)).^-1;
%figure; plot (lamb, ecu, 'r')



figure; plot(lamb,(potatt./min(potatt)), lamb, ecu) %ver como llega y la ecualización a la vez
ecu1=flipud(ecu);   %Amplificador normalizado
ecu1dB=10*log10(ecu1);
%figure; plot (lamb, ecu1, 'g')
figure; plot (lamb, ecu1dB, 'g')
frecu=flipud(fre);            %Inversamente proporcional a la longitud de onda

figure; plot (frecu, ecu1, 'g')
frecu=(frecu - frecu(1));
frecu=frecu./max(frecu);

%Apagado de Canales
% Hay 74 puntos.   6 canales tendrán 7 puntos y 4 canales tendrá 8
%ecu(1-7) canal 1, ecu(8-14) canal 2, ecu(1-21) canal 3...
%ecu(22-29) canal 4, ecu(30-37) canal 5, ecu(38-45) canal 6
%ecu(46-53)  canal 7,  ecu(54-60)  canal 8, ecu(61-67) canal 9
%ecu(68-74)  canal 10
canal=[1 1 0 1 1 1 1 1 1 1];  %Se apagan canales 3 y 8
for i=1:74
    if i<=7
        ecu1(i)=ecu1(i)*canal(1);
    elseif i>7 & i<=14
        ecu1(i)=ecu1(i)*canal(2);
    elseif i>14 & i<=21
        ecu1(i)=ecu1(i)*canal(3);
    elseif i>21 & i<=29
        ecu1(i)=ecu1(i)*canal(4);
    elseif i>29 & i<=37
        ecu1(i)=ecu1(i)*canal(5);
    elseif i>37 & i<=45
        ecu1(i)=ecu1(i)*canal(6);
    elseif i>45 & i<=53
        ecu1(i)=ecu1(i)*canal(7);
    elseif i>53 & i<=60
        ecu1(i)=ecu1(i)*canal(8);
    elseif i>60 & i<=67
        ecu1(i)=ecu1(i)*canal(9);    
    else i>67 & i<=74
        ecu1(i)=ecu1(i)*canal(10);
    end    
end
a=ecu1.^0.5;
figure; plot(frecu,20*log10(a),'k');   %ASE


b=firls(30,frecu,a);
den=[1 zeros(1,20)];
[H,F]=freqz(b,den,74);

%hold
plot(F./pi,20*log10(abs(H)), 'y')   %Resp a la frec





