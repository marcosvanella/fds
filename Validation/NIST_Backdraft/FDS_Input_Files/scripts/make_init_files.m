% Routine to generate initial zone .fds file to test changes in Temp and
% composition in flamability at Scaled compartment spark locations.
% Propane 25kW, FFt 210s, 240s, 270s, 300s.
% Linearly fitted curves in height from Ryan Falkentein-Smith:
%
% FFT   X_F_Function      X_O2_Function      X_CO2_Function      X_CO_Function      X_H2O_Function      X_N2_Function      T (K)
% 210   0.00007*z+0.0416  -0.0006*z+0.0991   0.0003*z+0.0517     0.00005*z+0.0058   0.0006*z+0.0835     -0.0002*z+0.7247   0.9031*z+445.47
% 240   0.00006*z+0.0462  -0.0007*z+0.1114   0.0003*z+0.0463     0.00009*z+0.0035   0.0007*z+0.0717     -0.0003*z+0.7264   0.7744*z+436.93
% 270   -0.0002*z+0.0662  -0.0012*z+0.1538   0.0006*z+0.0224     0.0001*z+0.0013    0.0011*z+0.0313     -0.0001*z+0.7258   0.4743*z+420.85
% 300   -0.0003*z+0.0618  -0.0009*z+0.1328   0.0005*z+0.0305     0.00007*z+0.0027   0.0009*z+0.0454     0.0003*z+0.7044    0.4279*z+425.08

close all
clear all
clc

basedir='/Users/mnv/Documents/FIREMODELS_FORK/fds/Validation/NIST_Backdraft/FDS_Input_Files/scripts/'

FFT =[210 240 270 300];
TVEC={'TVAR','TMIN','TMAX'};
ZVEC=linspace(0.005,0.995,100);

for ift=1:length(FFT)
    filename_TMIN=['Propane25kW_' num2str(FFT(ift)) 's_Tmin_init.fds'];
    filename_TMAX=['Propane25kW_' num2str(FFT(ift)) 's_Tmax_init.fds'];
    filename_TVAR=['Propane25kW_' num2str(FFT(ift)) 's_Tvar_init.fds'];
    [fid_Tmin]=fopen([basedir filename_TMIN],'w');
    fprintf(fid_Tmin,'%s\n',['# Initialization file for Propane 25kW, FFT=' num2str(FFT(ift)) 's, minimum temperature across compartment height.']);
    [fid_Tmax]=fopen([basedir filename_TMAX],'w');
    fprintf(fid_Tmax,'%s\n',['# Initialization file for Propane 25kW, FFT=' num2str(FFT(ift)) 's, maximum temperature across compartment height.']);
    [fid_Tvar]=fopen([basedir filename_TVAR],'w');
    fprintf(fid_Tvar,'%s\n',['# Initialization file for Propane 25kW, FFT=' num2str(FFT(ift)) 's, varying temperature across compartment height.']);

    for iz=1:length(ZVEC)
        z = ZVEC(iz);
        ! Get temperatures and composition for &INIT XB:
        switch FFT(ift)
            case 210
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX]=get_TComp_FFT210(z);
            case 240
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX]=get_TComp_FFT240(z);
            case 270
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX]=get_TComp_FFT270(z);
            case 300
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX]=get_TComp_FFT300(z);
        end
        ! Write fid_Tmin line for this z:
        initL_Tmin = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                      'VOLUME_FRACTION(1)=' num2str(X_F)   ', SPEC_ID(1)=''PROPANE'''          ', ' ...
                      'VOLUME_FRACTION(2)=' num2str(X_O2)  ', SPEC_ID(2)=''OXYGEN'''           ', ' ...
                      'VOLUME_FRACTION(3)=' num2str(X_CO2) ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                      'VOLUME_FRACTION(4)=' num2str(X_CO)  ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                      'VOLUME_FRACTION(5)=' num2str(X_H2O) ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                      'TEMPERATURE=' num2str(TMIN) '/'];
        fprintf(fid_Tmin,'%s\n',initL_Tmin);

        ! Write fid_Tmax line for this z:
        initL_Tmax = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                      'VOLUME_FRACTION(1)=' num2str(X_F)   ', SPEC_ID(1)=''PROPANE'''          ', ' ...
                      'VOLUME_FRACTION(2)=' num2str(X_O2)  ', SPEC_ID(2)=''OXYGEN'''           ', ' ...
                      'VOLUME_FRACTION(3)=' num2str(X_CO2) ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                      'VOLUME_FRACTION(4)=' num2str(X_CO)  ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                      'VOLUME_FRACTION(5)=' num2str(X_H2O) ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                      'TEMPERATURE=' num2str(TMAX) '/'];
        fprintf(fid_Tmax,'%s\n',initL_Tmax);

        ! Write fid_Tvar line for this z:
        initL_Tvar = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                      'VOLUME_FRACTION(1)=' num2str(X_F)   ', SPEC_ID(1)=''PROPANE'''          ', ' ...
                      'VOLUME_FRACTION(2)=' num2str(X_O2)  ', SPEC_ID(2)=''OXYGEN'''           ', ' ...
                      'VOLUME_FRACTION(3)=' num2str(X_CO2) ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                      'VOLUME_FRACTION(4)=' num2str(X_CO)  ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                      'VOLUME_FRACTION(5)=' num2str(X_H2O) ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                      'TEMPERATURE=' num2str(T) '/'];
        fprintf(fid_Tvar,'%s\n',initL_Tvar);
        X_FV(iz)  = X_F;
        X_O2V(iz) = X_O2;
        X_CO2V(iz)= X_CO2;
        X_COV(iz) = X_CO;
        TMPV(iz)  = T;
    end
    fclose(fid_Tmin); fclose(fid_Tmax); fclose(fid_Tvar);
    figure
   
    subplot(1,2,1)
    hold on
    plot(X_FV,ZVEC,'--k')
    plot(X_O2V,ZVEC,'b')
    plot(X_CO2V,ZVEC,'g')
    plot(X_COV,ZVEC,'r')
    xlabel('Vol Fraction X')
    ylabel('Height [m]')
    title(['Propane 25kW, FFT=' num2str(FFT(ift)) 'sec.'])
    legend('Propane','O_2','CO_2','CO')
    grid on 
    box on
    subplot(1,2,2)
    plot(TMPV,ZVEC,'k')
    xlabel('Temperature [C]')
    ylabel('Height [m]')
    title(['Propane 25kW, FFT=' num2str(FFT(ift)) 'sec.'])
    grid on 
    box on
    pause

end

return

%% Composition and Temperature definition functions:

function [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX]=get_TComp_FFT210(y)
% FFT   X_F_Function      X_O2_Function      X_CO2_Function      X_CO_Function      X_H2O_Function      X_N2_Function      T (K)
% 210   0.00007*z+0.0416  -0.0006*z+0.0991   0.0003*z+0.0517     0.00005*z+0.0058   0.0006*z+0.0835     -0.0002*z+0.7247   0.9031*z+445.47

z=100*y; % from m to cm. 

X_F = 0.00007*z+0.0416;
X_O2=-0.0006*z+0.0991;
X_CO2=0.0003*z+0.0517;
X_CO=0.00005*z+0.0058;
X_H2O=0.0006*z+0.0835;
X_N2=-0.0002*z+0.7247;
T   =(0.9031*z+445.47)-273.15; % In deg C.
TMIN=(0.9031*0+445.47)-273.15; % In deg C at the floor   z=0 m.
TMAX=(0.9031*1+445.47)-273.15; % In deg C at the Ceiling z=1 m.
end

function [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX]=get_TComp_FFT240(y)          
% FFT   X_F_Function      X_O2_Function      X_CO2_Function      X_CO_Function      X_H2O_Function      X_N2_Function      T (K)
% 240   0.00006*z+0.0462  -0.0007*z+0.1114   0.0003*z+0.0463     0.00009*z+0.0035   0.0007*z+0.0717     -0.0003*z+0.7264   0.7744*z+436.93

z=100*y; % from m to cm. 

X_F = 0.00006*z+0.0462;
X_O2=-0.0007*z+0.1114;
X_CO2=0.0003*z+0.0463;
X_CO=0.00009*z+0.0035;
X_H2O=0.0007*z+0.0717;
X_N2=-0.0003*z+0.7264;
T   =(0.7744*z+436.93)-273.15; % In deg C.
TMIN=(0.7744*0+436.93)-273.15; % In deg C at the floor   z=0 m.
TMAX=(0.7744*1+436.93)-273.15; % In deg C at the Ceiling z=1 m.
end

function [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX]=get_TComp_FFT270(y)          
% FFT   X_F_Function      X_O2_Function      X_CO2_Function      X_CO_Function      X_H2O_Function      X_N2_Function      T (K)
% 270   -0.0002*z+0.0662  -0.0012*z+0.1538   0.0006*z+0.0224     0.0001*z+0.0013    0.0011*z+0.0313     -0.0001*z+0.7258   0.4743*z+420.85

z=100*y; % from m to cm. 

X_F = -0.0002*z+0.0662;
X_O2=-0.0012*z+0.1538;
X_CO2=0.0006*z+0.0224;
X_CO=0.0001*z+0.0013;
X_H2O=0.0011*z+0.0313;
X_N2=-0.0001*z+0.7258;
T   =(0.4743*z+420.85)-273.15; % In deg C.
TMIN=(0.4743*0+420.85)-273.15; % In deg C at the floor   z=0 m.
TMAX=(0.4743*1+420.85)-273.15; % In deg C at the Ceiling z=1 m.
end

function [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX]=get_TComp_FFT300(y)          
% FFT   X_F_Function      X_O2_Function      X_CO2_Function      X_CO_Function      X_H2O_Function      X_N2_Function      T (K)
% 300   -0.0003*z+0.0618  -0.0009*z+0.1328   0.0005*z+0.0305     0.00007*z+0.0027   0.0009*z+0.0454     0.0003*z+0.7044    0.4279*z+425.08

z=100*y; % from m to cm. 

X_F = -0.0003*z+0.0618;
X_O2=-0.0009*z+0.1328;
X_CO2=0.0005*z+0.0305;
X_CO=0.00007*z+0.0027;
X_H2O=0.0009*z+0.0454;
X_N2=0.0003*z+0.7044;
T   =(0.4279*z+425.08)-273.15; % In deg C.
TMIN=(0.4279*0+425.08)-273.15; % In deg C at the floor   z=0 m.
TMAX=(0.4279*1+425.08)-273.15; % In deg C at the Ceiling z=1 m.
end