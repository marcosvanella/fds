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

basedir='/Users/mnv/Documents/FIREMODELS_FORK/fds/Validation/NIST_Backdraft/FDS_Input_Files/'

% First Cases with varying temperature:
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
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT210(z);
            case 240
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT240(z);
            case 270
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT270(z);
            case 300
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT300(z);
        end
        ! Write fid_Tmin line for this z:
        initL_Tmin = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                      'VOLUME_FRACTION(1)=' num2str(X_F)   ', SPEC_ID(1)=''PROPANE'''          ', ' ...
                      'VOLUME_FRACTION(2)=' num2str(X_O2)  ', SPEC_ID(2)=''OXYGEN COMP'''           ', ' ...
                      'VOLUME_FRACTION(3)=' num2str(X_CO2) ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                      'VOLUME_FRACTION(4)=' num2str(X_CO)  ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                      'VOLUME_FRACTION(5)=' num2str(X_H2O) ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                      'TEMPERATURE=' num2str(TMIN) '/'];
        fprintf(fid_Tmin,'%s\n',initL_Tmin);

        ! Write fid_Tmax line for this z:
        initL_Tmax = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                      'VOLUME_FRACTION(1)=' num2str(X_F)   ', SPEC_ID(1)=''PROPANE'''          ', ' ...
                      'VOLUME_FRACTION(2)=' num2str(X_O2)  ', SPEC_ID(2)=''OXYGEN COMP'''           ', ' ...
                      'VOLUME_FRACTION(3)=' num2str(X_CO2) ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                      'VOLUME_FRACTION(4)=' num2str(X_CO)  ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                      'VOLUME_FRACTION(5)=' num2str(X_H2O) ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                      'TEMPERATURE=' num2str(TMAX) '/'];
        fprintf(fid_Tmax,'%s\n',initL_Tmax);

        ! Write fid_Tvar line for this z:
        initL_Tvar = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                      'VOLUME_FRACTION(1)=' num2str(X_F)   ', SPEC_ID(1)=''PROPANE'''          ', ' ...
                      'VOLUME_FRACTION(2)=' num2str(X_O2)  ', SPEC_ID(2)=''OXYGEN COMP'''           ', ' ...
                      'VOLUME_FRACTION(3)=' num2str(X_CO2) ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                      'VOLUME_FRACTION(4)=' num2str(X_CO)  ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                      'VOLUME_FRACTION(5)=' num2str(X_H2O) ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                      'TEMPERATURE=' num2str(T) '/'];
        fprintf(fid_Tvar,'%s\n',initL_Tvar);
        X_FV(iz)  = X_F;
        X_O2V(iz) = X_O2;
        X_CO2V(iz)= X_CO2;
        X_COV(iz) = X_CO;
        X_H2OV(iz)= X_H2O;
        TMPV(iz)  = T;
    end
    fprintf(fid_Tmin,'%s\n',"! Initialize the whole domain to dry air:");
    fprintf(fid_Tmin,'%s\n',"&INIT XB=-1.54,3.03,-1.11,1.11,-1.00,3.84, MASS_FRACTION(1)=0.23, SPEC_ID(1)='OXYGEN' /");
    fprintf(fid_Tmax,'%s\n',"! Initialize the whole domain to dry air:");
    fprintf(fid_Tmax,'%s\n',"&INIT XB=-1.54,3.03,-1.11,1.11,-1.00,3.84, MASS_FRACTION(1)=0.23, SPEC_ID(1)='OXYGEN' /");
    fprintf(fid_Tvar,'%s\n',"! Initialize the whole domain to dry air:");
    fprintf(fid_Tvar,'%s\n',"&INIT XB=-1.54,3.03,-1.11,1.11,-1.00,3.84, MASS_FRACTION(1)=0.23, SPEC_ID(1)='OXYGEN' /");
    fclose(fid_Tmin); fclose(fid_Tmax); fclose(fid_Tvar);
    figure
   
    subplot(1,2,1)
    hold on
    plot(X_FV,ZVEC,'--k')
    plot(X_O2V,ZVEC,'b')
    plot(X_CO2V,ZVEC,'g')
    plot(X_COV,ZVEC,'r')
    plot(X_H2OV,ZVEC,'m')
    xlabel('Vol Fraction X')
    ylabel('Height [m]')
    %title(['Propane 25kW, FFT=' num2str(FFT(ift)) 'sec.'])
    legend('C_3H_8','O_2','CO_2','CO','H_2O')
    grid on 
    box on
    subplot(1,2,2)
    plot(TMPV,ZVEC,'k')
    xlabel('Temperature [^oC]')
    %title(['Propane 25kW, FFT=' num2str(FFT(ift)) 'sec.'])
    grid on 
    box on
    pause

end 


% Cases with min and max X_F for 300sec FFT:
FFT =[300];
TVEC={'TVAR'};
for ift=1:length(FFT)
    filename_FMIN =['Propane25kW_' num2str(FFT(ift)) 's_Fmin_init.fds'];
    filename_FMEAN=['Propane25kW_' num2str(FFT(ift)) 's_Fmean_init.fds'];
    filename_FMAX =['Propane25kW_' num2str(FFT(ift)) 's_Fmax_init.fds'];
    [fid_Fmin]=fopen([basedir filename_FMIN],'w');
    fprintf(fid_Fmin,'%s\n',['# Initialization file for Propane 25kW, FFT=' num2str(FFT(ift)) 's, minimum X_F across compartment height.']);
    [fid_Fmean]=fopen([basedir filename_FMEAN],'w');
    fprintf(fid_Fmean,'%s\n',['# Initialization file for Propane 25kW, FFT=' num2str(FFT(ift)) 's, mean X_F across compartment height.']);
    [fid_Fmax]=fopen([basedir filename_FMAX],'w');
    fprintf(fid_Fmax,'%s\n',['# Initialization file for Propane 25kW, FFT=' num2str(FFT(ift)) 's, maximum X_F across compartment height.']);
    for iz=1:length(ZVEC)
        z = ZVEC(iz);
        ! Get temperatures and composition for &INIT XB:
        switch FFT(ift)
            case 210
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT210(z);
            case 240
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT240(z);
            case 270
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT270(z);
            case 300
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT300(z);
        end
        ! Write fid_Fmin line for this z:
        initL_Fmin = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                      'VOLUME_FRACTION(1)=' num2str(X_FMIN) ', SPEC_ID(1)=''PROPANE'''         ', ' ...
                      'VOLUME_FRACTION(2)=' num2str(X_O2)  ', SPEC_ID(2)=''OXYGEN COMP'''      ', ' ...
                      'VOLUME_FRACTION(3)=' num2str(X_CO2) ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                      'VOLUME_FRACTION(4)=' num2str(X_CO)  ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                      'VOLUME_FRACTION(5)=' num2str(X_H2O) ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                      'TEMPERATURE=' num2str(T) '/'];
        fprintf(fid_Fmin,'%s\n',initL_Fmin);

        ! Write fid_Fmean line for this z:
        initL_Fmean = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                       'VOLUME_FRACTION(1)=' num2str(0.5*(X_FMIN+X_FMAX)) ', SPEC_ID(1)=''PROPANE'''         ', ' ...
                       'VOLUME_FRACTION(2)=' num2str(X_O2)  ', SPEC_ID(2)=''OXYGEN COMP'''      ', ' ...
                       'VOLUME_FRACTION(3)=' num2str(X_CO2) ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                       'VOLUME_FRACTION(4)=' num2str(X_CO)  ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                       'VOLUME_FRACTION(5)=' num2str(X_H2O) ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                       'TEMPERATURE=' num2str(T) '/'];
        fprintf(fid_Fmean,'%s\n',initL_Fmean);


        ! Write fid_Fmax line for this z:
        initL_Fmax = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                      'VOLUME_FRACTION(1)=' num2str(X_FMAX) ', SPEC_ID(1)=''PROPANE'''         ', ' ...
                      'VOLUME_FRACTION(2)=' num2str(X_O2)  ', SPEC_ID(2)=''OXYGEN COMP'''      ', ' ...
                      'VOLUME_FRACTION(3)=' num2str(X_CO2) ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                      'VOLUME_FRACTION(4)=' num2str(X_CO)  ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                      'VOLUME_FRACTION(5)=' num2str(X_H2O) ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                      'TEMPERATURE=' num2str(T) '/'];
        fprintf(fid_Fmax,'%s\n',initL_Fmax);

        X_O2V(iz) = X_O2;
        X_CO2V(iz)= X_CO2;
        X_COV(iz) = X_CO;
        TMPV(iz)  = T;
    end
    fprintf(fid_Fmin,'%s\n',"! Initialize the whole domain to dry air:");
    fprintf(fid_Fmin,'%s\n',"&INIT XB=-1.54,3.03,-1.11,1.11,-1.00,3.84, MASS_FRACTION(1)=0.23, SPEC_ID(1)='OXYGEN' /");
    fprintf(fid_Fmean,'%s\n',"! Initialize the whole domain to dry air:");
    fprintf(fid_Fmean,'%s\n',"&INIT XB=-1.54,3.03,-1.11,1.11,-1.00,3.84, MASS_FRACTION(1)=0.23, SPEC_ID(1)='OXYGEN' /");
    fprintf(fid_Fmax,'%s\n',"! Initialize the whole domain to dry air:");
    fprintf(fid_Fmax,'%s\n',"&INIT XB=-1.54,3.03,-1.11,1.11,-1.00,3.84, MASS_FRACTION(1)=0.23, SPEC_ID(1)='OXYGEN' /");
    fclose(fid_Fmin); fclose(fid_Fmean); fclose(fid_Fmax);

    figure
    subplot(1,2,1)
    hold on
    plot(X_FMIN*ones(1,length(ZVEC)),ZVEC,'--k')
    plot(0.5*(X_FMIN+X_FMAX)*ones(1,length(ZVEC)),ZVEC,'-xk')
    plot(X_FMAX*ones(1,length(ZVEC)),ZVEC,'k')
    plot(X_O2V,ZVEC,'b')
    plot(X_CO2V,ZVEC,'g')
    plot(X_COV,ZVEC,'r')
    xlabel('Vol Fraction X')
    ylabel('Height [m]')
    title(['Propane 25kW, FFT=' num2str(FFT(ift)) 'sec.'])
    legend('Propane FMIN','Propane FMEAN','Propane FMAX','O_2','CO_2','CO')
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


% Cases with min and max X_O2 for 300sec FFT:
for ift=1:length(FFT)
    filename_O2MIN =['Propane25kW_' num2str(FFT(ift)) 's_O2min_init.fds'];
    filename_O2MEAN=['Propane25kW_' num2str(FFT(ift)) 's_O2mean_init.fds'];
    filename_O2MAX =['Propane25kW_' num2str(FFT(ift)) 's_O2max_init.fds'];
    [fid_O2min]=fopen([basedir filename_O2MIN],'w');
    fprintf(fid_O2min,'%s\n',['# Initialization file for Propane 25kW, FFT=' num2str(FFT(ift)) 's, minimum X_O2 across compartment height.']);
    [fid_O2mean]=fopen([basedir filename_O2MEAN],'w');
    fprintf(fid_O2mean,'%s\n',['# Initialization file for Propane 25kW, FFT=' num2str(FFT(ift)) 's, mean X_O2 across compartment height.']);
    [fid_O2max]=fopen([basedir filename_O2MAX],'w');
    fprintf(fid_O2max,'%s\n',['# Initialization file for Propane 25kW, FFT=' num2str(FFT(ift)) 's, maximum X_O2 across compartment height.']);
    for iz=1:length(ZVEC)
        z = ZVEC(iz);
        ! Get temperatures and composition for &INIT XB:
        switch FFT(ift)
            case 210
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT210(z);
            case 240
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT240(z);
            case 270
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT270(z);
            case 300
                [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT300(z);
        end
        ! Write fid_O2min line for this z:
        initL_O2min = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                      'VOLUME_FRACTION(1)=' num2str(X_F)     ', SPEC_ID(1)=''PROPANE'''         ', ' ...
                      'VOLUME_FRACTION(2)=' num2str(X_O2MIN) ', SPEC_ID(2)=''OXYGEN COMP'''      ', ' ...
                      'VOLUME_FRACTION(3)=' num2str(X_CO2)   ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                      'VOLUME_FRACTION(4)=' num2str(X_CO)    ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                      'VOLUME_FRACTION(5)=' num2str(X_H2O)   ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                      'TEMPERATURE=' num2str(T) '/'];
        fprintf(fid_O2min,'%s\n',initL_O2min);

        ! Write fid_O2mean line for this z:
        initL_O2mean = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                        'VOLUME_FRACTION(1)=' num2str(X_F)     ', SPEC_ID(1)=''PROPANE'''         ', ' ...
                        'VOLUME_FRACTION(2)=' num2str(0.5*(X_O2MIN+X_O2MAX)) ', SPEC_ID(2)=''OXYGEN COMP'''      ', ' ...
                        'VOLUME_FRACTION(3)=' num2str(X_CO2)   ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                        'VOLUME_FRACTION(4)=' num2str(X_CO)    ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                        'VOLUME_FRACTION(5)=' num2str(X_H2O)   ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                        'TEMPERATURE=' num2str(T) '/'];
        fprintf(fid_O2mean,'%s\n',initL_O2mean);

        ! Write fid_O2max line for this z:
        initL_O2max = ['&INIT XB=-1.52,0.,-0.495,0.495,' num2str(z-0.005) ',' num2str(z+0.005)  ', ' ...
                      'VOLUME_FRACTION(1)=' num2str(X_F)     ', SPEC_ID(1)=''PROPANE'''         ', ' ...
                      'VOLUME_FRACTION(2)=' num2str(X_O2MAX) ', SPEC_ID(2)=''OXYGEN COMP'''      ', ' ...
                      'VOLUME_FRACTION(3)=' num2str(X_CO2)   ', SPEC_ID(3)=''CARBON DIOXIDE'''   ', ' ...
                      'VOLUME_FRACTION(4)=' num2str(X_CO)    ', SPEC_ID(4)=''CARBON MONOXIDE'''  ', ' ...
                      'VOLUME_FRACTION(5)=' num2str(X_H2O)   ', SPEC_ID(5)=''WATER VAPOR'''      ', ' ...
                      'TEMPERATURE=' num2str(T) '/'];
        fprintf(fid_O2max,'%s\n',initL_O2max);

        X_FV(iz) = X_F;
        X_CO2V(iz)= X_CO2;
        X_COV(iz) = X_CO;
        TMPV(iz)  = T;
    end
    fprintf(fid_O2min,'%s\n',"! Initialize the whole domain to dry air:");
    fprintf(fid_O2min,'%s\n',"&INIT XB=-1.54,3.03,-1.11,1.11,-1.00,3.84, MASS_FRACTION(1)=0.23, SPEC_ID(1)='OXYGEN' /");
    fprintf(fid_O2mean,'%s\n',"! Initialize the whole domain to dry air:");
    fprintf(fid_O2mean,'%s\n',"&INIT XB=-1.54,3.03,-1.11,1.11,-1.00,3.84, MASS_FRACTION(1)=0.23, SPEC_ID(1)='OXYGEN' /");
    fprintf(fid_O2max,'%s\n',"! Initialize the whole domain to dry air:");
    fprintf(fid_O2max,'%s\n',"&INIT XB=-1.54,3.03,-1.11,1.11,-1.00,3.84, MASS_FRACTION(1)=0.23, SPEC_ID(1)='OXYGEN' /");
    fclose(fid_O2min); fclose(fid_O2mean); fclose(fid_O2max);

    figure
    subplot(1,2,1)
    hold on
    plot(X_O2MIN*ones(1,length(ZVEC)),ZVEC,'--k')
    plot(X_O2MAX*ones(1,length(ZVEC)),ZVEC,'k')
    plot(X_FV,ZVEC,'b')
    plot(X_CO2V,ZVEC,'g')
    plot(X_COV,ZVEC,'r')
    xlabel('Vol Fraction X')
    ylabel('Height [m]')
    title(['Propane 25kW, FFT=' num2str(FFT(ift)) 'sec.'])
    legend('O2 MIN','O2 MAX','Propane','CO_2','CO')
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

function [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT210(y)
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
X_FMIN = 0.00007*0+0.0416;
X_FMAX = 0.00007*100+0.0416;
X_O2MIN=-0.0006*100+0.0991;
X_O2MAX=-0.0006*0+0.0991;
TMIN=(0.9031*0+445.47)-273.15; % In deg C at the floor   z=0 m.
TMAX=(0.9031*100+445.47)-273.15; % In deg C at the Ceiling z=1 m.
end

function [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT240(y)          
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
X_FMIN = 0.00006*0+0.0462;
X_FMAX = 0.00006*100+0.0462;
X_O2MIN=-0.0007*100+0.1114;
X_O2MAX=-0.0007*0+0.1114;
TMIN=(0.7744*0+436.93)-273.15; % In deg C at the floor   z=0 m.
TMAX=(0.7744*100+436.93)-273.15; % In deg C at the Ceiling z=1 m.
end

function [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT270(y)          
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
X_FMAX = -0.0002*0+0.0662;
X_FMIN = -0.0002*100+0.0662;
X_O2MIN= -0.0012*100+0.1538;
X_O2MAX= -0.0012*0+0.1538;
TMIN=(0.4743*0+420.85)-273.15; % In deg C at the floor   z=0 m.
TMAX=(0.4743*100+420.85)-273.15; % In deg C at the Ceiling z=1 m.
end

function [X_F,X_O2,X_CO2,X_CO,X_H2O,X_N2,T,TMIN,TMAX,X_FMIN,X_FMAX,X_O2MIN,X_O2MAX]=get_TComp_FFT300(y)          
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
X_FMAX = -0.0003*0+0.0618;
X_FMIN = -0.0003*100+0.0618;
X_O2MIN= -0.0009*100+0.1328;
X_O2MAX= -0.0009*0+0.1328;
TMIN=(0.4279*0+425.08)-273.15; % In deg C at the floor   z=0 m.
TMAX=(0.4279*100+425.08)-273.15; % In deg C at the Ceiling z=1 m.
end