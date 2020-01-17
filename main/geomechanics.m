%%CO2 Injection Monitoring Geomechanics Modelling
%%Yohanes Nuwara
%%Institut Teknologi Bandung
 
clear,clc;
 
%%Read injection data
data_injeksi=struct('date',[],'cumvol',[],'pp_injected',[]);
fid=fopen('InjectionData_800tpd.dat', 'r');
i=0;
porosity=0.14; %input porosity
bulkvol=4.8*10^11; %input reservoir bulk volume
porevol=porosity*bulkvol;
compressibility=8.0906*10^-11; %input pore compressibility
while ~feof(fid)
 
    i=i+1;
    str=fscanf(fid,'%s',1);
    data_injeksi.date(i)=datenum(str,'dd/mm/yyyy');
    data_injeksi.cumvol(i)=fscanf(fid,'%f',1);
    data_injeksi.pp_injected(i)=(1/porevol)*(data_injeksi.cumvol(i)/compressibility)/10^6;
end
 
ndata=i;
 
fclose(fid);
 
%% Before Injection (in-situ state)
 
%input Sv, SH max, Sh min, pore pressure
 
%Sv hydrostatic=0.92 Psi/ft, Sv overpressure=0.92 Psi/ft
%SHmax hydrostatic=1.52 Psi/ft, SHmax overpressure=1.25 Psi/ft
%Shmin hydrostatic=0.75 Psi/ft, Shmin overpressure=0.77 Psi/ft
%Pp hydrostatic=0.43 Psi/ft, Pp overpressure=0.54 Psi/ft
 
input_depth='What is the depth of CO2 injection reservoir (in metre)?  ';
depth=input(input_depth); %injection depth range=2784 to 2874 m
Sv=0.92*3.28084*0.00689*depth; %input in Psi, converted to MPa
SHmax=1.52*3.28084*0.00689*depth; %input in Psi, converted to MPa
Shmin=0.75*3.28084*0.00689*depth; %input in Psi, converted to MPa
pore_pressure=0.43*3.28084*0.00689*depth; %input in Psi, converted to MPa
 
%convert to S1,S2,S3
sigma_matrix=[Sv SHmax Shmin];
ordered_sigma=sort(sigma_matrix,'descend');
sigma_1=ordered_sigma(1); %sigma 1 is the largest principal stress
sigma_2=ordered_sigma(2);
sigma_3=ordered_sigma(3);
 
tau_1  =(sigma_1-sigma_3)/2;
tau_2  =(sigma_1-sigma_2)/2;
tau_3  =(sigma_2-sigma_3)/2;
 
%Mohr-Coulomb envelope
friction_angle=28; %input friction angle
c_zero=6.7; %input cohesive strength
sigma_normal=0:1:sigma_1;
coeff_friction_angle=tand(friction_angle);
sigma_shear=c_zero+(coeff_friction_angle*(sigma_normal));
 
%Plotting the 3-D Mohr's Circle
angles=0:0.01:2*pi;
 
for i=1:30:ndata %time increment 30 days
   
Centre1=[(((sigma_1+sigma_3)/2)-pore_pressure) 0];
Centre2=[(((sigma_1+sigma_2)/2)-pore_pressure) 0];
Centre3=[(((sigma_2+sigma_3)/2)-pore_pressure) 0];
cirlce1=[Centre1(1)+tau_1*cos(angles') Centre1(2)+tau_1*sin(angles')];
cirlce2=[Centre2(1)+tau_2*cos(angles') Centre2(2)+tau_2*sin(angles')];
cirlce3=[Centre3(1)+tau_3*cos(angles') Centre3(2)+tau_3*sin(angles')];
 
subplot(1,2,1)
 
plot(sigma_shear,'r');
axis([0 80 -40 40]);
hold on
plot(cirlce1(:,1),cirlce1(:,2),'b',cirlce2(:,1),cirlce2(:,2),'g',...
    cirlce3(:,1),cirlce3(:,2),'r'),axis Equal;grid on;
title('Mohr-Coulomb before injection')
 
%Display annotation
if sigma_1>0
    text((sigma_1-pore_pressure)*1.01,0,'\sigma\prime_1','fontsize',15);
else
    text((sigma_1-pore_pressure)*0.95,0,'\sigma\prime_1','fontsize',15)
end
if sigma_2>0
    text((sigma_2-pore_pressure)*1.1,0,'\sigma\prime_2','fontsize',15)
else
    text((sigma_2-pore_pressure)*0.99,0,'\sigma\prime_2','fontsize',15)
end
if sigma_3>0
    text((sigma_3-pore_pressure)*1.1,0,'\sigma\prime_3','fontsize',15);
else
    text((sigma_3-pore_pressure)*0.99,0,'\sigma\prime_3','fontsize',15);
end
text(Centre1(1),tau_1*0.9,'\tau_1','fontsize',15);
text(Centre2(1),tau_2*0.9,'\tau_2','fontsize',15);
text(Centre3(1),tau_3*0.9,'\tau_3','fontsize',15);
xlabel('Effective Normal Stress, \sigma\prime','fontsize',15);
ylabel('Shear Stress, \tau','fontsize',15)
 
%% After injection 
 
%pore pressure increase
pore_pressure_new=pore_pressure+data_injeksi.pp_injected(i);
 
%Mohr-Coulomb envelope
sigma_shear=c_zero+(coeff_friction_angle*(sigma_normal));
 
%Plotting the 3-D Mohr's Circle after injection
angles=0:0.01:2*pi;
Centre1_new=[(((sigma_1+sigma_3)/2)-pore_pressure_new) 0];
Centre2_new=[(((sigma_1+sigma_2)/2)-pore_pressure_new) 0];
Centre3_new=[(((sigma_2+sigma_3)/2)-pore_pressure_new) 0];
cirlce1_new=[Centre1_new(1)+tau_1*cos(angles') Centre1_new(2)+tau_1*sin(angles')];
cirlce2_new=[Centre2_new(1)+tau_2*cos(angles') Centre2_new(2)+tau_2*sin(angles')];
cirlce3_new=[Centre3_new(1)+tau_3*cos(angles') Centre3_new(2)+tau_3*sin(angles')];
 
subplot(1,2,2)
plot(sigma_shear,'b');
axis([0 80 -40 40]);
hold on
plot(cirlce1_new(:,1),cirlce1_new(:,2),'b',cirlce2_new(:,1),cirlce2_new(:,2),'g',...
    cirlce3_new(:,1),cirlce3_new(:,2),'r'),axis Equal;grid on;
hold off
title('Mohr-Coulomb after injection')
 
%Display annotation
if sigma_1>0
    text((sigma_1-pore_pressure_new)*1.01,0,'\sigma\prime_1','fontsize',15);
else
    text((sigma_1-pore_pressure_new)*0.95,0,'\sigma\prime_1','fontsize',15)
end
if sigma_2>0
    text((sigma_2-pore_pressure_new)*1.1,0,'\sigma\prime_2','fontsize',15)
else
    text((sigma_2-pore_pressure_new)*0.99,0,'\sigma\prime_2','fontsize',15)
end
if sigma_3>0
    text((sigma_3-pore_pressure_new)*1.1,0,'\sigma\prime_3','fontsize',15);
else
    text((sigma_3-pore_pressure_new)*0.99,0,'\sigma\prime_3','fontsize',15);
end
text(Centre1_new(1),tau_1*0.9,'\tau_1','fontsize',15);
text(Centre2_new(1),tau_2*0.9,'\tau_2','fontsize',15);
text(Centre3_new(1),tau_3*0.9,'\tau_3','fontsize',15);
xlabel('Effective Normal Stress, \sigma\prime','fontsize',15);
ylabel('Shear Stress, \tau','fontsize',15)
 
%Critical injection
sigma_shear_new=c_zero+(coeff_friction_angle*(cirlce1_new(:,1)));
i1=(cirlce1_new(:,2)>=sigma_shear_new);
if sum(i1)>0 %sigma dan envelope adalah vektor jadi harus dikasih sum...
    stopdate=datetime(data_injeksi.date(i),'ConvertFrom','datenum');
    disp(stopdate);
    break
end
 
truedate=datetime(data_injeksi.date,'ConvertFrom','datenum');
subplot(1,2,2)
hold on
ket=['Injection date: ' str2mat(truedate(i)) ''];
text(4,40,ket)
hold off
 
%movie
pause(0.1);
 
end
%% Maximum allowable pore pressure change
 
critical_angle=90+friction_angle;
critical_normal_stress=((sigma_1+sigma_3)/2)+(((sigma_1-sigma_3)/2)*cosd(critical_angle));
critical_shear_stress=((sigma_1-sigma_3)/2)*sind(critical_angle);
pore_pressure_crit=critical_normal_stress+((c_zero-critical_shear_stress)/tand(friction_angle));
 
txta=['Initial pore pressure: ' num2str(pore_pressure) ' MPa'];
text(4,50,txta)
txtb=['Critical pore pressure: ' num2str(pore_pressure_crit) ' MPa'];
text(4,45,txtb)
txtc=['Stop injection at: ' str2mat(stopdate) ''];
text(4,40,txtc)
