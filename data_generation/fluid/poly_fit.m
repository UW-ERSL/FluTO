num = xlsread("fluidMicrostructureData_exp");

%data gen
vol_nums=linspace(0,1,7);
startp=1;
endp=7;
vol=num(startp:endp,2);
area=num(startp:endp,3);
c0=num(startp:endp,4);
c1=num(startp:endp,5);

%coeff gen
coeff_area = polyfit(vol,area,5)
area_poly = polyval(coeff_area,vol_nums);

coeff_c0 = polyfit(vol,c0,5)
c0_poly = polyval(coeff_c0,vol_nums);

coeff_c1 = polyfit(vol,c1,5)
c1_poly = polyval(coeff_c1,vol_nums);

writeDict([coeff_c0; coeff_c1; coeff_area]);

plot(vol,c0,'g');
hold on;
plot(vol_nums,c0_poly,'b');
%scatter(vol,c0,250,'r');
%hold on;
% scatter(vol,c0,'g');
% xlabel('Volume Fraction'), ylabel('C ');
% title(' C  vs Volume Frac. Squircle');
legend('actual','polyfit')
