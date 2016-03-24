function [gyy_123,gyy_12,gyy_1,gam3y_12,l1y,l12,l13,l23_1,l2y_1,l3y_12]=multicoher(g11,g22,g33,gyy,g1y,g2y,g3y,g12,g13,g23,f)
%here "d" is the pressure, x,y, z are the horizontals and the vertical. 
%[pzz_xyd,pzz_xy,pzz_x,tdz_xy,gamzd_xy]=multicoher(px,py,pd,pz,conj(czx),conj(czy),conj(czd),cxy,conj(cdx),conj(cdy),f);
%[gyy_123,gyy_12,gyy_1,l3y_12,gam3y_12]=multicoher(px,py,pd,pz,conj(czx),conj(czy),conj(czd),cxy,conj(cdx),conj(cdy),f);

%The model is an output channel y-  that is related to a series of inputs,
%1, 2, 3.  The signals from 1, 2, 3 (here x, y and presssure) are removed
%in turn from channel y. 
%autospectra are called gyy, g11, etc.
%cross spectra are called g1y 

%% figure(111)
%% clf
%these are the transfer functions between channels (1 and y);  (1 and 2);  (1 and 3). 


l1y=g1y./g11;
l12=g12./g11;
l13=g13./g11;
%%loglog(f,abs(l1y),'r');hold on;
%%loglog(f,abs(l12),'b');hold on;
%%oglog(f,abs(l13),'k');


%coherences between same
gam1y=abs(g1y).^2./(g11.*gyy);
gam12=abs(g12).^2./(g11.*g22);
gam13=abs(g13).^2./(g11.*g33);
%%semilogx(f,smooth(abs(gam1y),100),'r');hold on;
%%semilogx(f,smooth(abs(gam12),100),'b');hold on;
%%semilogx(f,smooth(abs(gam13),100),'k');
%%ylim([0 1.2]);

%this is predicting the channels y, 2 and 3  from 1
gyy_1=gyy.*(1-gam1y);
g22_1=g22.*(1-gam12);
g33_1=g33.*(1-gam13);

%removing the effect of channel 1 from the cross spectra
g2y_1=g2y-conj(l12).*g1y;
g3y_1=g3y-conj(l13).*g1y;
g23_1=g23-conj(l12).*g13;

%transfer function between 2 and 3 after removing the effect of channel 1
l23_1=g23_1./g22_1;
%transfer function between 2 and y after removing the effect of channel 1
l2y_1=g2y_1./g22_1;


%coherence between (2 and y), (3 and y) and (2 and 3) after removing effect
%of channel 1
gam2y_1=abs(g2y_1).^2./(g22_1.*gyy_1);
gam3y_1=abs(g3y_1).^2./(g33_1.*gyy_1);
gam23_1=abs(g23_1).^2./(g33_1.*g22_1);

%autospectra after removing effects of channels 1 and 3
gyy_12=gyy_1.*(1-gam2y_1);
g33_12=g33_1.*(1-gam23_1);
g3y_12=g3y_1-conj(l23_1).*g2y_1;

%coherence between 3 and y after removing effects of 2 and 3 
gam3y_12=abs(g3y_12).^2./(g33_12.*gyy_12);
%%semilogx(f,smooth(abs(gam3y_12),100),'r');hold on;
%%ylim([0 1.2]);

%autospectr for y after removing effects of 1, 2 and 3 
gyy_123=gyy_12.*(1-gam3y_12);

%transfer function between 3 and y after removing the effects of channels 1
%and 2.  This is what is used in compliance where "y" is the vertical and 1
%and 2 are the two horizontals- the transfer function after removing the
%tilt noise... 
l3y_12=g3y_12./g33_12;
