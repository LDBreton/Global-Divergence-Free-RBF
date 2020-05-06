function fbraux = fbrGram24(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM24
%    FBRAUX = FBRGRAM24(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Nov-2019 13:45:07

t2 = c.^2;
t3 = c.^3;
t4 = x1.*2.0;
t5 = x2.*2.0;
t6 = y1.*2.0;
t7 = y2.*2.0;
t8 = -x2;
t11 = -y2;
t13 = dtt.*muu.*6.0e+1;
t9 = -t5;
t10 = -t6;
t12 = -t7;
t14 = t8+x1;
t15 = t11+y1;
t16 = t14.^2;
t18 = t15.^2;
t32 = t7+t10+t14;
t33 = t6+t12+t14;
t17 = t16.^2;
t19 = t18.^2;
t20 = t16.*2.0;
t21 = t16.*5.0;
t23 = t18.*2.0;
t26 = t16.*2.3e+1;
t27 = t18.*9.0;
t28 = t18.*1.1e+1;
t29 = t18.*3.3e+1;
t34 = t16+t18;
t35 = t16.*t18.*4.0;
t22 = t17.*3.0;
t24 = t19.*7.0;
t30 = -t27;
t31 = -t29;
t36 = t20+t23;
t39 = t13+t16+t28;
t25 = -t22;
t37 = t21+t30;
t38 = t26+t31;
t41 = c.*t39;
t43 = -t33.*t36.*(t6+t12-x1+x2);
t40 = dtt.*muu.*t37.*5.0;
t42 = dtt.*muu.*t38.*5.0;
t44 = t24+t25+t35+t40;
t45 = t42+t43;
t48 = t3.*t34.*(-t22+t24+t35+t40);
t46 = t2.*t45.*2.0;
t47 = -t46;
t49 = t41+t47+t48+2.0;
fbraux = muu.*t2.*1.0./(c.*t34+1.0).^(1.1e+1./2.0).*(-NX1.^2.*t14.*t49+NY1.^2.*t14.*t49+NX1.*NY1.*t15.*(t2.*(t17.*7.0-t19.*3.0+t35+dtt.*muu.*(t16.*5.1e+1-t18.*5.0).*5.0)-c.*(t18.*3.0-t20+dtt.*muu.*3.0e+1)+t3.*t34.*(t34.*(t4+t9+t15).*(t4+t9-y1+y2)-dtt.*muu.*(t16.*6.0-t18).*5.0)-1.0).*4.0).*3.0;
