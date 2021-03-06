function fbraux = fbrGram24(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM24
%    FBRAUX = FBRGRAM24(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version mp('8.3').
%    07-Nov-2019 13:45:07

t2 = c.^2;
t3 = c.^3;
t4 = x1.*mp('2.0');
t5 = x2.*mp('2.0');
t6 = y1.*mp('2.0');
t7 = y2.*mp('2.0');
t8 = -x2;
t11 = -y2;
t13 = dtt.*muu.*mp('6.0e+1');
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
t20 = t16.*mp('2.0');
t21 = t16.*mp('5.0');
t23 = t18.*mp('2.0');
t26 = t16.*mp('2.3e+1');
t27 = t18.*mp('9.0');
t28 = t18.*mp('1.1e+1');
t29 = t18.*mp('3.3e+1');
t34 = t16+t18;
t35 = t16.*t18.*mp('4.0');
t22 = t17.*mp('3.0');
t24 = t19.*mp('7.0');
t30 = -t27;
t31 = -t29;
t36 = t20+t23;
t39 = t13+t16+t28;
t25 = -t22;
t37 = t21+t30;
t38 = t26+t31;
t41 = c.*t39;
t43 = -t33.*t36.*(t6+t12-x1+x2);
t40 = dtt.*muu.*t37.*mp('5.0');
t42 = dtt.*muu.*t38.*mp('5.0');
t44 = t24+t25+t35+t40;
t45 = t42+t43;
t48 = t3.*t34.*(-t22+t24+t35+t40);
t46 = t2.*t45.*mp('2.0');
t47 = -t46;
t49 = t41+t47+t48+mp('2.0');
fbraux = muu.*t2.*mp('1.0')./(c.*t34+mp('1.0')).^(mp('1.1e+1')./mp('2.0')).*(-NX1.^2.*t14.*t49+NY1.^2.*t14.*t49+NX1.*NY1.*t15.*(t2.*(t17.*mp('7.0')-t19.*mp('3.0')+t35+dtt.*muu.*(t16.*mp('5.1e+1')-t18.*mp('5.0')).*mp('5.0'))-c.*(t18.*mp('3.0')-t20+dtt.*muu.*mp('3.0e+1'))+t3.*t34.*(t34.*(t4+t9+t15).*(t4+t9-y1+y2)-dtt.*muu.*(t16.*mp('6.0')-t18).*mp('5.0'))-mp('1.0')).*mp('4.0')).*mp('3.0');
