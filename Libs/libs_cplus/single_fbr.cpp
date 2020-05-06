#include <fstream>
#include <stdio.h>
#include <math.h>
#include <Eigen/MPRealSupport>

using namespace mpfr;
using namespace Eigen;
using namespace std;  

mpreal fbrGram11 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(c*(x1*x1)+c*(x2*x2)-c*(y1*y1)*2.0-c*(y2*y2)*2.0-c*x1*x2*2.0+c*y1*y2*4.0+1.0);

return t0;
}

mpreal fbrGram12 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = (c*c)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(3.0/4.0);

return t0;
}

mpreal fbrGram13 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(c*(x1*x1)+c*(x2*x2)-c*(y1*y1)*2.0-c*(y2*y2)*2.0-c*x1*x2*2.0+c*y1*y2*4.0+1.0);

return t0;
}

mpreal fbrGram14 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = (c*c)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(3.0/4.0);

return t0;
}

mpreal fbrGram15 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = x1*x1;
  t5 = c*t4;
  t6 = x2*x2;
  t7 = c*t6;
  t8 = y1*y1;
  t9 = y2*y2;
  t10 = c*c;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(t5+t7-c*t8*2.0-c*t9*2.0-c*x1*x2*2.0+c*y1*y2*4.0+1.0)+dt*mus*t10*1.0/pow(t5+t7+c*t8+c*t9-c*x1*x2*2.0-c*y1*y2*2.0+1.0,9.0/2.0)*(c*t4*3.0+c*t6*3.0-c*t8*2.7E1-c*t9*2.7E1-(t4*t4)*t10-(t6*t6)*t10+(t8*t8)*t10*4.0+(t9*t9)*t10*4.0-c*x1*x2*6.0+c*y1*y2*5.4E1-t4*t6*t10*6.0+t4*t8*t10*3.0+t4*t9*t10*3.0+t6*t8*t10*3.0+t6*t9*t10*3.0+t8*t9*t10*2.4E1+t4*t10*x1*x2*4.0+t6*t10*x1*x2*4.0-t8*t10*x1*x2*6.0-t9*t10*x1*x2*6.0-t4*t10*y1*y2*6.0-t6*t10*y1*y2*6.0-t8*t10*y1*y2*1.6E1-t9*t10*y1*y2*1.6E1+t10*x1*x2*y1*y2*1.2E1+4.0)*3.0;

return t0;
}

mpreal fbrGram16 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = c*c;
  t5 = t2*t2;
  t6 = t3*t3;
  t7 = t5+t6;
  t8 = c*t7;
  t9 = t8+1.0;
  t0 = t4*1.0/pow(t9,5.0/2.0)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*(3.0/4.0)-c*dt*mus*t2*t3*t4*1.0/pow(t9,9.0/2.0)*(c*(x1*x1)+c*(x2*x2)+c*(y1*y1)+c*(y2*y2)-c*x1*x2*2.0-c*y1*y2*2.0-6.0)*1.5E1;

return t0;
}

mpreal fbrGram21 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = (c*c)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(3.0/4.0);

return t0;
}

mpreal fbrGram22 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(c*(x1*x1)*-2.0-c*(x2*x2)*2.0+c*(y1*y1)+c*(y2*y2)+c*x1*x2*4.0-c*y1*y2*2.0+1.0);

return t0;
}

mpreal fbrGram23 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = (c*c)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(3.0/4.0);

return t0;
}

mpreal fbrGram24 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(c*(x1*x1)*-2.0-c*(x2*x2)*2.0+c*(y1*y1)+c*(y2*y2)+c*x1*x2*4.0-c*y1*y2*2.0+1.0);

return t0;
}

mpreal fbrGram25 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = c*c;
  t5 = t2*t2;
  t6 = t3*t3;
  t7 = t5+t6;
  t8 = c*t7;
  t9 = t8+1.0;
  t0 = t4*1.0/pow(t9,5.0/2.0)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*(3.0/4.0)-c*dt*mus*t2*t3*t4*1.0/pow(t9,9.0/2.0)*(c*(x1*x1)+c*(x2*x2)+c*(y1*y1)+c*(y2*y2)-c*x1*x2*2.0-c*y1*y2*2.0-6.0)*1.5E1;

return t0;
}

mpreal fbrGram26 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = x1*x1;
  t5 = x2*x2;
  t6 = y1*y1;
  t7 = c*t6;
  t8 = y2*y2;
  t9 = c*t8;
  t10 = c*c;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(t7+t9-c*t4*2.0-c*t5*2.0+c*x1*x2*4.0-c*y1*y2*2.0+1.0)+dt*mus*t10*1.0/pow(t7+t9+c*t4+c*t5-c*x1*x2*2.0-c*y1*y2*2.0+1.0,9.0/2.0)*(c*t4*-2.7E1-c*t5*2.7E1+c*t6*3.0+c*t8*3.0+(t4*t4)*t10*4.0+(t5*t5)*t10*4.0-(t6*t6)*t10-(t8*t8)*t10+c*x1*x2*5.4E1-c*y1*y2*6.0+t4*t5*t10*2.4E1+t4*t6*t10*3.0+t5*t6*t10*3.0+t4*t8*t10*3.0+t5*t8*t10*3.0-t6*t8*t10*6.0-t4*t10*x1*x2*1.6E1-t5*t10*x1*x2*1.6E1-t6*t10*x1*x2*6.0-t8*t10*x1*x2*6.0-t4*t10*y1*y2*6.0-t5*t10*y1*y2*6.0+t6*t10*y1*y2*4.0+t8*t10*y1*y2*4.0+t10*x1*x2*y1*y2*1.2E1+4.0)*3.0;

return t0;
}

mpreal fbrGram31 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(c*(x1*x1)+c*(x2*x2)-c*(y1*y1)*2.0-c*(y2*y2)*2.0-c*x1*x2*2.0+c*y1*y2*4.0+1.0);

return t0;
}

mpreal fbrGram32 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = (c*c)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(3.0/4.0);

return t0;
}

mpreal fbrGram33 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(c*(x1*x1)+c*(x2*x2)-c*(y1*y1)*2.0-c*(y2*y2)*2.0-c*x1*x2*2.0+c*y1*y2*4.0+1.0);

return t0;
}

mpreal fbrGram34 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = (c*c)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(3.0/4.0);

return t0;
}

mpreal fbrGram35 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = x1*x1;
  t5 = c*t4;
  t6 = x2*x2;
  t7 = c*t6;
  t8 = y1*y1;
  t9 = y2*y2;
  t10 = c*c;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(t5+t7-c*t8*2.0-c*t9*2.0-c*x1*x2*2.0+c*y1*y2*4.0+1.0)+dt*mus*t10*1.0/pow(t5+t7+c*t8+c*t9-c*x1*x2*2.0-c*y1*y2*2.0+1.0,9.0/2.0)*(c*t4*3.0+c*t6*3.0-c*t8*2.7E1-c*t9*2.7E1-(t4*t4)*t10-(t6*t6)*t10+(t8*t8)*t10*4.0+(t9*t9)*t10*4.0-c*x1*x2*6.0+c*y1*y2*5.4E1-t4*t6*t10*6.0+t4*t8*t10*3.0+t4*t9*t10*3.0+t6*t8*t10*3.0+t6*t9*t10*3.0+t8*t9*t10*2.4E1+t4*t10*x1*x2*4.0+t6*t10*x1*x2*4.0-t8*t10*x1*x2*6.0-t9*t10*x1*x2*6.0-t4*t10*y1*y2*6.0-t6*t10*y1*y2*6.0-t8*t10*y1*y2*1.6E1-t9*t10*y1*y2*1.6E1+t10*x1*x2*y1*y2*1.2E1+4.0)*3.0;

return t0;
}

mpreal fbrGram36 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = c*c;
  t5 = t2*t2;
  t6 = t3*t3;
  t7 = t5+t6;
  t8 = c*t7;
  t9 = t8+1.0;
  t0 = t4*1.0/pow(t9,5.0/2.0)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*(3.0/4.0)-c*dt*mus*t2*t3*t4*1.0/pow(t9,9.0/2.0)*(c*(x1*x1)+c*(x2*x2)+c*(y1*y1)+c*(y2*y2)-c*x1*x2*2.0-c*y1*y2*2.0-6.0)*1.5E1;

return t0;
}

mpreal fbrGram41 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = (c*c)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(3.0/4.0);

return t0;
}

mpreal fbrGram42 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(c*(x1*x1)*-2.0-c*(x2*x2)*2.0+c*(y1*y1)+c*(y2*y2)+c*x1*x2*4.0-c*y1*y2*2.0+1.0);

return t0;
}

mpreal fbrGram43 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = (c*c)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(3.0/4.0);

return t0;
}

mpreal fbrGram44 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(c*(x1*x1)*-2.0-c*(x2*x2)*2.0+c*(y1*y1)+c*(y2*y2)+c*x1*x2*4.0-c*y1*y2*2.0+1.0);

return t0;
}

mpreal fbrGram45 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = c*c;
  t5 = t2*t2;
  t6 = t3*t3;
  t7 = t5+t6;
  t8 = c*t7;
  t9 = t8+1.0;
  t0 = t4*1.0/pow(t9,5.0/2.0)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*(3.0/4.0)-c*dt*mus*t2*t3*t4*1.0/pow(t9,9.0/2.0)*(c*(x1*x1)+c*(x2*x2)+c*(y1*y1)+c*(y2*y2)-c*x1*x2*2.0-c*y1*y2*2.0-6.0)*1.5E1;

return t0;
}

mpreal fbrGram46 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = x1*x1;
  t5 = x2*x2;
  t6 = y1*y1;
  t7 = c*t6;
  t8 = y2*y2;
  t9 = c*t8;
  t10 = c*c;
  t0 = c*1.0/pow(c*(t2*t2+t3*t3)+1.0,5.0/2.0)*(t7+t9-c*t4*2.0-c*t5*2.0+c*x1*x2*4.0-c*y1*y2*2.0+1.0)+dt*mus*t10*1.0/pow(t7+t9+c*t4+c*t5-c*x1*x2*2.0-c*y1*y2*2.0+1.0,9.0/2.0)*(c*t4*-2.7E1-c*t5*2.7E1+c*t6*3.0+c*t8*3.0+(t4*t4)*t10*4.0+(t5*t5)*t10*4.0-(t6*t6)*t10-(t8*t8)*t10+c*x1*x2*5.4E1-c*y1*y2*6.0+t4*t5*t10*2.4E1+t4*t6*t10*3.0+t5*t6*t10*3.0+t4*t8*t10*3.0+t5*t8*t10*3.0-t6*t8*t10*6.0-t4*t10*x1*x2*1.6E1-t5*t10*x1*x2*1.6E1-t6*t10*x1*x2*6.0-t8*t10*x1*x2*6.0-t4*t10*y1*y2*6.0-t5*t10*y1*y2*6.0+t6*t10*y1*y2*4.0+t8*t10*y1*y2*4.0+t10*x1*x2*y1*y2*1.2E1+4.0)*3.0;

return t0;
}

mpreal fbrGram51 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = t2*t2;
  t5 = t3*t3;
  t6 = t4+t5;
  t7 = c*t6;
  t8 = t7+1.0;
  t9 = c*c;
  t10 = 1.0/pow(t8,5.0/2.0);
  t11 = 1.0/pow(t8,7.0/2.0);
  t12 = t9*t9;
  t13 = 1.0/pow(t8,9.0/2.0);
  t0 = c*1.0/pow(t8,3.0/2.0)+dt*mus*(t9*t10*1.2E1+(t5*t5)*t12*t13*1.05E2-c*t4*t9*t11*1.5E1-c*t5*t9*t11*7.5E1+t4*t5*t12*t13*1.05E2-c*t9*t11*(y1*2.0-y2*2.0)*(y1*8.0-y2*8.0)*(1.5E1/8.0))-t5*t9*t10*3.0;

return t0;
}

mpreal fbrGram52 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = c*c;
  t5 = t2*t2;
  t6 = t3*t3;
  t7 = t5+t6;
  t8 = c*t7;
  t9 = t8+1.0;
  t0 = t4*1.0/pow(t9,5.0/2.0)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*(3.0/4.0)-c*dt*mus*t2*t3*t4*1.0/pow(t9,9.0/2.0)*(c*(x1*x1)+c*(x2*x2)+c*(y1*y1)+c*(y2*y2)-c*x1*x2*2.0-c*y1*y2*2.0-6.0)*1.5E1;

return t0;
}

mpreal fbrGram53 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = t2*t2;
  t5 = t3*t3;
  t6 = t4+t5;
  t7 = c*t6;
  t8 = t7+1.0;
  t9 = c*c;
  t10 = 1.0/pow(t8,5.0/2.0);
  t11 = 1.0/pow(t8,7.0/2.0);
  t12 = t9*t9;
  t13 = 1.0/pow(t8,9.0/2.0);
  t0 = c*1.0/pow(t8,3.0/2.0)+dt*mus*(t9*t10*1.2E1+(t5*t5)*t12*t13*1.05E2-c*t4*t9*t11*1.5E1-c*t5*t9*t11*7.5E1+t4*t5*t12*t13*1.05E2-c*t9*t11*(y1*2.0-y2*2.0)*(y1*8.0-y2*8.0)*(1.5E1/8.0))-t5*t9*t10*3.0;

return t0;
}

mpreal fbrGram54 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = c*c;
  t5 = t2*t2;
  t6 = t3*t3;
  t7 = t5+t6;
  t8 = c*t7;
  t9 = t8+1.0;
  t0 = t4*1.0/pow(t9,5.0/2.0)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*(3.0/4.0)-c*dt*mus*t2*t3*t4*1.0/pow(t9,9.0/2.0)*(c*(x1*x1)+c*(x2*x2)+c*(y1*y1)+c*(y2*y2)-c*x1*x2*2.0-c*y1*y2*2.0-6.0)*1.5E1;

return t0;
}

mpreal fbrGram55 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = t2*t2;
  t5 = t3*t3;
  t6 = t4+t5;
  t7 = c*t6;
  t8 = t7+1.0;
  t9 = c*c;
  t10 = 1.0/pow(t8,7.0/2.0);
  t11 = t9*t9;
  t12 = 1.0/pow(t8,9.0/2.0);
  t13 = t5*t5;
  t14 = 1.0/pow(t8,1.1E1/2.0);
  t15 = 1.0/pow(t8,1.3E1/2.0);
  t16 = y1*2.0;
  t21 = y2*2.0;
  t17 = t16-t21;
  t18 = y1*8.0;
  t22 = y2*8.0;
  t19 = t18-t22;
  t20 = t4*t4;
  t23 = 1.0/pow(t8,5.0/2.0);
  t24 = t9*t23*1.2E1;
  t25 = t11*t12*t13*1.05E2;
  t26 = t4*t5*t11*t12*1.05E2;
  t0 = c*1.0/pow(t8,3.0/2.0)+dt*mus*(t24+t25+t26+dt*mus*(c*t9*t10*2.7E2-t4*t11*t12*3.15E2-t5*t11*t12*5.145E3+c*t11*t13*t14*5.67E3-t11*t12*t17*t19*(1.05E2/8.0)+c*t4*t5*t11*t14*5.67E3-t4*t9*t11*t13*t15*1.0395E4-t5*t9*t11*t13*t15*1.0395E4+c*t3*t5*t11*t14*t17*3.78E3+c*t3*t5*t11*t14*t19*(9.45E2/4.0))+dt*mus*(c*t9*t10*9.0E1-t4*t11*t12*9.45E2-t5*t11*t12*7.35E2+c*t11*t13*t14*9.45E2+c*t11*t14*t20*9.45E2-t11*t12*t17*t19*(1.05E2/8.0)+c*t4*t5*t11*t14*9.45E3-t4*t9*t11*t13*t15*1.0395E4-t5*t9*t11*t15*t20*1.0395E4+c*t4*t11*t14*t17*t19*(9.45E2/8.0))-c*t4*t9*t10*1.5E1-c*t5*t9*t10*1.05E2)-t5*t9*t23*3.0+dt*mus*(t24+t25+t26-c*t4*t9*t10*1.5E1-c*t5*t9*t10*7.5E1-c*t9*t10*t17*t19*(1.5E1/8.0))+c*(dt*dt)*t23*(c*(x1*x1)*-2.0-c*(x2*x2)*2.0+c*(y1*y1)+c*(y2*y2)+c*x1*x2*4.0-c*y1*y2*2.0+1.0);

return t0;
}

mpreal fbrGram56 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = c*c;
  t3 = x1-x2;
  t4 = y1-y2;
  t5 = x1*2.0;
  t15 = x2*2.0;
  t6 = t5-t15;
  t7 = t3*t3;
  t8 = t4*t4;
  t9 = t7+t8;
  t10 = c*t9;
  t11 = t10+1.0;
  t12 = 1.0/pow(t11,7.0/2.0);
  t13 = y1*2.0;
  t17 = y2*2.0;
  t14 = t13-t17;
  t16 = t2*t2;
  t18 = 1.0/pow(t11,9.0/2.0);
  t19 = 1.0/pow(t11,1.1E1/2.0);
  t20 = x1*x1;
  t21 = c*t20;
  t22 = x2*x2;
  t23 = c*t22;
  t24 = y1*y1;
  t25 = c*t24;
  t26 = y2*y2;
  t27 = c*t26;
  t29 = c*x1*x2*2.0;
  t30 = c*y1*y2*2.0;
  t28 = t21+t23+t25+t27-t29-t30-6.0;
  t31 = 1.0/pow(t11,1.3E1/2.0);
  t32 = 1.0/pow(t11,5.0/2.0);
  t0 = -dt*mus*(c*t2*t12*t14*(x1*8.0-x2*8.0)*(-1.5E1/8.0)-c*t2*t6*t12*(y1*8.0-y2*8.0)*(1.5E1/8.0)-c*t2*t6*t12*t14*(1.5E1/2.0)+t4*t6*t8*t16*t18*(1.05E2/2.0)+t3*t7*t14*t16*t18*(1.05E2/2.0)-dt*mus*t3*t4*t16*t18*1.8E2+dt*mus*t3*t4*t16*t19*t28*2.7E2+dt*mus*t4*t6*t16*t19*t28*1.35E2+dt*mus*t3*t14*t16*t19*t28*1.35E2+c*dt*mus*t4*t6*t7*t16*t19*2.7E2+c*dt*mus*t3*t8*t14*t16*t19*2.7E2-c*dt*mus*t3*t4*t7*t16*t28*t31*1.485E3-c*dt*mus*t3*t4*t8*t16*t28*t31*1.485E3)+t2*t6*t14*t32*(3.0/4.0)-(dt*dt)*t2*t3*t14*t32*(3.0/2.0)-c*dt*mus*t2*t3*t4*t18*t28*1.5E1;

return t0;
}

mpreal fbrGram61 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = c*c;
  t5 = t2*t2;
  t6 = t3*t3;
  t7 = t5+t6;
  t8 = c*t7;
  t9 = t8+1.0;
  t0 = t4*1.0/pow(t9,5.0/2.0)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*(3.0/4.0)-c*dt*mus*t2*t3*t4*1.0/pow(t9,9.0/2.0)*(c*(x1*x1)+c*(x2*x2)+c*(y1*y1)+c*(y2*y2)-c*x1*x2*2.0-c*y1*y2*2.0-6.0)*1.5E1;

return t0;
}

mpreal fbrGram62 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = t2*t2;
  t5 = t3*t3;
  t6 = t4+t5;
  t7 = c*t6;
  t8 = t7+1.0;
  t9 = c*c;
  t10 = 1.0/pow(t8,5.0/2.0);
  t11 = 1.0/pow(t8,7.0/2.0);
  t12 = t9*t9;
  t13 = 1.0/pow(t8,9.0/2.0);
  t0 = c*1.0/pow(t8,3.0/2.0)+dt*mus*(t9*t10*1.2E1+(t4*t4)*t12*t13*1.05E2-c*t4*t9*t11*7.5E1-c*t5*t9*t11*1.5E1+t4*t5*t12*t13*1.05E2-c*t9*t11*(x1*2.0-x2*2.0)*(x1*8.0-x2*8.0)*(1.5E1/8.0))-t4*t9*t10*3.0;

return t0;
}

mpreal fbrGram63 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = c*c;
  t5 = t2*t2;
  t6 = t3*t3;
  t7 = t5+t6;
  t8 = c*t7;
  t9 = t8+1.0;
  t0 = t4*1.0/pow(t9,5.0/2.0)*(x1*2.0-x2*2.0)*(y1*2.0-y2*2.0)*(3.0/4.0)-c*dt*mus*t2*t3*t4*1.0/pow(t9,9.0/2.0)*(c*(x1*x1)+c*(x2*x2)+c*(y1*y1)+c*(y2*y2)-c*x1*x2*2.0-c*y1*y2*2.0-6.0)*1.5E1;

return t0;
}

mpreal fbrGram64 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = t2*t2;
  t5 = t3*t3;
  t6 = t4+t5;
  t7 = c*t6;
  t8 = t7+1.0;
  t9 = c*c;
  t10 = 1.0/pow(t8,5.0/2.0);
  t11 = 1.0/pow(t8,7.0/2.0);
  t12 = t9*t9;
  t13 = 1.0/pow(t8,9.0/2.0);
  t0 = c*1.0/pow(t8,3.0/2.0)+dt*mus*(t9*t10*1.2E1+(t4*t4)*t12*t13*1.05E2-c*t4*t9*t11*7.5E1-c*t5*t9*t11*1.5E1+t4*t5*t12*t13*1.05E2-c*t9*t11*(x1*2.0-x2*2.0)*(x1*8.0-x2*8.0)*(1.5E1/8.0))-t4*t9*t10*3.0;

return t0;
}

mpreal fbrGram65 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = c*c;
  t3 = x1-x2;
  t4 = y1-y2;
  t5 = x1*2.0;
  t15 = x2*2.0;
  t6 = t5-t15;
  t7 = t3*t3;
  t8 = t4*t4;
  t9 = t7+t8;
  t10 = c*t9;
  t11 = t10+1.0;
  t12 = 1.0/pow(t11,7.0/2.0);
  t13 = y1*2.0;
  t17 = y2*2.0;
  t14 = t13-t17;
  t16 = t2*t2;
  t18 = 1.0/pow(t11,9.0/2.0);
  t19 = 1.0/pow(t11,1.1E1/2.0);
  t20 = x1*x1;
  t21 = c*t20;
  t22 = x2*x2;
  t23 = c*t22;
  t24 = y1*y1;
  t25 = c*t24;
  t26 = y2*y2;
  t27 = c*t26;
  t29 = c*x1*x2*2.0;
  t30 = c*y1*y2*2.0;
  t28 = t21+t23+t25+t27-t29-t30-6.0;
  t31 = 1.0/pow(t11,1.3E1/2.0);
  t32 = 1.0/pow(t11,5.0/2.0);
  t0 = -dt*mus*(c*t2*t12*t14*(x1*8.0-x2*8.0)*(-1.5E1/8.0)-c*t2*t6*t12*(y1*8.0-y2*8.0)*(1.5E1/8.0)-c*t2*t6*t12*t14*(1.5E1/2.0)+t4*t6*t8*t16*t18*(1.05E2/2.0)+t3*t7*t14*t16*t18*(1.05E2/2.0)-dt*mus*t3*t4*t16*t18*1.8E2+dt*mus*t3*t4*t16*t19*t28*2.7E2+dt*mus*t4*t6*t16*t19*t28*1.35E2+dt*mus*t3*t14*t16*t19*t28*1.35E2+c*dt*mus*t4*t6*t7*t16*t19*2.7E2+c*dt*mus*t3*t8*t14*t16*t19*2.7E2-c*dt*mus*t3*t4*t7*t16*t28*t31*1.485E3-c*dt*mus*t3*t4*t8*t16*t28*t31*1.485E3)+t2*t6*t14*t32*(3.0/4.0)-(dt*dt)*t2*t4*t6*t32*(3.0/2.0)-c*dt*mus*t2*t3*t4*t18*t28*1.5E1;

return t0;
}

mpreal fbrGram66 (mpreal x1,mpreal y1, mpreal x2, mpreal y2, mpreal c,mpreal mus, mpreal dt){
mpreal t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50;
  t2 = x1-x2;
  t3 = y1-y2;
  t4 = t2*t2;
  t5 = t3*t3;
  t6 = t4+t5;
  t7 = c*t6;
  t8 = t7+1.0;
  t9 = c*c;
  t10 = 1.0/pow(t8,5.0/2.0);
  t11 = 1.0/pow(t8,7.0/2.0);
  t12 = t9*t9;
  t13 = 1.0/pow(t8,9.0/2.0);
  t14 = t9*t10*1.2E1;
  t15 = t4*t4;
  t16 = t12*t13*t15*1.05E2;
  t17 = x1*2.0;
  t21 = x2*2.0;
  t18 = t17-t21;
  t19 = x1*8.0;
  t23 = x2*8.0;
  t20 = t19-t23;
  t22 = 1.0/pow(t8,1.1E1/2.0);
  t24 = 1.0/pow(t8,1.3E1/2.0);
  t25 = t5*t5;
  t26 = t4*t5*t12*t13*1.05E2;
  t0 = c*1.0/pow(t8,3.0/2.0)+dt*mus*(t14+t16+t26+dt*mus*(c*t9*t11*2.7E2-t4*t12*t13*5.145E3-t5*t12*t13*3.15E2+c*t12*t15*t22*5.67E3-t12*t13*t18*t20*(1.05E2/8.0)+c*t4*t5*t12*t22*5.67E3-t4*t9*t12*t15*t24*1.0395E4-t5*t9*t12*t15*t24*1.0395E4+c*t2*t4*t12*t18*t22*3.78E3+c*t2*t4*t12*t20*t22*(9.45E2/4.0))+dt*mus*(c*t9*t11*9.0E1-t4*t12*t13*7.35E2-t5*t12*t13*9.45E2+c*t12*t15*t22*9.45E2+c*t12*t22*t25*9.45E2-t12*t13*t18*t20*(1.05E2/8.0)+c*t4*t5*t12*t22*9.45E3-t5*t9*t12*t15*t24*1.0395E4-t4*t9*t12*t24*t25*1.0395E4+c*t5*t12*t18*t20*t22*(9.45E2/8.0))-c*t4*t9*t11*1.05E2-c*t5*t9*t11*1.5E1)-t4*t9*t10*3.0+dt*mus*(t14+t16+t26-c*t4*t9*t11*7.5E1-c*t5*t9*t11*1.5E1-c*t9*t11*t18*t20*(1.5E1/8.0))+c*(dt*dt)*t10*(c*(x1*x1)+c*(x2*x2)-c*(y1*y1)*2.0-c*(y2*y2)*2.0-c*x1*x2*2.0+c*y1*y2*4.0+1.0);

return t0;
}
