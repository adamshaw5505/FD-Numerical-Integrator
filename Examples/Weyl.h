#pragma once
#define NNN  32
#define NX NNN
#define NY NNN
#define NZ NNN
#define NF 25
#define NG 3

// #define KOEPS Base(-0.0001)

#include "../Headers/Definitions.h"
#include "../Headers/Evolution.h"

#define F_DEF(u,I,DX,DY,DZ) F[6]=DY(u,21+I,dy)-DZ(u,20+I,dz)+Base(0.5)*u[14+I]*u[18+I]-Base(0.5)*u[17+I]*u[15+I]-Base(0.5)*u[15+I]*u[17+I]+Base(0.5)*u[18+I]*u[14+I]; \
F[7]=-DX(u,21+I,dx)+DZ(u,19+I,dz)-Base(0.5)*u[13+I]*u[18+I]+Base(0.5)*u[16+I]*u[15+I]+Base(0.5)*u[15+I]*u[16+I]-Base(0.5)*u[18+I]*u[13+I];\
F[8]=DX(u,20+I,dx)-DY(u,19+I,dy)+Base(0.5)*u[13+I]*u[17+I]-Base(0.5)*u[16+I]*u[14+I]-Base(0.5)*u[14+I]*u[16+I]+Base(0.5)*u[17+I]*u[13+I];\
F[3]=DY(u,18+I,dy)-DZ(u,17+I,dz)-Base(0.5)*u[14+I]*u[21+I]+Base(0.5)*u[20+I]*u[15+I]+Base(0.5)*u[15+I]*u[20+I]-Base(0.5)*u[21+I]*u[14+I];\
F[4]=-DX(u,18+I,dx)+DZ(u,16+I,dz)+Base(0.5)*u[13+I]*u[21+I]-Base(0.5)*u[19+I]*u[15+I]-Base(0.5)*u[15+I]*u[19+I]+Base(0.5)*u[21+I]*u[13+I];\
F[5]=DX(u,17+I,dx)-DY(u,16+I,dy)-Base(0.5)*u[13+I]*u[20+I]+Base(0.5)*u[19+I]*u[14+I]+Base(0.5)*u[14+I]*u[19+I]-Base(0.5)*u[20+I]*u[13+I];\
F[0]=DY(u,15+I,dy)-DZ(u,14+I,dz)+Base(0.5)*u[17+I]*u[21+I]-Base(0.5)*u[20+I]*u[18+I]-Base(0.5)*u[18+I]*u[20+I]+Base(0.5)*u[21+I]*u[17+I];\
F[1]=-DX(u,15+I,dx)+DZ(u,13+I,dz)-Base(0.5)*u[16+I]*u[21+I]+Base(0.5)*u[19+I]*u[18+I]+Base(0.5)*u[18+I]*u[19+I]-Base(0.5)*u[21+I]*u[16+I];\
F[2]=DX(u,14+I,dx)-DY(u,13+I,dy)+Base(0.5)*u[16+I]*u[20+I]-Base(0.5)*u[19+I]*u[17+I]-Base(0.5)*u[17+I]*u[19+I]+Base(0.5)*u[20+I]*u[16+I];

#define E_DEF(u,I) E[0]=F[0]*u[4+I]*u[8+I]-F[0]*u[5+I]*u[7+I]-F[3]*u[3+I]*u[8+I]+F[3]*u[5+I]*u[6+I]+F[6]*u[3+I]*u[7+I]-F[6]*u[4+I]*u[6+I];\
E[1]=-F[0]*u[1+I]*u[8+I]+F[0]*u[2+I]*u[7+I]+F[3]*u[0+I]*u[8+I]-F[3]*u[2+I]*u[6+I]-F[6]*u[0+I]*u[7+I]+F[6]*u[1+I]*u[6+I];\
E[2]=F[0]*u[1+I]*u[5+I]-F[0]*u[2+I]*u[4+I]-F[3]*u[0+I]*u[5+I]+F[3]*u[2+I]*u[3+I]+F[6]*u[0+I]*u[4+I]-F[6]*u[1+I]*u[3+I];\
E[3]=F[1]*u[4+I]*u[8+I]-F[1]*u[5+I]*u[7+I]-F[4]*u[3+I]*u[8+I]+F[4]*u[5+I]*u[6+I]+F[7]*u[3+I]*u[7+I]-F[7]*u[4+I]*u[6+I];\
E[4]=-F[1]*u[1+I]*u[8+I]+F[1]*u[2+I]*u[7+I]+F[4]*u[0+I]*u[8+I]-F[4]*u[2+I]*u[6+I]-F[7]*u[0+I]*u[7+I]+F[7]*u[1+I]*u[6+I];\
E[5]=F[1]*u[1+I]*u[5+I]-F[1]*u[2+I]*u[4+I]-F[4]*u[0+I]*u[5+I]+F[4]*u[2+I]*u[3+I]+F[7]*u[0+I]*u[4+I]-F[7]*u[1+I]*u[3+I];\
E[6]=F[2]*u[4+I]*u[8+I]-F[2]*u[5+I]*u[7+I]-F[5]*u[3+I]*u[8+I]+F[5]*u[5+I]*u[6+I]+F[8]*u[3+I]*u[7+I]-F[8]*u[4+I]*u[6+I];\
E[7]=-F[2]*u[1+I]*u[8+I]+F[2]*u[2+I]*u[7+I]+F[5]*u[0+I]*u[8+I]-F[5]*u[2+I]*u[6+I]-F[8]*u[0+I]*u[7+I]+F[8]*u[1+I]*u[6+I];\
E[8]=F[2]*u[1+I]*u[5+I]-F[2]*u[2+I]*u[4+I]-F[5]*u[0+I]*u[5+I]+F[5]*u[2+I]*u[3+I]+F[8]*u[0+I]*u[4+I]-F[8]*u[1+I]*u[3+I];\
Field det_Psi = Base(1.0/3)*u[0+I]*u[0+I]*u[0+I]+Base(1.0/3)*u[0+I]*u[3+I]*u[1+I]+Base(1.0/3)*u[0+I]*u[6+I]*u[2+I]+Base(1.0/3)*u[3+I]*u[1+I]*u[0+I]+Base(1.0/3)*u[3+I]*u[4+I]*u[1+I]+Base(1.0/3)*u[3+I]*u[7+I]*u[2+I]+Base(1.0/3)*u[6+I]*u[2+I]*u[0+I]+Base(1.0/3)*u[6+I]*u[5+I]*u[1+I]+Base(1.0/3)*u[6+I]*u[8+I]*u[2+I]+Base(1.0/3)*u[1+I]*u[0+I]*u[3+I]+Base(1.0/3)*u[1+I]*u[3+I]*u[4+I]+Base(1.0/3)*u[1+I]*u[6+I]*u[5+I]+Base(1.0/3)*u[4+I]*u[1+I]*u[3+I]+Base(1.0/3)*u[4+I]*u[4+I]*u[4+I]+Base(1.0/3)*u[4+I]*u[7+I]*u[5+I]+Base(1.0/3)*u[7+I]*u[2+I]*u[3+I]+Base(1.0/3)*u[7+I]*u[5+I]*u[4+I]+Base(1.0/3)*u[7+I]*u[8+I]*u[5+I]+Base(1.0/3)*u[2+I]*u[0+I]*u[6+I]+Base(1.0/3)*u[2+I]*u[3+I]*u[7+I]+Base(1.0/3)*u[2+I]*u[6+I]*u[8+I]+Base(1.0/3)*u[5+I]*u[1+I]*u[6+I]+Base(1.0/3)*u[5+I]*u[4+I]*u[7+I]+Base(1.0/3)*u[5+I]*u[7+I]*u[8+I]+Base(1.0/3)*u[8+I]*u[2+I]*u[6+I]+Base(1.0/3)*u[8+I]*u[5+I]*u[7+I]+Base(1.0/3)*u[8+I]*u[8+I]*u[8+I];\
E[0] = E[0]/det_Psi;\
E[1] = E[1]/det_Psi;\
E[2] = E[2]/det_Psi;\
E[3] = E[3]/det_Psi;\
E[4] = E[4]/det_Psi;\
E[5] = E[5]/det_Psi;\
E[6] = E[6]/det_Psi;\
E[7] = E[7]/det_Psi;\
E[8] = E[8]/det_Psi;


#define WEYL_RHS(u,u_rhs,I,DX,DY,DZ) u_rhs[8+I]=u[12+I]*DZ(u,8+I,dz)+u[11+I]*DY(u,8+I,dy)+u[10+I]*DX(u,8+I,dx)+IM*u[9+I]*E[2]*DZ(u,7+I,dz)+IM*u[9+I]*E[1]*DY(u,7+I,dy)+IM*u[9+I]*E[0]*DX(u,7+I,dx)-IM*u[9+I]*E[5]*DZ(u,6+I,dz)-IM*u[9+I]*E[4]*DY(u,6+I,dy)-IM*u[9+I]*E[3]*DX(u,6+I,dx)+u[13+I]*IM*u[9+I]*E[0]*u[4+I]+u[14+I]*IM*u[9+I]*E[1]*u[4+I]+u[15+I]*IM*u[9+I]*E[2]*u[4+I]-u[16+I]*IM*u[9+I]*E[0]*u[1+I]-u[17+I]*IM*u[9+I]*E[1]*u[1+I]-u[18+I]*IM*u[9+I]*E[2]*u[1+I]-u[13+I]*IM*u[9+I]*E[3]*u[3+I]-u[14+I]*IM*u[9+I]*E[4]*u[3+I]-u[15+I]*IM*u[9+I]*E[5]*u[3+I]+u[16+I]*IM*u[9+I]*E[3]*u[0+I]+u[17+I]*IM*u[9+I]*E[4]*u[0+I]+u[18+I]*IM*u[9+I]*E[5]*u[0+I]-u[13+I]*IM*u[9+I]*E[0]*u[8+I]-u[14+I]*IM*u[9+I]*E[1]*u[8+I]-u[15+I]*IM*u[9+I]*E[2]*u[8+I]+u[19+I]*IM*u[9+I]*E[0]*u[6+I]+u[20+I]*IM*u[9+I]*E[1]*u[6+I]+u[21+I]*IM*u[9+I]*E[2]*u[6+I]-u[16+I]*IM*u[9+I]*E[3]*u[8+I]-u[17+I]*IM*u[9+I]*E[4]*u[8+I]-u[18+I]*IM*u[9+I]*E[5]*u[8+I]+u[19+I]*IM*u[9+I]*E[3]*u[7+I]+u[20+I]*IM*u[9+I]*E[4]*u[7+I]+u[21+I]*IM*u[9+I]*E[5]*u[7+I]-u[2+I]*u[5+I]*u[22+I]+u[2+I]*u[2+I]*u[23+I]; \
u_rhs[5+I]=u[12+I]*DZ(u,5+I,dz)+u[11+I]*DY(u,5+I,dy)+u[10+I]*DX(u,5+I,dx)+IM*u[9+I]*E[2]*DZ(u,4+I,dz)+IM*u[9+I]*E[1]*DY(u,4+I,dy)+IM*u[9+I]*E[0]*DX(u,4+I,dx)-IM*u[9+I]*E[5]*DZ(u,3+I,dz)-IM*u[9+I]*E[4]*DY(u,3+I,dy)-IM*u[9+I]*E[3]*DX(u,3+I,dx)-u[13+I]*IM*u[9+I]*E[0]*u[7+I]-u[14+I]*IM*u[9+I]*E[1]*u[7+I]-u[15+I]*IM*u[9+I]*E[2]*u[7+I]+u[19+I]*IM*u[9+I]*E[0]*u[1+I]+u[20+I]*IM*u[9+I]*E[1]*u[1+I]+u[21+I]*IM*u[9+I]*E[2]*u[1+I]+u[13+I]*IM*u[9+I]*E[3]*u[6+I]+u[14+I]*IM*u[9+I]*E[4]*u[6+I]+u[15+I]*IM*u[9+I]*E[5]*u[6+I]-u[19+I]*IM*u[9+I]*E[3]*u[0+I]-u[20+I]*IM*u[9+I]*E[4]*u[0+I]-u[21+I]*IM*u[9+I]*E[5]*u[0+I]-u[13+I]*IM*u[9+I]*E[0]*u[5+I]-u[14+I]*IM*u[9+I]*E[1]*u[5+I]-u[15+I]*IM*u[9+I]*E[2]*u[5+I]+u[19+I]*IM*u[9+I]*E[0]*u[3+I]+u[20+I]*IM*u[9+I]*E[1]*u[3+I]+u[21+I]*IM*u[9+I]*E[2]*u[3+I]-u[16+I]*IM*u[9+I]*E[3]*u[5+I]-u[17+I]*IM*u[9+I]*E[4]*u[5+I]-u[18+I]*IM*u[9+I]*E[5]*u[5+I]+u[19+I]*IM*u[9+I]*E[3]*u[4+I]+u[20+I]*IM*u[9+I]*E[4]*u[4+I]+u[21+I]*IM*u[9+I]*E[5]*u[4+I]+u[8+I]*u[22+I]-u[2+I]*u[24+I]-u[4+I]*u[22+I]+u[1+I]*u[23+I];\
u_rhs[2+I]=u[12+I]*DZ(u,2+I,dz)+u[11+I]*DY(u,2+I,dy)+u[10+I]*DX(u,2+I,dx)+IM*u[9+I]*E[2]*DZ(u,1+I,dz)+IM*u[9+I]*E[1]*DY(u,1+I,dy)+IM*u[9+I]*E[0]*DX(u,1+I,dx)-IM*u[9+I]*E[5]*DZ(u,0+I,dz)-IM*u[9+I]*E[4]*DY(u,0+I,dy)-IM*u[9+I]*E[3]*DX(u,0+I,dx)+u[16+I]*IM*u[9+I]*E[0]*u[7+I]+u[17+I]*IM*u[9+I]*E[1]*u[7+I]+u[18+I]*IM*u[9+I]*E[2]*u[7+I]-u[19+I]*IM*u[9+I]*E[0]*u[4+I]-u[20+I]*IM*u[9+I]*E[1]*u[4+I]-u[21+I]*IM*u[9+I]*E[2]*u[4+I]-u[16+I]*IM*u[9+I]*E[3]*u[6+I]-u[17+I]*IM*u[9+I]*E[4]*u[6+I]-u[18+I]*IM*u[9+I]*E[5]*u[6+I]+u[19+I]*IM*u[9+I]*E[3]*u[3+I]+u[20+I]*IM*u[9+I]*E[4]*u[3+I]+u[21+I]*IM*u[9+I]*E[5]*u[3+I]-u[13+I]*IM*u[9+I]*E[0]*u[2+I]-u[14+I]*IM*u[9+I]*E[1]*u[2+I]-u[15+I]*IM*u[9+I]*E[2]*u[2+I]+u[19+I]*IM*u[9+I]*E[0]*u[0+I]+u[20+I]*IM*u[9+I]*E[1]*u[0+I]+u[21+I]*IM*u[9+I]*E[2]*u[0+I]-u[16+I]*IM*u[9+I]*E[3]*u[2+I]-u[17+I]*IM*u[9+I]*E[4]*u[2+I]-u[18+I]*IM*u[9+I]*E[5]*u[2+I]+u[19+I]*IM*u[9+I]*E[3]*u[1+I]+u[20+I]*IM*u[9+I]*E[4]*u[1+I]+u[21+I]*IM*u[9+I]*E[5]*u[1+I]-u[8+I]*u[23+I]+u[5+I]*u[24+I]-u[3+I]*u[22+I]+u[0+I]*u[23+I];\
u_rhs[7+I]=u[12+I]*DZ(u,7+I,dz)+u[11+I]*DY(u,7+I,dy)+u[10+I]*DX(u,7+I,dx)-IM*u[9+I]*E[2]*DZ(u,8+I,dz)-IM*u[9+I]*E[1]*DY(u,8+I,dy)-IM*u[9+I]*E[0]*DX(u,8+I,dx)+IM*u[9+I]*E[8]*DZ(u,6+I,dz)+IM*u[9+I]*E[7]*DY(u,6+I,dy)+IM*u[9+I]*E[6]*DX(u,6+I,dx)-u[13+I]*IM*u[9+I]*E[0]*u[5+I]-u[14+I]*IM*u[9+I]*E[1]*u[5+I]-u[15+I]*IM*u[9+I]*E[2]*u[5+I]+u[16+I]*IM*u[9+I]*E[0]*u[2+I]+u[17+I]*IM*u[9+I]*E[1]*u[2+I]+u[18+I]*IM*u[9+I]*E[2]*u[2+I]+u[13+I]*IM*u[9+I]*E[6]*u[3+I]+u[14+I]*IM*u[9+I]*E[7]*u[3+I]+u[15+I]*IM*u[9+I]*E[8]*u[3+I]-u[16+I]*IM*u[9+I]*E[6]*u[0+I]-u[17+I]*IM*u[9+I]*E[7]*u[0+I]-u[18+I]*IM*u[9+I]*E[8]*u[0+I]-u[13+I]*IM*u[9+I]*E[0]*u[7+I]-u[14+I]*IM*u[9+I]*E[1]*u[7+I]-u[15+I]*IM*u[9+I]*E[2]*u[7+I]+u[16+I]*IM*u[9+I]*E[0]*u[6+I]+u[17+I]*IM*u[9+I]*E[1]*u[6+I]+u[18+I]*IM*u[9+I]*E[2]*u[6+I]+u[16+I]*IM*u[9+I]*E[6]*u[8+I]+u[17+I]*IM*u[9+I]*E[7]*u[8+I]+u[18+I]*IM*u[9+I]*E[8]*u[8+I]-u[19+I]*IM*u[9+I]*E[6]*u[7+I]-u[20+I]*IM*u[9+I]*E[7]*u[7+I]-u[21+I]*IM*u[9+I]*E[8]*u[7+I]-u[4+I]*u[22+I]+u[1+I]*u[23+I]+u[8+I]*u[22+I]-u[2+I]*u[24+I];\
u_rhs[4+I]=u[12+I]*DZ(u,4+I,dz)+u[11+I]*DY(u,4+I,dy)+u[10+I]*DX(u,4+I,dx)-IM*u[9+I]*E[2]*DZ(u,5+I,dz)-IM*u[9+I]*E[1]*DY(u,5+I,dy)-IM*u[9+I]*E[0]*DX(u,5+I,dx)+IM*u[9+I]*E[8]*DZ(u,3+I,dz)+IM*u[9+I]*E[7]*DY(u,3+I,dy)+IM*u[9+I]*E[6]*DX(u,3+I,dx)+u[13+I]*IM*u[9+I]*E[0]*u[8+I]+u[14+I]*IM*u[9+I]*E[1]*u[8+I]+u[15+I]*IM*u[9+I]*E[2]*u[8+I]-u[19+I]*IM*u[9+I]*E[0]*u[2+I]-u[20+I]*IM*u[9+I]*E[1]*u[2+I]-u[21+I]*IM*u[9+I]*E[2]*u[2+I]-u[13+I]*IM*u[9+I]*E[6]*u[6+I]-u[14+I]*IM*u[9+I]*E[7]*u[6+I]-u[15+I]*IM*u[9+I]*E[8]*u[6+I]+u[19+I]*IM*u[9+I]*E[6]*u[0+I]+u[20+I]*IM*u[9+I]*E[7]*u[0+I]+u[21+I]*IM*u[9+I]*E[8]*u[0+I]-u[13+I]*IM*u[9+I]*E[0]*u[4+I]-u[14+I]*IM*u[9+I]*E[1]*u[4+I]-u[15+I]*IM*u[9+I]*E[2]*u[4+I]+u[16+I]*IM*u[9+I]*E[0]*u[3+I]+u[17+I]*IM*u[9+I]*E[1]*u[3+I]+u[18+I]*IM*u[9+I]*E[2]*u[3+I]+u[16+I]*IM*u[9+I]*E[6]*u[5+I]+u[17+I]*IM*u[9+I]*E[7]*u[5+I]+u[18+I]*IM*u[9+I]*E[8]*u[5+I]-u[19+I]*IM*u[9+I]*E[6]*u[4+I]-u[20+I]*IM*u[9+I]*E[7]*u[4+I]-u[21+I]*IM*u[9+I]*E[8]*u[4+I]+u[2+I]*u[7+I]*u[22+I]-u[2+I]*u[1+I]*u[24+I];\
u_rhs[1+I]=u[12+I]*DZ(u,1+I,dz)+u[11+I]*DY(u,1+I,dy)+u[10+I]*DX(u,1+I,dx)-IM*u[9+I]*E[2]*DZ(u,2+I,dz)-IM*u[9+I]*E[1]*DY(u,2+I,dy)-IM*u[9+I]*E[0]*DX(u,2+I,dx)+IM*u[9+I]*E[8]*DZ(u,0+I,dz)+IM*u[9+I]*E[7]*DY(u,0+I,dy)+IM*u[9+I]*E[6]*DX(u,0+I,dx)-u[16+I]*IM*u[9+I]*E[0]*u[8+I]-u[17+I]*IM*u[9+I]*E[1]*u[8+I]-u[18+I]*IM*u[9+I]*E[2]*u[8+I]+u[19+I]*IM*u[9+I]*E[0]*u[5+I]+u[20+I]*IM*u[9+I]*E[1]*u[5+I]+u[21+I]*IM*u[9+I]*E[2]*u[5+I]+u[16+I]*IM*u[9+I]*E[6]*u[6+I]+u[17+I]*IM*u[9+I]*E[7]*u[6+I]+u[18+I]*IM*u[9+I]*E[8]*u[6+I]-u[19+I]*IM*u[9+I]*E[6]*u[3+I]-u[20+I]*IM*u[9+I]*E[7]*u[3+I]-u[21+I]*IM*u[9+I]*E[8]*u[3+I]-u[13+I]*IM*u[9+I]*E[0]*u[1+I]-u[14+I]*IM*u[9+I]*E[1]*u[1+I]-u[15+I]*IM*u[9+I]*E[2]*u[1+I]+u[16+I]*IM*u[9+I]*E[0]*u[0+I]+u[17+I]*IM*u[9+I]*E[1]*u[0+I]+u[18+I]*IM*u[9+I]*E[2]*u[0+I]+u[16+I]*IM*u[9+I]*E[6]*u[2+I]+u[17+I]*IM*u[9+I]*E[7]*u[2+I]+u[18+I]*IM*u[9+I]*E[8]*u[2+I]-u[19+I]*IM*u[9+I]*E[6]*u[1+I]-u[20+I]*IM*u[9+I]*E[7]*u[1+I]-u[21+I]*IM*u[9+I]*E[8]*u[1+I]-u[7+I]*u[23+I]+u[4+I]*u[24+I]+u[6+I]*u[22+I]-u[0+I]*u[24+I];\
u_rhs[6+I]=u[12+I]*DZ(u,6+I,dz)+u[11+I]*DY(u,6+I,dy)+u[10+I]*DX(u,6+I,dx)+IM*u[9+I]*E[5]*DZ(u,8+I,dz)+IM*u[9+I]*E[4]*DY(u,8+I,dy)+IM*u[9+I]*E[3]*DX(u,8+I,dx)-IM*u[9+I]*E[8]*DZ(u,7+I,dz)-IM*u[9+I]*E[7]*DY(u,7+I,dy)-IM*u[9+I]*E[6]*DX(u,7+I,dx)+u[13+I]*IM*u[9+I]*E[3]*u[5+I]+u[14+I]*IM*u[9+I]*E[4]*u[5+I]+u[15+I]*IM*u[9+I]*E[5]*u[5+I]-u[16+I]*IM*u[9+I]*E[3]*u[2+I]-u[17+I]*IM*u[9+I]*E[4]*u[2+I]-u[18+I]*IM*u[9+I]*E[5]*u[2+I]-u[13+I]*IM*u[9+I]*E[6]*u[4+I]-u[14+I]*IM*u[9+I]*E[7]*u[4+I]-u[15+I]*IM*u[9+I]*E[8]*u[4+I]+u[16+I]*IM*u[9+I]*E[6]*u[1+I]+u[17+I]*IM*u[9+I]*E[7]*u[1+I]+u[18+I]*IM*u[9+I]*E[8]*u[1+I]+u[13+I]*IM*u[9+I]*E[3]*u[7+I]+u[14+I]*IM*u[9+I]*E[4]*u[7+I]+u[15+I]*IM*u[9+I]*E[5]*u[7+I]-u[16+I]*IM*u[9+I]*E[3]*u[6+I]-u[17+I]*IM*u[9+I]*E[4]*u[6+I]-u[18+I]*IM*u[9+I]*E[5]*u[6+I]+u[13+I]*IM*u[9+I]*E[6]*u[8+I]+u[14+I]*IM*u[9+I]*E[7]*u[8+I]+u[15+I]*IM*u[9+I]*E[8]*u[8+I]-u[19+I]*IM*u[9+I]*E[6]*u[6+I]-u[20+I]*IM*u[9+I]*E[7]*u[6+I]-u[21+I]*IM*u[9+I]*E[8]*u[6+I]-u[3+I]*u[22+I]+u[0+I]*u[23+I]-u[8+I]*u[23+I]+u[5+I]*u[24+I];\
u_rhs[3+I]=u[12+I]*DZ(u,3+I,dz)+u[11+I]*DY(u,3+I,dy)+u[10+I]*DX(u,3+I,dx)+IM*u[9+I]*E[5]*DZ(u,5+I,dz)+IM*u[9+I]*E[4]*DY(u,5+I,dy)+IM*u[9+I]*E[3]*DX(u,5+I,dx)-IM*u[9+I]*E[8]*DZ(u,4+I,dz)-IM*u[9+I]*E[7]*DY(u,4+I,dy)-IM*u[9+I]*E[6]*DX(u,4+I,dx)-u[13+I]*IM*u[9+I]*E[3]*u[8+I]-u[14+I]*IM*u[9+I]*E[4]*u[8+I]-u[15+I]*IM*u[9+I]*E[5]*u[8+I]+u[19+I]*IM*u[9+I]*E[3]*u[2+I]+u[20+I]*IM*u[9+I]*E[4]*u[2+I]+u[21+I]*IM*u[9+I]*E[5]*u[2+I]+u[13+I]*IM*u[9+I]*E[6]*u[7+I]+u[14+I]*IM*u[9+I]*E[7]*u[7+I]+u[15+I]*IM*u[9+I]*E[8]*u[7+I]-u[19+I]*IM*u[9+I]*E[6]*u[1+I]-u[20+I]*IM*u[9+I]*E[7]*u[1+I]-u[21+I]*IM*u[9+I]*E[8]*u[1+I]+u[13+I]*IM*u[9+I]*E[3]*u[4+I]+u[14+I]*IM*u[9+I]*E[4]*u[4+I]+u[15+I]*IM*u[9+I]*E[5]*u[4+I]-u[16+I]*IM*u[9+I]*E[3]*u[3+I]-u[17+I]*IM*u[9+I]*E[4]*u[3+I]-u[18+I]*IM*u[9+I]*E[5]*u[3+I]+u[13+I]*IM*u[9+I]*E[6]*u[5+I]+u[14+I]*IM*u[9+I]*E[7]*u[5+I]+u[15+I]*IM*u[9+I]*E[8]*u[5+I]-u[19+I]*IM*u[9+I]*E[6]*u[3+I]-u[20+I]*IM*u[9+I]*E[7]*u[3+I]-u[21+I]*IM*u[9+I]*E[8]*u[3+I]+u[6+I]*u[22+I]-u[0+I]*u[24+I]-u[7+I]*u[23+I]+u[4+I]*u[24+I];\
u_rhs[0+I]=u[12+I]*DZ(u,0+I,dz)+u[11+I]*DY(u,0+I,dy)+u[10+I]*DX(u,0+I,dx)+IM*u[9+I]*E[5]*DZ(u,2+I,dz)+IM*u[9+I]*E[4]*DY(u,2+I,dy)+IM*u[9+I]*E[3]*DX(u,2+I,dx)-IM*u[9+I]*E[8]*DZ(u,1+I,dz)-IM*u[9+I]*E[7]*DY(u,1+I,dy)-IM*u[9+I]*E[6]*DX(u,1+I,dx)+u[16+I]*IM*u[9+I]*E[3]*u[8+I]+u[17+I]*IM*u[9+I]*E[4]*u[8+I]+u[18+I]*IM*u[9+I]*E[5]*u[8+I]-u[19+I]*IM*u[9+I]*E[3]*u[5+I]-u[20+I]*IM*u[9+I]*E[4]*u[5+I]-u[21+I]*IM*u[9+I]*E[5]*u[5+I]-u[16+I]*IM*u[9+I]*E[6]*u[7+I]-u[17+I]*IM*u[9+I]*E[7]*u[7+I]-u[18+I]*IM*u[9+I]*E[8]*u[7+I]+u[19+I]*IM*u[9+I]*E[6]*u[4+I]+u[20+I]*IM*u[9+I]*E[7]*u[4+I]+u[21+I]*IM*u[9+I]*E[8]*u[4+I]+u[13+I]*IM*u[9+I]*E[3]*u[1+I]+u[14+I]*IM*u[9+I]*E[4]*u[1+I]+u[15+I]*IM*u[9+I]*E[5]*u[1+I]-u[16+I]*IM*u[9+I]*E[3]*u[0+I]-u[17+I]*IM*u[9+I]*E[4]*u[0+I]-u[18+I]*IM*u[9+I]*E[5]*u[0+I]+u[13+I]*IM*u[9+I]*E[6]*u[2+I]+u[14+I]*IM*u[9+I]*E[7]*u[2+I]+u[15+I]*IM*u[9+I]*E[8]*u[2+I]-u[19+I]*IM*u[9+I]*E[6]*u[0+I]-u[20+I]*IM*u[9+I]*E[7]*u[0+I]-u[21+I]*IM*u[9+I]*E[8]*u[0+I]-u[2+I]*u[6+I]*u[23+I]+u[2+I]*u[3+I]*u[24+I];\
u_rhs[21+I]=DZ(u,24+I,dz)+u[15+I]*u[23+I]-u[18+I]*u[22+I]+Base(0.5)*IM*u[9+I]*u[6+I]*E[3]*E[7]-Base(0.5)*IM*u[9+I]*u[6+I]*E[4]*E[6]-Base(0.5)*IM*u[9+I]*u[6+I]*E[6]*E[4]+Base(0.5)*IM*u[9+I]*u[6+I]*E[7]*E[3]-Base(0.5)*IM*u[9+I]*u[7+I]*E[0]*E[7]+Base(0.5)*IM*u[9+I]*u[7+I]*E[1]*E[6]+Base(0.5)*IM*u[9+I]*u[7+I]*E[6]*E[1]-Base(0.5)*IM*u[9+I]*u[7+I]*E[7]*E[0]+Base(0.5)*IM*u[9+I]*u[8+I]*E[0]*E[4]-Base(0.5)*IM*u[9+I]*u[8+I]*E[1]*E[3]-Base(0.5)*IM*u[9+I]*u[8+I]*E[3]*E[1]+Base(0.5)*IM*u[9+I]*u[8+I]*E[4]*E[0]-u[10+I]*F[7]+u[11+I]*F[6];\
u_rhs[20+I]=DY(u,24+I,dy)+u[14+I]*u[23+I]-u[17+I]*u[22+I]-Base(0.5)*IM*u[9+I]*u[6+I]*E[3]*E[8]+Base(0.5)*IM*u[9+I]*u[6+I]*E[5]*E[6]+Base(0.5)*IM*u[9+I]*u[6+I]*E[6]*E[5]-Base(0.5)*IM*u[9+I]*u[6+I]*E[8]*E[3]+Base(0.5)*IM*u[9+I]*u[7+I]*E[0]*E[8]-Base(0.5)*IM*u[9+I]*u[7+I]*E[2]*E[6]-Base(0.5)*IM*u[9+I]*u[7+I]*E[6]*E[2]+Base(0.5)*IM*u[9+I]*u[7+I]*E[8]*E[0]-Base(0.5)*IM*u[9+I]*u[8+I]*E[0]*E[5]+Base(0.5)*IM*u[9+I]*u[8+I]*E[2]*E[3]+Base(0.5)*IM*u[9+I]*u[8+I]*E[3]*E[2]-Base(0.5)*IM*u[9+I]*u[8+I]*E[5]*E[0]+u[10+I]*F[8]-u[12+I]*F[6];\
u_rhs[19+I]=DX(u,24+I,dx)+u[13+I]*u[23+I]-u[16+I]*u[22+I]+Base(0.5)*IM*u[9+I]*u[6+I]*E[4]*E[8]-Base(0.5)*IM*u[9+I]*u[6+I]*E[5]*E[7]-Base(0.5)*IM*u[9+I]*u[6+I]*E[7]*E[5]+Base(0.5)*IM*u[9+I]*u[6+I]*E[8]*E[4]-Base(0.5)*IM*u[9+I]*u[7+I]*E[1]*E[8]+Base(0.5)*IM*u[9+I]*u[7+I]*E[2]*E[7]+Base(0.5)*IM*u[9+I]*u[7+I]*E[7]*E[2]-Base(0.5)*IM*u[9+I]*u[7+I]*E[8]*E[1]+Base(0.5)*IM*u[9+I]*u[8+I]*E[1]*E[5]-Base(0.5)*IM*u[9+I]*u[8+I]*E[2]*E[4]-Base(0.5)*IM*u[9+I]*u[8+I]*E[4]*E[2]+Base(0.5)*IM*u[9+I]*u[8+I]*E[5]*E[1]-u[11+I]*F[8]+u[12+I]*F[7];\
u_rhs[18+I]=DZ(u,23+I,dz)-u[15+I]*u[24+I]+u[21+I]*u[22+I]+Base(0.5)*IM*u[9+I]*u[3+I]*E[3]*E[7]-Base(0.5)*IM*u[9+I]*u[3+I]*E[4]*E[6]-Base(0.5)*IM*u[9+I]*u[3+I]*E[6]*E[4]+Base(0.5)*IM*u[9+I]*u[3+I]*E[7]*E[3]-Base(0.5)*IM*u[9+I]*u[4+I]*E[0]*E[7]+Base(0.5)*IM*u[9+I]*u[4+I]*E[1]*E[6]+Base(0.5)*IM*u[9+I]*u[4+I]*E[6]*E[1]-Base(0.5)*IM*u[9+I]*u[4+I]*E[7]*E[0]+Base(0.5)*IM*u[9+I]*u[5+I]*E[0]*E[4]-Base(0.5)*IM*u[9+I]*u[5+I]*E[1]*E[3]-Base(0.5)*IM*u[9+I]*u[5+I]*E[3]*E[1]+Base(0.5)*IM*u[9+I]*u[5+I]*E[4]*E[0]-u[10+I]*F[4]+u[11+I]*F[3];\
u_rhs[17+I]=DY(u,23+I,dy)-u[14+I]*u[24+I]+u[20+I]*u[22+I]-Base(0.5)*IM*u[9+I]*u[3+I]*E[3]*E[8]+Base(0.5)*IM*u[9+I]*u[3+I]*E[5]*E[6]+Base(0.5)*IM*u[9+I]*u[3+I]*E[6]*E[5]-Base(0.5)*IM*u[9+I]*u[3+I]*E[8]*E[3]+Base(0.5)*IM*u[9+I]*u[4+I]*E[0]*E[8]-Base(0.5)*IM*u[9+I]*u[4+I]*E[2]*E[6]-Base(0.5)*IM*u[9+I]*u[4+I]*E[6]*E[2]+Base(0.5)*IM*u[9+I]*u[4+I]*E[8]*E[0]-Base(0.5)*IM*u[9+I]*u[5+I]*E[0]*E[5]+Base(0.5)*IM*u[9+I]*u[5+I]*E[2]*E[3]+Base(0.5)*IM*u[9+I]*u[5+I]*E[3]*E[2]-Base(0.5)*IM*u[9+I]*u[5+I]*E[5]*E[0]+u[10+I]*F[5]-u[12+I]*F[3];\
u_rhs[16+I]=DX(u,23+I,dx)-u[13+I]*u[24+I]+u[19+I]*u[22+I]+Base(0.5)*IM*u[9+I]*u[3+I]*E[4]*E[8]-Base(0.5)*IM*u[9+I]*u[3+I]*E[5]*E[7]-Base(0.5)*IM*u[9+I]*u[3+I]*E[7]*E[5]+Base(0.5)*IM*u[9+I]*u[3+I]*E[8]*E[4]-Base(0.5)*IM*u[9+I]*u[4+I]*E[1]*E[8]+Base(0.5)*IM*u[9+I]*u[4+I]*E[2]*E[7]+Base(0.5)*IM*u[9+I]*u[4+I]*E[7]*E[2]-Base(0.5)*IM*u[9+I]*u[4+I]*E[8]*E[1]+Base(0.5)*IM*u[9+I]*u[5+I]*E[1]*E[5]-Base(0.5)*IM*u[9+I]*u[5+I]*E[2]*E[4]-Base(0.5)*IM*u[9+I]*u[5+I]*E[4]*E[2]+Base(0.5)*IM*u[9+I]*u[5+I]*E[5]*E[1]-u[11+I]*F[5]+u[12+I]*F[4];\
u_rhs[15+I]=DZ(u,22+I,dz)+u[18+I]*u[24+I]-u[21+I]*u[23+I]+Base(0.5)*IM*u[9+I]*u[0+I]*E[3]*E[7]-Base(0.5)*IM*u[9+I]*u[0+I]*E[4]*E[6]-Base(0.5)*IM*u[9+I]*u[0+I]*E[6]*E[4]+Base(0.5)*IM*u[9+I]*u[0+I]*E[7]*E[3]-Base(0.5)*IM*u[9+I]*u[1+I]*E[0]*E[7]+Base(0.5)*IM*u[9+I]*u[1+I]*E[1]*E[6]+Base(0.5)*IM*u[9+I]*u[1+I]*E[6]*E[1]-Base(0.5)*IM*u[9+I]*u[1+I]*E[7]*E[0]+Base(0.5)*IM*u[9+I]*u[2+I]*E[0]*E[4]-Base(0.5)*IM*u[9+I]*u[2+I]*E[1]*E[3]-Base(0.5)*IM*u[9+I]*u[2+I]*E[3]*E[1]+Base(0.5)*IM*u[9+I]*u[2+I]*E[4]*E[0]-u[10+I]*F[1]+u[11+I]*F[0];\
u_rhs[14+I]=DY(u,22+I,dy)+u[17+I]*u[24+I]-u[20+I]*u[23+I]-Base(0.5)*IM*u[9+I]*u[0+I]*E[3]*E[8]+Base(0.5)*IM*u[9+I]*u[0+I]*E[5]*E[6]+Base(0.5)*IM*u[9+I]*u[0+I]*E[6]*E[5]-Base(0.5)*IM*u[9+I]*u[0+I]*E[8]*E[3]+Base(0.5)*IM*u[9+I]*u[1+I]*E[0]*E[8]-Base(0.5)*IM*u[9+I]*u[1+I]*E[2]*E[6]-Base(0.5)*IM*u[9+I]*u[1+I]*E[6]*E[2]+Base(0.5)*IM*u[9+I]*u[1+I]*E[8]*E[0]-Base(0.5)*IM*u[9+I]*u[2+I]*E[0]*E[5]+Base(0.5)*IM*u[9+I]*u[2+I]*E[2]*E[3]+Base(0.5)*IM*u[9+I]*u[2+I]*E[3]*E[2]-Base(0.5)*IM*u[9+I]*u[2+I]*E[5]*E[0]+u[10+I]*F[2]-u[12+I]*F[0];\
u_rhs[13+I]=DX(u,22+I,dx)+u[16+I]*u[24+I]-u[19+I]*u[23+I]+Base(0.5)*IM*u[9+I]*u[0+I]*E[4]*E[8]-Base(0.5)*IM*u[9+I]*u[0+I]*E[5]*E[7]-Base(0.5)*IM*u[9+I]*u[0+I]*E[7]*E[5]+Base(0.5)*IM*u[9+I]*u[0+I]*E[8]*E[4]-Base(0.5)*IM*u[9+I]*u[1+I]*E[1]*E[8]+Base(0.5)*IM*u[9+I]*u[1+I]*E[2]*E[7]+Base(0.5)*IM*u[9+I]*u[1+I]*E[7]*E[2]-Base(0.5)*IM*u[9+I]*u[1+I]*E[8]*E[1]+Base(0.5)*IM*u[9+I]*u[2+I]*E[1]*E[5]-Base(0.5)*IM*u[9+I]*u[2+I]*E[2]*E[4]-Base(0.5)*IM*u[9+I]*u[2+I]*E[4]*E[2]+Base(0.5)*IM*u[9+I]*u[2+I]*E[5]*E[1]-u[11+I]*F[2]+u[12+I]*F[1];

#define WEYL_CONSTRAINT(G,DX,DY,DZ) G[2]=E[8]*DZ(u,8+I,dz)+E[7]*DY(u,8+I,dy)+E[6]*DX(u,8+I,dx)+E[5]*DZ(u,7+I,dz)+E[4]*DY(u,7+I,dy)+E[3]*DX(u,7+I,dx)+E[2]*DZ(u,6+I,dz)+E[1]*DY(u,6+I,dy)+E[0]*DX(u,6+I,dx)+u[13+I]*E[0]*u[3+I]+u[13+I]*E[3]*u[4+I]+u[13+I]*E[6]*u[5+I]+u[14+I]*E[1]*u[3+I]+u[14+I]*E[4]*u[4+I]+u[14+I]*E[7]*u[5+I]+u[15+I]*E[2]*u[3+I]+u[15+I]*E[5]*u[4+I]+u[15+I]*E[8]*u[5+I]-u[16+I]*E[0]*u[0+I]-u[16+I]*E[3]*u[1+I]-u[16+I]*E[6]*u[2+I]-u[17+I]*E[1]*u[0+I]-u[17+I]*E[4]*u[1+I]-u[17+I]*E[7]*u[2+I]-u[18+I]*E[2]*u[0+I]-u[18+I]*E[5]*u[1+I]-u[18+I]*E[8]*u[2+I]+u[13+I]*E[6]*u[7+I]+u[14+I]*E[7]*u[7+I]+u[15+I]*E[8]*u[7+I]-u[13+I]*E[3]*u[8+I]-u[14+I]*E[4]*u[8+I]-u[15+I]*E[5]*u[8+I]-u[16+I]*E[6]*u[6+I]-u[17+I]*E[7]*u[6+I]-u[18+I]*E[8]*u[6+I]+u[16+I]*E[0]*u[8+I]+u[17+I]*E[1]*u[8+I]+u[18+I]*E[2]*u[8+I]+u[19+I]*E[3]*u[6+I]+u[20+I]*E[4]*u[6+I]+u[21+I]*E[5]*u[6+I]-u[19+I]*E[0]*u[7+I]-u[20+I]*E[1]*u[7+I]-u[21+I]*E[2]*u[7+I];\
G[1]=E[8]*DZ(u,5+I,dz)+E[7]*DY(u,5+I,dy)+E[6]*DX(u,5+I,dx)+E[5]*DZ(u,4+I,dz)+E[4]*DY(u,4+I,dy)+E[3]*DX(u,4+I,dx)+E[2]*DZ(u,3+I,dz)+E[1]*DY(u,3+I,dy)+E[0]*DX(u,3+I,dx)-u[13+I]*E[0]*u[6+I]-u[13+I]*E[3]*u[7+I]-u[13+I]*E[6]*u[8+I]-u[14+I]*E[1]*u[6+I]-u[14+I]*E[4]*u[7+I]-u[14+I]*E[7]*u[8+I]-u[15+I]*E[2]*u[6+I]-u[15+I]*E[5]*u[7+I]-u[15+I]*E[8]*u[8+I]+u[19+I]*E[0]*u[0+I]+u[19+I]*E[3]*u[1+I]+u[19+I]*E[6]*u[2+I]+u[20+I]*E[1]*u[0+I]+u[20+I]*E[4]*u[1+I]+u[20+I]*E[7]*u[2+I]+u[21+I]*E[2]*u[0+I]+u[21+I]*E[5]*u[1+I]+u[21+I]*E[8]*u[2+I]+u[13+I]*E[6]*u[4+I]+u[14+I]*E[7]*u[4+I]+u[15+I]*E[8]*u[4+I]-u[13+I]*E[3]*u[5+I]-u[14+I]*E[4]*u[5+I]-u[15+I]*E[5]*u[5+I]-u[16+I]*E[6]*u[3+I]-u[17+I]*E[7]*u[3+I]-u[18+I]*E[8]*u[3+I]+u[16+I]*E[0]*u[5+I]+u[17+I]*E[1]*u[5+I]+u[18+I]*E[2]*u[5+I]+u[19+I]*E[3]*u[3+I]+u[20+I]*E[4]*u[3+I]+u[21+I]*E[5]*u[3+I]-u[19+I]*E[0]*u[4+I]-u[20+I]*E[1]*u[4+I]-u[21+I]*E[2]*u[4+I];\
G[0]=E[8]*DZ(u,2+I,dz)+E[7]*DY(u,2+I,dy)+E[6]*DX(u,2+I,dx)+E[5]*DZ(u,1+I,dz)+E[4]*DY(u,1+I,dy)+E[3]*DX(u,1+I,dx)+E[2]*DZ(u,0+I,dz)+E[1]*DY(u,0+I,dy)+E[0]*DX(u,0+I,dx)+u[16+I]*E[0]*u[6+I]+u[16+I]*E[3]*u[7+I]+u[16+I]*E[6]*u[8+I]+u[17+I]*E[1]*u[6+I]+u[17+I]*E[4]*u[7+I]+u[17+I]*E[7]*u[8+I]+u[18+I]*E[2]*u[6+I]+u[18+I]*E[5]*u[7+I]+u[18+I]*E[8]*u[8+I]-u[19+I]*E[0]*u[3+I]-u[19+I]*E[3]*u[4+I]-u[19+I]*E[6]*u[5+I]-u[20+I]*E[1]*u[3+I]-u[20+I]*E[4]*u[4+I]-u[20+I]*E[7]*u[5+I]-u[21+I]*E[2]*u[3+I]-u[21+I]*E[5]*u[4+I]-u[21+I]*E[8]*u[5+I]+u[13+I]*E[6]*u[1+I]+u[14+I]*E[7]*u[1+I]+u[15+I]*E[8]*u[1+I]-u[13+I]*E[3]*u[2+I]-u[14+I]*E[4]*u[2+I]-u[15+I]*E[5]*u[2+I]-u[16+I]*E[6]*u[0+I]-u[17+I]*E[7]*u[0+I]-u[18+I]*E[8]*u[0+I]+u[16+I]*E[0]*u[2+I]+u[17+I]*E[1]*u[2+I]+u[18+I]*E[2]*u[2+I]+u[19+I]*E[3]*u[0+I]+u[20+I]*E[4]*u[0+I]+u[21+I]*E[5]*u[0+I]-u[19+I]*E[0]*u[1+I]-u[20+I]*E[1]*u[1+I]-u[21+I]*E[2]*u[1+I];

#define CANONICAL_CONSTRAINT(G) FOR(n,0,NG) {G[0] += u[n+I]-u_0[n+I];}

#define INDXMAT(i,j) (3*(i-1) + (j-1))
#define INDVEC(i) (i-1)

#define SYM_PSI(u,I) u[0+I+INDXMAT(1,2)] = u[0+I+INDXMAT(1,2)]+u[0+I+INDXMAT(2,1)]/Base(2.0);\
					 u[0+I+INDXMAT(2,1)] = u[0+I+INDXMAT(1,2)];\
					 u[0+I+INDXMAT(2,3)] = u[0+I+INDXMAT(2,3)]+u[0+I+INDXMAT(3,2)]/Base(2.0);\
					 u[0+I+INDXMAT(3,2)] = u[0+I+INDXMAT(2,3)];\
					 u[0+I+INDXMAT(1,3)] = u[0+I+INDXMAT(1,3)]+u[0+I+INDXMAT(3,1)]/Base(2.0);\
					 u[0+I+INDXMAT(1,3)] = u[0+I+INDXMAT(3,1)]; \
					 Field trace = (u[0+I+INDXMAT(1,1)] + u[0+I+INDXMAT(2,2)] + u[0+I+INDXMAT(3,3)])/Base(3.0);\
					 u[0+INDXMAT(1,1)] = u[0+INDXMAT(1,1)] - trace;\
					 u[0+INDXMAT(2,2)] = u[0+INDXMAT(2,2)] - trace;\
					 u[0+INDXMAT(3,3)] = u[0+INDXMAT(3,3)] - trace;

#define CDWX CDX
#define CDWY CDY
#define CDWZ CDZ


class Weyl_System : public SysParams {


	public:

	Base M = Base(1.0);
	Field* u_0;
	
	Weyl_System(Base LX0, Base LX1, Base LY0, Base LY1, Base LZ0, Base LZ1, Base M) : SysParams(LX0, LX1, LY0, LY1, LZ0, LZ1), M(M) {
		Name = "WO";
		u_0 = new Field[NF*NX*NY*NZ];
	}

	void Rhs_Single(Field* u, Field* u_rhs, int I) {
		Field* F = new Field[9];
		Field* E = new Field[9];
		F_DEF(u,I,CDWX,CDWY,CDWZ);
		E_DEF(u,I);

		WEYL_RHS(u,u_rhs,I,CDWX,CDWY,CDWZ);

		// SYM_PSI(u_rhs,I);

		delete[] F;
		delete[] E;
	}

	void Constraint_Single(Field* u, int I, Field* G_ret) {
		#if NG == 1
			CANONICAL_CONSTRAINT(G_ret)
		#else

		Field* F = new Field[9];
		Field* E = new Field[9];
		F_DEF(u,I,CDWX,CDWY,CDWZ);
		E_DEF(u,I);
		WEYL_CONSTRAINT(G_ret,CDWX,CDWY,CDWZ);

		delete[] F;
		delete[] E;

		#endif
	}

	void Init_Isotropic_Spherical(Field*& u) {
		#ifdef KOEPS
		Name += "IS-N-" + std::to_string(NNN) + "-E-" +std::to_string(int(1.0/std::abs(KOEPS)));
		#else
		Name += "IS-N-" + std::to_string(NNN);
		#endif

		std::fstream log_file;
		log_file.open("./Logs/"+Name+".log", std::fstream::in | std::fstream::out | std::fstream::trunc);
		log_file.close();
		log_file.open(("Data/"+Name+"-Solution"+".dat").c_str(), std::fstream::trunc);
		log_file.close();
		FOR(i,0,NX) FOR(j,0,NY) FOR(k,0,NZ) {
			Base r     = X[i];
			Base theta = Y[j];
			Base phi   = Z[k];

			FOR(n,0,NF) {u[n+INDX3D(i,j,k)] = 0;}

			// \Psi^{ij}
			u[0+INDXMAT(1,1)+INDX3D(i,j,k)] = -Base(4)*M*std::pow(r,3)/std::pow((r+M),6);
			u[0+INDXMAT(2,2)+INDX3D(i,j,k)] =  Base(2)*M*std::pow(r,3)/std::pow((r+M),6);
			u[0+INDXMAT(3,3)+INDX3D(i,j,k)] =  Base(2)*M*std::pow(r,3)/std::pow((r+M),6);

			// \utilde{N}
			u[9+INDX3D(i,j,k)] = (r-M)*std::pow(r,4)/(std::pow(r,3)*std::sin(theta));

			// N^a

			// A^i_a
			u[13+INDXMAT(2,3)+INDX3D(i,j,k)] = (M-r)/(M+r)*std::sin(theta);
			u[13+INDXMAT(3,2)+INDX3D(i,j,k)] = (r-M)/(r+M);
			u[13+INDXMAT(1,3)+INDX3D(i,j,k)] = std::cos(theta);


			// A^i_0
			u[22+INDVEC(1)+INDX3D(i,j,k)]    = -Base(2.0)*IM*M*r*r/std::pow(M+r,5);

			FOR(n,0,NF) {u_0[n+INDX3D(i,j,k)] = u[n+INDX3D(i,j,k)];}
		}
	}

	void Init_Isotropic_Euclidian(Field*& u) {
		#ifdef KOEPS
		Name += "IE-N-" + std::to_string(NNN) + "-E-" + std::to_string(int(1.0/std::abs(KOEPS)));
		#else
		Name += "IE-N-" + std::to_string(NNN);
		#endif

		std::fstream log_file;
		log_file.open("./Logs/"+Name+".log", std::fstream::in | std::fstream::out | std::fstream::trunc);
		log_file.close();
		log_file.open(("Data/"+Name+"-Solution"+".dat").c_str(), std::fstream::trunc);
		log_file.close();

		FOR(i,0,NX) FOR(j,0,NY) FOR(k,0,NZ) {
			Base x = X[i];
			Base y = Y[j];
			Base z = Z[k];
			Base r = std::sqrt(x*x+y*y+z*z);

			FOR(n,0,NF) {u[n+INDX3D(i,j,k)] = 0;}
			
			u[0+INDXMAT(1,1)+INDX3D(i,j,k)] = (2*M*r/std::pow(r+M,6))*(r*r-3*x*x);
			u[0+INDXMAT(1,2)+INDX3D(i,j,k)] = (2*M*r/std::pow(r+M,6))*(-3*x*y);
			u[0+INDXMAT(1,3)+INDX3D(i,j,k)] = (2*M*r/std::pow(r+M,6))*(-3*x*z);

			u[0+INDXMAT(2,1)+INDX3D(i,j,k)] = (2*M*r/std::pow(r+M,6))*(-3*y*x);
			u[0+INDXMAT(2,2)+INDX3D(i,j,k)] = (2*M*r/std::pow(r+M,6))*(r*r-3*y*y);
			u[0+INDXMAT(2,3)+INDX3D(i,j,k)] = (2*M*r/std::pow(r+M,6))*(-3*y*z);

			u[0+INDXMAT(3,1)+INDX3D(i,j,k)] = (2*M*r/std::pow(r+M,6))*(-3*z*x);
			u[0+INDXMAT(3,2)+INDX3D(i,j,k)] = (2*M*r/std::pow(r+M,6))*(-3*z*y);
			u[0+INDXMAT(3,3)+INDX3D(i,j,k)] = (2*M*r/std::pow(r+M,6))*(r*r-3*z*z);

			u[9+INDX3D(i,j,k)] = std::pow(r,6)*(r-M)/std::pow(r+M,7);

			u[13+INDXMAT(1,2)+INDX3D(i,j,k)] =  2*M*z/(r*r*(r+M));
			u[13+INDXMAT(1,3)+INDX3D(i,j,k)] = -2*M*y/(r*r*(r+M));
			u[13+INDXMAT(2,1)+INDX3D(i,j,k)] = -2*M*z/(r*r*(r+M));
			u[13+INDXMAT(2,3)+INDX3D(i,j,k)] =  2*M*x/(r*r*(r+M));
			u[13+INDXMAT(3,1)+INDX3D(i,j,k)] =  2*M*y/(r*r*(r+M));
			u[13+INDXMAT(3,2)+INDX3D(i,j,k)] = -2*M*x/(r*r*(r+M));

			u[22+INDVEC(1)+INDX3D(i,j,k)] = -Base(2)*IM*M*r*x/std::pow(r+M,4);
			u[22+INDVEC(2)+INDX3D(i,j,k)] = -Base(2)*IM*M*r*y/std::pow(r+M,4);
			u[22+INDVEC(3)+INDX3D(i,j,k)] = -Base(2)*IM*M*r*z/std::pow(r+M,4);

			FOR(n,0,NF) {u_0[n+INDX3D(i,j,k)] = u[n+INDX3D(i,j,k)];}
		}
	}

	Base Tau_Function(Field* u) {
		return dt;
		Field* F = new Field[9];
		Field* E = new Field[9];
		Field N_avg = Base(0.0);
		Field detE = Base(0.0);

		Base detE_max = Base(0.0);
		#pragma omp parallel for
		#if NX > 1
		FOR(i,1,NX-1) {
		#else
			int i = 0;
		#endif

		#if NY > 1
			FOR(j,1,NY-1) {
		#else
			int j = 0;
		#endif
		
		#if NZ > 1
				FOR(k,1,NZ-1) {
		#else
				int k = 0;
		#endif

			F_DEF(u,INDX3D(i,j,k),CDWX,CDWY,CDWY);
			E_DEF(u,INDX3D(i,j,k));
			
			detE = E[0]*E[4]*E[8]-E[0]*E[5]*E[7]-E[1]*E[3]*E[8]+E[1]*E[5]*E[6]+E[2]*E[3]*E[7]-E[2]*E[4]*E[6];	
			N_avg += u[9+INDX3D(i,j,k)]*std::sqrt(detE);
		
		#if NZ > 1
				}
		#endif

		#if NY > 1
			}
		#endif

		#if NX > 1
		}
		#endif

		N_avg /= Field(NX*NY*NZ);

		return std::abs(std::real(N_avg*dt));

	}

	void Apply_BCs(Field* u, Field* u_rhs) {
	}

	void Log_Solution(Field* u, Base Tau) {
		std::fstream log_file;
		log_file.open("./Logs/"+Name+".log", std::fstream::in | std::fstream::out | std::fstream::app);

		if (!log_file) {
			log_file.open("./Logs/"+Name+".log", std::fstream::in | std::fstream::out | std::fstream::trunc);
		}

		log_file << "\n\n\n########################################################################################\n";
		log_file << "Tau:" << Tau;
		FOR(i,0,NX) FOR(j,0,NY) FOR(k,0,NZ) {
			log_file << "r: " + std::to_string(X[i]) + "\n";

			log_file << "Psi: \n";
			FOR(x,1,4) {
				FOR(y,1,4) {
					log_file << u[0+INDXMAT(x,y)+INDX3D(i,j,k)] << ", "; 
				}
				log_file << "\n";
			}
			log_file << "\n";

			log_file << "A: \n";
			FOR(x,1,4) {
				FOR(y,1,4) {
					log_file << u[13+INDXMAT(x,y)+INDX3D(i,j,k)] << ", ";
				}
				log_file << "\n";
			}
			log_file << "\n";
		}
		log_file << "########################################################################################\n";
		log_file.close();
	}


};
