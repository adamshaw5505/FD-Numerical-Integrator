#pragma once
#define NNN  (2<<2)
#define NX NNN
#define NY NNN
#define NZ NNN
#define NF  4
#define NG  1

#define BETA Base(0.0)

// #define KOEPS Base(0)

#include "../Headers/Definitions.h"
#include "../Headers/Evolution.h"

#define MAXWELL_MOD_RHS(u,u_rhs,I,DX,DY,DZ,dx,dy,dz) u_rhs[0+I] = IM*DY(u,2+I,dy) - IM*DZ(u,1+I,dz) + DX(u,3+I,dx); \
										  u_rhs[1+I] = IM*DZ(u,0+I,dz) - IM*DX(u,2+I,dx) + DY(u,3+I,dy); \
										  u_rhs[2+I] = IM*DX(u,1+I,dx) - IM*DY(u,0+I,dy) + DZ(u,3+I,dz); \
										  u_rhs[3+I] = DX(u,0+I,dx) + DY(u,1+I,dy) + DZ(u,2+I,dz) - BETA * u[3+I];

class Maxwell_System : public SysParams {
	public:
	Maxwell_System(Base LX0, Base LX1, Base LY0, Base LY1, Base LZ0, Base LZ1) : SysParams(LX0, LX1, LY0, LY1, LZ0, LZ1) {
		Name = "MM";
	}

	void Apply_Periodic_BCs(Field* u, Field* u_rhs) {
		#pragma omp parallel for	
		FOR(i,1,NX-1) FOR(j,1,NY-1) {
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(i,j,0),   CDX,CDY,CDZ0_P,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(i,j,NZ-1),CDX,CDY,CDZN_P,dx,dy,dz)
		}
		#pragma omp parallel for
		FOR(j,1,NY-1) FOR(k,1,NZ-1) {
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(0,j,k),   CDX0_P,CDY,CDZ,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(NX-1,j,k),CDXN_P,CDY,CDZ,dx,dy,dz)
		}
		#pragma omp parallel for
		FOR(i,1,NX-1) FOR(k,1,NZ-1) {
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(i,0,k),   CDX,CDY0_P,CDZ,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(i,NY-1,k),CDX,CDYN_P,CDZ,dx,dy,dz)
		}
		#pragma omp parallel for
		FOR(i,1,NX-1) {
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(i,0,0),      CDX,CDY0_P,CDZ0_P,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(i,0,NZ-1),   CDX,CDY0_P,CDZN_P,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(i,NY-1,0),   CDX,CDYN_P,CDZ0_P,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(i,NY-1,NZ-1),CDX,CDYN_P,CDZN_P,dx,dy,dz)
		}
		#pragma omp parallel for
		FOR(j,1,NY-1) {
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(0,j,0),      CDX0_P,CDY,CDZ0_P,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(0,j,NZ-1),   CDX0_P,CDY,CDZN_P,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(NX-1,j,0),   CDXN_P,CDY,CDZ0_P,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(NY-1,j,NZ-1),CDXN_P,CDY,CDZN_P,dx,dy,dz)
		}
		#pragma omp parallel for
		FOR(k,1,NZ-1) {
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(0,0,k),      CDX0_P,CDY0_P,CDZ,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(0,NY-1,k),   CDX0_P,CDYN_P,CDZ,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(NX-1,0,k),   CDXN_P,CDY0_P,CDZ,dx,dy,dz)
			MAXWELL_MOD_RHS(u,u_rhs,INDX3D(NX-1,NY-1,k),CDXN_P,CDYN_P,CDZ,dx,dy,dz)
		}

		MAXWELL_MOD_RHS(u,u_rhs,INDX3D(0,0,0),         CDX0_P,CDY0_P,CDZ0_P,dx,dy,dz)
		MAXWELL_MOD_RHS(u,u_rhs,INDX3D(0,0,NZ-1),      CDX0_P,CDY0_P,CDZN_P,dx,dy,dz)
		MAXWELL_MOD_RHS(u,u_rhs,INDX3D(0,NY-1,0),      CDX0_P,CDYN_P,CDZ0_P,dx,dy,dz)
		MAXWELL_MOD_RHS(u,u_rhs,INDX3D(NX-1,0,0),      CDXN_P,CDY0_P,CDZ0_P,dx,dy,dz)
		MAXWELL_MOD_RHS(u,u_rhs,INDX3D(0,NY-1,NZ-1),   CDX0_P,CDYN_P,CDZN_P,dx,dy,dz)
		MAXWELL_MOD_RHS(u,u_rhs,INDX3D(NX-1,0,NZ-1),   CDXN_P,CDY0_P,CDZN_P,dx,dy,dz)
		MAXWELL_MOD_RHS(u,u_rhs,INDX3D(NX-1,NY-1,0),   CDXN_P,CDYN_P,CDZ0_P,dx,dy,dz)
		MAXWELL_MOD_RHS(u,u_rhs,INDX3D(NZ-1,NY-1,NZ-1),CDXN_P,CDYN_P,CDZN_P,dx,dy,dz)

	}

	void Apply_BCs(Field* u, Field* u_rhs) {
		Apply_Periodic_BCs(u,u_rhs);
	}

	void Rhs_Single(Field* u, Field* u_rhs, int I) {
		// u_rhs[0+I] = IM*CDY(u,2+I,dy) - IM*CDZ(u,1+I,dz);
		// u_rhs[1+I] = IM*CDZ(u,0+I,dz) - IM*CDX(u,2+I,dx);
		// u_rhs[2+I] = IM*CDX(u,1+I,dx) - IM*CDY(u,0+I,dy);
		MAXWELL_MOD_RHS(u,u_rhs,I,CDX,CDY,CDZ,dx,dy,dz)
	}

	void Constraint_Single(Field* u, int I, Field* G_ret) {
		G_ret[0] = CDX(u,0+I,dx) + CDY(u,1+I,dy) + CDZ(u,2+I,dz);
	}

	void Init_Gaussian(Field*& u) {
		Name += "G-N-" + std::to_string(NNN) + "-B-" + std::to_string(BETA);
		Base A = Base(100);
		Base B = Base(1);
		FOR(i,0,NX) FOR(j,0,NY) FOR(k,0,NZ) {
			Base x = X[i];
			Base y = Y[j];
			Base z = Z[k];
			FOR(n,0,NF) {u[n+INDX3D(i,j,k)] = 0;}


			u[0 + INDX3D(i,j,k)] = -IM*A*y*std::exp(-B*(x*x+y*y));
			u[1 + INDX3D(i,j,k)] = IM*A*y*std::exp(-B*(x*x+y*y));
		}
	}

	void Init_Periodic_At_T(Field*& u, Base t) {
		Name += "P-N-" + std::to_string(NNN) + "-B-" + std::to_string(BETA);
		Base s3 = std::sqrt(3);

		FOR(i,0,NX) FOR(j,0,NY) FOR(k,0,NZ) {
			Base x = X[i];
			Base y = Y[j];
			Base z = Z[k];
			FOR(n,0,NF) {u[n+INDX3D(i,j,k)] = 0;}
			u[0 + INDX3D(i,j,k)] = std::exp(IM*(x+y+z)/s3-IM*t);
			u[1 + INDX3D(i,j,k)] = -(Base(1.0)+IM*s3)*Base(0.5)*std::exp(IM*(x+y+z)/s3-IM*t);
			u[2 + INDX3D(i,j,k)] = (-Base(1.0)+IM*s3)*Base(0.5)*std::exp(IM*(x+y+z)/s3-IM*t);
			u[3 + INDX3D(i,j,k)] = 0;
		}
	}

	void Init_Radial(Field*& u) {
		Name += "R-N-" + std::to_string(NNN) + "-B-" + std::to_string(BETA);

		FOR(i,0,NX) FOR(j,0,NY) FOR(k,0,NZ) {
			Base x = X[i];
			Base y = Y[j];
			Base z = Z[k];
			Base r = std::sqrt(x*x+y*y+z*z);
			Base r3 = r*r*r;
			FOR(n,0,NF) {u[n+INDX3D(i,j,k)] = 0;}
			u[0 + INDX3D(i,j,k)] = IM*x/r3;
			u[1 + INDX3D(i,j,k)] = IM*y/r3;
			u[2 + INDX3D(i,j,k)] = IM*z/r3;
		}		
	}
};
