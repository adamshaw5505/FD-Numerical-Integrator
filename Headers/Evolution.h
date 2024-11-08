#pragma once

#ifndef NG
#define NG 0
#endif

#ifndef NX
#define NX 0
#endif

#ifndef NY
#define NY 0
#endif

#ifndef NZ
#define NZ 0
#endif


#ifndef Base
#define Base float
#endif

#ifndef Field
#define Field float
#endif

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif

#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif

void Calculate_Constraint(Field* u, SysParams* sys) {
	int I = 0;
	std::vector<Base> G_ret(NG);
	Field G_temp[NG] = {Field(0)};
	Field G_value[NG] = {Field(0)};

	//#pragma omp parallel for
	#if NX != 1
	FOR(i,1,NX-1) {
	#else
	int i = 0;
	#endif
	
	#if NY != 1
		FOR(j,1,NY-1) {
	#else
		int j = 0;
	#endif

	#if NZ != 1
		FOR(k,1,NZ-1) {
	#else
		int k = 0;
	#endif

		sys->Constraint_Single(u,INDX3D(i,j,k),G_temp);
		FOR(n,0,NG) {
			G_value[n] += std::pow(std::abs(G_temp[n]),2)*sys->dx*sys->dy*sys->dz;
		}

	#if NZ != 1
			}
	#endif

	#if NY != 1
		}
	#endif

	#if NX != 1
	}
	#endif
	Base NN = Base(NX*NY*NZ);
	FOR(n,0,NG) {
			G_ret[n]  = std::real(std::sqrt(G_value[n]/NN));
	}

	sys->G.push_back(G_ret);
}

void Rhs(Field* u, Field* u_rhs, SysParams* sys) {

	#pragma omp parallel for
	#if NX != 1
	FOR(i,1,NX-1) {
	#else
	int i = 0;
	#endif
	
	#if NY != 1
		FOR(j,1,NY-1) {
	#else
		int j = 0;
	#endif

	#if NZ != 1
		FOR(k,1,NZ-1) {
	#else
		int k = 0;
	#endif

		sys->Rhs_Single(u,u_rhs,INDX3D(i,j,k));

	#if NZ != 1
			}
	#endif

	#if NY != 1
		}
	#endif

	#if NX != 1
	}
	#endif
}

void ICN_Step(Field* u,Field* u_new, SysParams* sys) {
	Field* u_half = new Field[NF*NX*NY*NZ];
	Field* u_rhs  = new Field[NF*NX*NY*NZ];

	Rhs(u,u_rhs,sys);

	#pragma omp parallel for
	FOR(I,0,NF*NX*NY*NZ) {
		u_new[I] = u[I] + sys->dt*u_rhs[I];
		u_half[I] = (u_new[I]+u[I])/Base(2.0);
		u_rhs[I] = Base(0.0);
	}

	Rhs(u_half,u_rhs,sys);

	#pragma omp parallel for
	FOR(I,0,NF*NX*NY*NZ) {
		u_new[I] = u[I] + sys->dt*u_rhs[I];
		u_half[I] = (u_new[I]+u[I])/Base(2.0);
		u_rhs[I] = Base(0.0);
	}

	Rhs(u_half,u_rhs,sys);

	#pragma omp parallel for
	FOR(I,0,NF*NX*NY*NZ) {
		u_new[I] = u[I] + sys->dt*u_rhs[I];
	}

	delete[] u_half;
	delete[] u_rhs;
	
}

void RK4_Step(Field* u, Field* u_new, SysParams* sys) {
	Field* k1 = new Field[NF*NX*NY*NZ];
	Field* k2 = new Field[NF*NX*NY*NZ];
	Field* k3 = new Field[NF*NX*NY*NZ];
	Field* k4 = new Field[NF*NX*NY*NZ];
	Field* ko = new Field[NF*NX*NY*NZ];

	Rhs(u,k1,sys);

	#pragma omp parallel for
	FOR(I,0,NF*NX*NY*NZ) {
		u_new[I] = u[I] + sys->dt*k1[I]/Base(2.0);
	}
	
	Rhs(u_new,k2,sys);

	
	#pragma omp parallel for
	FOR(I,0,NF*NX*NY*NZ) {
		u_new[I] = u[I] + sys->dt*k2[I]/Base(2.0);
	}

	Rhs(u_new,k3,sys);

	#pragma omp parallel for
	FOR(I,0,NF*NX*NY*NZ) {
		u_new[I] = u[I] + sys->dt*k3[I];
	}

	Rhs(u_new,k4,sys);


	#ifdef KOEPS
	sys->Apply_KO(u,ko);
	#endif
	
	FOR(I,0,NF*NX*NY*NZ) {
		#ifdef KOEPS
		u_new[I] = u[I] + sys->dt/Base(6.0) * (k1[I] + Base(2)*k2[I] + Base(2)*k3[I] + k4[I] ) + ko[I];
		#else
		u_new[I] = u[I] + sys->dt/Base(6.0) * (k1[I] + Base(2)*k2[I] + Base(2)*k3[I] + k4[I] );
		#endif
		
	}	

	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
	delete[] ko;
}

void Integrate(Field* u, Field*& u_new, SysParams* sys, Base t0, Base t1) {
	printf("Integrating program: %s \n From tau: %f to tau %f\n",sys->Name.c_str(),t0,t1);
	u_new = new Field[NF*NX*NY*NZ];
	Field* u_old = new Field[NF*NX*NY*NZ];

	FOR(I,0,NF*NX*NY*NZ) {
		u_old[I] = u[I];
		u_new[I] = Field(0.0);
	}

	Base t = t0;
	sys->Tau.push_back(t0);
	Calculate_Constraint(u,sys);
	// sys->Write_X_To_Solution_Data();
	// sys->Write_Solution_Data(u);
	Base f = Base(0.0);
	Base df = Base(0.01);
	Base temp = Base(0.0);
	while(temp < Base(1.0)) {
		std::cout << "#";
		temp += df;
	}
	std::cout << std::endl;

	timestamp_t time_start = get_time();
	while(t < t1) {
		// RK4_Step(u,u_new,sys);
		ICN_Step(u,u_new,sys);
		sys->Apply_BCs(u,u_rhs);

		t += sys->Tau_Function(u_new);
		sys->Tau.push_back(t);

		Calculate_Constraint(u_new,sys);

		sys->Log_Solution(u_new,t);
		// sys->Write_Solution_Data(u_new);

		
		// RK4_Step(u_new,u,sys);
		ICN_Step(u_new,u,sys);

		t += sys->Tau_Function(u);
		sys->Tau.push_back(t);

		Calculate_Constraint(u,sys);

		// sys->Log_Solution(u,t);
		// sys->Write_Solution_Data(u);

		
		while((t-t0)/t1 > f+df && t < t1) {
			std::cout << "#" << std::flush;
			f+=df;
		}
	}
	// sys->Write_Tau_To_Solution_Data();
	std::cout << "#" << std::endl;
	timestamp_t time_end = get_time();
	std::cout << "Time Taken: " << (time_end-time_start)/1000000.0L << std::endl;

	sys->Write_Constraint_Data(t1);
}
