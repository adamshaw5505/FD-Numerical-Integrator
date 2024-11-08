#pragma once
#define NNN  (2<<7)
#define NX NNN
#define NY 1
#define NZ 1
#define NF 5
#define NG 2
// #define NG 5


#include "../Headers/Definitions.h"
#include "../Headers/Evolution.h"

#define CDWX CDX
#define CDWY DZERO
#define CDWZ DZERO

#define WEYL_RHS(u,u_rhs,I,DX,DY,DZ) u_rhs[0+I] = Base(3.0)*u[2+I]*(u[3+I]*u[3+I])*u[1+I];\
                                     u_rhs[1+I] = 0;\
                                     u_rhs[2+I] = Base(2.0)*(u[1+I]/u[0+I])*((u[2+I]*u[2+I])*(u[3+I]*u[3+I]) - DX(u,3+I,dx)*DX(u,3+I,dx)) - DX(u,4+I,dx);\
                                     u_rhs[3+I] = u[2+I]*u[3+I]*u[1+I]*(u[3+I]*u[3+I] - Base(1.0))/(Base(2.0)*u[0+I]);\
                                     u_rhs[4+I] = 0;

#define WEYL_CONSTRAINT(G,DX,DY,DZ) G[0] = Base(2.0)*u[4+I]*u[3+I]*u[0+I]+(u[3+I]*u[3+I]-Base(1.0))*u[1+I]*DX(u,3+I,dx);\
                                    G[1] = -Base(6.0)*u[3+I]*u[0+I]*DX(u,3+I,dx)+(u[3+I]*u[3+I]-Base(1.0))*DX(u,0+I,dx);

#define CANONICAL_CONSTRAINT(G) Field temp = 0;\
								FOR(n,0,NF) {temp += (u[n+I] - u_0[n+I])*(u[n+I] - u_0[n+I]);}\
								G[0] = std::sqrt(temp);

#define CANONICAL_CONSTRAINT_N(G) FOR(n,0,NF) {G[n] = (u[n+I] - u_0[n+I])*(u[n+I] - u_0[n+I]);}


class Weyl_Spherical_System : public SysParams {


	public:

	Field* u_0 = new Field[NF*NX*NY*NZ];
	Base M = Base(1.0);
	
	Weyl_Spherical_System(Base LX0, Base LX1, Base M) : SysParams(LX0, LX1, 0.0, 0.0, 0.0, 0.0), M(M) {
		Name = "WS";
	}

	void Rhs_Single(Field* u, Field* u_rhs, int I) {

		WEYL_RHS(u,u_rhs,I,CDWX,CDWY,CDWZ);
	}

	void Constraint_Single(Field* u, int I, Field* G_ret) {
		#if NG == 1
			CANONICAL_CONSTRAINT(G_ret);
		#elif NG == NF
			CANONICAL_CONSTRAINT_N(G_ret)
		#else
			WEYL_CONSTRAINT(G_ret,CDWX,CDWY,CDWZ);
		#endif
	}

	void Init_Isotropic_Schwarzschild(Field*& u) {
        #ifdef KOEPS
		Name += "IS-N-" + std::to_string(NNN) + "-E-" +std::to_string(std::abs(KOEPS));
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

			u[0+INDX3D(i,j,k)] = Base(2.0)*M*std::pow(r,3)/std::pow((r+M),6); // h - Checked
			u[1+INDX3D(i,j,k)] = (r-M)*std::pow(r,4)/(std::pow(r+M,7)); // dlap - Checked
            u[2+INDX3D(i,j,k)] = 0.0; // B - Checked
			u[3+INDX3D(i,j,k)] = (r-M)/(r+M); // C - Checked
			u[4+INDX3D(i,j,k)] = -Base(2.0)*M*r*r/std::pow(r+M,4); // A - Checked

			// FOR(n,0,NF) {u_0[n+INDX3D(i,j,k)] = u[n+INDX3D(i,j,k)];}
	}
	}

	Base Tau_Function(Field* u) {
		return dt;

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
			FOR(f,0,NF) {
					log_file << u[f+INDX3D(i,j,k)] << ", "; 
			}
			log_file << "\n";
		}
		log_file << "########################################################################################\n";
		log_file.close();
	}


};
