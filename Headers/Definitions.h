#pragma once
#include <iostream>
#include <fstream>
#include <complex>
#include <ctime>
#include <sys/time.h>
#include <cmath>
#include <vector>

typedef unsigned long long timestamp_t;

static timestamp_t get_time() {
    struct timeval now;
    gettimeofday(&now,NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

#define PI 3.1415962
#define Base double
#define Complex std::complex<Base>
#define Field Complex
#define IM Complex(0,1)

#define FOR(II,NI,NF) for(int II = NI; II < NF; II++)
#define INDX3D(i,j,k) ((i)*NF + (j)*NF*NX + (k)*NF*NX*NY)

#define DZERO(u,I,dd) Base(0.0)

#define CDX(u,I,dx) (u[I + NF] - u[I - NF])/(Base(2.0)*dx)
#define CDY(u,I,dy) (u[I + NF*NX] - u[I - NF*NX])/(Base(2.0)*dy)
#define CDZ(u,I,dz) (u[I + NF*NX*NY] - u[I - NF*NX*NY])/(Base(2.0)*dz)

#define CDX0_P(u,I,dx) (u[I+NF]-u[I+(NX-1)*NF])/(Base(2.0)*dx)
#define CDY0_P(u,I,dy) (u[I+NF*NX]-u[I+(NY-1)*NX*NF])/(Base(2.0)*dy)
#define CDZ0_P(u,I,dz) (u[I+NF*NX*NY]-u[I+(NZ-1)*NY*NX*NF])/(Base(2.0)*dz)

#define CDXN_P(u,I,dx) (u[I-(NX-1)*NF]-u[I-NF])/(Base(2.0)*dx)
#define CDYN_P(u,I,dy) (u[I-(NY-1)*NX*NF]-u[I-NF*NX])/(Base(2.0)*dy)
#define CDZN_P(u,I,dz) (u[I-(NZ-1)*NY*NX*NF]-u[I-NF*NX*NY])/(Base(2.0)*dz)


#define CDKOX(u,I,dx) (u[I+INDX3D(2,0,0)]-Base(4.0)*u[I+INDX3D(1,0,0)] + Base(6)*u[I+INDX3D(0,0,0)] + Base(4)*u[I+INDX3D(-1,0,0)] + u[I+INDX3D(-2,0,0)])
#define CDKOY(u,I,dy) (u[I+INDX3D(0,2,0)]-Base(4.0)*u[I+INDX3D(0,1,0)] + Base(6)*u[I+INDX3D(0,0,0)] + Base(4)*u[I+INDX3D(0,-1,0)] + u[I+INDX3D(0,-2,0)])
#define CDKOZ(u,I,dz) (u[I+INDX3D(0,0,2)]-Base(4.0)*u[I+INDX3D(0,0,1)] + Base(6)*u[I+INDX3D(0,0,0)] + Base(4)*u[I+INDX3D(0,0,-1)] + u[I+INDX3D(0,0,-2)])


class SysParams {
	public:

		Base dt, dx, dy, dz;
		Base* X;
		Base* Y;
		Base* Z;

		std::vector<Base> Tau;

		std::vector<std::vector<Base>> G;
		
		std::string Name = "";
		
		void Write_Constraint_Data(Base t1) {
			std::ofstream file;
			file.open(("Data/"+Name+"-Const"+".dat").c_str(), std::fstream::trunc);
			int N = Tau.size();
			if (G.size() < Tau.size()) {N = G.size();}

			FOR(i,0,N) {
				if (Tau[i] > t1 || std::isnan(Tau[i])) {continue;}
				file << Tau[i];
				FOR(n,0,NG) {
					file <<  ",  " << G[i][n];
				}
				file << "\n";
			}
			file.close();
		};

		void Write_Solution_Data(Field* u) {
			std::ofstream file;
			// if(Tau.size() == 0) {
			// 	file.open(("Data/"+Name+"-Sol-0.0"+".dat").c_str(), std::fstream::trunc);
			// } else {
			// 		file.open(("Data/"+Name+"-Sol-"+std::to_string(Tau[Tau.size()-1])+".dat").c_str(), std::fstream::trunc);
			// }

			file.open(("Data/"+Name+"-Sol-0.0"+".dat").c_str(), std::fstream::app);
			if (!file) {
				file.open(("Data/"+Name+"-Sol-0.0"+".dat").c_str(), std::fstream::trunc);
			}

			FOR(I,0,NF*NX*NY*NZ) {
				file  << std::real(u[I]) << ", ";
				file << std::imag(u[I]) << ", ";
			}
			file << "\n";
			file.close();
		}

		void Write_X_To_Solution_Data() {
			std::ofstream file;
			file.open(("Data/"+Name+"-Sol-0.0"+".dat").c_str(), std::fstream::app);
			if (!file) {
				file.open(("Data/"+Name+"-Sol-0.0"+".dat").c_str(), std::fstream::trunc);
			}

			#if NX > 1
				FOR(i,0,NX) {
					file << X[i] << ", ";
				}
				file << "\n";
			#endif

			file.close();
		}

		void Write_Tau_To_Solution_Data() {
			std::ofstream file;
			file.open(("Data/"+Name+"-Sol-0.0"+".dat").c_str(), std::fstream::app);
			if (!file) {
				file.open(("Data/"+Name+"-Sol-0.0"+".dat").c_str(), std::fstream::trunc);
			}

			FOR(i,0,Tau.size()) {
				file << Tau[i] << ", ";
			}
			file << "\n";
			file.close();
		}

		Base Convergence_Factor(Field* u, Field* u_exact) {
			Base u_rms = Base(0);

			FOR(I,0,NF*NX*NY*NZ) {
				u_rms += std::pow(std::abs(u[I]-u_exact[I]),2);
			}
			u_rms /= NX*NY*NZ*NF;

			u_rms = std::sqrt(u_rms);

			return u_rms;
		}

		virtual void Apply_BCs(Field* u, Field* u_new) {};
		virtual void Rhs_Single(Field* u, Field* u_rhs,int I) {};
		virtual void Constraint_Single(Field* u, int I, Field* G_return) {};
		virtual Base Tau_Function(Field* u) {return dt;};
		virtual void Log_Solution(Field* u, Base tau) {};

		SysParams(Base LX0, Base LX1, Base LY0, Base LY1, Base LZ0, Base LZ1) {
			dt = 0.0;

			#if NX == 1
				dx = 1.0;
			#else
				dx = (LX1-LX0)/(NX-1);
				dt += dx*dx;
			#endif
			
			#if NY == 1
				dy = 1.0;
			#else
				dy = (LY1-LY0)/(NY-1);
				dt += dy*dy;
			#endif

			#if NZ == 1
				dz = 1.0;
			#else
				dz = (LZ1-LZ0)/(NZ-1);
				dt += dz*dz;
			#endif

			dt = Base(0.1)*std::sqrt(dt);

			X = new Base[NX];
			Y = new Base[NY];
			Z = new Base[NZ];

			FOR(i,0,NX) FOR(j,0,NY) FOR(k,0,NZ) {
				X[i] = LX0 + i*dx;
				Y[j] = LY0 + j*dy;
				Z[k] = LZ0 + k*dz;
			}
		}

		void Apply_KO(Field* u, Field* u_rhs) {
			#ifdef KOEPS
			if(KOEPS == Base(0)) {return;}
			#if NX == 1
				int i = 0;
			#else
				FOR(i,2,NX-2) {
			#endif

			#if NY == 1
				int j = 0;
			#else
				FOR(j,2,NY-2) {
			#endif

			#if NZ == 1
				int k = 0;
			#else
				FOR(k,2,NZ-2) {
			#endif

				FOR(f,0,NF) {
					#if NX > 2
					u_rhs[f+INDX3D(i,j,k)]	+= KOEPS*dt/dx * CDKOX(u,f+INDX3D(i,j,k),dx);
					#endif
					#if NY > 2
					u_rhs[f+INDX3D(i,j,k)]	+= KOEPS*dt/dy * CDKOY(u,f+INDX3D(i,j,k),dy);
					#endif
					#if NZ > 2
					u_rhs[f+INDX3D(i,j,k)]	+= KOEPS*dt/dz * CDKOZ(u,f+INDX3D(i,j,k),dz);
					#endif
				}
			#if NZ > 1
				}
			#endif
				
			#if NY > 1
				}
			#endif

			#if NX > 1
				}
			#endif
			#endif
			
		}

		~SysParams() {
			delete[] X;
			delete[] Y;
			delete[] Z;
		}
};
