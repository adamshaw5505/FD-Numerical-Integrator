#pragma once
#define NNN 1000
#define NX NNN
#define NY 1
#define NZ 1
#define NF 10
#define NG 1
#define KOEPS Base(0.000)

#include "../Headers/Definitions.h"
#include "../Headers/Evolution.h"

#define CDWX CDX
#define CDWY DZERO
#define CDWZ DZERO

#define h   u[0+I]
#define Nr  u[1+I]
#define g   u[2+I]
#define p   u[3+I]
#define chi u[4+I]
#define f   u[5+I]
#define B   u[6+I]
#define C   u[7+I]
#define D   u[8+I]
#define E   u[9+I]

#define DXh   CDWX(u,0+I,dx)
#define DXNr  CDWX(u,1+I,dx)
#define DXg   CDWX(u,2+I,dx)
#define DXp   CDWX(u,3+I,dx)
#define DXchi CDWX(u,4+I,dx)
#define DXf   CDWX(u,5+I,dx)
#define DXB   CDWX(u,6+I,dx)
#define DXC   CDWX(u,7+I,dx)
#define DXD   CDWX(u,8+I,dx)
#define DXE   CDWX(u,9+I,dx)

#define ASH_RHS(u,u_rhs,I,DX,DY,DZ) u_rhs[0+I] = -h*DXNr + Base(2.0)*h*(D*g-E*p); \
                                    u_rhs[1+I] = -h*DXh + Base(2.0)*h*(D*p-E*g); \
                                    u_rhs[2+I] = -h*DXp +D*(g*g+p*p)-Base(2.0)*E*p*g+B*h*g-p*(B*Nr-C); \
                                    u_rhs[3+I] = -h*DXg -E*(g*g+p*p)+Base(2.0)*D*p*g+B*h*p-g*(B*Nr-C); \
                                    u_rhs[4+I] = h*DXf - DXh + Base(2.0)*E*g-Base(2.0)*D*p; \
                                    u_rhs[5+I] = h*DXchi - DXNr + Base(2.0)*chi*(E*g-D*p); \
                                    u_rhs[6+I] = -Nr*DXB + DXC + (g*g-p*p)*(D*D-E*E+Base(1.0))/h; \
                                    u_rhs[7+I] = (h*h-Nr*Nr)*DXB + Nr*DXC + Nr*(g*g-p*p)*(D*D-E*E+Base(1.0))/h; \
                                    u_rhs[8+I] = -h*DXE + B*D*h + B*E*Nr - E*C; \
                                    u_rhs[9+I] = -h*DXD + B*E*h + B*D*Nr - D*C;


#define ASH_CONSTRAINT(G,DX,DY,DZ) G[0] = DXh - h*DXf + Base(2.0)*(D*p-E*g);

class Ash_Spherical_System : public SysParams {
    public:

    Field* u_0 = new Field[NF*NX*NY*NZ];
    Base M = Base(1.0);
    Ash_Spherical_System(Base LX0, Base LX1, Base M) : SysParams(LX0, LX1, 0.0, 0.0, 0.0, 0.0), M(M) {
        Name = "ASH";
    }

    void Rhs_Single(Field* u, Field* u_rhs, int I) {
		ASH_RHS(u,u_rhs,I,CDWX,CDWY,CDWZ);
	}

    void Constraint_Single(Field* u, int I, Field* G_ret) {
        ASH_CONSTRAINT(G_ret,CDWX,CDWY,CDWZ);
    }

    void Init_Isotropic_Minkowski(Field*& u) {
        Name += "IS-N-" + std::to_string(NNN);

        std::fstream log_file;
		log_file.open("./Logs/"+Name+".log", std::fstream::in | std::fstream::out | std::fstream::trunc);
		log_file.close();
		log_file.open(("Data/"+Name+"-Solution"+".dat").c_str(), std::fstream::trunc);
		log_file.close();

        FOR(i,0,NX) FOR(j,0,NY) FOR(k,0,NZ) {
            Base r     = X[i];
            Base theta = Y[j];
            Base phi   = Z[k];

            FOR(n,0,NF) {u[n+INDX3D(i,j,k)] = Base(0.0);}
            int I = INDX3D(i,j,k);

            h = Base(1.0);
            g = Base(1.0)/r;
            f = -Base(2.0)*std::log(r);
            E = Base(1.0);

        }
    }

    Base Tau_Function(Field* u) {
        return dt;
    }

	void Log_Solution(Field* u, Base Tau) {
		std::fstream log_file;
		log_file.open("./Logs/"+Name+".log", std::fstream::in | std::fstream::out | std::fstream::app);

		if (!log_file) {
			log_file.open("./Logs/"+Name+".log", std::fstream::in | std::fstream::out | std::fstream::trunc);
		}
		log_file << "Tau:" << Tau << " ";
        int I = 0;
		FOR(i,0,NX) FOR(j,0,NY) FOR(k,0,NZ) {
            I = INDX3D(i, j, k);
            log_file << chi << ", ";
		}
        log_file << "\n";
		log_file.close();
	}
};