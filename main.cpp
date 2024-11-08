
// #include "Examples/Weyl.h"
// void Test() {
// 	Field* u = new Field[NF*NX*NY*NZ];
// 	Field* u_new;

// 	Base L = Base(5.0);
// 	Weyl_System sys(-L,L,-L,L,-L,L,Base(1.0));
// 	sys.Init_Isotropic_Euclidian(u);
// 	Integrate(u,u_new,&sys,Base(0.0),Base(10.0));

// 	delete[] u;
// }


// #include "Examples/Weyl_Mod.h"
// void Run() {
// 	Field* u = new Field[NF*NX*NY*NZ];
// 	Field* u_new;

// 	Base L = Base(5.0);
// 	Weyl_System sys(-L,L,-L,L,-L,L,Base(1.0),Base(00.0));
// 	sys.Init_Isotropic_Euclidian(u);
// 	Integrate(u,u_new,&sys,Base(0 .0),Base(10.0));
// }

// #include "Examples/Ash_Spherical_Mod.h"
// void Run() {
// 	Field* u = new Field[NF*NX*NY*NZ];
// 	Field* u_new;
// 	Base M = Base(1.0);
// 	Base R0 = Base(10.0)*M;
// 	Base R1 = R0 + Base(10.0);
// 	Ash_Spherical_System sys(R0,R1,1.0);

// 	sys.Init_Isotropic_Minkowski(u);
// 	Integrate(u,u_new,&sys,Base(0.0),Base(6.0));

// 	delete[] u;
// }

// #include "Examples/Weyl_Spherical.h"
// // #include "Examples/Weyl_Spherical_Mod.h"
// void Run() {
// 	Field* u = new Field[NF*NX*NY*NZ];
// 	Field* u_new;

// 	Base R0 = Base(00.1);
// 	Base R1 = Base(10.1);
// 	Weyl_Spherical_System sys(R0,R1,1.0);

// 	sys.Init_Isotropic_Schwarzschild(u);

// 	Integrate(u,u_new,&sys,Base(0.0),Base(100.0));

// 	delete[] u;
// }

// #include "Examples/Weyl_Simple.h"
// void Run() {
// 	Field* u = new Field[NF*NX*NY*NZ];
// 	Field* u_new;

// 	Base R0 = Base(00.1);
// 	Base R1 = Base(10.1);
// 	Weyl_Spherical_System sys(R0,R1,1.0);

// 	sys.Init_Isotropic_Schwarzschild(u);
// 	Integrate(u,u_new,&sys,Base(0.0),Base(100.0));

// 	delete[] u;
// }


#include "Examples/Maxwell_Mod.h"
void Run() {
	Field* u = new Field[NF*NX*NY*NZ];
	Field* u_new;

	Base L = Base(5.0);
	Maxwell_System sys(-L,L,-L,L,-L,L);
	// sys.Init_Periodic_At_T(u,Base(0.0));
	sys.Init_Gaussian(u);

	Integrate(u,u_new,&sys,Base(0.0),Base(100.0));
	
	delete[] u;
}

int main() {

	Run();
	
	return 0;
}