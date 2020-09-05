//*************************************************************************************************
//          name:       Belen.c		
//
//   description:       Este programa implementa el "algoritmo de Belen", un nuevo metodo de Monte Carlo para
//			sistemas de partículas con dos niveles de energía que simula una dinámica con "tiempo físico"
//
//            by:       Rogelio Díaz-Méndez  
//
// starting date:       February 9, 2014
//
//         place:       Belén (Habana Vieja), Cuba
//
//*************************************************************************************************

//*************************************************************************************************
// Hace la dinámica del modelo de Ising ferromagnético
//*************************************************************************************************

//*************************************************************************************************
//
//   bifurcation: 	From now on this is just the function of the Belén step in its generalized
//			formalism of minimal changes
//
//			Strasbourg, June 28, 2014
//
//*************************************************************************************************

//*************************************************************************************************
// 			Strasbourg, February 27, 2015
//
// Includes a more general function relaxing the "single spin flip" approach, so considering the 
// flip of a fraction of the spins, in average. This fraction is given as input.
//
//*************************************************************************************************


// includes
//*************************************
#include <stdlib.h>
#include <math.h>


/**********************************************************************/
/***** THIS IS COMPLETELY DEVOTED TO THE RANDOM NUMEBR GENERATOR ******/
/**********************************************************************/
#define FNORM   (2.3283064365e-10)
#define RANDOM  ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM)
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)         

unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;

unsigned rand4init(void)
{
  unsigned long long y;

  y = (myrand*16807LL);
  myrand = (y&0x7fffffff) + (y>>31);
  if (myrand&0x80000000)
    myrand = (myrand&0x7fffffff) + 1;
  return myrand;
}

void Init_Random(void)
{
  int i;

  ip = 128;
  ip1 = ip - 24;
  ip2 = ip - 55;
  ip3 = ip - 61;

  for (i=ip3; i<ip; i++)
    ira[i] = rand4init();
}

/**********************************************************************/
/*********************************** END ******************************/
/**********************************************************************/


// Functions:
//*******************
int BL_Step(double, int, double *, int *, double *);
int BL_Step_P(double, int, double *, int *, double *, int);



// paso del algoritmo de Belen
// 
// inputs:  
//		BL_beta: inverse temperature
//		BL_nmc: number of posible changes 
//		*BL_mc: array of minimal changes, in terms of energy, i. e. array with the energy-differences associated to the BL_nmc minimal changes
//		*BL_prop: number of the change proposed (spin that will flip), memory addres of a single integer
//		BL_time: time step
//
// return:	number of accepted minimal changes, before selecting only one	 	 
//
//*********************************************************
int BL_step(double BL_beta, int BL_nmc, double *BL_mc, int *BL_prop, double *BL_time){

	int i, flips;
	double sump2e;
	double dice, lastDice; 
	

	// calcula y guarda las probabilidades de equilibrio de todos los cambios mínimos, y también la suma
	sump2e=0;
	for (i=0;i<BL_nmc;i++){
		BL_mc[i]=1/(1+exp(-BL_beta*(-BL_mc[i])));
		sump2e+=BL_mc[i];
	}

	if (sump2e<=1.0) sump2e=1/(1-exp(-1)); // ajuste en caso de que las probabilidades sean muy pequeñas 

	flips=0;	// pone a cero los flips
	lastDice=0;	// pone a cero el dado	
	for (i=0;i<BL_nmc;i++) { 	// para cada cambio minimo
		if (FRANDOM < (BL_mc[i]/sump2e)) {    // acepta con cierta probabilidad 
			dice=FRANDOM;		// tira el dado 
			if (dice>lastDice){	// si ganó
				*BL_prop = i;	// actualiza
				lastDice = dice; // y recuerda el resultado del dado
			} 
			flips++;	// cuenta todos los aceptados
		}
	}

	*BL_time = -log(1-(1/sump2e));	// halla el tiempo fisico

	if (flips==0) *BL_prop=-1;	// si nadie se aceptó devuelve el mc -1

	return (flips);
}






// paso del algoritmo de Belen soft (general)
// 
// inputs:  
//		BL_beta: inverse temperature
//		BL_nmc: number of posible changes 
//		*BL_mc: array of minimal changes, in terms of energy, i. e. array with the energy-differences associated to the BL_nmc minimal changes
//		*BL_sp: array of numbers of the changes proposed  (spins that will flip), array of size BL_nmc
//		*BL_time: time step to be updated
//		BL_nn: normalization condition (average number of proposed)
//
// return:	number of accepted minimal changes, all of them are allowed 	 
//
//*********************************************************
int BL_step_P(double BL_beta, int BL_nmc, double *BL_mc, int *BL_prop, double *BL_time, int BL_nn){

	int i, flips;
	double sump2e;
	//double dice, lastDice; 
	

	// calcula y guarda las probabilidades de equilibrio de todos los cambios mínimos, y también la suma
	sump2e=0;
	for (i=0;i<BL_nmc;i++){
		BL_mc[i]=1/(1+exp(-BL_beta*(-BL_mc[i])));
		sump2e+=BL_mc[i];
	}

	if (sump2e<=BL_nn) sump2e=1/(1-exp(-1)); // ajuste en caso de que las probabilidades sean muy pequeñas (BL_nn adjust)

	
	for (i=0;i<BL_nmc;i++) BL_prop[i]=-1;  // erase propositions


	flips=0;	// pone a cero los flips
	//lastDice=0;	// pone a cero el dado	
	for (i=0;i<BL_nmc;i++) { 	// para cada cambio minimo
		if (FRANDOM < (BL_nn*BL_mc[i]/sump2e)) {    // acepta con cierta probabilidad  (BL_nn adjust)
			//dice=FRANDOM;		// tira el dado 
			//if (dice>lastDice){	// si ganó
			//	*BL_prop = i;	// actualiza
			//	lastDice = dice; // y recuerda el resultado del dado
			//}
			BL_prop[flips]=i;	// adiciona el cambio i a los propuestos 
			flips++;	// cuenta todos los aceptados
		}
	}

	*BL_time = -log(1-(1/sump2e));	// halla el tiempo fisico

	//if (flips==0) *BL_prop=-1;	// si nadie se aceptó devuelve el mc -1

	return (flips);
}





