//*************************************************************************************************
//          name:       bising.c		
//
//   description:       Este programa implementa el "algoritmo de Belen", un novedoso metodo de Monte Carlo para
//			sistemas de partículas con dos niveles de energía que simula una dinámica con "tiempo físico"
//
//        author:       Dr. Rogelio Díaz-Méndez  
//
// starting date:       February 9, 2014
//
//         place:       Belén (Habana Vieja), Cuba
//
//*************************************************************************************************

//*************************************************************************************************
// Hace la dinámica del modelo de Ising ferromagnético
//*************************************************************************************************



//*******************************************************************
//  bifurcation:        bperc.c
//
//  description:        this program evolves the Belen algorithm on the 2D Ising model up to
//                      the formation of the first percolation cluster
//  
//           by:        Rogelio Diaz-Mendez
//
//         date:        January 26, 2015
//
//        place:        Strasbourg, France
//
//*******************************************************************
                
//*******************************************************************
//	Boundary conditions should be free
//*******************************************************************

//*******************************************************************
//  bifurcation:        bpising2d.c
//
//  description:        now it continues until the final arrest, tracking the final configurations, to be 
//			to be passed to a CNN 
//  
//           by:        Rogelio Diaz-Mendez
//
//         date:        August 3, 2017 
//
//        place:        Stockholm, Sweden
//
//*******************************************************************




// inclusions
//*************************************

#include <stdio.h>
#include <math.h>
#include "belen.c"



// definiciones
//*************************************

#define	pi  acos(-1)
#define	N   l*l
#define Z   4
#define delta 1


// functions pre-declaration
//*************************************

void allocate_memory();
void free_memory();
void set_nearest_neighbours();
void init_S(int);
double set_exchange_energy();
double set_delta_energy(int); 

//int testperc();

int find_percolation(int*, int);

double update_nx();
double update_ny();
double update_m();



// variables globales
//*************************************

int l;			// tamanho lineal del sistema
int *S, *neighbours;	// espines y arreglo de vecinos
double energy, *E2;	// energia y arreglo de energias de los segundos niveles
double T, beta;		// temperatura e inverso
unsigned myrand_back;	// copia de la semilla
double	energy, m;	// energia y magnetizacion
//int flips;		// contador de flips en el paso de Belen




// main program	 
//*************************************

int main(int argc, char *argv[]) {


	//time_t tbeg, tend;
	int i;
	int conf;//, nsmax;	// configuracion inicial, numero del paso y numero total de pasos
	double t, dt; 	// tiempo total y paso de tiempo (variable)
	int sp, acc;

	double tave=0;

	//double temporal;
	//char filename[30];

	int nt, samples, ss;
	//double tstep;

	int perc;

	//FILE *ale;



	FILE *opf;

	char *dir,ord[300];

	double nx,ny;
	double sfrac=0;

	int cp;
	double mp,mf;





 	// recording starting time
        //tbeg = time(NULL);

        if(argc!=7) {
                printf("usage: %s <l> <conf> <T> <samples> <dir> <seed>\n", argv[0]);
                exit(1);
        }

        l=(int)atoi(argv[1]);		// tamanho lineal del sistema 
        conf=(int)atoi(argv[2]);   	// configuracion inicial
	T=(double)atof(argv[3]);	// Temperatura	
	samples=(int)atoi(argv[4]);	//  
	dir=argv[5];
	myrand=(unsigned)atoi(argv[6]);	// seed of the random numbers generator

	myrand_back=myrand;	// copiando la semilla original
	Init_Random();		// inicializando el generador

	beta = ((T == 0.0) ? 99999999999.0 : 1.0/T);	// correccion para temperatura cero

	allocate_memory();	// da memoria a los arreglos

	set_nearest_neighbours();	// llena el arreglo de vecinos




	sprintf(ord, "touch %s; rm -r %s; mkdir %s",dir,dir,dir);
	system(ord);


	for (ss=0;ss<samples;ss++){

		init_S(conf); 	// conf=-9 aleatorio

		energy = set_exchange_energy(); // actualizando la energia

		// add the snapshot to the file "firsts.dat"
		sprintf(ord, "%s/firsts.dat",dir); 
		opf = fopen (ord, "a"); 
		for (i=0;i<N;i++){
			if (i%l==0) fprintf(opf,"\n");
			fprintf(opf,"%i\t", S[i]);
		}
		fprintf(opf,"\n\n"); 
		fclose(opf);


		t=0;
		nt=1;
		//for (ns=1; ns<=nsmax; ns++) {	// ciclo por el total de pasos
		perc=0;
		//while (nt<nsmax){
		while (perc==0){

			for (i=0;i<N;i++) E2[i]=set_delta_energy(i); // actualizando las diferencias de energía de los cambios mínimos
			acc=BL_step(beta,N,E2,&sp,&dt);  // hace el paso de Belén

			if (acc!=0) { 		// si alguno se aceptó
				energy+=set_delta_energy(sp);
				S[sp]=-S[sp];
			}

			t+=dt; // siempre incrementar el tiempo

			perc=(find_percolation(S,1)||find_percolation(S,-1));			

			nt++;
		}

		//printf("%f\t%i\n",t,nt-1);

		tave+=t;


		// add the snapshot to the file
		sprintf(ord, "%s/snapshots.dat",dir); 
		opf = fopen (ord, "a"); 
		for (i=0;i<N;i++){
			if (i%l==0) fprintf(opf,"\n");
			fprintf(opf,"%i\t", S[i]);
		}
		fprintf(opf,"\n\n"); 
		fclose(opf);

                // saving percolation measures 
                cp=(find_percolation(S,1)?1:-1);
                mp=update_m();


		// loop until one of the two final configurations is reached
		nx=update_nx();
		ny=update_ny();
		while (nx!=0&&ny!=0){

			for (i=0;i<N;i++) E2[i]=set_delta_energy(i); // actualizando las diferencias de energía de los cambios mínimos
			acc=BL_step(beta,N,E2,&sp,&dt);  // hace el paso de Belén

			if (acc!=0) { 		// si alguno se aceptó
				energy+=set_delta_energy(sp);
				S[sp]=-S[sp];
			}

			t+=dt; // siempre incrementar el tiempo
			nt++;

			nx=update_nx();
			ny=update_ny();
		}		

		

                // add the label of the configuration to the file
                sprintf(ord, "%s/targets.dat",dir);
                opf = fopen (ord, "a");
                if (nx==ny) // and equals to zero, of course
                        if (S[0]==1) fprintf(opf,"1\t0\t0\n");
                        else fprintf(opf,"0\t1\t0\n");

                else
                        fprintf(opf,"0\t0\t1\n");
                fclose(opf);



		

		
		// add the snapshot to the file "lasts.dat"
		sprintf(ord, "%s/lasts.dat",dir); 
		opf = fopen (ord, "a"); 
		for (i=0;i<N;i++){
			if (i%l==0) fprintf(opf,"\n");
			fprintf(opf,"%i\t", S[i]);
		}
		fprintf(opf,"\n\n"); 
		fclose(opf);

               // add magnetization measures to the file
                mf=update_m();
                sprintf(ord, "%s/magnetizations.dat",dir);
                opf = fopen (ord, "a");
                fprintf(opf,"%i\t%f\t%f\n",cp,mp,mf);
                fclose(opf);


		if (nx!=ny) sfrac++;


	}	


	sfrac=sfrac/samples;

	// write the stripes fraction
	sprintf(ord, "%s/readme.dat",dir); 
	opf = fopen (ord, "a");
	fprintf(opf,"stripes fraction = %f\n",sfrac); 
	fclose(opf);




	
	//printf("%f\n",tave/samples);

	free_memory();	// libera la memoria

	return 1;
}



//*******************
double update_nx(){
	int i;
	double nx_=0;

	for (i=0;i<N;i++) nx_+=abs(S[i]*S[neighbours[Z*i]]-1)/2;
	
	return (nx_/N);
}

//*******************
double update_ny(){
	int i;
	double ny_=0;

	for (i=0;i<N;i++) ny_+=abs(S[i]*S[neighbours[(Z*i)+3]]-1)/2;
	
	return (ny_/N);
}


//*******************
double update_m(){
        int i;
        double m_=0;

        for (i=0;i<N;i++) m_+=(double)S[i];

        return (m_/((double)N));
}










//****************************************************************
// Devuelve 1 si encuentra la percolación de sitios con target FP_target
// en la matriz FP_mat
//
int find_percolation(int *FP_mat, int FP_target){

        int i, j;
        int pepe = 0;

        int label[N];

        // horizontal percolation
        for (i=0;i<N;i++) label[i]=0; // clean
        for (i=0;i<l;i++) if (FP_mat[i*l]==FP_target) label[i*l]=1; // updating first column
        for (j=1;j<l;j++){ // resto de las columnas
                // updating desde la izquierda
                for (i=0;i<l;i++)
                        if ((label[(i*l)+(j-1)]==1) && (FP_mat[i*l+j]==FP_target)) {
                                label[i*l+j]=1;
                                if (j==l-1) pepe=1;
                        }
                // updating de arriba hacia abajo
                for (i=1;i<l;i++)
                        if ((label[((i-1)*l)+(j)]==1) && (FP_mat[i*l+j]==FP_target)) {
                                label[i*l+j]=1;
                                //if (j==l-1) pepe=1;
                        }
                // updating de abajo hacia arriba
                for (i=l-2;i>=0;i--)
                        if ((label[((i+1)*l)+(j)]==1) && (FP_mat[i*l+j]==FP_target)) {
                                label[i*l+j]=1;
                                //if (j==l-1) pepe=1;
                        }


        }

        // imprimiendo la matriz label      
        //for(i=0;i<N;i++){
        //        printf("%i\t",label[i]);
        //        if (div(i+1,l).rem==0) printf("\n");
        //}
        //printf("\n");



        // vertical percolation
        for (i=0;i<N;i++) label[i]=0; // clean
        for (i=0;i<l;i++) if (FP_mat[i]==FP_target) label[i]=1; // updating first line
        for (j=1;j<l;j++){ // resto de las lineas
                // updating desde arriba
                for (i=0;i<l;i++)
                        if ((label[((j-1)*l)+i]==1) && (FP_mat[j*l+i]==FP_target)) {
                                label[j*l+i]=1;
                                if (j==l-1) pepe=1;
                        }
                // updating de izquierda a derecha
                for (i=1;i<l;i++)
                        if ((label[(j*l)+(i-1)]==1) && (FP_mat[j*l+i]==FP_target)) {
                                label[j*l+i]=1;
                                //if (j==l-1) pepe=1;
                        }
                // updating de abajo hacia arriba
                for (i=l-2;i>=0;i--)
                        if ((label[(j*l)+(i+1)]==1) && (FP_mat[j*l+i]==FP_target)) {
                                label[j*l+i]=1;
                                //if (j==l-1) pepe=1;
                        }


        }

        // imprimiendo la matriz label      
        //for(i=0;i<N;i++){
        //        printf("%i\t",label[i]);
        //        if (div(i+1,l).rem==0) printf("\n");
        //}
        //printf("\n");

        return (pepe);

}






/*
// comprueba si hay percolación en una salida de hk
// *********************************************************
int testperc(){

	int i, j, jj;
	int up=0,down=0,left=0,right=0;
	int percolation=0;

	for (i=1;i<N;i++)//clusters[0];i++)//if (clusters[3+i]>l)  // VERY SLOW!!!
		{
			
			up=0,down=0,left=0,right=0;
			for (j=0;j<N;j++){
				jj=3+N+j;
				if (clusters[jj]==i){
					if (j<l) up=1;
					if (j>=N-l) down=1;
					if (j%l==0) left=1;
					if ((j+1)%l==0) right=1;
				}
			}
			if ((up && down) || (left && right)) percolation = 1;
			//exit?
		}

	return (percolation);
}
*/




// diferencia de energia al variar el espin n en un dS 
//*********************************************************
double set_delta_energy(int n) {


  	double ene=0;
  	double sum;

	

    	//exchange   OBC
    	sum=( (neighbours[Z*n]<0?0:S[neighbours[Z*n]]) + (neighbours[(Z*n)+1]<0?0:S[neighbours[(Z*n)+1]]) + (neighbours[(Z*n)+2]<0?0:S[neighbours[(Z*n)+2]]) + (neighbours[(Z*n)+3]<0?0:S[neighbours[(Z*n)+3]]) );

	ene=2*(delta)*sum*S[n];

  	return (ene);
}






// devuelve la energia de intercambio ie. corto alcance.
//*********************************************************
double set_exchange_energy() {

    int i, j, kdown=0, krigth=0;
    double energy_=0;

    for (j=0; j<l; j++)
    {
                for (i=0; i<l; i++)
                {
                        // The interactions at this poits are betwen:
                        // kdown+i             and         neighbours[(Z*(kdown+i))+3]
                        // krigth+(i*l)        and         neighbours[Z*(krigth+(i*l))]

                        // Calculation of the hamiltonian:  OBC
                        energy_+=-delta*S[kdown+i]* (neighbours[(Z*(kdown+i))+3]<0?0:S[neighbours[(Z*(kdown+i))+3]]);
                        energy_+=-delta*S[krigth+(i*l)]* (neighbours[Z*(krigth+(i*l))]<0?0:S[neighbours[Z*(krigth+(i*l))]]);
                }
                kdown=neighbours[(Z*kdown)+3];
                krigth=neighbours[Z*krigth];
    }
    return(energy_);
}




// inicializa el valor de los espines segun el parámetro hh
//*********************************************************
void init_S(int hh) {
    
	int i, j;

        // ferro
    	if (hh==0) for (i=0;i<l*l;i++) S[i]=1;

    	// stripes
    	if (hh>0) {
		for (i=0;i<l*l;i++) S[i]=1;
              	for (i=0;i<l;i+=2*hh)
                	for (j=0;(j<hh*l)&&((i*l+j)<(l*l));j++) S[i*l+j]=-1;
    	}

	// antiferro AF
    	if (hh==-1) {
        	for (i=0;i<l*l;i++) S[i]=1;
                for (i=0;i<l;i++)
                	for(j=((div(i,2).rem==0) ? 0 : 1);j<l;j+=2) S[(i*l)+j]=-1;
   	}

    	// aleatorio
    	if (hh==-9) for (i=0;i<l*l;i++) S[i]=pm1;

}





// Genera el arreglo de los vecinos cercanos
//*********************************************************
void set_nearest_neighbours() {
   
	int i;

	int tipo=2;

	if (tipo==1){ // PBC
		for (i = 0; i < N; i++) // I know this is not the best method but I don't care
   		{                       // becouse 1-This do not slow the program and 2-This 
                	          	// is mine.
        		// rigth               
        		if (div(i+1, l).rem == 0) neighbours[(Z * i)] = i - (l -1);
	        	else neighbours[(Z * i)] = i + 1;
        		// left
	        	if (div(i, l).rem == 0) neighbours[(Z * i) + 1] = i + (l -1);
        		else neighbours[(Z * i) + 1] = i -1;
        		// up
	        	if (i < l)  neighbours[(Z * i) + 2] = i + ((l -1) * l);
        		else neighbours[(Z * i) + 2] = i - l;
        		// down
	        	if ((i + l) >= N) neighbours[(Z * i) + 3] = i -((l - 1) * l);
        		else neighbours[(Z * i) + 3] = i + l;
    		}
	}
	else 
	if (tipo==2){ // OBC
                for (i = 0; i < N; i++) 
                {                       
                        // rigth               
                        if (div(i+1, l).rem == 0) neighbours[(Z * i)] = -1;
                        else neighbours[(Z * i)] = i + 1;
                        // left
                        if (div(i, l).rem == 0) neighbours[(Z * i) + 1] = -1;
                        else neighbours[(Z * i) + 1] = i -1;
                        // up
                        if (i < l)  neighbours[(Z * i) + 2] = -1;
                        else neighbours[(Z * i) + 2] = i - l;
                        // down
                        if ((i + l) >= N) neighbours[(Z * i) + 3] = -1;
                        else neighbours[(Z * i) + 3] = i + l;
                }
        }

}




// function that allocates the memory for the arrays
//*********************************************************
void allocate_memory(){

  	S = (int *) malloc(N*sizeof(int));
 	neighbours=(int *) malloc(N*Z*sizeof(int));
  	E2 = (double *) malloc(N*sizeof(double));	
}



// function that free the memory of the arrays
//*********************************************************
void free_memory(){

	free(S);
	free(neighbours);
	free(E2);
}












