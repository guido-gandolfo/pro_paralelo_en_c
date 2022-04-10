#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

//Variables globales 


/*R = (PromB * (A*C)) + (PromA * (B*D))

/*   Función para calcular el tiempo   */
double dwalltime(){
    double sec;
    struct timeval tv;
    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}

void paraleloRoot(int, int, int);

void paralelo(int, int, int);

int main(int argc, char* argv[2]){
    
    MPI_Init(&argc, &argv);

    int numeroProcesadores,idProceso, thilos;
    MPI_Comm_size(MPI_COMM_WORLD, &numeroProcesadores);
    MPI_Comm_rank(MPI_COMM_WORLD, &idProceso);

    int M,i,j,k,promA,promB;       
  
    //Comprobación de los parametros ingresados por teclado
    if (argc < 3){
             printf("\n Falta un argumento: \n");
             printf("\n M dimensión de la matriz \n");
             printf("\n T cantidad de hilos \n");
             MPI_Abort(MPI_COMM_WORLD,1);
        } 

    M = atoi(argv[1]);
    thilos = atoi(argv[2]);

   omp_set_num_threads(thilos);

    if(idProceso == 0) {  
	//Llama a la funcion para el proceso 0
        paraleloRoot(numeroProcesadores, thilos, M);
    }
    else{
	//Llama a la funcion para los demás procesos
        paralelo(numeroProcesadores, thilos, M);
    }
 
    MPI_Finalize(); 
    
    return 0;
} 

//Funcion para el proceso 0
void paraleloRoot(int numProc, int hilos, int N){

    //Variables locales
    int i,j,k,promA,promB;
    double *A, *AA, *B, *BB, *C, *D, *AC, *BD, *Resultado, *R;  
    double tiempo, tiempoTotal=0;
    register double acumulador = 0, acumulador2=0;
    double *suma, *vectorSuma;
    double check = 1;

    //Alocación de memoria
    A=(double*)malloc(sizeof(double)*N*N);
    B=(double*)malloc(sizeof(double)*N*N); 
    C=(double*)malloc(sizeof(double)*N*N);
    D=(double*)malloc(sizeof(double)*N*N); 
    R=(double*)malloc(sizeof(double)*N*N);
    AA=(double*)malloc(sizeof(double)*N*N);			
    BB=(double*)malloc(sizeof(double)*N*N); 	
    AC=(double*)malloc(sizeof(double)*N*N);  
    BD=(double*)malloc(sizeof(double)*N*N);
    suma=(double*)malloc(sizeof(double)*2); 
    vectorSuma=(double*)malloc(sizeof(double)*2);  
    Resultado=(double*)malloc(sizeof(double)*N*N);

    //Inicialización de matrices
    for(i=0;i<N;i++){
   		for(j=0;j<N;j++){
			A[i*N+j]= 1; //Por filas
            B[i*N+j]= 1; //Por filas
            C[i+j*N]= 1; //Por columnas
            D[i+j*N]= 1; //Por columnas
        }
    }

    //Comienza a contabilizar el tiempo
    tiempo = dwalltime();
    
    MPI_Scatter(A,N*N/numProc,MPI_DOUBLE,AA,N*N/numProc,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(B,N*N/numProc,MPI_DOUBLE,BB,N*N/numProc,MPI_DOUBLE,0,MPI_COMM_WORLD);

    MPI_Bcast(C,N*N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(D,N*N,MPI_DOUBLE,0,MPI_COMM_WORLD);

    suma[0] = 0;
    suma[1] = 0;

     //Multiplicación de matrices A*C y B*D
     #pragma omp parallel 
    {
        //Multiplicacion de matriz A*C
        #pragma omp for private(i,j,k,acumulador) reduction(+:suma[0]) nowait
        for(i=0;i<N/numProc;i++){         
            for(j=0;j<N;j++){
                acumulador = 0;
                for(k=0;k<N;k++){
                    acumulador += (AA[i*N+k]*C[k+N*j]);  
                       suma[0] += AA[i*N+k]; 
                }
                AC[i*N+j] = acumulador;  
            }   
        }  

        //Multiplicacion de matriz B*D
        #pragma omp for private(i,j,k,acumulador2) reduction(+:suma[1])
        for(i=0;i<N/numProc;i++){         
            for(j=0;j<N;j++){
                acumulador2 = 0; 
               
                for(k=0;k<N;k++){
                    acumulador2 += (BB[i*N+k]*D[k+N*j]); 
                    suma[1] += BB[i*N+k];
                    
                }
                BD[i*N+j] = acumulador2;
            }  
        }    
    }   

    MPI_Allreduce(suma,vectorSuma,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  

    //Promedio de A y B
    promA = (vectorSuma[0]/(N*N));
    promB = (vectorSuma[1]/(N*N));

    //Aplica el promedio a matriz AC y la suma de la matriz resultante
    #pragma omp parallel for private(i,j)
     
        for(i=0;i<N/numProc;i++){         
            for(j=0;j<N;j++){
               
                    //Primer termino = PromB * AC
                  AC[i*N+j] = AC[i*N+j]*promB; 
                    //Segundo termino = PromA * BD
                  BD[i*N+j] = BD[i*N+j]*promA;
                 //Resultado= Primer termino + Segundo termino
                 Resultado[i*N+j] = (AC[i*N+j] + BD[i*N+j]); 
            }  
        
        } 

     MPI_Gather( Resultado,N*N/numProc,MPI_DOUBLE,R,N*N/numProc,MPI_DOUBLE,0,MPI_COMM_WORLD);    

    tiempoTotal += dwalltime()- tiempo;  
    
    //Validación de resultado obtenido
    printf("Tiempo total: %.6f \n",tiempoTotal);
        //Verifica el resultad
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
	            check=check&&(R[i*N+j]==(double)2*N*N);
            }
        }   

        if(check){
            printf("Resultado Correcto\n");
            printf("\n");
        }else{
            printf("Resultado Incorrecto\n");
            printf("\n");
        }
    
    //Liberación de memoria
    free(A); 
    free(B);
    free(C);
    free(D);
    free(AC);
    free(BD);
    free(suma); 
    free(vectorSuma);
    free(Resultado);  

}

void paralelo(int numProc, int hilos, int N){

    //Variables locales
    int i,j,k,promA,promB;
    double *A, *AA, *B, *BB, *C, *D, *AC, *BD, *Resultado, *R;  
    register double acumulador = 0, acumulador2=0;
    double *suma, *vectorSuma;

    //Alocación de memoria
    A=(double*)malloc(sizeof(double)*N*N/numProc);
    B=(double*)malloc(sizeof(double)*N*N/numProc); 
    C=(double*)malloc(sizeof(double)*N*N);
    D=(double*)malloc(sizeof(double)*N*N); 
    R=(double*)malloc(sizeof(double)*N*N);
    AA=(double*)malloc(sizeof(double)*N*N);			
    BB=(double*)malloc(sizeof(double)*N*N); 	
    AC=(double*)malloc(sizeof(double)*N*N);  
    BD=(double*)malloc(sizeof(double)*N*N);
    suma=(double*)malloc(sizeof(double)*2); 
    vectorSuma=(double*)malloc(sizeof(double)*2);  
    Resultado=(double*)malloc(sizeof(double)*N*N);

    MPI_Scatter(A,N*N/numProc,MPI_DOUBLE,AA,N*N/numProc,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(B,N*N/numProc,MPI_DOUBLE,BB,N*N/numProc,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
    MPI_Bcast(C,N*N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(D,N*N,MPI_DOUBLE,0,MPI_COMM_WORLD);

    suma[0] = 0;
    suma[1] = 0;

     //Multiplicación de matrices A*C y B*D
    #pragma omp parallel 
    {
        //Multiplicacion de matriz A*C
        #pragma omp for private(i,j,k,acumulador) reduction(+:suma[0]) nowait
        for(i=0;i<N/numProc;i++){         
            for(j=0;j<N;j++){
                acumulador = 0;
                for(k=0;k<N;k++){
                    acumulador += (AA[i*N+k]*C[k+N*j]);  
                       suma[0] += AA[i*N+k]; 
                }
                AC[i*N+j] = acumulador;  
            }   
        }  

        //Multiplicacion de matriz B*D
        #pragma omp for private(i,j,k,acumulador2) reduction(+:suma[1])
        for(i=0;i<N/numProc;i++){         
            for(j=0;j<N;j++){
                acumulador2 = 0; 
               
                for(k=0;k<N;k++){
                    acumulador2 += (BB[i*N+k]*D[k+N*j]); 
                    suma[1] += BB[i*N+k];
                    
                }
                BD[i*N+j] = acumulador2;
            }  
        }    
    }   
    
    MPI_Allreduce(suma,vectorSuma,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  

    //Promedio de A y B
    promA = (vectorSuma[0]/(N*N));
    promB = (vectorSuma[1]/(N*N));

    //Aplica el promedio a matriz AC y la suma de la matriz resultante
    #pragma omp parallel for private(i,j)
     
        for(i=0;i<N/numProc;i++){         
            for(j=0;j<N;j++){
               
                    //Primer termino = PromB * AC
                  AC[i*N+j] = AC[i*N+j]*promB; 
                    //Segundo termino = PromA * BD
                  BD[i*N+j] = BD[i*N+j]*promA;
                 //Resultado = Primer termino + Segundo termino
                 Resultado[i*N+j] = (AC[i*N+j] + BD[i*N+j]); 
            }  
        
        } 

    MPI_Gather( Resultado,N*N/numProc,MPI_DOUBLE,R,N*N/numProc,MPI_DOUBLE,0,MPI_COMM_WORLD);    

    //Liberación de memoria
    free(A); 
    free(B);
    free(C);
    free(D);
    free(AC);
    free(BD);
    free(suma); 
    free(vectorSuma);
    free(Resultado);    

}
