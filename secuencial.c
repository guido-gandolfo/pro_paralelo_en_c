
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

void secuencial(int);

int main(int argc, char* argv[2]){
  
    int M,i,j,k,promA,promB;       
  

    if (argc < 2){

        printf("\n Falta un argumento:: N dimensión de la matriz \n");
           
    } 

    M = atoi(argv[1]);

    secuencial(M);
    
    return 0;
} 

void secuencial(int N){
    int i,j,k,promA,promB;
    double *A, *B, *C, *D, *AC, *BD, *Resultado;  
    double tiempo, tiempoTotal=0;
    register double acumulador = 0, acumulador2=0;
    double *suma;
    double check = 1;

    A=(double*)malloc(sizeof(double)*N*N);
    B=(double*)malloc(sizeof(double)*N*N); 
    C=(double*)malloc(sizeof(double)*N*N);
    D=(double*)malloc(sizeof(double)*N*N); 

    AC=(double*)malloc(sizeof(double)*N*N);  
    BD=(double*)malloc(sizeof(double)*N*N);
    suma=(double*)malloc(sizeof(double)*2); 
     
    Resultado=(double*)malloc(sizeof(double)*N*N);

    for(i=0;i<N;i++){
   		for(j=0;j<N;j++){
			A[i*N+j]= 1; //Por filas
            B[i*N+j]= 1; //Por filas
            C[i+j*N]= 1; //Por columnas
            D[i+j*N]= 1; //Por columnas
        }
    }

    tiempo = dwalltime();

    suma[0] = 0;
    

        //Multiplicacion de matriz A*C
       
        for(i=0;i<N;i++){         
            for(j=0;j<N;j++){
                acumulador = 0;
                for(k=0;k<N;k++){
                    acumulador += (A[i*N+k]*C[k+N*j]);  
                       suma[0] += A[i*N+k]; 
                }
                AC[i*N+j] = acumulador;  
            }   
        }  

        //Multiplicacion de matriz B*D
       suma[1] = 0;
        for(i=0;i<N;i++){         
            for(j=0;j<N;j++){
                acumulador2 = 0; 
               
                for(k=0;k<N;k++){
                    acumulador2 += (B[i*N+k]*D[k+N*j]); 
                    suma[1] += B[i*N+k];
                    
                }
                BD[i*N+j] = acumulador2;
            }  
        }    
     

    //Promedio de A y B
    promA = (suma[0]/(N*N));
    promB = (suma[1]/(N*N));

       //Aplica el promedio a matriz AC
  
        for(i=0;i<N;i++){         
            for(j=0;j<N;j++){
               
                    //Primer termino = PromB * AC
                  AC[i*N+j] = AC[i*N+j]*promB; 
                    //Segundo termino = PromA * BD
                  BD[i*N+j] = BD[i*N+j]*promA;
                 //Resultado= Primer termino + Segundo termino
                 Resultado[i*N+j] = (AC[i*N+j] + BD[i*N+j]); 
            }  
        
        } 

  
    tiempoTotal += dwalltime()- tiempo;  

    printf("Tiempo total: %.6f \n",tiempoTotal);
        //Verifica el resultad
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
	            check=check&&(Resultado[i*N+j]==(double)2*N*N);
            }
        }   

        if(check){
            printf("Resultado Correcto\n");
            printf("\n");
        }else{
            printf("Resultado Incorrecto\n");
            printf("\n");
        }
    
    free(A); 
    free(B);
    free(C);
    free(D);
    free(AC);
    free(BD);
    free(suma); 
    free(Resultado);  

}