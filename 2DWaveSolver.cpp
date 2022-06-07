#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>

void printMatrix(double **u, double M, double N, double T, float dt, int *fileCounter){
	char fileName[20];
	char fileNumber[10];
	sprintf(fileNumber, "%d", *fileCounter);
	
	strcpy(fileName, "results/data");
	strncat(fileName, fileNumber, strlen(fileNumber));
	strncat(fileName, ".txt", 4);
	
	FILE *fptr;
	fptr = fopen(fileName, "w");
	bool saveToFile = true;
	
	*fileCounter = *fileCounter + 1;
	
	printf(" --- TIME: %f of T: %f --- \n", dt,  T);
	for(int k = 0; k<= M+1; k++){
		for(int l = 0; l<=N+1; l++){
			//printf("%f\t", u[k][l]);
		}
		//printf("\n");
	}
	
	if(saveToFile){
		for(int k = 0; k<= M+1; k++){
			for(int l = 0; l<=N+1; l++){
				fprintf(fptr, "%f\t", u[k][l]);
			}
			fprintf(fptr, "\n");
		}
		fprintf(fptr, "\n");
	}
	
	
	fclose(fptr);
}

int main (int nargs, char** args)
{
  int M, N, L, i, j;
  double T, dx, dy, dt, x, y, t;
  double **up, **u, **um, **tmp, *up_data, *u_data, *um_data;
  int counter = 0;

  M = 48; N = 48; L = 100;
  
  T = 1.0; dx = 1./(M+1); dy = 1./(N+1); dt = T/L;
  
  
  
  printf(" --- dx = %f, dy = %f, dt = %f --- \n\n", dx, dy, dt);

  /* data allocation */
  up_data = (double*)malloc((M+2)*(N+2)*sizeof(double));
  u_data = (double*)malloc((M+2)*(N+2)*sizeof(double));
  um_data = (double*)malloc((M+2)*(N+2)*sizeof(double));

  up = (double**)malloc((N+2)*sizeof(double*));
  u = (double**)malloc((N+2)*sizeof(double*));
  um = (double**)malloc((N+2)*sizeof(double*));

  for (j=0; j<=N+1; j++) {
    up[j] = &(up_data[j*(M+2)]);
    u[j] = &(u_data[j*(M+2)]);
    um[j] = &(um_data[j*(M+2)]);
  }

  /* Enforce initial condition 1 */
  for (j=0; j<=N+1; j++) {
    y = j*dy;
    for (i=0; i<=M+1; i++) {
      x = i*dx;
      um[j][i] = sin(2*M_PI*(x+y));
      um[j][i] = 0;
    }
  }
  um[24][24] = 10;
  printMatrix(um, M, N, T, 0, &counter);
  

  /* Compute u at t=dt (enforcing initial condition 2) */
  for (j=1; j<=N; j++) {
    y = j*dy;
    for (i=1; i<=M; i++) {
      x = i*dx;
      u[j][i] = um[j][i]+((dt*dt)/(4*dx*dx))*(um[j][i-1]-2*um[j][i]+um[j][i+1])
	+((dt*dt)/(4*dy*dy))*(um[j-1][i]-2*um[j][i]+um[j+1][i]);
	//+dt*2*M_PI*cos(2*M_PI*(x+y));
    }
  }

  /* enforcing boundary conditions for t = dt */
  t = dt;
  
  for (j=0; j<=N+1; j++) {
    y = j*dy;
    // x = 0
    u[j][0] = sin(2*M_PI*(y+t));   
	u[j][0] = 0;  
    // x = 1
    u[j][M+1] = sin(2*M_PI*(1+y+t)); 
    u[j][M+1] = 0;
  }
  for (i=1; i<=M; i++) {
    x = i*dx;
    // y = 0
    u[0][i] = sin(2*M_PI*(x+t));   
	u[0][i] = 0;     
    // y = 1
    u[N+1][i] = sin(2*M_PI*(x+1+t)); 
    u[N+1][i] = 0;
  }
  
	printMatrix(u, M, N, T, t, &counter);


  /* main computation starting from t=2*dt */
  while (t<T) {
    t += dt;

    for (j=1; j<=N; j++)
      for (i=1; i<=M; i++)   /* interior points */
	up[j][i] = 2*u[j][i]-um[j][i]
	  +((dt*dt)/(2*dx*dx))*(u[j][i-1]-2*u[j][i]+u[j][i+1])
	  +((dt*dt)/(2*dy*dy))*(u[j-1][i]-2*u[j][i]+u[j+1][i]);

    /* enforcing boundary conditions */
    
    for (j=0; j<=N+1; j++) {
      y = j*dy;
      // x = 0
      up[j][0] = sin(2*M_PI*(y+t));    
	  up[j][0] = 0;    
      // x = 1
      up[j][M+1] = sin(2*M_PI*(1+y+t)); 
      up[j][M+1] = 0;
    }
    for (i=1; i<=M; i++) {
      x = i*dx;
      // y = 0 
      up[0][i] = sin(2*M_PI*(x+t));     
      up[0][i] = 0;
      // y = 1
      up[N+1][i] = sin(2*M_PI*(x+1+t)); 
      up[N+1][i] = 0; 
    }
	
    /* data shuffle */
    tmp = um;
    um = u;
    u = up;
    up = tmp;
    
    printMatrix(u, M, N, T, t, &counter);
  }

	

  /* clean up */
  free (up_data); free (up);
  free (u_data); free (u);
  free (um_data); free (um);

  return 0;
}
