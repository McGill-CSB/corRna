#include <string.h>
#include <stdio.h>
#include "nr.h"
#include "nrutil.h"

int sizeOfLaplacianMatrix =0;
char s[200];
int nrot;


/*this function takes as a parameter the Shapiro representation of
  the RNA second structure and return the Laplacian matrix that
  matches the Shapiro representation*/
float **makeLaplacianMatrix(const char *shapiro){
  float **matrix;
  int sizeOfMatrix , sizeOfShapiro;
  int i,k;
  int numOfParens=0;
  int node1=0, node2=0;
  int isNode=0;
  int isS=0;
  int count;
  char newShapiro[200];
  int j=0;
  
  sizeOfShapiro = strlen(shapiro);
  sizeOfMatrix = 0;
   
  for (i=0; i<sizeOfShapiro; i++){
    if ((((int) shapiro[i])>64 && ((int) shapiro[i])<91
         && shapiro[i]!='S' && shapiro[i]!='E')&&
        !((shapiro[i]=='I') && (shapiro[i+1]=='2') && (shapiro[i+2]==')'))&&
        !((shapiro[i]=='B') && (shapiro[i+1]=='1') && (shapiro[i+2]==')')))
      sizeOfMatrix++;
  }
  
  matrix = (float**) malloc(sizeOfMatrix*sizeof(float*));
  
  for (i=0; i<sizeOfMatrix; i++){
    matrix[i] = (float*) calloc(sizeOfMatrix,sizeof(float));
  }

  count = -1;
  for (i=0; i<sizeOfShapiro; i++){
    if (shapiro[i] == 'S') isS=1;
    else if ((shapiro[i]=='I') && (shapiro[i+1]=='2') && (shapiro[i+2]==')')){
      newShapiro[j++] = 'I';
    }
    else if ((shapiro[i]=='B') && (shapiro[i+1]=='1') && (shapiro[i+2]==')')){
      newShapiro[j++] = 'B';
    }
    else if (((int) shapiro[i])>64 && ((int) shapiro[i])<91){
      count++;
      isNode = 1;
    }
    else if (((int) shapiro[i])>47 && ((int) shapiro[i])<58){
    }
    else if (isNode!=0){
      isNode = 0;
      newShapiro[j++] = (char) ((count / 10)+48);
      newShapiro[j++] = (char) ((count % 10)+48);
      i--;
    }
    else if (isS!=0){
      isS = 0;
      newShapiro[j] = 'S';
      j++;
      i--;
    }
    else{
      newShapiro[j] = shapiro[i];
      j++;
    }
  }

  newShapiro[j] = '\0';
  
  for (i=0; i<j; i++){
    if (((int) newShapiro[i])>47 && ((int) newShapiro[i])<58){
      node1 = 10*(((int) newShapiro[i]) - 48) + (((int) newShapiro[i+1]) - 48);
      i = i+5;
      numOfParens = 0;
      for (k=i; k<j; k++){
        if (newShapiro[k] == '(') numOfParens ++;
        if (newShapiro[k] == ')') numOfParens --;
        if ((newShapiro[k] == 'I' || newShapiro[k] == 'B') && numOfParens <=0) {
          numOfParens=0;
          k = k+3;
        }
        if (((int) newShapiro[k])>47 &&
            ((int) newShapiro[k])<58 && numOfParens==0){
        node2 =
            10*(((int) newShapiro[k]) - 48) + (((int) newShapiro[k+1]) - 48);
          matrix[node1][node1]++;
          matrix[node2][node2]++;
          matrix[node1][node2]=-1;
          matrix[node2][node1]=-1;
          break;
        }
      }
      i--;
    }
  }
  
  strcpy(s,newShapiro);
  sizeOfLaplacianMatrix = sizeOfMatrix;
  return matrix;
}


/*this function takes as a parameter the Laplacian matrix
  and returns the array of its eigenvalues*/
float eig2Value(float **matrix1){
  int i,j;
  float *d, *r, **v, **e;
  float matrix2[sizeOfLaplacianMatrix][sizeOfLaplacianMatrix];
  
  for (i=0; i<sizeOfLaplacianMatrix; i++)
    for (j=0; j<sizeOfLaplacianMatrix; j++)
      matrix2[i][j] = matrix1[i][j];


  d = vector(1,sizeOfLaplacianMatrix);
  r = vector(1,sizeOfLaplacianMatrix);
  v = matrix(1,sizeOfLaplacianMatrix,1,sizeOfLaplacianMatrix);
  e = convert_matrix(&matrix2[0][0],1,sizeOfLaplacianMatrix,
                     1,sizeOfLaplacianMatrix);


  jacobi(e,sizeOfLaplacianMatrix,d,v,&nrot);
  eigsrt(d,v,sizeOfLaplacianMatrix);
  return d[sizeOfLaplacianMatrix-1];
}



int main(){
  char input[200];
  float eigVal2;
  int i,j;
  scanf("%s",input);
  float **matrix = makeLaplacianMatrix(input);
  
 
  eigVal2 = eig2Value(matrix);

  printf("%12.6f  %d  %s\n",eigVal2,sizeOfLaplacianMatrix,input);

}
