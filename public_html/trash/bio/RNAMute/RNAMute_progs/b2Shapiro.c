#include "RNAstruct.h"
#include <stdio.h>
#include <string.h>

int main(int argc,char *argv[]){
  char inString[1000];
  char deltaG[15];
  char *shapiro = (char*) malloc(1000*sizeof(char));
  int i;
  FILE *result;

   
  result = fopen("result","a"); 
  scanf("%s",inString);
  scanf("%s",inString);
  fprintf(result," %s ",inString);
  scanf("%s",deltaG);
  fprintf(result,"%s",deltaG);
  if (deltaG[0] == '(' && strlen(deltaG) == 1){
  	scanf("%s",deltaG);
  	fprintf(result,"%s",deltaG);
  }
  
  shapiro = b2Shapiro(inString);
  printf("%s\n",shapiro);
  fclose(result);
  free(shapiro);
}
