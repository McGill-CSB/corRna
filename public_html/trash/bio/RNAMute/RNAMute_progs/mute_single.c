
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define MAX 1000 /* max sequence length */

/*****************MAIN PROGRAM***********************/
int main()

{
   int i;
   char RNA[1000];
   char RNA_flex[1000];
   FILE *outfile;
   FILE *mutfile;
   FILE *our_sequence;

  
   scanf("%s",RNA);
   strcpy (RNA_flex, RNA);
   remove("result");
   system("rm -r -f myDir");
   system("rm -r -f  htmlDir");
   system("rm  -f our_sequence");
   system ("mkdir psPictures" );
   /******************************/
   system("rm -r -f jpgPictures");
   system("rm -r -f  htmlDir");
   system ("mkdir jpgPictures" );
   /******************************/
   system("date +\"DATE: %m/%d/%y  TIME: %H:%M:%S\n\" ");

   outfile = fopen("seq","w");
   our_sequence=fopen("our_sequence","w");
   fprintf(outfile,"%s", RNA_flex);
   fprintf(our_sequence,"%s", RNA_flex);
   fclose(outfile);
   /*handles the wild type before creating the mutations*/
   system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");/*puts in a temporary file(foo_file) the Shapiro Representation*/
   system ("bin/calcEig2 <foo_file.ct >> result");
   system("convert -resize 350x350 rna.ps jpgPictures/wild_typ.jpg"); 
   system("mv rna.ps psPictures/wild_typ.ps");

 char  buffer[200], s[] = "convert -resize 350X350 rna.ps jpgPictures/", c[] =".jpg";/*initializes a buffer to create a jpeg picture*/
 char  buffer2[200], s2[] = "mv rna.ps psPictures/", c2[] =".ps";/*initializes a buffer to create a ps picture*/
   for (i = 0; i < strlen(RNA); i = i+1)/* passes over the RNA sequence and for each nucleotide makes three different  mutations*/ 
     {
	  if (RNA_flex[i]=='A'|| RNA_flex[i]=='a') {/*if the current nucleotide is A*/
		  int j;
		  int j2;
		 RNA_flex[i]='U';/*first mutation*/
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'U');
         fclose (mutfile);
         system ("cat mut >>result");/*writes the mut name in the result file*/
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");/*creates the picture of the mutation*/
         system ("bin/calcEig2 <foo_file.ct >> result");/*calculates the 2nd eigen-value and write it
														  in the result file in the line mute name */
		 /*give a name to the jpeg picture file*/
         j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-U");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
		 /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-U");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);
         
         RNA_flex[i]='C';/*second mutation*/
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'C');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		 /*give a name to the jpeg picture file*/
         j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-C");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
		 /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-C");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);
         
         RNA_flex[i]='G';/*third mutation*/
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'G');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		/*give a name to the jpeg picture file*/
         j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-G");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
		 /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-G");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);

         strcpy (RNA_flex, RNA);                      
       }
       if (RNA_flex[i]=='U'|| RNA_flex[i]=='T'|| RNA_flex[i]=='u' || RNA_flex[i]== 't') {/*if the current nucleotide is U*/
		int j;
		int j2;
         RNA_flex[i]='A';
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'A');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		 /*give a name to the jpeg picture file*/
         j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-A");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
		 /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-A");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);
         
         RNA_flex[i]='C';
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'C');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		 /*give a name to the jpeg picture file*/
         j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-C");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
		 /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-C");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);
         
         RNA_flex[i]='G';
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'G');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		 /*give a name to the jpeg picture file*/
         j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-G");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
		 /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-G");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);

         strcpy (RNA_flex, RNA);            
       }
       if (RNA_flex[i]=='C'|| RNA_flex[i]=='c') {/*if the current nucleotide is C*/
		   int j;
		   int j2;

         RNA_flex[i]='U';
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'U');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		 /*give a name to the jpeg picture file*/
		 j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-U");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
		 /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-U");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);
         
         RNA_flex[i]='A';
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'A');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		 /*give a name to the jpeg picture file*/
         j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-A");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
         /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-A");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);

         RNA_flex[i]='G';
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'G');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		 /*give a name to the jpeg picture file*/
         j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-G");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
		 /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-G");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);

         strcpy (RNA_flex, RNA); 
       }
       if (RNA_flex[i]=='G'|| RNA_flex[i]=='g') {/*if the current nucleotide is G*/
		   int j;
		   int j2;

         RNA_flex[i]='U';
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'U');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		  /*give a name to the jpeg picture file*/
		 j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-U");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
         /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-U");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);

         RNA_flex[i]='C';
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'C');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		 /*give a name to the jpeg picture file*/
         j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-C");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
		 /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-C");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);
         
         RNA_flex[i]='A';
         outfile = fopen("seq","w");
         fprintf (outfile,"%s", RNA_flex);
         fclose (outfile);
         mutfile = fopen("mut","w");
         fprintf (mutfile,"%d-",i+1);
         fprintf (mutfile,"%c",'A');
         fclose (mutfile);
         system ("cat mut >>result");
         system ("bin/RNAfold -noLP <seq | bin/b2Shapiro > foo_file.ct");
         system ("bin/calcEig2 <foo_file.ct >> result");
		/*give a name to the jpeg picture file*/
         j  = sprintf( buffer, "%s", s );
         j += sprintf( buffer + j, "%d", i+1);
		 j += sprintf( buffer + j, "%s", "-A");
         j += sprintf( buffer + j, "%s", c );
		 system ( buffer);
		 /*give a name to the ps picture file*/
         j2  = sprintf( buffer2 , "%s", s2 );
         j2 += sprintf( buffer2 + j2, "%d", i+1);
		 j2 += sprintf( buffer2 + j2, "%s", "-A");
         j2 += sprintf( buffer2 + j2, "%s", c2 );
		 system ( buffer2);

         strcpy (RNA_flex, RNA);        
       }
     }
   system("date +\"DATE: %m/%d/%y  TIME: %H:%M:%S\n\" ");
   
   return;
}







