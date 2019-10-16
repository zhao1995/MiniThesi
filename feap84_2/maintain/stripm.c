#include <stdio.h>
#define CR '\015'

main(argc,argv)
int argc;
char *argv[];
{

        FILE *fin, *fout, *fopen();
         int c;
        char cp[100];
        int  i;

        if(argc < 2) printf("ERROR -- Need to specify a file\n");

        else {

          for(i = 1; i <= argc-1; i++) {
            fin = fopen( argv[i], "r" );
            fout = fopen( "temp", "w" );

            while( (c = getc(fin)) != EOF) if (c != CR)  putc(c,fout);

            fclose(fin); fclose(fout);

            sprintf(cp,"cp temp %s",argv[i]);
            system(cp); system("rm temp");
          }


        }

}
