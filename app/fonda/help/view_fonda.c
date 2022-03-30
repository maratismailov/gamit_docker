/*
      FONDA_VIEW  --- a guid book for FONDA learner  
*/
#include <stdio.h>
#include <ctype.h>

main ()
{
   int opt,pac,loop,vers,inf;
   char name[80],command[80];

/* set up links */
   sprintf(name,"%s",getenv("FONDA_H"));
   sprintf(command,"ln -s %shelp FOHELP",name);
   system(command);
   sprintf(command,"ln -s %sexample EXAMPLE",name);
   system(command);
   sprintf(command,"ln -s %sbin UTIL",name);
   system(command);

/* general introduction to FONDA  */
   system("more FOHELP/fonda.help");

/* choose description mode  */
   printf("\n Please choose text mode:\n");
   printf("   1 = terse,  2 = verbose,  3 = movie(not valid now)\n");
   scanf("%s",name);
   vers = atoi(name);

   loop = 0;
   while (loop<=20)
   {
   static char line[60] = " ----------------------------------------------------------";
   static char lin1[60] = " ----------------------- options --------------------------";
   static char lin2[60] = " ----------------------- example --------------------------";
   char c[1];
   loop = loop+1;
   printf("\n%s\n",lin1);
   printf(" Please type name or pick up number to see the packages:\n");
   printf(" package : 1=maked  2=solvem  3=diagno  4=fmodel  5=disply  6=design  7=utility\n");
   printf(" other   : 8=comment 9=mail  10=exit\n");
   printf("%s\n",line);
   scanf("%s",name);
   pac = atoi(name);
   printf("%s\n",line);
/* show structure of MAKED */
   if(strcmp("maked",name)==0||strcmp("MAKED",name)==0||pac==1)
   {
      static char *option[] = {"more EXAMPLE/maked/maked.ngs.drv",
         "more EXAMPLE/maked/maked.ngs.in ",
         "more EXAMPLE/maked/maked.ngs.map",
         "more EXAMPLE/maked/maked.ngs.out",
         "more EXAMPLE/maked/vmodel.ngs ",
         "more EXAMPLE/maked/ngsdat.pri ",
         "more EXAMPLE/maked/quake.mdl "};
      if (vers == 1) 
         system("head -13 FOHELP/maked.help");
      else
         system("more FOHELP/maked.help");
      opt = 1;
      while(opt<=7 && opt>0)
      {
         printf("\n%s\n",lin1);
         printf(" Please pick up a number to see the examples:\n");
         printf(" 1=driving file,  2=input file,    3=mapping file\n");
         printf(" 4=output file,   5=vmodel file,   6=priori file\n");
         printf(" 7=quake mdlfil,  8=other package, 9=exit\n");
         printf("%s\n",line);
         scanf("%s",c);
         opt = atoi(c);
         if(opt<=7) printf("%s\n",lin2);
         if(opt<=7 && opt>0) system(option[opt-1]);
      }
      if(opt==8) continue;
      if(opt>=9 || opt<=0) break;
   }
/* show structure of SOLVEM */
   if(strcmp("solvem",name)==0||strcmp("SOLVEM",name)==0||pac==2)
   {
      static char *option[] = {"more EXAMPLE/solvem/solvem_slt.drv",
         "more EXAMPLE/solvem/solvem.ngs.dat",
         "more EXAMPLE/solvem/solvem_slt.map",
         "more EXAMPLE/solvem/solvem_slt.out",
         "more EXAMPLE/solvem/solvem_slt.res",
         "more EXAMPLE/solvem/solvem_slt.str",
         "more EXAMPLE/solvem/solvem_slt.evt",
         "more EXAMPLE/solvem/quake.list.salto",
         "more EXAMPLE/solvem/outer_salto.list"};
      if (vers == 1) 
         system("head -16 FOHELP/solvem.help");
      else
         system("more FOHELP/solvem.help");
      opt = 1;
      while(opt<=9 && opt>0)
      {
         printf("\n%s\n",lin1);
         printf(" Please pick up a number to see the examples:\n");
         printf(" 1=driving file,  2=input file,    3=mapping file\n");
         printf(" 4=solution file, 5=residual file, 6=strain file\n");
         printf(" 7=event file,    8=earthquake file\n");
         printf(" 9=outer coordinate list file\n");
         printf(" 10=other package, 11=exit\n");
         printf("%s\n",line);
         scanf("%s",c);
         opt = atoi(c);
         if(opt<=9) printf("%s\n",lin2);
         if(opt<=9 && opt>0) system(option[opt-1]);
      }
      if(opt==10) continue;
      if(opt>=11 || opt<=0) break;
   }
/* show structure of DIAGNO */
   if(strcmp("diagno",name)==0||strcmp("DIAGNO",name)==0||pac==3)
   {
      if (vers == 1) 
         system("head -10 FOHELP/diagno.help");
      else
         system("more FOHELP/diagno.help");
   }
/* show structure of FMODEL */
   if(strcmp("fmodel",name)==0||strcmp("FMODEL",name)==0||pac==4)
   {
      if (vers == 1) 
         system("head -10 FOHELP/fmodel.help");
      else
         system("more FOHELP/fmodel.help");
   }
/* show structure of DISPLY */
   if(strcmp("disply",name)==0||strcmp("DISPLY",name)==0||pac==5)
   {
      if (vers == 1) 
         system("head -10 FOHELP/disply.help");
      else
         system("more FOHELP/disply.help");
      continue;
   }
/* show structure of DESIGN */
   if(strcmp("design",name)==0||strcmp("DESIGN",name)==0||pac==6)
   {
      if (vers == 1) 
         system("head -10 FOHELP/design.help");
      else
         system("more FOHELP/design.help");
      continue;
   }
/* show structure of UTILITY */
   if(strcmp("utility",name)==0||strcmp("UTILITY",name)==0||pac==7)
   {
      static char *option[] = {"UTIL/align_frame",
         "UTIL/compare_coor",
         "UTIL/get_line",
         "UTIL/get_v_rel",
         "UTIL/net_update",
         "UTIL/plot_line",
         "more UTIL/history.note"};
      if (vers == 1) 
         system("head -7 FOHELP/utility.help");
      else
         system("more FOHELP/utility.help");
      opt = 1;
      while(opt<=7 && opt>0)
      {
         printf("\n%s\n",lin1);
         printf(" Please pick up a number to see the examples:\n");
         printf(" 1=align_frame,   2=compare_coor,  3=get_line\n");
         printf(" 4=get_v_rel,     5=net_update,    6=plot_line\n");
         printf(" 7=history,       8=other package, 9=exit\n");
         printf("%s\n",line);
         scanf("%s",c);
         opt = atoi(c);
         if(opt<=7) printf("%s\n",lin2);
         if(opt<=7 && opt>0) system(option[opt-1]);
      }
      if(opt==8) continue;
      if(opt>=9 || opt<=0) break;
   }

/* read or write comments  */
   if(strcmp("comment",name)==0||strcmp("COMMENT",name)==0||pac == 8)
   {
         system("vi FOHELP/comment.help");
      continue;
   }

/* mail to supervisor  */
   if(strcmp("mail",name)==0||strcmp("MAIL",name)==0||pac == 9)
   {
         system("mail dong@freia.jpl.nasa.gov king@chandler.mit.edu tvd@gemini.gsfc.nasa.gov andrea@cobra.jpl.nasa.gov");
         system("sleep 20");
      continue;
   }

/* exit or retry */
   if(strcmp("exit",name)==0 || pac<=0 || pac>=10) break;

   }

/* remove links */
   system("rm -f FOHELP EXAMPLE UTIL");

   exit(0);
}
