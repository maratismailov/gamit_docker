/*  read modified USGS net file and create FONDA format site-table

     unit:
         x, y, z : km
        vx,vy,vz : mm|year
        time     : year
*/

   read_usgs_new_net(ifil1,ifil2,ifil3,nnet,netfmt)
{
# include <stdio.h>
   int ld,md,i,nnet,netid,lm,mm;
   char netid[3],*netfmt[49],name[7],id1,id2;
   char string[8];
   FILE *ifil1,*ifil2,*ifil3;
   double ve,vn,vu,ht,sla,slo,s1,s2;
   
   fprintf(ifil3,"%s\n"," Network distribution");
   ve = 0.0d0;
   vn = 0.0d0;
   vu = 0.0d0;
/*
   default format
*/
   if (strcmp(netfmt[0],'*') == 0) netfmt = "(%s%s%s%d%d%d%f\n)";
   fscanf(ifil1,"%d%d%s%s%s\n",&nobs,&nnet,netid,string,id1);
   fscanf(ifil1,"%f%d%f\n",&s1,&ld,&s2);
   for(i = 0; i < nnet; ++i) {
/*    shift format line by line  */
      fscanf(ifil1,netfmt,ntid,name,id1,&ld,&lm,&sla,&md,&mm,&slo,id2,&ht);
      s1 = dble(ld)+dble(lm)/60.0d0+sla/3600.0d0
c     westward longitude
      s2 = (dble(md)+dble(mm)/60.0d0+slo/3600.0d0)
      id1 = 'N'
      if (ld.lt.0) then
         id1 = 'S'
         ld = -ld
         s1 = -s1
      endif
      id2 = 'E'
      if (md.lt.0) then
         id2 = 'W'
         md = -md
         s2 = -s2
      endif
         write (ifil2,30) netid,name,
     .         id1,ld,lm,sla,id2,md,mm,slo,ht,ve,vn,vu
         write (ifil3,40) s2,s1,ve,vn,name
c         if(i.lt.nnet) read (ifil1,'(2x)')
 20   continue
c
c10   format (a4,2x,a8,3x,a1,2i3,f9.5,i4,i3,f9.5,3x,a1,f8.2)
 30   format (a4,1x,a8,4x,a1,i2,1x,i2,1x,f8.5,
     .        1x,a1,i3,1x,i2,1x,f8.5,f13.4,3f8.2)
 40   format (1x,2f15.8,2f10.4,3x,a8)
c
      return
      end
