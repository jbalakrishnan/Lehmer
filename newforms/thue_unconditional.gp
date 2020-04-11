read("F2m.gp");
dl = [[[6],7],[[10],11],[[6,12],13],[[16],17],[[18],19],[[10,22],23],[[6, 28],29],[[30],31],[[18, 36],37]];
y = 1
default(parisizemax, 10^9)
{
for(i= 1, #dl,
   ell = dl[i][2];
   print("ell = ", ell);
   print("--------------------------------------------");
   for(j=1, #dl[i][1],
       dmin1 = dl[i][1][j];
       print("for d - 1 =  ", dmin1);
       T = getwalltime();
       tnf = thueinit(eval(A[dmin1/2][2]),1);
       print("solutions for ell = ", ell);
       print (thue(tnf, ell));
       print("solutions for -ell = -", ell);
       print(thue(tnf, -ell));
       print(getwalltime() - T);
   print("++++++++++++++++++++++++++++++++++++++++++++")));
   
}