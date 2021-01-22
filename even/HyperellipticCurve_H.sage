import itertools
from pprint import pprint
import copy

gp.parisizemax=512000000000000;

def valuation(P,I):
     F = list(I.factor());
     for i in range(0,len(F)):
         if not F[i][0].is_coprime(P):
             return F[i][1];
     return 0;

def group_by_conjugates(primes,taus):
     already_viewed = [];
     conj_pairs = [];
     for i in range(0,len(primes)):
         if i not in already_viewed:
             p = primes[i];
             already_viewed.append(i);
             found = False;
             for j in range(i+1,len(primes)):
                 q = primes[j].apply_morphism(taus[1]);
                 if (j not in already_viewed) and (not q.is_coprime(p)) and (not found):
                     found = True;
                     already_viewed.append(j);
                     conj_pairs.append([p,primes[j]]);
             if not found:
                 conj_pairs.append([p,p]);
     return conj_pairs;


    


[cstart,dstart, lstart, sgn] = [5,11,29^2, -1];
[C,D, n] = [cstart,-sgn*lstart,dstart] #H-curve: x^2 + D = C*y^n. EX: x^2+29^2=5*y^11
print("Eqn: x^2 +",D,"=",C,"y^",n)
print("(c,d,l,+/-) = (",cstart,dstart,lstart,sgn,")")


d = squarefree_part(D);
q = sqrt(abs(D/d));


var('x')

K.<a>=NumberField(x^2+d);
taus=K.embeddings(K);
O = K.maximal_order(); 
H, h = K.class_group(), K.class_number();
U = UnitGroup(K);

U_reps = [];
for u in U.gens():
     i=0;
     while i<n:
         if u**i not in U_reps:
             U_reps.append(u**i);
             U_reps.append(U.gens()[0]*u**i); 
         i += 1;
print("Quotient group--U/U^n:",U_reps) 
#print()



discK = K.discriminant();
w = K.ideal(1).basis()[1]; 
wbar = taus[1]((K.ideal(1).basis()[1]));


TDC =2*D*C;
TDC_ideal = K.ideal(TDC);
S_2DC = list(TDC_ideal.prime_factors());
print("Prime ideals dividing (2DC):")
for x in S_2DC: print(x)
print()


conj_primes = group_by_conjugates(S_2DC,taus); # prime factors of (2DC)
print("list of the primes", conj_primes);

F= list(factor(TDC));
Fnew = sorted(set([f[0] for f in F for _ in range(f[1])])); 
print("Prime ideals dividing (2DC)",Fnew)
f=len(Fnew);

#print(conj_primes)


SP = [];
for i in range(f):
    setname = [];
    for k in (0..n-1):
        for l in (0..n-1): 
            if min(k,l)<=valuation(conj_primes[i][0],K.ideal(4*q*a)) and mod(k+l-valuation(conj_primes[i][0],K.ideal(5)),n)==0:
                if conj_primes[i][0].is_coprime(conj_primes[i][1]): setname.append(conj_primes[i][0]^(k)*conj_primes[i][1]^(l));
                elif k==l: setname.append(conj_primes[i][0]^(k));
    SP.append(setname)
#    if not setname == []: 
#        SP.append(setname);
if [] in SP:
    SP.remove([])
print("SP",SP)


a_plus_ideals = [];

if len(SP) == 1:
    for p in list(itertools.product(*[SP[0]])):
        a_plus_ideals.append(p)
if len(SP) == 2:
    for p,k in list(itertools.product(*[SP[0],SP[1]])):
        a_plus_ideals.append(p*k)
if len(SP) == 3:
    for p,k,r in list(itertools.product(*[SP[0],SP[1],SP[2]])):
        a_plus_ideals.append(p*k*r)
if len(SP) == 4:
    for p,k,r in list(itertools.product(*[SP[0],SP[1],SP[2]])):
        a_plus_ideals.append(p*k*r)


    
print("aplus",a_plus_ideals)


H_elts = [i.ideal() for i in list(H)];

gammas = [1];
for a_ideal in a_plus_ideals:
     for g in H_elts:
         ideal = a_ideal*g^(-n);
         if ideal.is_principal():
             gammas.append(ideal.gens_reduced()[0]);

Gamma = [];
for g in gammas:
     for u in U_reps:
         Gamma.append([taus[0](K(u)*g), taus[1](K(u)*g)]);
#print()
#print(Gamma)
print()


R.<X,Y> = K['X','Y'];

thues = [];
for r in Gamma:
     thues.append(1/a*(r[0]*(X+Y*w)**n-r[1]*(X+Y*wbar)**n)); 
        
print(thues) #All Thue equations derived from the H-curve.
SOLNS = [];


count = 0;
for r,f in zip(Gamma,thues):
     print(count+1,"of",len(thues))
     count = count + 1;
     fnew = simplify(f.denominator()*expand(f(Y=1)));
     th=gp.thueinit(fnew,0); #grh assuming is 0 not grh assuming is 1
     #print(2*q*f.denominator(),th[1])
     th_solns = gp.thue(th,2*q*f.denominator());
     #print(th_solns)
     #print()
     for s in th_solns:
         SOLNS.append([s,r]);
#print()
#print('SOLUTION',SOLNS)

FINAL_SOLNS = [];
for e in SOLNS:
     x = 1/2 * (e[1][0]*(pari(e[0][1])+pari(e[0][2])*w)**n +
e[1][1]*(pari(e[0][1])+pari(e[0][2])*wbar)**n);
     FINAL_SOLNS.append(x);
print("Final Solutions",FINAL_SOLNS)
