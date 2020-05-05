import itertools
from pprint import pprint
import copy
import sys
sys.setrecursionlimit(1000000000);
gp.set_default('parisizemax',99000000000);

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

#Construction of the Equation x^2 + D = 5y^n
[dstart, lstart, sgn] = [ , , ]; #Y^2 = 5X^(2d) \pm 4*l
[D, n] = [-sgn*4*lstart,dstart] #x^2 + D = 5y^n
print("Eqn: x^2 +",D,"= 5y^",n)
#print()

d = squarefree_part(D);
q = sqrt(abs(D/d));

#Establish the Number Field K = Q(sqrt(-d))
var('x')
K.<a>=NumberField(x^2+d);
taus=K.embeddings(K);
O = K.maximal_order();
H, h = K.class_group(), K.class_number();
U = UnitGroup(K);

U_reps = [];
for u in U.gens():
    for i in (0..n-1):
        if u==U.gens()[0] and sgn==-1 and is_odd(n)==true and i==n-1:
            U_reps.append(U.gens()[0]^2);
        elif u not in U_reps and sgn==-1 and is_even(n)==true and i==n-1:
            U_reps.append(u);
            U_reps.append(u^2);
        elif u**i not in U_reps and sgn==1 and is_odd(n)==true:
            U_reps.append(U.gens()[1]**i);
        elif u**i not in U_reps and sgn==1 and is_even(n)==true:
            U_reps.append(u**i);
            U_reps.append(U.gens()[0]*u**i);


discK = K.discriminant();
w = K.ideal(1).basis()[1];
wbar = taus[1]((K.ideal(1).basis()[1]));

#Establish the set of primes S(2DC) = S(10D)
TwoDC = 10*D;
TwoDC_ideal = K.ideal(TwoDC);
S_2DC = list(TwoDC_ideal.prime_factors());
#print("Primes dividing <2DC>:")


conj_primes = group_by_conjugates(S_2DC,taus);



S2=[]
#prime 2
for k in (0..n-1):
     for l in (0..n-1):
         if min(k,l)<=valuation(conj_primes[0][0],K.ideal(4*q*a)) and mod(k+l-valuation(conj_primes[0][0],K.ideal(5)),n)==0:
            if conj_primes[0][0].is_coprime(conj_primes[0][1]): S2.append(conj_primes[0][0]^(k)*conj_primes[0][1]^(l));
            elif k==l: S2.append(conj_primes[0][0]^(k));
S5=[];
#prime 5
for k in (0..n-1):
     for l in (0..n-1):
         if min(k,l)<=valuation(conj_primes[1][0],K.ideal(4*q*a)) and mod(k+l-valuation(conj_primes[1][0],K.ideal(5)),n)==0:
            if conj_primes[1][0].is_coprime(conj_primes[1][1]): S5.append(conj_primes[1][0]^(k)*conj_primes[1][1]^(l));
            elif k==l: S5.append(conj_primes[1][0]^(k));

Slstart=[];
#prime \ell
for k in (0..n-1):
     for l in (0..n-1):
         if min(k,l)<=valuation(conj_primes[2][0],K.ideal(4*q*a)) and mod(k+l-valuation(conj_primes[2][0],K.ideal(5)),n)==0:
            if conj_primes[2][0].is_coprime(conj_primes[2][1]): Slstart.append(conj_primes[2][0]^(k)*conj_primes[2][1]^(l));
            elif k==l: Slstart.append(conj_primes[2][0]^(k));




a_plus_ideals = [];
for p,k,r in list(itertools.product(*[S2,S5,Slstart])):
    a_plus_ideals.append(p*k*r)
    

H_elts = [i.ideal() for i in list(H)];

gammas = [];
for a_ideal in a_plus_ideals:
    for g in H_elts:
        ideal = a_ideal*g^(-n);
        if ideal.is_principal():
            gammas.append(ideal.gens_reduced()[0]);

Gamma = [];
for g in gammas:
    for u in U_reps:
        Gamma.append([taus[0](K(u)*g), taus[1](K(u)*g)]);
print(Gamma)


R.<X,Y> = K['X','Y'];

thues = [];
for R in Gamma:
    thues.append(1/a*(R[0]*(X+Y*w)**n-R[1]*(X+Y*wbar)**n));
print(thues)
SOLNS = [];
for R,f in zip(Gamma,thues):
    th=gp.thueinit(simplify(f.denominator()*expand(f(Y=1))),0);## flag=0 under GRH assumotion; flag=1 without GRH assumotion ##
    th_solns = gp.thue(th,2*q*f.denominator());
    #print(th_solns)
    #print()
    for s in th_solns:
        SOLNS.append([s,R]);
#print()
#print('SOLUTION',SOLNS)




FINAL_SOLNS = [];
for e in SOLNS:
   x = 1/2 * (e[1][0]*(pari(e[0][1])+pari(e[0][2])*w)**n + e[1][1]*(pari(e[0][1])+pari(e[0][2])*wbar)**n);
   FINAL_SOLNS.append(x);
print(FINAL_SOLNS)
