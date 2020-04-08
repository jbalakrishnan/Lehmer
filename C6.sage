#Computing integral points on Y^2 = X^11 - 691

K.<a> = QuadraticField(-691) #K = Q(\sqrt{-691}), a=\sqrt{-691}
taus = K.embeddings(K) #Two complex embeddings
O = K.maximal_order() #Ring of integers
I = K.class_group() #Ideal class group
AI = I.gen().ideal() #Generator of I

g = AI^5 #Principal idea g_i
gbar = (AI^5).apply_morphism(taus[1]) #Conjugate of principal ideal g_i

gneg11 = g^(-11) #Fractional ideal (495538587136045198/277555756156289135105907917022705078125*a + 10386281270974935631/277555756156289135105907917022705078125)
gbarneg11 = gbar^(-11) #Fractional ideal (-495538587136045198/277555756156289135105907917022705078125*a + 10386281270974935631/277555756156289135105907917022705078125)

#Elements in \Gamma
r1 = gneg11.gens_reduced()[0]
r2 = gbarneg11.gens_reduced()[0]
d = r1.denominator()
if r1*r2 != 1/d:
    r2 = -r2

K.<x,y> = K[] #The ring of polynomials in x,y over K

#Thue equation of degree 11
g(x,y)=1/a*(r1*(x+y*(1/2+1/2*a))^11-r2*(x+y*(1/2-1/2*a))^11)*d

#Calling PARI/GP to use the Thue equation solver
th=gp.thueinit(simplify(expand(g(x,1))),1)
gp.thue(th,2*d)