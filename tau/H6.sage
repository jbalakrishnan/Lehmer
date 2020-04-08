#Computing integral points on Y^2 = 5X^22 - 4*691

K.<a> = QuadraticField(-691) #K = Q(\sqrt{-691}), a=\sqrt{-691}
taus = K.embeddings(K) #Two complex embeddings
O = K.maximal_order() #Ring of integers
I = K.class_group() #Ideal class group
AI = I.gen().ideal() #Generator of I

g = AI^5 #Principal idea g_i
gbar = (AI^5).apply_morphism(taus[1]) #Conjugate of principal ideal g_i

gneg22 = g^(-22) #Fractional ideal (-10293606293232974813990340335913299876/77037197775489434122239117703397092741524065928615527809597551822662353515625*a + 61806078876679686569163192573518341803/77037197775489434122239117703397092741524065928615527809597551822662353515625)

gbarneg22 = gbar^(-22) #Fractional ideal (10293606293232974813990340335913299876/77037197775489434122239117703397092741524065928615527809597551822662353515625*a + 61806078876679686569163192573518341803/77037197775489434122239117703397092741524065928615527809597551822662353515625)

#Elements in \Gamma
r1 = gneg22.gens_reduced()[0]
r2 = gbarneg22.gens_reduced()[0]
d = r1.denominator()
if r1*r2 != 1/d:
    r2 = -r2

K.<x,y> = K[] #The ring of polynomials in x,y over K

#Thue equation of degree 22
g(x,y)=1/a*(r1*(x+y*(1/2+1/2*a))^22-r2*(x+y*(1/2-1/2*a))^22)*d

#Calling PARI/GP to use the Thue equation solver
th=gp.thueinit(simplify(expand(g(x,1))),1)
gp.thue(th,4*d)