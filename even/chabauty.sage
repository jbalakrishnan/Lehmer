def prune(D):
    newlist=[]
    p = D[-1][0].parent().order()
    GFp = GF(p)
    XF = D[0].parent()
    for x in D:
        if XF(x[0],-x[1]) not in newlist:
            newlist.append(x)
    return newlist

def find_Q(C,P,p):
    i = 0
    while True:
        try:
            Q = C.lift_x(ZZ(P[0]) + i*p)
            break
        except ValueError:
            i += 1
    try:
        if ZZ(Q[1][0]) != ZZ(P[1][0]):
            Q = C(Q[0],-Q[1])
    except TypeError:
        pass
    return Q

def effective_chabauty(X, gens, p, prec):
    K = Qp(p,prec)
    XK = X.change_ring(K)
    F = GF(p)
    XF = X.change_ring(F)
    disks = XF.rational_points()
    disks.remove(XF(0,1,0))
    disks = prune(disks)
    g = X.genus()
    m = []
    S.<x> = K['x']
    if len(gens)>0:
        v = vector(2*g*[0])
        for poly in gens:
            if (poly[0]).degree() == 1:
                Rs = [XK(poly[0].roots()[0][0],poly[1])]
            else:
                r  = [r[0] for r in S(poly[0]).roots()]
                Rs = [XK(rt,poly[1](rt)) for rt in r]
            for Ri in Rs:
                if Ri[1].valuation() < 0:
                    #Ri is in the disk at infinity
                    xinf,yinf = XK.local_coordinates_at_infinity(prec=prec)
                    t = xinf.parent().variable_name()
                    L = K[[t]]
                    L.set_default_prec(prec)
                    dx = xinf.derivative()
                    omegas = vector([dx*xinf^i/(2*yinf) for i in range(g)])
                    integrals = [i.integral() for i in omegas]
                    v = v + vector([i(Ri[0]^g/Ri[1]) for i in integrals] + g*[0]) #padding with 0s
                elif Ri[1].valuation() > 0:
                    #Ri is in Weierstrass disk
                    W = XK.find_char_zero_weier_point(Ri)
                    v = v + XK.tiny_integrals_on_basis(W,Ri)
                else:
                    v = v + XK.coleman_integrals_on_basis(XK(0,1,0),Ri)
            m = m + v[0:g].list()
        M = matrix(len(gens),g, m)
        M = M.transpose()
        ann = M.kernel().basis()
    else:
        ann = identity_matrix(g).rows()
    for locP in disks:
        print "working in residue disk corresponding to ", locP
        Q = find_Q(XK,locP,p)
        if (Q[1])%p == 0:
            Q = XK.find_char_zero_weier_point(Q)
        print Q
        x,y = XK.local_coord(Q,prec=prec)
        t = x.parent().variable_name()
        S = K[[t]]
        S.set_default_prec(prec)
        t = S.gens()[0]
        dx = x.derivative()
        ann_local = [vector(S,v.list())*vector([x^i*dx/(2*y) for i in range(g)]) for v in ann]
        if (Q[1])%p == 0:
            #weierstrass case
            I = [v.integral()(p*t) for v in ann_local] #no global constant of integration since weierstrass disk
        else:
            try:
                I = [ann[i]*(XK.coleman_integrals_on_basis(XK(0,1,0),Q)[0:g]) + (ann_local[i]).integral()(p*t) for i in range(len(ann))]
            except TypeError:
                I = [ann[i]*(XK.coleman_integrals_on_basis(XK(0,1,0),Q)[0:g]) + (ann_local[i]).power_series().integral()(p*t) for i in range(len(ann))]
        for i in I:
            coeffval = min(c.valuation() for c in i.list())
            i = i/p^coeffval
            i = i.truncate(prec)
            i = Qp(p,prec-5)['t'](i)
            rI = i.roots()
            print [rt for rt in rI if (rt[0]).valuation() > -1]
        print 20*"-"

##############################################

R.<x> = QQ[]
print 20*"="
X = HyperellipticCurve(x^11-4*41^2*5^10)  #genus 5 rank 0
print X
effective_chabauty(X, [], 3, 50)

