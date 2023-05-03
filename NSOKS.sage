def NSOKS(n,r,maxsq=-1): ### This is the same as sors, but faster. We only do recursion on the multiplicities of each square.
    if n==0:
        return [[0,r]]
    n0=RDF(n)
    M=floor(sqrt(n0))
    if maxsq!=-1: M=min(M,maxsq)
    L=ceil(sqrt(n0/r))
    Lsqr=[]
    for s in range(L,M+1):
        s2=s^2
        ns=floor(n0/s2)
        for i in range(1,min(r+1,ns+1)):
            nnew=n-i*s2
            if i==r:
                Lsqr.append([[s,r]])
            else:
                rem=sorsm(nnew,r-i,maxsq=s-1)
                if rem==[]:
                    continue
                for lsqr in rem:
                    Lsqr.append([[s,i]]+[lsqr])
    return Lsqr