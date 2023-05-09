
"""


def IsBegWithPos(L): #Returns True if L[0]>0. Otherwise return False.
    x=1/0
    for x in L:
        if x>0:
            return(True)
        if x<0:
            return(False)
    return(False)





def normAboveAgainstMatrix(AboveMatrix,R2):#AboveMatrix is a given IW partial matrix. R2 is a new row that we want to add to the matrix. We normalize R2[i]<0 if the column above it in AboveMatrix is zero.

    for i in range(len(R2)):
        if AboveMatrix.column(i).norm()==0:
            R2[i]=-abs(R2[i])
    return R2

def normAbove_old(myMAT):#MyMat is a given IW partial matrix. We normalize each column of MyMat to begin with a negative entry.
    for i in range(myMAT.ncols()):
        for j in range(myMAT.nrows()):
            if myMAT[j,i]!=0:
                myMAT[j,:]*=-sign(myMAT[j,i])
                break
    return myMAT

"""

def dpAgainstMatrix(AboveMatrix,L2): #We compute the norm of the product of L2 (a candidate ro given in a list form) against the rows of AboveMatrix.
    return((AboveMatrix*vector(L2)).norm())








def AcceptOrRejectRow(AboveMat,Row,Mon,base): #This function returns True if AboveMat+Row is minimal. This means acceptance. Otherwise reject.
    NewMat=AboveMat.stack(matrix(Row))
    if NewMat==FindMinInHadClass(NewMat,Mon,base):
        return True
    else:
        return False



def OrderMat(SPMat,base):#This procedure reorders the columns of SPMat (a partial IW) by increasing lex order. We interpret each column as represention of a number in base 'base'. This entry should be sufficiently large to accomodate all possible entries in the matrix. 

    nR=SPMat.nrows();
    baseVec = vector([base^(nR-j) for j in range(nR)])  
    SPMat=normAbove(SPMat)
    SPMat=sorted(SPMat.columns(),key=lambda vec:baseVec*vec)
    return matrix(SPMat).transpose()  


def normAbove(myMAT):#MyMat is a given IW partial matrix. We normalize each column of MyMat to begin with a negative entry.
    
    for i in range(myMAT.ncols()):
        for j in range(myMAT.nrows()):
            if myMAT[j,i]!=0:
                myMAT[:,i]*=-sign(myMAT[j,i])
                break
    return myMAT


def NormalizeByColumns(Mat,base): #This procedure normalizes a matrix by colums. That is, it returns the minimal matrix which is column equivalent to Mat.
    MatN=normAbove(Mat)
    MatN=OrderMat(Mat,base)
    return MatN

def FindMinInHadClass(A,Mon,base): #This returns the minimal matrix in the class of A, by row lex ordering.

    m=A.nrows(); n=A.ncols()
    CodeMin=n*[10^10]
    for M in Mon[m]:
        MA=M*A
        MAnormC=NormalizeByColumns(MA,base)
        MAnormCode=list(MAnormC*vector([base^(n-j) for j in range(n)]))
        #print(MAnormC)
        #print(MAnormCode)
        if MAnormCode<CodeMin:
            CodeMin=MAnormCode
            MinMat=MAnormC
    return MinMat


def AllMonMatUpTo(n): #This gives the full list of all monomial matrices up to order n.

    Mon={}
    for m in range(1,n+1):
        H=HadamardSpace(m)
        P=[pi.matrix() for pi in SymmetricGroup(m)]
        Lm=[-diagonal_matrix(h)*p for h in H for p in P]   
        Mon[m]=Lm
    return Mon


def goodPerm(L,s):#Given two lists L,s of the same length, check that there is no -1 in s corresponding to 0 in L.

    for i in range(1,len(L)):
        if (L[i]==0) and (s[i]==-1) : return False
    return True

def HadamardSpace(d): #Return a list of all Hadamard (+-1)-vectors of size d that begin with -1.
    SpaceSize = 2^(d-1)
    L=[ZZ(n).digits(2) for n in range(SpaceSize)]
    L1 = [x+[0]*(d-len(x)-1) for x in L]
    return([[-1]+[2*y-1 for y in x] for x in L1])

def AllVectorsWeight(n,w): #This returns a list of all vectors from sors, including permutation and signing, that start with a negative entry. This is the full list of rows from which the rows of the weighing are taken.
    Pos=sors(w,n,sqrt(w))
    AboveMat=matrix(n*[0])
    AllVectors=[]
    count=0
    for Typ in Pos:
        MTyp=OrthogonalSignedPermutationsWithFirstNonZeroPositiveAgainstMatrix(AboveMat,Typ,(0,0),floor(sqrt(w)))
        for x in MTyp:
            for j in range(len(x)):
                if x[j]>0:
                    break
                elif x[j]<0:
                    AllVectors.append(vector(x))
                    count+=1
                    break
    AllVectors.sort()
    return tuple(AllVectors)
        
    
def OrthogonalToAboveMatrix(AboveMatrix,FullRowList,HighestIndex,IndexSubset,Mon,base,MinimizeRowLimit=4):
    """
    This function computes a list of vectors that can be appended to the matrix AboveMat, plus a list of indices (see explanation below). The rows should be: 
    a) Orthogonal to AboveMat.
    b) The augmented matrix must be minimal in the col lex ordering.

    Input: AboveMat=The matrix we want to enlarge. Assumed to be minimal.
           FullRowList=The full list of all possible row vectors.
           HighestIndex=The index of the last row of AboveMatrix.
           IndexSubset=A subset of indices that is attached to AboveMat. The orthogonal vectors indices must belong to this list.
           Mon=A list of monomial matrices required for minimization.
           base=The code invariant base.

    The function loops over all row vectors whose index is taken from IndexSubset, not less than HighestIndex+1. For each choice it is first checked that it is orthogonal to the last row of AboveMat. If yes, it is then checked that the augmnted matrix is minimal. The function returns the list of all such vectors, together with the list of all indices of all orthogonal vectors (not just the minimal ones). This list will serve as the new IndexSubset in the next application of this function.
    """
    RowIndexList=[]
    OrthIndexList=[]
    for ii in IndexSubset:
        if ii<=HighestIndex:
            continue
        Row=FullRowList[ii]
        if AboveMatrix.nrows()==0 or Row*AboveMatrix[-1]==0:
            OrthIndexList.append(ii)
            if AboveMatrix.nrows()>=MinimizeRowLimit:
                RowIndexList.append(ii)
                continue
            NewMat=AboveMatrix.stack(matrix(Row))
            NewMatMin=FindMinInHadClass(NewMat,Mon,base)
            if NewMat==NewMatMin:
                RowIndexList.append(ii)
    return RowIndexList,OrthIndexList


def ExhaustiveListIW(n,w,MinimizeRowLimit=4):
    WeightVectors=AllVectorsWeight(n,w)
    Mon=AllMonMatUpTo(MinimizeRowLimit)
    base=1+2*floor(sqrt(w))
    MatList=IterateInitials([],[],[],WeightVectors,Mon,base,MinimizeRowLimit,n)
    return MatList



def IterateInitials(MatList,HighestIndexList,IndexSubsetList,WeightVectors,Mon,base,MinimizeRowLimit,n,m=0):
    LRows=len(WeightVectors)
    if m==0:
        AboveMatrix=matrix(0,n)
        HighestIndex=-1
        IndexSubset=list(range(len(WeightVectors)))
        Short,Long=OrthogonalToAboveMatrix(AboveMatrix,WeightVectors,HighestIndex,IndexSubset,Mon,base,MinimizeRowLimit)
        MatList=[matrix(WeightVectors[s]) for s in Short]
        HighestIndexList=Short
        IndexSubsetList=[[i for i in range(s,LRows)] for s in Short]
        m=1
    if m==n:
        return MatList
    else:
        print("m=",m)
        print(len(MatList))
        NewMatList=[]
        NewHighestIndexList=[]
        NewIndexSubsetList=[]
        for ii in range(len(MatList)):
            AboveMatrix=MatList[ii]
            HighestIndex=HighestIndexList[ii]
            IndexSubset=IndexSubsetList[ii]
            Short,Long=OrthogonalToAboveMatrix(AboveMatrix,WeightVectors,HighestIndex,IndexSubset,Mon,base,MinimizeRowLimit)
            for s in Short:
                row=WeightVectors[s]
                NewAboveMatrix=AboveMatrix.stack(matrix(row))
                NewMatList.append(NewAboveMatrix)
                high=s
                index_sub=[jj for jj in Long if jj>high]
                NewHighestIndexList.append(high)
                NewIndexSubsetList.append(index_sub)
        #print([NewMatList[randrange(len(NewMatList))] for k in range(4)])
        MatList=IterateInitials(NewMatList,NewHighestIndexList,NewIndexSubsetList,WeightVectors,Mon,base,MinimizeRowLimit,n,m+1)
        return MatList




        
def GreedyMinimize(MM,base):
    """
    This is a documentation.
    """
    n=MM.ncols()
    Sn=SymmetricGroup(n)
    pi=Sn.random_element()
    M=copy(MM)
    M.permute_rows(pi)
    M=NormalizeByColumns(M,base)
    for i in range(n):
        CodeMin=list(M*vector([base^(n-j) for j in range(n)]))
        Nplus=M
        for j in range(i+1,n):
            TestN=copy(M)
            TestN[i]=M[j]
            TestN[j]=M[i]
            T1=NormalizeByColumns(TestN,base)
            TCode=list(T1*vector([base^(n-j) for j in range(n)]))
            if TCode<CodeMin:
                #print('here1')
                CodeMin=TCode
                Nplus=T1
            TestN=copy(M)
            TestN[i]=-M[j]
            TestN[j]=M[i]
            T2=NormalizeByColumns(TestN,base)
            #print('T2=')
            #print(T2)
            TCode=list(T2*vector([base^(n-j) for j in range(n)]))
            if TCode<CodeMin:
                #print('here2')
                CodeMin=TCode
                Nplus=T2
            #print('Nplus=')
            #print(Nplus)
            #print(i,j)
        M=Nplus
        #print(M)
        #print()
    return M
            
            
            



def OrthogonalSignedPermutationsWithFirstNonZeroPositiveAgainstMatrix(AboveMat,L,Lindex,MaxV):
    """
    This precedure produces a list of all vectors of type L (i.e. Hadamard equivalent to L) that are orthogonal to AboveMat.
    Input:
          AboveMat = a partial IW matrix.
          L = The given type of the orthogonal row.
          Lindex = The index of L in the list of 'sors'.
          maxV = An upper bound of all modulii of all matrix entries.

    Output: A list of pairs: ( [index in sors, index in permutations, index in signings], new orthogonal row)
    """
    P = Permutations(L).list()
    d = len(L)
    S = HadamardSpace(d)
    ln=AboveMat.nrows()+1
    SL=[diagonal_matrix(v)*p.matrix() for v in HadamardSpace(ln) for p in SymmetricGroup(ln)]
    SL = HadamardSpace(AboveMat.nrows()+1)
    SP = []
    s = []
    sd = [] 
    sd1 = []
    sd2 = []
    codeBase = 1+2*MaxV
    nR=AboveMat.nrows()
    cB=(vector([codeBase^(nR-j) for j in range(nR)]))*AboveMat
    for pi in range(Lindex[1],len(P)):
        p = P[pi]
        #s0=Lindex[2] if pi==Lindex[1] else 0
        s0=0
        for si in range(s0,len(S)):
            sn=S[si]
            if (goodPerm(p,sn)): #Take only the signing at nonzero entries of the vector. This saves a lot of time.
                #print(p,sn)
                r =[p[i]*sn[i] for i in range(d)]
                sd1.append(r)
                
                
                if true:
                    if dpAgainstMatrix(AboveMat,r)==0: # Is the vector orthogonal?
                        #r=normAboveAgainstMatrix(AboveMat,r)
                        NewMat=matrix(nR+1,d,AboveMat.list()+r)
                        signedPermutaionsM = [diagonal_matrix(v)*p.matrix() for v in HadamardSpace(nR+1) for p in SymmetricGroup(nR+1)]
                        #for sPM in signedPermutaionsM:
                        #    OrederedSnewMat=orderMat(sPM*NewMat) #This is for the new algorithm (TODO: continue)
                        if not(r in SP):
                            SP.append(r)
                            R2i=vector(r)
                            M1 = cB+R2i #This computes the new code invariant of the (potential) augmented matrix.
                            M1s = sorted(M1)
                            if not(M1s in s):
                                M1m = cB-R2i #This is the same with -the row.
                                M1sm = sorted(M1m)
                                if not(M1sm in s): #Only consider row if it contributes new code invariant.
                                    s.append(M1s)
                                    sd.append([[Lindex[0],pi,si],R2i])   
                                     
    return(sd1) #Originally return(sd)
"""
def LEcodeVector(cv1,cv2): #This asks if two vectors cv1<=cv2 under the lex order (entries are ordered by integer ordering).
    for i in range(len(cv1)):
        if cv1[i]<cv2[i]: return true
        if cv1[i]>cv2[i]: return false
    return true

"""
def sors(n,r,maxsq):  ### Find all ways to express n=\sum_i s_i^2 as a sum of r squares, with maxsq>=s_1>=s_2>=...>=s_r.
    #n0=RDF(n). This is NOT our best algorithm (TODO: look for the more recent 'nsoks').
    n0=n
    if r==1:
        s=sqrt(n0)
        fs=floor(s)
        if s==fs and s<=maxsq:
            return [[s]]
        else:
            return false
    else:
        sors2=[]
        for x in range(floor(sqrt(n0/r)),min(floor(sqrt(n0)),maxsq)+1):
            m=n-x^2
            #print x
            sors1=sors(m,r-1,x)
            if sors1:
                sors2+=[[x]+sors1[i] for i in range(len(sors1))   ]
           # print "sors2=",sors2
        if len(sors2)>0:
            return sors2
        else:
            return false

"""
def OurMain3(WMWeight, WMOrderDiv,lines): 
    
    Return a list of (partial) IW matrices of a given size.

    Input: WMWeight = the desired weight.
           WMOrderDiv = the row size.
           lines = the column size.

    Output: A list of IW matrices with the given parameters. The list is an exhaustive list of all possible solution up to Hadamard equivalence.
            It is possible though that the same Hadamard type will repeat more than once. 
    
    MaxSquare = WMWeight^(1/2)
    A = sors(WMWeight, WMOrderDiv, MaxSquare) #This gives the list of all possible rows up to sign and ordering.
    AllSoks=[list(map(lambda x:-x , A[len(A)-i-1])) for i in range(len(A)) ] #This is just to represent each row by the minimal lex representative.
    print("Nsoks: ",AllSoks)
    Sols= [[matrix(1,WMOrderDiv,AllSoks[i]),[i,0,0]] for i in range(len(AllSoks))] #Sols is the list of pairs: A partial matrix and the index of its last row in Allsoks. We start with just 1 row per matrix, and update as we go through the main loop.
    print(Sols)
    for line in range(lines-1): #This is the main loop, which adds one row at a time to each matrix. 
        nSols=[]
        countSol=0
        countScan=0
        for Sol in Sols:
            countScan+=1
            AboveMat=Sol[0]
            R1Index=Sol[1]
            for R2i in range(R1Index[0],len(AllSoks)): # We want to add a new row, but one that appears later (or same) in Allsoks. This is because we can allow our final matrix to be monotonely increasing in Allsoks.
                R2Index=R1Index if R2i==R1Index[0] else [R2i,0,0]
                R2 = vector(AllSoks[R2i])
                cpF = OrthogonalSignedPermutationsWithFirstNonZeroPositiveAgainstMatrix(AboveMat,R2,R2Index,MaxSquare) #This gives a list of all rows which are Hadamard equivalent to R2 and orthogonal to AboveMat.
                countSol+=len(cpF)
                if len(cpF)>0:
                    nSols.extend([AboveMat,x] for x in cpF) #This is a list of pairs of AboveMat and its new row x.
            print("countScan,countSol  = ",countScan,countSol)
        Sols=[] #We now update Sols.
        for nSol in nSols:
            NL=matrix(WMOrderDiv,1,list(nSol[1][1])) #This practically adds the row to the matrix.
            AM=((nSol[0].transpose()).augment(NL)).transpose()
            Sols.append([AM,nSol[1][0]])
        print("# of sols in line",line, "  =",len(Sols),[Sols[i] for i in range(min(len(Sols),3))])
    return(Sols)

from gc import collect
collect()
#%time res=OurMain3(8,8,8)



def AcceptOrRejectNewRow_Old(AboveMat,Row,base):
    
    This procedure tests whether a new row, added to a minimal matrix keeps the matrix minimal. This is the main rejection criterion used in the main algorithm. If the result is not minimal then necessarily then AboveMat is assumed to have been treated in the past and can reject this row. The function returns True if Row is accepted and False if Row is rejected.

    The algorithm inserts the vector 'row' at some position (we loop on that position).
    Accordingly we divide the matrix into three parts: Top,Row and Bottom.

    We insert Row and -Row and then colum-minimize Top+Row.

    1)If Row became smaller than the original vector of AboveMatrix, reject and finish.
    2)If Row became larger or equal, then move Row one place down (continue to the next loop).
    3)If Row reached all the way down, then accept.
    
    
    Remarks: 
        a) In the case that Row equals (in case 2), it may still be possible to reject Row, if the minimization of the part below Row ends with a smaller matrix. We feel that this case occurs seldom and we can skip this test.
        b) As a consequence, we might generate extra solutions with initials not being minimal. We feel that this is not too much of an excess, but this should be tested.
    
    m=AboveMat.nrows(); n=AboveMat.ncols()
    CodeVec=vector([base^(n-j) for j in range(n)])
    for ii in range(m-1):
        TopPlusRow= AboveMat[:ii].stack(matrix(Row))
        NormTopPlusRow=NormalizeByColumns(TopPlusRow,base)
        Code=CodeVec*NormTopPlusRow[-1]
        OldCode=CodeVec*AboveMat[ii]
        if Code<OldCode:
            print(ii,NormTopPlusRow[-1],AboveMat[ii])
            return False  #Reject in permuted Row is smaller than original vector.
    return True

"""