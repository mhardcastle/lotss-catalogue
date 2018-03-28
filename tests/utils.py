def test_duplicate(t,colname,description,exclude=[]):
    print 'Constructing list...'
    count=0
    idlist=[]
    for i,r in enumerate(t):
        value=r[colname].rstrip()
        if value and value not in exclude:
            idlist.append((value,i))
    print 'Sorting list...'
    idlist.sort()
    print 'Looking for duplicates...'
    
    i=0
    while i<len(idlist)-1:
        j=0
        while i+j+1<len(idlist) and idlist[i][0]==idlist[i+j+1][0]:
            j+=1
        # now j counts how many duplicates there were
        if j>0:
            value=idlist[i][0]
            print 'Duplicate',description,'with name',value,':'
            for k in range(i,i+j+1):
                ix=idlist[k][1]
                print '          ',ix,t[ix]['Source_Name'],t[ix]['ID_flag']
            count+=1
        i+=j+1
    return count

def banner(s):
    print '------------------------------------------------------------------------------'
    print s
    print '------------------------------------------------------------------------------'
