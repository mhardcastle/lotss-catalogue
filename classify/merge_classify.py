lines=[]
for dir in ['/data/lofar/mjh/hetdex_v4/toclassify', '/home/jcroston', '/data/lofar/wwilliams/hetdex/source_sorter/toclassify_171124']:

    filename=dir+'/large_faint_toclassify.txt'
    lines.append(open(filename,'r').readlines())

outfile=open('toclassify_merged.txt','w')
    
assert(len(lines[0])==len(lines[1]) and len(lines[1])==len(lines[2]))
    
for i in range(len(lines[0])):
    for j in range(3):
        bits=lines[j][i].split(',')
        if bits[0]:
            outfile.write(lines[j][i])
            break
    else:
        print i
        raise RuntimeError('We screwed up!')
