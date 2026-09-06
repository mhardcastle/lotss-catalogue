lines=open('green.txt').readlines()
region=open('green.reg','w')
region.write("""# Region file format: DS9 version 4.1
global color=yellow dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
""")
for l in lines[4:-1]:
    bits=l.split()
    size=bits[5]
    size=size.replace('?','')
    if 'x' in size:
        sbits=size.split('x')
        size=float(sbits[0])+float(sbits[1])*2
    else:
        size=float(size)
    size/=2.0
    region.write(f'circle({bits[0]}:{bits[1]}:{bits[2]},{bits[3]}:{bits[4]}:00,{size}\')\n')

region.close()
