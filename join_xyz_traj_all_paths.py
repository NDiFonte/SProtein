nin=input("Input traj filename\n")
fin=open(nin,'r').read()
lines=fin.splitlines()
i=0
fout=open("joint_"+nin,'w')
while i<len(lines):
    line=lines[i]
    if not line and i<len(lines)-1  and "FRAME" not in lines[i+1]:
        i+=1
        continue
    fout.write(line+"\n")
    i+=1
