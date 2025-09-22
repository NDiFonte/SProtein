inp=open("traj_sp.xyz",'r').read()
out=open("VMD_path_traj.xyz",'w')
lines=inp.splitlines()
#Number of frames
N=11
gen_frame_zero=True

i=0
j=0
k=0
while i<len(lines):
    if "FRAME" in lines[i]:
        i+=1
    if len(lines[i])>0:
        j+=1
    else:
        k=max(j,k)
        j=0
    i+=1
i=0

if gen_frame_zero:
    out.write(str(k)+"\n")
    out.write("FRAME 0"+"\n")
    j=0
    while j<k:
        out.write("  O  0.000     0.000     0.000"+"\n")
        j+=1

while i<len(lines):
    if "FRAME" in lines[i]:
        out.write(str(k)+"\n")
        out.write(lines[i]+"\n")
        i+=1
    elif len(lines[i])>0:
        j=0
        while j<k:
            if i==len(lines):
                if j<k:
                    while j<k:
                        out.write("  Se  0.000     0.000     0.000"+"\n")
                        j+=1
                    break
            if len(lines[i])>0 and not "FRAME" in lines[i]:
                out.write("  Se  "+lines[i]+"\n")
                j+=1
                i+=1
            else:
                if j<k:
                    out.write("  Se  0.000     0.000     0.000"+"\n")
                    j+=1
            if len(lines[i])==0:
                i+=1
    else:
        j=0
        while j<k:
            out.write("  O  0.000     0.000     0.000"+"\n")
            j+=1
        i+=1
