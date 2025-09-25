import statistics
import matplotlib.pyplot as plt


def analyze_md_file(filename):
    lengths = []
    fin=open(filename,'r').read()
    lines=fin.splitlines()
    f=0
    totframes=0
    i=0
    prov_path=0
    while i<len(lines):
        line=lines[i]
        if "FRAME" in line:
            totframes+=1
        elif not line:
            if prov_path==0:
                i=i
            else:
                lengths.append(prov_path-1)#-1 to remove start node
                prov_path=0
                f+=1
        else:
            prov_path+=1
        i+=1


    perc=f/totframes*100
    avg_len=sum(lengths)/len(lengths)
    std_len=statistics.stdev(lengths)
    print(f"Total frames: {totframes}")
    print(f"Frames with path: {f} ({perc:.2f}%)")
    print(f"Average path length: {avg_len:.2f}")
    print(f"Standard deviation: {std_len:.2f}")

    # Distribution plot
    plt.hist(lengths, bins=range(0, max(lengths)+2), align="left", rwidth=0.8)
    plt.xlabel("Path length")
    plt.ylabel("Frequency")
    plt.title("Distribution of path lengths")
    plt.show()

    fout=open("lengths_"+filename,'w')
    for l in lengths:
        fout.write(str(l)+"\n")

# Usage
nin=input("input file name\n")
analyze_md_file(nin)

