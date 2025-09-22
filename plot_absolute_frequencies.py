from collections import Counter
import matplotlib.pyplot as plt

n =int(input("How many residues would you like to display?\n"))  # Change this to display more or fewer amino acids

#OPENS NEIGHBORS FILE AND READS IT
with open('neighbors_traj_path_nodes.dat', 'r') as f:
    lines = f.readlines()

# Count amino acids in columns 4 and beyond
amino_acid_counter = Counter()

for line in lines:
    if line.strip().startswith("FRAME") or not line.strip():
        continue
    parts = line.strip().split()
    amino_acids = parts[3:]  # Columns 4 and onward, depending on the number of neighbors
    amino_acid_counter.update(amino_acids)

# Get top n most common
top_amino_acids = amino_acid_counter.most_common(n)
labels, counts = zip(*top_amino_acids)

# Plotting
plt.figure(figsize=(10, 6))
plt.bar(labels, counts, color='blue')
plt.xticks(rotation=45,fontsize=16)
plt.ylabel("Frequency",fontsize=18)
plt.title(f"Top {n} Most Frequent Amino Acids",fontsize=24)
plt.tight_layout()
plt.show()

