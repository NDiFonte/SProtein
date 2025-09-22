from md_functions_water_network import *
from collections import defaultdict
from scipy.spatial import KDTree
np.set_printoptions(threshold=np.inf)
import plotly.graph_objects as go
#mdanalysis imported as mda

def extract_data_sp(U, kk, list_molecules, list_nodes):
    '''
    Extracts the coordinates of the atoms and the values of the box for the frame kk
    '''
    
    for tr in U.trajectory[kk:kk+1]:
    
        box = np.array(tr.dimensions)
        box[0] = box[0]/10                              # we convert in nm
        box[1] = box[1]/10
        box[2] = box[2]/10
        
        box_matrix = mda.lib.mdamath.triclinic_vectors(box)  # a triclinic matrix representation
        diag_box = np.diag(box_matrix) 
                      
        coord = np.array([[None, None, None]])
        for ii in range(len(list_molecules)):
            ml = U.select_atoms('resname '+str(list_molecules[ii]))
            for jj in range(len(list_nodes[ii])):
                at = ml.select_atoms('name '+str(list_nodes[ii][jj]))
                pos = at.positions
                coord = np.vstack((coord, pos))
            
    coord = np.delete(coord, 0, 0)
    coord = coord * (1/10)                            # we convert in nm
    coord = np.asarray(coord, dtype = "float")
    
    return {'coord': coord, 'box': box, 'diag_box': diag_box}

def extract_nodes(U,reference_atom,ox_atom,selections,cutoff_distance=25.0):
    gprotein=U.select_atoms("protein")
    protein_selection_string=""
    for selection in selections:
        protein_selection_string=protein_selection_string+f"({selection}) or "
    protein_selection_string=protein_selection_string[:-3]
    ref_atom=reference_atom
    if selections !=[]:
        graph_nodes = U.select_atoms(f"({ref_atom}) or ({protein_selection_string}) or (name {ox_atom} and around {cutoff_distance} group gprotein)",gprotein=gprotein)#group with reference atom, key atoms in the protein and water oxygens within cutoff from the protein
    else:
        graph_nodes = U.select_atoms(f"({ref_atom}) or (name {ox_atom} and around {cutoff_distance} group gprotein)",gprotein=gprotein)
    coords=graph_nodes.positions
    return graph_nodes,coords

def extract_H_coordinates(U,all_hydrogens,h1_atom,h2_atom,cutoff_distance):
    gprotein=U.select_atoms("protein")
    protein_selection_string=""
    for selection in all_hydrogens:
        protein_selection_string=protein_selection_string+f"({selection}) or "
    protein_selection_string=protein_selection_string[:-3]
    if all_hydrogens !=[]:
        graph_nodes = U.select_atoms(f"({protein_selection_string}) or ((name {h1_atom} {h2_atom}) and around {cutoff_distance} group gprotein)",gprotein=gprotein)#group with reference atom, key atoms in the protein and water oxygens within cutoff from the protein
    else:
        graph_nodes = U.select_atoms(f"((name {h1_atom} {h2_atom}) and around {cutoff_distance} group gprotein)",gprotein=gprotein)
    coords=graph_nodes.positions
    return graph_nodes,coords


def display_atoms_for_residue(residue):
    # Mostra tutti gli atomi per un residuo di esempio
    print(f"\nAtoms of {residue.resname}:")
    for i, atom in enumerate(residue.atoms):
        print(f"{i + 1}: {atom.name}")

def display_residues(residue_list):
    print("Available residues in the topology file")
    for i, residue in enumerate(residue_list):
        print(f"{i + 1}: {residue}")

def select_residues_and_atoms(U):
    fresidues =set([res.resname for res in U.residues])
    fresidues=list(fresidues)
    fresidues.sort()
    selections = [] 
    display_residues(fresidues)
    str_residues =input("Select residues to consider (separated by a space) ").split()
    residues=[]
    for resid in str_residues:
        resid=int(resid)-1
        residues.append(fresidues[resid])
    for residue_name in residues:
        residue = U.select_atoms(f'resname {residue_name}').residues[0]
        display_atoms_for_residue(residue)
        str_atom_names = input(f"Select atoms for {residue_name}:  ").split()
        atom_names=[]
        fatom_names = [atom.name for residue in U.residues if residue.resname == residue_name for atom in residue.atoms]
        fatom_names=list(fatom_names)
        for atom in str_atom_names:
            atom=int(atom)-1
            atom_names.append(fatom_names[atom])
        if atom_names:
            selection_string = f"resname {residue_name} and name " + " ".join(atom_names)
            selections.append(selection_string)
    return selections


def select_solvent(U):
    fresidues =set([res.resname for res in U.residues])
    fresidues=list(fresidues)
    fresidues.sort()
    selections = [] 
    display_residues(fresidues)
    str_residue =input("Select water residue  ")
    wat_residue=int(str_residue)-1
    wat_residue=fresidues[wat_residue]
    residue= U.select_atoms(f'resname {wat_residue}').residues[0]
    H1_atom=""
    for atom in residue.atoms:
        atom=str(atom)        
        if " O" in atom:
            elem=atom.split(": ")
            ns=elem[1].split(" of")
            ox_atom=ns[0]
        if "H" in atom and " O" not in atom:
            if H1_atom=="":
                elem=atom.split(": ")
                ns=elem[1].split(" of")
                H1_atom=ns[0]
            else:
                elem=atom.split(": ")
                ns=elem[1].split(" of")
                H2_atom=ns[0]
    return ox_atom,H1_atom,H2_atom,wat_residue

def select_ref_atom(U):
    fresidues =set([res.resname for res in U.residues])
    fresidues=list(fresidues)
    fresidues.sort()

    selections = [] 
    display_residues(fresidues)
    str_residue =input("Select reference atom residue  ")
    ref_residue=int(str_residue)-1
    ref_residue=fresidues[ref_residue]
    residue= U.select_atoms(f'resname {ref_residue}').residues[0]
    display_atoms_for_residue(residue)
    ref_atom_num = input(f"Select atoms for {ref_residue} ")
    i=1
    for atom in residue.atoms:
        if i==int(ref_atom_num):
            atom=str(atom)
            elem=atom.split(": ")
            ns=elem[1].split(" of")
            ref_atom=ns[0]
        i=i+1
    residue_dict = defaultdict(list)
    for res in U.residues:
        residue_dict[res.resname].append(res.resid)
    print(f"Select {ref_residue}: {', '.join(map(str, residue_dict[ref_residue]))}")
    res_number=int(input())
    ref_atom_string=f"resname {ref_residue} and resid {res_number} and name {ref_atom}"
    return ref_atom_string

def find_hydrogens_from_selection(U,selection,hdist):
    sel_atoms = U.select_atoms(selection)
    #print(sel_atoms)
    hydrogen_selections = []
    for atom in sel_atoms:
        # Find hydrogen atoms within a cutoff (e.g., 1.3 Å) around the atom
        nearby_hydrogens = U.select_atoms(f"(around {hdist} (index {atom.index})) and name H*")
        for hydrogen in nearby_hydrogens:
            if hydrogen.resid == atom.resid and hydrogen.resname == atom.resname:
                hydrogen_sel = f"resname {hydrogen.resname} and name {hydrogen.name}"
                hydrogen_selections.append(hydrogen_sel)
    return list(set(hydrogen_selections))

def filter_sequences(sequences):
    # Sort sequences by length in descending order
    sequences = sorted(sequences, key=len, reverse=True)
    filtered = []
    
    for seq in sequences:
        # Check if seq is a suffix of any sequence already in filtered
        if not any(seq == other[-len(seq):] for other in filtered):
            filtered.append(seq)
    
    # Restore the original order (optional)
    filtered = sorted(filtered, key=sequences.index)
    return filtered


###########################################################################################################END OF SELECTION FUNCTIONS    

def get_graph_directed_sp(U,coord_nodes,nodes,hydrogens,coord_H,box, boundary,dist):

    A = new_oriented_adjacency_matrix_sp(dist, coord_nodes,nodes,hydrogens, coord_H, box, boundary)

    
    #G = nx.from_scipy_sparse_matrix(A, create_using=nx.DiGraph)
    
    #n_edges = G.number_of_edges()
    return A

def get_graph_undirected_sp(U,coord,box, boundary, dist):
    
    if int(boundary) == 1:
        A = get_adjacency_matrix_sp(dist, coord)
    elif int(boundary) == 2:
        A = get_adjacency_matrix_pbc_sp(dist, coord, box)
    else:
        print('Value not valid')
        sys.exit(1)
    
    G = nx.from_scipy_sparse_matrix(A)
    
    n_edges = G.number_of_edges()
        
    return G,n_edges,A

def get_adjacency_matrix_sp(dist, coord): #NON VA MODIFICATA :)

    '''
    Returns the adjacency matrix without periodic boundary conditions  
    '''

    coord = np.asmatrix(coord)    
    D = distance_matrix(coord, coord)
    
    D1 = D <= dist #boolean matrix whose entry {i,j} is True if the distance between node v_i and node v_j is less or equal than dist
    
    A = D1*1
    A = A - scipy.sparse.diags(A.diagonal())
    
    return csr_matrix(A)
    
    
def get_adjacency_matrix_pbc_sp(dist, coord, box):  #NON VA MODIFICATA :)

    '''
    Returns the adjacency matrix with periodic boundary conditions
    '''
    
    coord = np.squeeze(np.array(coord))
    
    N = np.shape(coord)[0]
    
    M = np.zeros((N,N))
    
    triu = np.triu_indices_from(M, k=1)
    
    self_distance_pbc = mda.lib.distances.self_distance_array(coord, box = box)
    
    M[triu] = self_distance_pbc
    M.T[triu] = self_distance_pbc
    
    D1 = M <= dist #boolean matrix whose entry {i,j} is True if the distance between node v_i and node v_j is less or equal than dist
    
    A = csr_matrix(D1*1)
    A = A - scipy.sparse.diags(A.diagonal())
    return csr_matrix(A)



def find_shortest_path(adjacency, src, dst, num_nodes):
    """
    Finds the shortest path in a directed graph represented as a CSR matrix.

    Parameters:
        adjacency (scipy.sparse.csr_matrix): The adjacency matrix in CSR format.
        src (int): The source node index.
        dst (int): The destination node index.
        num_nodes (int): Total number of nodes in the graph.

    Returns:
        path (list): The shortest path from src to dst as a list of node indices.
        sp_length (int): The length of the shortest path, or -1 if no path exists.
    """
    # Initialize BFS variables
    queue = deque()
    visited = [-1] * num_nodes  # -1 means unvisited
    sp_length = 0

    # Start BFS
    queue.append(src)
    visited[src] = src  # Mark the source as visited with itself as the parent

    while queue:
        current = queue.popleft()

        # If destination is reached, exit the loop
        if current == dst:
            break

        # Traverse neighbors using CSR row slicing
        neighbors = adjacency.indices[adjacency.indptr[current] : adjacency.indptr[current + 1]]
        for neighbor in neighbors:
            if visited[neighbor] == -1:  # Unvisited
                queue.append(neighbor)
                visited[neighbor] = current  # Record the parent node

    # Reconstruct the shortest path
    path = []
    i = dst
    while i != src:
        if visited[i] == -1:
            return [], -1  # No path found
        path.append(i)
        sp_length += 1
        i = visited[i]

    path.append(src)
    path.reverse()

    return path, sp_length


def check_if_protein(nodes,index):
    protein_resnames = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "HSD", "HSE", "HSP", "MSE", "SEP", "TPO", "PTR",
        "ACE", "NH2", "CYX", "HID", "HIE", "HIP", "ASH", "GLH", "LYN"
    }
    is_protein=False
    resid=str(nodes[index].residue)
    l=resid.split(",")
    ll=l[0].split()
    resid=ll[1]
    if resid in protein_resnames:
        is_protein=True
    return is_protein

def is_same_residue(node,hydrogen):
    return node.residue == hydrogen.residue

def set_row_to_zero(A, i):
    A = A.tolil()  # Convert to LIL format to efficiently modify rows
    A[i, :] = 0  # Set the row to zero
    return A.tocsr()  # Convert back to CSR format


def new_oriented_adjacency_matrix_sp(dist, coord_nodes,nodes,hydrogens,coords_H, box, boundary):#fa la matrice di adiacenza solo sulla dist e poi mette 0 se non è buono l'angolo
    if int(boundary) == 1:
        A_dist = get_adjacency_matrix_sp(dist, coord_nodes)
    elif int(boundary) == 2:
        A_dist = get_adjacency_matrix_pbc_sp(dist, coord_nodes, box)

    n = np.shape(A_dist)[0]
    dim = 3
    kdtree = KDTree(coords_H)
    dist_dummy, indices = kdtree.query(coord_nodes, k=3)
    lower = ((30)*np.pi)/180
    
 
    for i in range(n):
        l=-1
        hpartners=[]
        nodes_1 = np.squeeze(np.array(coord_nodes[i,:]))
        for idx in indices[i]:
            if is_same_residue(nodes[i],hydrogens[idx]):
                l=idx
                hpartners.append(idx)
        if (l==-1) and (not check_if_protein(nodes,i)):
            A_dist=set_row_to_zero(A_dist,i)
            continue
        hydr_coords=[]
        for l in hpartners:
            if int(boundary) == 1:
                hyd = np.squeeze(np.array(coords_H[l,:]))
                
            elif int(boundary) == 2:
                hyd = np.squeeze(np.array(coords_H[l,:]))
                for d in range(dim):
                    dist = abs(nodes_1[d] - hyd[d])
                    if dist > box[d]/2:
                        if nodes_1[d] >= hyd[d]:
                            hyd[d] = box[d] + hyd[d]
                        else:
                            hyd[d] = hyd[d] - box[d]
            hydr_coords.append(hyd)
        vect = A_dist[i,:]
        _,J,_ = find(vect) #LISTA DEI VICINI
        for j in J:
            if check_if_protein(nodes,i) and check_if_protein(nodes,j) and is_same_residue(nodes[i],nodes[j]):
                continue
            else:
                # I select the coordinates of O2 among all the periodic coordinates. I chose the one with the shortest distance from O.
                nodes_2 = np.squeeze(np.array(coord_nodes[j,:]))
                for d in range(dim):
                    dist = abs(nodes_1[d] - nodes_2[d])
                    if dist > box[d]/2:
                        if nodes_1[d] >= nodes_2[d]:
                            nodes_2[d] = box[d] + nodes_2[d]
                        else:
                            nodes_2[d] = nodes_2[d] - box[d]
                zzz=0
                for hyd in hydr_coords:     
                    deg=angle_between(hyd, nodes_1, nodes_2)
                    if deg<=lower:
                        break
                    zzz=zzz+1
                if zzz>=len(hydr_coords):
                    A_dist[i,j] = 0
    #A_dist.eliminate_zeros()
    return A_dist

def filter_redundant_paths(paths):
    """
    Removes redundant paths from a list of paths. A path is considered redundant if it contains 
    the same destination node as a shorter path already in the list.

    Args:
        paths (list of list): A list of paths, where each path is a list of nodes.

    Returns:
        list of list: A list of non-redundant paths.
    """
    # Sort paths by length (shorter paths first)
    paths = sorted(paths, key=len)
    
    non_redundant = []
    destination_set = set()  # Keeps track of arrival points already seen
    
    for path in paths:
        # Check the destination node (last node in the path)
        if path==[]:
            continue
        destination = path[-1]
        if destination not in destination_set:
            non_redundant.append(path)
            destination_set.add(destination)  # Mark this destination as covered
    
    return non_redundant,destination_set


def truncate_paths(paths, destination_set):
    """
    Truncates each path in a list of paths at the first occurrence of any node in the destination set.

    Args:
        paths (list of list): A list of paths, where each path is a list of nodes.
        destination_set (set): A set of destination nodes to check against.

    Returns:
        list of list: A list of truncated paths.
    """
    truncated_paths =set()
    
    for path in paths:
        truncated_path = []
        for node in path:
            truncated_path.append(node)
            if node in destination_set:
                break
        truncated_paths.add(tuple(truncated_path))
    
    return [list(path) for path in truncated_paths]


def filter_shortest_path(paths):
    """
    Selects the shortest path (by length) from a list of paths. 
    If there are multiple shortest paths, one is selected randomly.

    Args:
        paths (list of list): A list of paths.

    Returns:
        list: The shortest path.
    """
    if not paths:
        return []
    
    # Find the minimum length of all paths
    min_length = min(len(path) for path in paths)
    
    # Collect all paths with the minimum length
    shortest_candidates = [path for path in paths if len(path) == min_length]
    
    # Randomly select one from the shortest candidates
    return random.choice(shortest_candidates)


################################################################################################## END OF GRAPH THEORY FUNCTIONS

def create_scatter_plot_from_coords(coords):
    x_coords = [p[0] for p in coords]
    y_coords = [p[1] for p in coords]
    z_coords = [p[2] for p in coords]

    # Create a scatter plot
    fig = go.Figure(
        data=[
            go.Scatter3d(
                x=x_coords,
                y=y_coords,
                z=z_coords,
                mode='markers',
                marker=dict(
                    size=1,
                    color='blue',
                    opacity=1
                )
            )
        ]
    )

    # Update layout for better visualization
    fig.update_layout(
        scene=dict(
            xaxis_title='X-axis',
            yaxis_title='Y-axis',
            zaxis_title='Z-axis'
        ),
        title='3D Scatter Plot of Nodes'
    )
    fig.show()
    return


def create_scatter_plot_from_mcoords(*datasets):
    """
    Create a 3D scatter plot with multiple datasets.

    Parameters:
    - datasets: A variable number of datasets, where each dataset is a list of coordinates [(x, y, z), ...]

    Each dataset will be plotted with a unique color.
    """
    colors = ['blue', 'red', 'green', 'purple', 'orange']  # Define a list of colors
    fig = go.Figure()

    for i, coords in enumerate(datasets):
        x_coords = [p[0] for p in coords]
        y_coords = [p[1] for p in coords]
        z_coords = [p[2] for p in coords]

        # Add a trace for each dataset
        fig.add_trace(
            go.Scatter3d(
                x=x_coords,
                y=y_coords,
                z=z_coords,
                mode='markers',
                marker=dict(
                    size=3,
                    color=colors[i % len(colors)],  # Cycle through colors if more datasets
                    opacity=0.8
                ),
                name=f'Dataset {i + 1}'  # Name for legend
            )
        )

    # Update layout for better visualization
    fig.update_layout(
        scene=dict(
            xaxis_title='X-axis',
            yaxis_title='Y-axis',
            zaxis_title='Z-axis'
        ),
        title='3D Scatter Plot of Multiple Datasets',
        legend=dict(
            title="Datasets"
        )
    )

    fig.show()


def get_sorted_edges(sparse_adj_matrix):
    # Extract the data, rows, and columns
    rows, cols = sparse_adj_matrix.nonzero()
    data = sparse_adj_matrix.data
    
    # Combine rows, cols, and data into a single list
    edges = list(zip(rows, cols, data))
    
    # Sort the edges by (row, col)
    sorted_edges = sorted(edges, key=lambda x: (x[0], x[1]))
    
    # Unpack the sorted edges
    sorted_rows = [edge[0] for edge in sorted_edges]
    sorted_cols = [edge[1] for edge in sorted_edges]
    sorted_data = [edge[2] for edge in sorted_edges]
    
    # Create a new CSR matrix from the sorted data
    sorted_matrix = csr_matrix(
        (sorted_data, (sorted_rows, sorted_cols)),
        shape=sparse_adj_matrix.shape
    )
    
    return sorted_matrix



def plot_graph_directed(coordinates, sparse_adj_matrix, additional_points=None, additional_colors=None):


    # Ensure coordinates are a NumPy array
    coordinates = np.array(coordinates)

    # Extract edges from the sparse adjacency matrix
    upper_edge_x = []
    upper_edge_y = []
    upper_edge_z = []

    lower_edge_x = []
    lower_edge_y = []
    lower_edge_z = []

    # Iterate over non-zero elements in the sparse matrix
    for row in range(sparse_adj_matrix.shape[0]):
        neighbors = sparse_adj_matrix.indices[
            sparse_adj_matrix.indptr[row]:sparse_adj_matrix.indptr[row + 1]
        ]
        for col in neighbors:
            x0, y0, z0 = coordinates[row]
            x1, y1, z1 = coordinates[col]
            if row < col:  # Upper triangle
                upper_edge_x.extend([x0, x1, None])  # None to break the line
                upper_edge_y.extend([y0, y1, None])
                upper_edge_z.extend([z0, z1, None])
            elif row > col:  # Lower triangle
                lower_edge_x.extend([x0, x1, None])
                lower_edge_y.extend([y0, y1, None])
                lower_edge_z.extend([z0, z1, None])

    # Create edge traces for upper and lower triangle edges
    upper_edge_trace = go.Scatter3d(
        x=upper_edge_x, y=upper_edge_y, z=upper_edge_z,
        mode='lines',
        line=dict(width=2, color='red'),
        hoverinfo='none',
        opacity=1
    )

    lower_edge_trace = go.Scatter3d(
        x=lower_edge_x, y=lower_edge_y, z=lower_edge_z,
        mode='lines',
        line=dict(width=5, color='blue'),
        hoverinfo='none',
        opacity=0.5
    )

    # Create node trace
    node_x = coordinates[:, 0]
    node_y = coordinates[:, 1]
    node_z = coordinates[:, 2]

    node_trace = go.Scatter3d(
        x=node_x, y=node_y, z=node_z,
        mode='markers',
        marker=dict(size=2, color='black'),
        text=[f'Node {i}, coord {coord}' for i,coord in zip(range(len(coordinates)),coordinates)],
        hoverinfo='text'
    )

    # Prepare additional points trace (if provided)
    additional_trace = None
    if additional_points is not None:
        additional_points = np.array(additional_points)
        additional_colors = np.array(additional_colors) if additional_colors is not None else ['red'] * len(additional_points)
        
        additional_x = additional_points[:, 0]
        additional_y = additional_points[:, 1]
        additional_z = additional_points[:, 2]

        additional_trace = go.Scatter3d(
            x=additional_x, y=additional_y, z=additional_z,
            mode='markers',
            marker=dict(size=7, color=additional_colors),
            text=[f'Node {i}, coord {coord}' for i,coord in zip(range(len(coordinates)),coordinates)],
            hoverinfo='text'
        )

    # Create the figure
    data = [upper_edge_trace, lower_edge_trace, node_trace]
    
    if additional_trace is not None:
        data.append(additional_trace)
    
    fig = go.Figure(data=data)
    fig.update_layout(
        showlegend=False,
        scene=dict(
            xaxis=dict(showbackground=False),
            yaxis=dict(showbackground=False),
            zaxis=dict(showbackground=False),
        ),
        title="3D Graph Representation",
    )

    # Show the figure
    fig.show()


def plot_graph_undirected(coordinates, sparse_adj_matrix, additional_points=None, additional_colors=None):
    """
    Plot a 3D graph using Plotly from node coordinates and a sparse adjacency matrix in CSR format, 
    with the option to plot additional points with custom colors.

    Parameters:
        coordinates (array-like): Coordinates of the nodes, shape (N, 3).
        sparse_adj_matrix (scipy.sparse.csr_matrix): Sparse adjacency matrix in CSR format, shape (N, N).
        additional_points (array-like, optional): Coordinates of additional points, shape (M, 3). Default is None.
        additional_colors (array-like, optional): Colors for the additional points, shape (M,). Default is None.
    """
    # Ensure coordinates are a NumPy array
    coordinates = np.array(coordinates)

    # Extract edges from the sparse adjacency matrix
    edge_x = []
    edge_y = []
    edge_z = []

    # Iterate over non-zero elements in the sparse matrix
    for row in range(sparse_adj_matrix.shape[0]):
        neighbors = sparse_adj_matrix.indices[
            sparse_adj_matrix.indptr[row]:sparse_adj_matrix.indptr[row + 1]
        ]
        for col in neighbors:
            x0, y0, z0 = coordinates[row]
            x1, y1, z1 = coordinates[col]
            edge_x.extend([x0, x1, None])  # None to break the line
            edge_y.extend([y0, y1, None])
            edge_z.extend([z0, z1, None])

    # Create edge trace
    edge_trace = go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z,
        mode='lines',
        line=dict(width=2, color='black'),
        hoverinfo='none'
    )

    # Create node trace
    node_x = coordinates[:, 0]
    node_y = coordinates[:, 1]
    node_z = coordinates[:, 2]

    node_trace = go.Scatter3d(
        x=node_x, y=node_y, z=node_z,
        mode='markers',
        marker=dict(size=2, color='blue'),
        text=[f'Node {i}, coord {coord}' for i,coord in zip(range(len(coordinates)),coordinates)],
        hoverinfo='text'
    )

    # Prepare additional points trace (if provided)
    additional_trace = None
    if additional_points is not None:
        additional_points = np.array(additional_points)
        additional_colors = np.array(additional_colors) if additional_colors is not None else ['red'] * len(additional_points)
        
        additional_x = additional_points[:, 0]
        additional_y = additional_points[:, 1]
        additional_z = additional_points[:, 2]

        additional_trace = go.Scatter3d(
            x=additional_x, y=additional_y, z=additional_z,
            mode='markers',
            marker=dict(size=7, color=additional_colors),
            text=[f'Point {i}' for i in range(len(additional_points))],
            hoverinfo='text'
        )

    # Create the figure
    data = [edge_trace, node_trace]
    
    if additional_trace is not None:
        data.append(additional_trace)
    
    fig = go.Figure(data=data)
    fig.update_layout(
        showlegend=False,
        scene=dict(
            xaxis=dict(showbackground=False),
            yaxis=dict(showbackground=False),
            zaxis=dict(showbackground=False),
        ),
        title="3D Graph Representation",
    )

    # Show the figure
    fig.show()

