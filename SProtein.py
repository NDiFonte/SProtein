
from md_link_to_main import *
from md_sp_link_to_main import *
 
if __name__ == '__main__': 

    #SETUP PARAMETERS
    topol="structure.gro" #STRUCTURE FILE NAME
    traj="traj.xtc" #TRAJECTORY FILE NAME
    U=mda.Universe(topol,traj) #CREATING MDA UNIVERSE OBJECT
    selections=[] # SPECIFIC SELECTIONS TO BE INCLUDED IN THE GRAPH, FOR EXAMPLE AMINOACIDS (e.g. ['resname GLU and name O2', resname ASP and name O1 O2'])
    h_dist=1.3 #DISTANCE TO DETERMINE HYDROGEN CLOSE TO SELECTION
    boundary=2 #BOUNDARY=1 NO PBC, BOUNDARY=2 PBC 
    cutoff_distance_from_prot=5 #CUTOFF DISTANCE FROM PROTEIN IN ANGSTROM (same residue as within X of protein) 
    dist=4 # CUTOFF DISTANCE TO CONNECT NODES WITH AN EDGE
    donor=True #SPECIFY IF SOURCE NODE IS A DONOR OR ACCEPTOR (ONLY MEANINGFUL IF WORKING WITH DIRECTED GRAPHS)
    directed=False #TOGGLE DIRECTED/UNDIRECTED GRAPHS
    deg_threshold=7 #DEGREE THRESHOLD TO IDENTIFY BULK WATER NODES (number of WATER neighbors within dist)
    ref_atom="resname PAL and resid 651 and name O2 and index 8562" #REFERENCE ATOM AS SOURCE NODE (e.g."resname PAL and resid 651 and name O2 and index 8562")
    all_hydrogens = []
    

    for selection in selections:
        hydrogens =find_hydrogens_from_selection(U, selection,h_dist)
        all_hydrogens.extend(hydrogens)
    ox_atom,h1_atom,h2_atom,wat_residue=select_solvent(U)
    str_ox_atom=f'resname {wat_residue} and name {ox_atom}'
    bulk_atom=""
    fout3=open("traj_sp.xyz",'w')
    fout4=open("traj_path_nodes.dat",'w')
    fout5=open("traj_all_paths.xyz",'w')
    fout6=open("traj_all_paths_nodes.dat",'w')
    f=1
    print("\n")
    print("Calculating shortest path...")
    for timestep in U.trajectory:
        fout3.write(f"FRAME {f}"+"\n")
        fout4.write(f"FRAME {f}"+"\n")
        fout5.write(f"FRAME {f}"+"\n")
        fout6.write(f"FRAME {f}"+"\n")
        box=timestep.dimensions
        nodes,coords=extract_nodes(U,ref_atom,ox_atom,selections,cutoff_distance_from_prot)
        r_atom=U.select_atoms(ref_atom)
        start_coords=r_atom.positions[0]
        start_node = np.where((coords == start_coords).all(axis=1))[0][0]
        Hydrogens,H_coords=extract_H_coordinates(U,all_hydrogens,h1_atom,h2_atom,cutoff_distance_from_prot+3)
        wat_nodes=nodes.select_atoms(f"name {ox_atom}")
        wat_coords=wat_nodes.positions

        ########GOT ALL OF THE ELEMENTS TO START WORKING WITH GRAPH THEORY

        Graph_water, _,adj_matrix_water =get_graph_undirected_sp(U,wat_coords,box,boundary,dist) #BUILDING WATER GRAPH
        deg_water=Graph_water.degree
        bulk_wat_nodes=[index_deg_pair[0] for index_deg_pair in deg_water if index_deg_pair[1]>=deg_threshold]
        bulk_wat_coords=[wat_coords[index_deg_pair[0]] for index_deg_pair in deg_water if index_deg_pair[1]>=deg_threshold] #FOUND BULK WATER NODES
        #ADJACENCY MATRIX MANIPULATION
        if not directed:
            _,_,Adj_matrix=get_graph_undirected_sp(U,coords,box, boundary,dist)
            Adj_matrix=get_sorted_edges(Adj_matrix)
            Adj_matrix.eliminate_zeros()
        else:
            Adj_matrix=get_graph_directed_sp(U,coords,nodes,Hydrogens,H_coords,box, boundary,dist)
            if donor:
                Adj_matrix=Adj_matrix.transpose()
                Adj_matrix.eliminate_zeros()
                Adj_matrix=get_sorted_edges(Adj_matrix)
            else:
                Adj_matrix.eliminate_zeros()
                Adj_matrix=get_sorted_edges(Adj_matrix)
        #END OF ADJACENCY MATRIX OPERATIONS, GLOBAL GRAPH BUILT
        nnodes=len(nodes)
        paths=[]
        for wat_node in bulk_wat_coords: #ITERATING OVER BULK WATER NODES
            end_coords=wat_node
            end_node = np.where((coords == end_coords).all(axis=1))[0][0]            
            shortest_path,sp_length=find_shortest_path(Adj_matrix,start_node,end_node,nnodes)#FINDING SHORTEST PATH
            paths.append(shortest_path)
        unique_paths,destinations=filter_redundant_paths(paths)
        unique_paths=truncate_paths(unique_paths,destinations)
        for path in unique_paths:#WRITING ALL PATHS
            for node in path:
                t=str(nodes[node])
                ggg=t.split()
                index=ggg[1].strip(":")
                resn=ggg[8].strip(",")
                resid=ggg[10]
                fout6.write(index+" "+resn+" "+resid+"\n")
                s=str(coords[node])
                s=s.replace("[","")
                s=s.replace("]","")
                fout5.write(s+"\n")
            fout5.write("\n")
            fout6.write("\n")
        shortest_path=filter_shortest_path(unique_paths)#FINDING SHORTEST PATH OVERALL
        for node in shortest_path:#WRITING SHORTEST PATH
            t=str(nodes[node])
            ggg=t.split()
            index=ggg[1].strip(":")
            resn=ggg[8].strip(",")
            resid=ggg[10]
            fout4.write(index+" "+resn+" "+resid+"\n")
            s=str(coords[node])
            s=s.replace("[","")
            s=s.replace("]","")
            fout3.write(s+"\n")               
        fout3.write("\n")
        fout4.write("\n")
        fout5.write("\n")
        fout6.write("\n")
        fout3.flush()
        fout4.flush()
        fout5.flush()
        fout6.flush()
        f=f+1
    fout3.close()
    fout4.close()
    fout5.close()
    fout6.close()

