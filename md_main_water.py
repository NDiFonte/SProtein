#!/usr/bin/env python3

# Chiara Faccio, Scuola Normale Superiore, Pisa
# chiara.faccio@sns.it
# March 2023

from md_link_to_main import *
from md_sp_link_to_main import *
 
if __name__ == '__main__':  
    
    # We choose the file.gro
    print('Hello, this code allows you to apply graph theory to Molecular Dynamics files. \nThe code uses MDAnalysis tools to access data in molecular dynamics trajectories. \nTo begin with, the topology and trajectory files need to be selected. The topology file can be in PSF, PDB, CRD, or GRO format. \nIf no trajectory files are provided, only a structure is loaded. \nThe supported trajectory file formats include single frames like \
 PDB, CRD, and GRO, and time-series data like DCD, XTC, TRR, and XYZ.\n ')
    
    #path, topol = select_topology_file()
    #print("Select trajectory")
    #answ = input('Do you want to upload a trajectory file? (yes/no) ')
    
    #if answ == 'yes':
    #    traj = select_trajectory_file()
    #    U = mda.Universe(topol, traj)
    #else:
    #    U = mda.Universe(topol)
    topol="Ann_Ca_PAL_FAD_WAT1_WAT2_restraints.gro"
    traj="centered_MD_0-150ns_dt100.xtc"
    U=mda.Universe(topol,traj)
    #print('What would you like to do?\n \
    #1) To find the shortest path for a possible proton transfer from protein to bulk\n \
    #2) To apply graph theory tools to water boxes')
    #choice=int(input())
    choice=1
    if choice==1:
        #selections=select_residues_and_atoms(U)
        selections=['resname WAT and name OW']
        h_dist=1.3
        boundary=2
        cutoff_distance_from_prot=7
        dist=3.5
        donor=True
        directed=False
        deg_threshold=5
        all_hydrogens = []
        for selection in selections:
            hydrogens =find_hydrogens_from_selection(U, selection,h_dist)
            all_hydrogens.extend(hydrogens)
        #print(all_hydrogens)
        #ox_atom,h1_atom,h2_atom,wat_residue=select_solvent(U)
        ox_atom,h1_atom,h2_atom,wat_residue="OW","HW1","HW2","SOL"
        str_ox_atom=f'resname {wat_residue} and name {ox_atom}'
        #ref_atom=select_ref_atom(U)
        ref_atom="resname PAL and resid 651 and name O2 and index 8562" 
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
            #print(all_hydrogens)
            #print(ref_atom)
            nodes,coords=extract_nodes(U,ref_atom,ox_atom,selections,cutoff_distance_from_prot)
            r_atom=U.select_atoms(ref_atom)
            start_coords=r_atom.positions[0]
            start_node = np.where((coords == start_coords).all(axis=1))[0][0]
            #print(start_node)
            Hydrogens,H_coords=extract_H_coordinates(U,all_hydrogens,h1_atom,h2_atom,cutoff_distance_from_prot+3)
            wat_nodes=nodes.select_atoms(f"name {ox_atom}")
            wat_coords=wat_nodes.positions

            ########GOT ALL OF THE ELEMENTS TO START WORKING WITH GRAPH THEORY

            Graph_water, _,adj_matrix_water =get_graph_undirected_sp(U,wat_coords,box,boundary,dist)
            deg_water=Graph_water.degree
            bulk_wat_nodes=[index_deg_pair[0] for index_deg_pair in deg_water if index_deg_pair[1]>=deg_threshold]
            bulk_wat_coords=[wat_coords[index_deg_pair[0]] for index_deg_pair in deg_water if index_deg_pair[1]>=deg_threshold]
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
            nnodes=len(nodes)
            paths=[]
            for wat_node in bulk_wat_coords:
                end_coords=wat_node
                end_node = np.where((coords == end_coords).all(axis=1))[0][0]            
                shortest_path,sp_length=find_shortest_path(Adj_matrix,start_node,end_node,nnodes)
                paths.append(shortest_path)
            unique_paths,destinations=filter_redundant_paths(paths)
            unique_paths=truncate_paths(unique_paths,destinations)
            for path in unique_paths:#WRITING ALL PATHS
                for node in path:
                    t=str(nodes[node])
                    fout6.write(t+"\n")
                    s=str(coords[node])
                    s=s.replace("[","")
                    s=s.replace("]","")
                    fout5.write(s+"\n")
                fout5.write("\n")
                fout6.write("\n")
            shortest_path=filter_shortest_path(unique_paths)
            for node in shortest_path:
                t=str(nodes[node])
                fout4.write(t+"\n")
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


    if choice==2:
        #fmt_und = 'Do you want undirected graph? (yes/no) (directed graphs are implemented only in the case of orthogonal boxes) '
        undirected = 'yes'#input(fmt_und)    
        if undirected.lower() == 'yes':
            par = parameters_graph_undirected(topol)
            dist = par['dist']                       # threshold for edges
            boundary = par['boundary']               # (1 = no boundary, 2 = pbc)
            list_molecules = par['list_molecules']   # which molecules consider
            list_nodes = par['list_nodes']           # list_nodes = for each molecule which atom to consider as a node
            nodes_to_save = par['nodes_to_save']     # nodes_to_save = the name of the molecule for which you would like to save the values of centrality measures of its nodes
            
            n_frames = len(U.trajectory)
            
            begin_loop = True
            
            while begin_loop:
            
                fmt = "What do you want to do?\n \
                1) to compute centrality measures \n \
                2) to compute other metrics \n \
                3) to plot a single box in 3D \n \
                4) to plot a single box in 2D \n \
                5) to save adjacency matrix in matlab format \n \
                6) to consider dynamic networks \n \
                7) to consider a weighted graph \n \
                8) to study the patches (unweighted, undirected graphs. It requires TC already computed) \n \
                9) exit \n "
                
                
                res = int(input(fmt))
                
                if res == 1:
                    fmt2 = 'Which centrality measures do you want to compute? \n \
                    1) node total communicability \n \
                    2) degree \n \
                    3) subgraph centrality \n \
                    4) closeness centrality \n \
                    5) betweenness centrality \n \
                    6) Katz centrality \n \
                    7) eigenvector centrality \n \
                    8) non-backtracking total communicability \n '
                    
                    res2 = input(fmt2)
                    if int(res2) == 1:
                    
                        res1_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path)
                                           
                    elif int(res2) == 2:
                        
                        res_1_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path)
                        
                    elif int(res2) == 3:
                        
                        res_1_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path)                 
                                                 
                    elif int(res2) == 4:
                        
                        res_1_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                                           
                    elif int(res2) == 5:
                    
                        res_1_5(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif int(res2) == 6:
                        
                        res_1_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif int(res2) == 7:
                        
                        res_1_7(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif int(res2) == 8:
                        
                        res_1_8(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                elif res == 2:
                    fmt = 'What do you want to do? \n \
                    1) to count edges and cycles of length 3, 4 and 5 \n \
                    2) to compute Estrada index \n \
                    3) to compute Watts-Strogatz clustering coefficient \n \
                    4) to compute transitivity index \n \
                    5) to compute bipartivity measures \n \
                    6) to compute the degree assortativity coefficient \n \
                    7) to compute eigenvalues \n \
                    8) sparsity pattern of a frame \n \
                    9) entropy of the graph with subgraph \n \
                    10) entropy of the graph with TC \n \
                    11) entropy of the graph with Laplacian (Von Neumann entropy) \n \
                    12) mean value of min/max max_min_eigenvalues (eigenvalues already calculated) \n \
                    13) mean values of the density of the graph \n \
                    14) mean values of the density of the boxes ({number of molecules}/{volume box}) \n \
                    15) energy of the graph (eigenvalues already calculated) \n \
                    16) ASPL and diameter \n \
                    17) algebraic connectivity \n'

                    res_metric = int(input(fmt))
                    
                    if res_metric == 1:
                        
                        res_2_1(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                                              
                    elif res_metric == 2:
                        
                        res_2_2(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                        
                    elif res_metric == 3:
                    
                        res_2_3(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                        
                    elif res_metric == 4:
                    
                        res_2_4(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                                            
                    elif res_metric == 5:    
                    
                        res_2_5(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                        
                    elif res_metric == 6:
                    
                        res_2_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 

                    elif res_metric == 7:

                        res_2_7(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                            
                    elif res_metric == 8:
                    
                        res_2_8(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                      
                    elif res_metric == 9:
                    
                        res_2_9(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif res_metric == 10:
                    
                        res_2_10(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif res_metric == 11:
                    
                        res_2_11(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                        
                    elif res_metric == 12:
                    
                        res_2_12(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                            
                    elif res_metric == 13:
                    
                        res_2_13(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                        
                    elif res_metric == 14: 
                    
                        res_2_14(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                                       
                    elif res_metric == 15:
                    
                        res_2_15(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                        
                    elif res_metric == 16:
                    
                        res_2_16(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                        
                    elif res_metric == 17:
                    
                        res_2_17(U, list_molecules, list_nodes, boundary, dist, n_frames, path)
                        
                    else:
                        print('Value not valid')
                        continue
                    
                    #endif              
                        
                elif res == 3:
                
                    res_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif res == 4:
                
                    res_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                elif res == 5:
                
                    res_5(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 

                elif res == 6:
                    
                    fmt_1 = 'What do you want to do?\n \
                    1) to compute dynamic communicability \n \
                    2) to compute aggregated degree \n'
                    
                    resp = int(input(fmt_1))
                    
                    if resp == 1:
                    
                        res_6_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif resp == 2:
                    
                        res_6_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif resp == 3:
                        
                        res_6_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                                         
                    
                elif res == 7:    
                
                    fmt_type = 'How are edge weights defined? \n \
                    1) w(i,j) = 1/d(i,j)  \n \
                    2) w(i,j) = 1/d^2(i,j)  \n \
                    3) w(i,j) = 1/d^3(i,j)  \n \
                    4) w(i,j) = 1/sqrt(d(i,j)) \n '
                    
                    weight_edges = int(input(fmt_type))  
                    if weight_edges == 1:
                        str_weight = '1fracd'
                    elif weight_edges ==2:
                        str_weight = '1fracd2'
                    elif weight_edges ==3:
                        str_weight = '1fracd3'
                    elif weight_edges ==4:
                        str_weight = '1fracsqrt(d)'
                    
                    
                    fmt_7 = 'What do you want to do? \n \
                    1) to compute degree e TC \n \
                    2) to compute Estrada index \n \
                    3) entropy of the graph with TC \n \
                    4) entropy of the graph with Laplacian (Von Neumann entropy) \n'
                    
                    res_7 = int(input(fmt_7))
                  
                    if res_7 == 1:
                    
                        res_7_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, str_weight, weight_edges) 
                            
                    elif res_7 == 2:
                    
                        res_7_2(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) 
                        
                    elif res_7 == 3:
                    
                        res_7_3(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) 
                        
                    elif res_7 == 4:
                    
                        res_7_4(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) 
                
                
                elif res == 8:
                
                    print('Insert the file with CM values to identify the patches ')
                    root = Tk()
                    root.withdraw()
                    file_txt_centrality = filedialog.askopenfilename()
                    root.destroy()
                
                    fmt_threshold = 'Insert the threshold  to identify the patches (molecule i belongs to a patch if CM(i) >= threshold): '
                    res_threshold = float(input(fmt_threshold))
                    
                    fmt_patches_8 = 'What do you want to do?\n \
                    0) to write information about the patches \n \
                    1) to analyze evolution patches with dimension >= 5 (calcola come evolve nel tempo TC patches, TC first shell in cui first shell è aggiornata per ogni t, TC degli altri nodi in cui altri nodi sono aggiornati per ogni t)\n \
                    2) to compute patches and first shell (per ogni t calcola i nodi che sono nei patches e nella first shell. La first shell è riportata solo se dim(patches) >= 5)\n \
                    3) to compute average lifetime patches \n \
                    4) to compute the NTC values of the nodes in a specific patch troughout the trajectory (requires NTC and patches already computed) \n \
                    5) to compute patches and first shell with a threshold on TC only \n \
                    6) to compute patches and first and second shell \n \
                    7) to compute TC of specific nodes from prompt \n \
                    8) to generate VMD input file \n \
                    9) to plot the cardinality of biggest clusters \n \
                    10)generate VMD input file for 1st/2nd biggest clusters \n \
                    11)to compute the percentage of the total HDL in clusters \n'
                    
                    resp_res_8_patches = int(input(fmt_patches_8))
                    
                    if resp_res_8_patches == 0:
                        res_8_write(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold) 
                    elif resp_res_8_patches == 1:
                        res_8(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold) 
                    elif resp_res_8_patches == 2:
                        res_8_info_patches(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold)
                    elif resp_res_8_patches == 3:
                        res_8_lifetime_patches(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold)
                    elif resp_res_8_patches == 4:
                        res_8_follow_node_TC(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, res_threshold)
                    elif resp_res_8_patches == 5:
                        res_8_patches_under_TC_constraint(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, res_threshold)                    
                    elif resp_res_8_patches == 6:
                        res_8_info_patches_second_shell(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold)
                    elif resp_res_8_patches == 7:
                        res_8_follow_node_TC_from_prompt(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, res_threshold)    
                    elif resp_res_8_patches == 8:
                        res_8_gen_vmd_input(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold)
                    elif resp_res_8_patches == 9:
                        res_8_plot_cardinality(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold)
                    elif resp_res_8_patches == 10:
                        res_8_gen_vmd_input_bigg_2bigg(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold)
                    elif resp_res_8_patches == 11:
                        res_8_HDLclu_vs_tot(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold)
                    else:
                        print('Value not valid')
                        continue
                
                elif res == 9:
                    sys.exit(1)
                    
                else:
                    print('Value not valid')
                    continue
                    
        else:

            par = parameters_graph_directed(topol)
            dist = par['dist']                       # threshold for edges
            boundary = par['boundary']               # (1 = no boundary, 2 = pbc)
            list_molecules = par['list_molecules']   # which molecules consider
            list_nodes = par['list_nodes']           # list_nodes = for each molecule which atom to consider as a node
            nodes_to_save = par['nodes_to_save']     # nodes_to_save = the name of the molecule for which you would like to save the values of centrality measures of its nodes
            
            n_frames = len(U.trajectory)
            
            begin_loop = True
            
            while begin_loop:
            
                fmt = 'What do you want to do?\n \
                1) to compute centrality measures \n \
                2) exit \n'
                
                res = int(input(fmt))
                
                if res == 1:
                    fmt2 = 'Which centrality measures do you want to compute? \n \
                    1) total communicability e^{beta A} \mathbf{1} and e^{beta A^T} \mathbf{1} \n \
                    2) degree A \mathbf{1} and A^T \mathbf{1}  \n \
                    3) eigenvector centrality left and right of A \n \
                    4) HITS \n \
                    5) gTC \n \
                    6) TC of the bipartization of A \n'
                    
                    res2 = input(fmt2)
                    
                    if int(res2) == 1:
                    
                        res_10_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif int(res2) == 2:
                    
                        res_10_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif int(res2) == 3:
                    
                        res_10_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif int(res2) == 4:
                    
                        res_10_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif int(res2) == 5:
                    
                        res_10_5(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                        
                    elif int(res2) == 6:
                    
                        res_10_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path)                     
                        
                    else:
                        print('Value not valid')
                        continue

                
                elif res == 2:
                    sys.exit(1)
                    
                else:
                    print('Value not valid')
                    continue
                        

        sys.exit(1)
    
    
    
