# This script contains the responses to the various "if"s of the "md_main_water.py" script

# Chiara Faccio, Scuola Normale Superiore, Pisa
# chiara.faccio@sns.it
# March 2023

from md_functions_water_network import *

def res1_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):
    
    fmt = 'Do you want normalized NTC? \n \
    1) no \n \
    2) yes, w.r.t. number of edges \n \
    3) yes, w.r.t. total network communicability (TNC) \n \
    4) yes, w.r.t. exp( beta * lambda_1) \n \
    5) yes, w.r.t. number of nodes in the graph \n \
    6) yes, w.r.t. the TNC of a complete graph \n \
    7) yes, w.r.t. the average TNC  of 100 random graphs with the same number of nodes and edges as the original graph  \n\
    8) yes, w.r.t the average degree \n'
    
    normalization = int(input(fmt))
    
    fmt = 'Value of beta (beta > 0) = '
    beta = float(input(fmt))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt = 'Value of beta (beta > 0) = '
        beta = float(input(fmt))
    
    
    with open(os.path.join(path, 'NTC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(normalization) + '_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        fobj.write('Type of normalization = {},  molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(normalization, list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    Tc_mean = []
    
    for aaa in range(n_frames):
    
        tc = compute_NTC(U, aaa, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, normalization, dist)
        Tc_mean.append(np.mean(tc))
                       
    print('mean value TC = ', np.mean(Tc_mean))  
    return
    
    
def res_1_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):

    deg_m = []

    with open(os.path.join(path, 'DEGREE'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for bbb in range(n_frames):
    
        deg = compute_degree(U, bbb, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        
        deg_m.append(np.mean(deg))
        
    print('mean degree = ', np.mean(deg_m))
    
    return
                    

def res_1_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):

    fmt = 'Value of beta (beta > 0) = '
    beta = float(input(fmt))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt = 'Value of beta (beta > 0) = '
        beta = float(input(fmt))
    
    with open(os.path.join(path, 'SUBGRAPH_beta_'+ str(beta) + '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for ccc in range(n_frames):
    
        _ = compute_subgraph(U, ccc, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    
    return
                    
 
def res_1_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    with open(os.path.join(path, 'CLOSENESS_boundary_' + str(boundary)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for ddd in range(n_frames):
    
        _ = compute_CLOSENESS(U, ddd, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    
    return
    

def res_1_5(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    with open(os.path.join(path, 'BETWEENNESS_boundary_' + str(boundary)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for eee in range(n_frames):
    
        _ = compute_BETWEENNESS(U, eee, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    
    return
    
    
def res_1_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    fmt = 'Value of alpha = '
    alpha = float(input(fmt))

    with open(os.path.join(path, 'KATZ_'+str(alpha)+ '_boundary_' + str(boundary)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for fff in range(n_frames):
    
        _ = compute_KATZ(U, fff, alpha, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    
    return
    
    
def res_1_7(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    with open(os.path.join(path, 'EIGENVECTOR'+ '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for ggg in range(n_frames):
    
        compute_eigenvector(U, ggg, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        
    return
    
 
def res_1_8(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    fmt = 'Value of beta (beta > 0) = '
    beta = float(input(fmt))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt = 'Value of beta (beta > 0) = '
        beta = float(input(fmt))
    
    
    with open(os.path.join(path, 'NBT_TC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for hhh in range(n_frames):
        compute_NBT_total_communicability(U, hhh, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        
    return
    
    
def res_2_1(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
    edges = []
    cycles3 = []
    cycles4 = []
    cycles5 = []

    for iii in range(n_frames):
        n_edges, F2, F5, F8 = count_cycles(U, iii, list_molecules, list_nodes, path, boundary, dist) # F2 = number cycles of length 3; F5 = number cycles of length 4; F8 = number cycles of length 5
            
        edges.append(n_edges)
        cycles3.append(F2)
        cycles4.append(F5)
        cycles5.append(F8)
        
    print('Average number of edges = {} \n'.format(np.mean(edges)))
    print('Average number of cycles of length 3 = {} \n'.format(np.mean(cycles3)))
    print('Average number of cycles of length 4 = {} \n'.format(np.mean(cycles4)))
    print('Average number of cycles of length 5 = {} \n'.format(np.mean(cycles5)))
    
    return
    

def res_2_2(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
    estrada_index_vect = []
    
    for jjj in range(n_frames):  
    
        estrada_index_values = estrada_index(U, jjj, list_molecules, list_nodes, path, boundary, dist)
        estrada_index_vect.append(estrada_index_values)
        
    print('Average number of estrada index = {} \n'.format(np.mean(estrada_index_vect)))
    
    return
    
    
def res_2_3(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :

    watts_strogatz_vect = []
    
    for kkk in range(n_frames):   
        watts_strogatz_values = watts_strogatz(U, kkk, list_molecules, list_nodes, path, boundary, dist)
        watts_strogatz_vect.append(watts_strogatz_values)
        
    print('Average number of  Watts-Strogatz coeff. clustering = {} \n'.format(np.mean(watts_strogatz_vect)))
    
    return
    
    
def res_2_4(U, list_molecules, list_nodes,  boundary, dist, n_frames, path) :
                
    transitivity_index_vect = []
    
    for lll in range(n_frames):   
        transitivity_index_values = transitivity_index(U, lll, list_molecules, list_nodes, path, boundary, dist)
        transitivity_index_vect.append(transitivity_index_values)
        
    print('Average number of transitivity index = {} \n'.format(np.mean(transitivity_index_vect)))
    
    return

    
def res_2_5(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
                    
    fmt_compute_eig = 'Have you already calculated the eigenvalues? (yes/no) (For each graph, the eigenvalues are written in a .txt file. All .txt files are inside a folder)  \n'

    compute_eig = input(fmt_compute_eig)
    
    if compute_eig == 'yes':
        
        G, _ = get_graph_G(U, 0, list_molecules, list_nodes, boundary, dist)
        n = len(G)
        bipartivity_for_files(path, n)
        
    elif compute_eig == 'no':
    
        bipartivity_vect = []
        path3 = os.path.join(path, 'eigenvalues' )
        try:
            os.mkdir(path3)
        except:
            print('Unable to create folder because it already exists')
            return
           
    
        for mmm in range(n_frames):   
            bipartivity_values = bipartivity(U, mmm, list_molecules, list_nodes,  path, path3, boundary, dist)
            bipartivity_vect.append(bipartivity_values)
        
        print('Average number of bipartivity valuex = {} \n'.format(np.mean(bipartivity_vect)))
            
    else: 
        print('Invalid answer')
        
    return

    
def res_2_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
                
    assortativity_vect = []
    path2_path1_vect = []
    path3_path2_vect = []
    
    for nnn in range(n_frames):  
        assortativity_values, path2_path1, path3_path2  = assortativity(U, nnn, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        assortativity_vect.append(assortativity_values)
        path2_path1_vect.append(path2_path1)
        path3_path2_vect.append(path3_path2)
        
    print('Average number of degree assortativity coefficient = {} \n'.format(np.mean(assortativity_vect)))
    print('Average number of P2/P1 = {} \n'.format(np.mean(path2_path1_vect)))
    print('Average number of P3/P2 = {} \n'.format(np.mean(path3_path2_vect)))
    
    return
    
    
def res_2_7(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :

    path3 = os.path.join(path, 'eigenvalues' )
    try:
        os.mkdir(path3)
    except:
        print('Unable to create folder because it already exists')
        return

    for ooo in range(n_frames):   
        eigenvalues(U, ooo, list_molecules, list_nodes, path3, boundary, dist)
        
    return
    
    
def res_2_8(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
    fmt_frame = 'Which frame do you want to see? (Python starts counting from 0) \n'

    frame = int(input(fmt_frame))

    sparsity_A(U, frame, list_molecules, list_nodes, boundary, dist, path)
    return
    
    
def res_2_9(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    

    entropy_vect = []
    fmt = 'Value of beta = '
    beta = float(input(fmt))
    
    for ppp in range(n_frames):   
    
        entropy_values = entropy_subgraph(U, ppp, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        entropy_vect.append(entropy_values)
        
    print('Average number of entropy (using subgraph) = {} \n'.format(np.mean(entropy_vect)))
        
    return
    
    
def res_2_10(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
                    
    entropy_vect = []
    fmt = 'Value of beta = '
    beta = float(input(fmt))
    
    for qqq in range(n_frames): 
    
        entropy_values = entropy_TC(U, qqq, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        entropy_vect.append(entropy_values)
        
    print('Average number of entropy (using TC) = {} \n'.format(np.mean(entropy_vect)))
    
    return
    
    
def res_2_11(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :

    entropy_vect = []
    
    for rrr in range(n_frames):
    
        entropy_values = entropy_von_neumann(U, rrr, list_molecules, list_nodes, path, boundary, dist)
        entropy_vect.append(entropy_values)
        
    print('Average number of entropy (Von Neumann entropy) = {} \n'.format(np.mean(entropy_vect)))
    return

    
def res_2_12(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :

    G, _ = get_graph_G(U, 0, list_molecules, list_nodes, boundary, dist)
    n = len(G)
    max_min_eigenvalues(n)  
    return


def res_2_13(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
            
    density_g = []

    for sss in range(n_frames):
        density_ = density_graph(U, sss, list_molecules, list_nodes, boundary, dist, path)
        density_g.append(density_)
        
    print('mean graph density ', np.mean(density_g))
    return


def res_2_14(U, list_molecules, list_nodes,  boundary, dist, n_frames, path) :
    density_v = []

    for ttt in range(n_frames):
        G, _ = get_graph_G(U, ttt, list_molecules, list_nodes, boundary, dist)
        n = len(G)
        density__ = density(U, ttt, list_molecules, list_nodes, n)
        density_v.append(density__)
        
    print('mean box density ', np.mean(density_v))
    
    with open(os.path.join(path, 'density.txt'), 'a+') as fobj:
        for ii in range(len(density_v)):
            fobj.write('{:f} \n'.format(density_v[ii]))
        fobj.write('\n')
        
    return


def res_2_15(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :

    G, _ = get_graph_G(U, 0, list_molecules, list_nodes, boundary, dist)
    n = len(G)
    energy(path, n)
    
    return
                    
    
def res_2_16(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
                
    ASPL_vect = []
    diam_vect = []
    isolated_vect = []
    
    for uuu in range(n_frames):
    
        ASPL, diam, n_isolated = ASPL_diameter_isolated(U, uuu, list_molecules, list_nodes, boundary, dist)
        ASPL_vect.append(ASPL)
        diam_vect.append(diam)
        isolated_vect.append(n_isolated)
        
        with open(os.path.join(path, 'ASPL_diameter_n_isolated.txt'), 'a+') as fobj:
            fobj.write('ASPL = {:f} \n'.format(ASPL))
            fobj.write('diameter = {:f} \n'.format(diam))
            fobj.write('n. isolated points = {:f} \n'.format(n_isolated))
            fobj.write('\n')
        
    print('Average number of ASPL = {} \n'.format(np.mean(ASPL_vect)))
    print('Average number of diam = {} \n'.format(np.mean(diam_vect)))
    print('Average number of isolated points = {} \n'.format(np.mean(isolated_vect)))
    
    return
    

def res_2_17(U, list_molecules, list_nodes, boundary, dist, n_frames, path):
                
    AG_vect = []
    for kkkk in range(n_frames):
    
        AG = algebraic_connectivity(U, kkkk, list_molecules, list_nodes, boundary, dist)
        AG_vect.append(AG)
        
        with open(os.path.join(path, 'algebraic_connectivity.txt'), 'a+') as fobj:
            fobj.write('AG = {:f} \n'.format(AG))
            fobj.write('\n')
        
    print('Average number of algebraic connectivity = {} \n'.format(np.mean(AG)))     
    
    return
    
    
def res_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
            
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())

    fmt_frame = 'Which frame do you want to plot? (Python starts counting from 0) \n'
    
    frame = int(input(fmt_frame))
    
    fmt = 'Do you want to color the nodes according to some measure of centrality already calculated? (yes / no) '
    
    res_c = input(fmt)
    
    if res_c.lower() == 'no':
        _, _, _, coord, _ , _= get_graph(U, frame, list_molecules, list_nodes, 1, dist, nodes_to_save)
        coord2 = coord[range_nodes,:]
        plot_network_3d(coord2, dist, path, 'plot 3D')
        
    elif res_c.lower() == 'yes':
        print('Choose the file.txt with the values: \n')
        
        root = Tk()
        root.withdraw()
        file_txt_centrality = filedialog.askopenfilename()
        root.destroy()
        
        fmt_color = 'Chromatic scale or two colors to identify LDL and HDL : (chro/bicolor)   '
        
        res_color = input(fmt_color)
        if res_color == 'bicolor':
            fmt_bicolor = 'Threshold to identify LDL/HDL phase : '
            try:
                threshold_color = float(input(fmt_bicolor))
                print('If the value of the node is <= threshold, then the node is in LDL phase (BLUE), otherwise it is in HDL phase (RED)')
            except ValueError:
                print('ERROR: value not valid')
                return
        
        color = []
        values = []
        
        n_nodes = len(range_nodes)
        
        
        with open(file_txt_centrality, 'r') as data:
            data.seek(0)
            for ff in range(((n_nodes)+1)*frame +2):
                data.readline()

            for i in range(n_nodes):
                cent = data.readline()
                 
                nn_for_read = 0
                while cent[nn_for_read] != '>':
                    nn_for_read +=1

                cent = float(cent[nn_for_read+1 : nn_for_read+25])
                values.append(cent)
                if res_color == 'chro':
                    color.append(cent)                    
                elif res_color == 'bicolor':
                    if cent <= threshold_color:
                        color.append(0)
                    else:
                        color.append(1)
        print(min(values))
        fmt = 'Do you want vmin and vmax by default? (yes/no) '
        res = input(fmt)
        if res.lower() == 'yes':
            vmin = None
            vmax = None
        elif res.lower() == 'no':
            fmt = 'vmin = '
            vmin = float(input(fmt))
            fmt = 'vmax = '
            vmax = float(input(fmt))
        else: 
            print('Value not valid')
            return
        
        fmt = 'What name do you want to save the plot with?  '
        name_img = input(fmt)                
        
        _, _, _, coord, _, _= get_graph(U, frame, list_molecules, list_nodes, 1, dist, nodes_to_save) 
        coord2 = coord[range_nodes,:]           
        plot_network_3d(coord2, dist, path, name_img, color = color, vmin = vmin, vmax =vmax)
        
    else:
        print('Value not valid')
        
    return


def res_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    fmt_frame = 'Which frame do you want to plot? (Python starts counting from 0) \n'
    
    frame = int(input(fmt_frame))
    
    fmt = 'Do you want to color the nodes according to some measure of centrality already calculated? (yes / no) '
    
    res_c2 = input(fmt)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    n_nodes = len(range_nodes)
    
    if res_c2.lower() == 'no':
    
        _, G, pos, _, _, _ = get_graph(U, frame, list_molecules, list_nodes, 1, dist, nodes_to_save)
        
        plot_network_2d(pos, G, range_nodes, path, 'plot 2D')
        
    elif res_c2.lower() == 'yes':
        print('Choose the file.txt with the values: \n')

        root = Tk()
        root.withdraw()
        file_txt_centrality = filedialog.askopenfilename()
        root.destroy()
        
        fmt_color = 'Chromatic scale or two colors to identify LDL and HDL : (chro/bicolor)   '
        
        res_color = input(fmt_color)
        if res_color == 'bicolor':
            fmt_bicolor = 'Threshold to identify LDL/HDL phase : '
            try:
                threshold_color = float(input(fmt_bicolor))
                print('If the value of the node is <= threshold, then the node is in LDL phase (BLUE), otherwise it is in HDL phase (RED)')
            except ValueError:
                print('ERROR: value not valid')
                return
        
        color = []
        
        
        with open(file_txt_centrality, 'r') as data:
            
            data.seek(0)
            for ff in range(((n_nodes)+1)*frame +2):
                data.readline()

            for i in range(n_nodes):
                cent = data.readline()
                
                nn_for_read = 0
                while cent[nn_for_read] != '>':
                    nn_for_read +=1

                cent = float(cent[nn_for_read+1 : nn_for_read+25])

                if res_color == 'chro':
                    color.append(cent)                    
                elif res_color == 'bicolor':
                    if cent <= threshold_color:
                        color.append(0)
                    else:
                        color.append(1)     
        
        fmt = 'Do you want vmin and vmax by default? (yes/no) '
        res = input(fmt)
        if res.lower() == 'yes':
            vmin = None
            vmax = None
        elif res.lower() == 'no':
            fmt = 'vmin = '
            vmin = float(input(fmt))
            fmt = 'vmax = '
            vmax = float(input(fmt))
        else: 
            print('Value not valid')
            return
        
        fmt = 'What name do you want to save the plot with?  '
        name_img = input(fmt)                
        
        _, G, pos, _, _ , _= get_graph(U, frame, list_molecules, list_nodes, 1, dist, nodes_to_save)
        
        plot_network_2d(pos, G, range_nodes, path, name_img, color = color, vmin = vmin, vmax =vmax)
    
    
        
    else:
        print('Value not valid')
        
    return
      
    
def res_5(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
                    
    new_path = os.path.join(path, 'matrices_matlab_format' )
    try:
        os.mkdir(new_path)
    except:
        print('Unable to create folder because it already exists')
        return
    
    for vvv in range(n_frames):  
        save_matrix_matlab_format(U, vvv, list_molecules, list_nodes, boundary, dist, new_path)
        
    return


def res_6_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
                
    fmt = 'Value of alpha (0 < alpha < 1/sigma, where sigma = max{rho(A^t), t = 0,1,...,T}) = '
    alpha = float(input(fmt))
    
    while alpha <= 0:
        print('ERROR: alpha is a positive parameter')
        
        fmt = 'Value of alpha = '
        alpha = float(input(fmt))
    
    katz_dynamic_graph(U, list_molecules, list_nodes, boundary, dist, n_frames, path, alpha, nodes_to_save)
    
    return
    
    
def res_6_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
                
    aggregated_degree(U, list_molecules, list_nodes, boundary, dist, n_frames, path, nodes_to_save)
    return


def res_7_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, str_weight, weight_edges) :
            
    fmt_b = 'Value of beta for TC (beta > 0) = '
    beta = float(input(fmt_b))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt_b = 'Value of beta (beta > 0) = '
        beta = float(input(fmt_b))
    
    with open(os.path.join(path, 'weighted_NTC_beta_'+ str(beta)+'_undirected_boundary_'+str(boundary)+'_dist_'+str(dist)+'_weight_'+str(str_weight)+'.txt'), 'a+') as fobj:
        fobj.write('\n')
        
    with open(os.path.join(path, 'weighted_DEGREE_undirected_boundary_'+str(boundary)+'_dist_'+str(dist)+'_weight_'+str(str_weight)+'.txt'), 'a+') as fobj:
        fobj.write('\n')
        
    
    for cccc in range(n_frames):
        compute_weighted_NTC_degree(U, cccc, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist, weight_edges= weight_edges)
        
    return
    
    
def res_7_2(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) :
                
    fmt_b = 'Value of beta = '
    beta = float(input(fmt_b))
    
    fmt_method = 'What method do you use to compute the Estrada index? (1 = eigenvalues, 2 = subgraph centrality values) '
    method = int(input(fmt_method))
    
    EE_vect = []
    
    for dddd in range(n_frames):  
        EE = compute_weighted_Estrada_index(U, dddd, beta, list_molecules, list_nodes, path, boundary, dist, method, weight_edges = weight_edges)
        EE_vect.append(EE)
        
    print('Average weighted EE = ', np.mean(EE_vect))
    
    return
    
    
def res_7_3(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) :
                    
    entropy_vect = []
    
    for eeee in range(n_frames): 
    
        entropy_values = weighted_entropy_TC(U, eeee, list_molecules, list_nodes, dist, path, boundary, weight_edges = weight_edges)
        entropy_vect.append(entropy_values)
        
    print('Average number of entropy (using TC) = {} \n'.format(np.mean(entropy_vect)))
        
    return
    
    
def res_7_4(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) :
    entropy_vect = []
    
    for ffff in range(n_frames):  
    
        entropy_values = weighted_entropy_von_neumann(U, ffff, list_molecules, list_nodes, dist, path, boundary, weight_edges = weight_edges)
        entropy_vect.append(entropy_values)
        
    print('Average number of entropy (Von Neumann entropy) = {} \n'.format(np.mean(entropy_vect)))
    
    return
    
def res_8_write(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold) :
  
    
    with open(os.path.join(path, 'patches_using_TC_larger_'+str(res_threshold)+'.txt'), 'a+') as fobj:
            fobj.write('For each frame there is: \n \
    - the number of patches with length greater than 1, \n \
    - the cardinality of each patch, \n \
    - the indexes of the nodes that form the patch (compared to the Python notation they have already been increased by one, in this way the nodes go from 1 to 710). \n \
    - the corresponding TC values for each node, \n \
    - the corresponding deg values for each node. \n ')
            fobj.write('\n')
    
    for qqq in range(n_frames):
    
        _, _, _, lst = nodes_patches(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
        
        _, G, _, _, _, _ = get_graph(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save)     
        dict_nodes = range_nodes_to_save(U, qqq, list_molecules, list_nodes, nodes_to_save)
        range_nodes = list(dict_nodes.keys())
        n_nodes = len(range_nodes)
        degree = G.degree
        
        list_set_nodes_WITHOUT_1 = [list(set(x)) for x in lst if len(x)>1]
        NC = len(list_set_nodes_WITHOUT_1)
        
        with open(os.path.join(path, 'patches_using_TC_larger_'+str(res_threshold)+'.txt'), 'a+') as fobj:
            fobj.write('FRAME {} => {} patches of length > 1 \n'.format(qqq, NC))
        
        TC = []
        with open(file_txt_centrality, 'r') as data:
            data.seek(0)
            for ff in range(((n_nodes)+1)*qqq +2):
                data.readline()

            for i in range(n_nodes):
                cent = data.readline()
                 
                nn_for_read = 0
                while cent[nn_for_read] != '>':
                    nn_for_read +=1

                cent = float(cent[nn_for_read+1 : nn_for_read+25])
                TC.append(cent) 
            
        for c in lst:
        
            values_deg = [degree[i] for i in c]
            values_tc = [TC[i] for i in c]
            index_plus_1 = [i+1 for i in c]

                
            with open(os.path.join(path, 'patches_using_TC_larger_'+str(res_threshold)+'.txt'), 'a+') as fobj:
                fobj.write('Cardinality = {:3d}, [{}] \n'.format(len(c), " ".join(map(str,index_plus_1))))
                fobj.write('              {} \n'.format(values_tc))   
                fobj.write('              {} \n'.format(values_deg))   
                fobj.write('\n')

    return
    
    
def res_8(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold):

    n_frames = 400

    print('Insert file with TC values')
    root = Tk()
    root.withdraw()
    file_TC = filedialog.askopenfilename()
    root.destroy()
    
    values_TC = []
    A = get_graph_A(U, 0, list_molecules, list_nodes, boundary, dist) 
    n_nodes = A.shape[0]
    
    with open(file_TC, 'r') as data:
        data.seek(0)
        data.readline()
        
        for frame in range(n_frames):
            data.readline()
            for i in range(n_nodes):
                cent = data.readline()
                 
                nn_for_read = 0
                while cent[nn_for_read] != '>':
                    nn_for_read +=1

                cent = float(cent[nn_for_read+1 : nn_for_read+25])
                values_TC.append(cent)
                
    matrix_patches = np.asmatrix(np.zeros((n_frames, n_frames)))
    matrix_first_shell = np.asmatrix(np.zeros((n_frames, n_frames)))
    matrix_other_nodes = np.asmatrix(np.zeros((n_frames, n_frames)))
            
    
    for qqq in range(n_frames):
        print(qqq)
    
        _, patches_m_larger_5, _, _ = nodes_patches(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
        patches = [item for sublist in patches_m_larger_5 for item in sublist]
        
        TC_frames_qqq = values_TC[n_nodes*qqq : n_nodes*(qqq+1)]
        
        if len(patches_m_larger_5) == 0:
        
            all_nodes = list(np.arange(n_nodes))
            TC_other_nodes = [TC_frames_qqq[i] for i in all_nodes]
            matrix_other_nodes[qqq,qqq] = np.mean(TC_other_nodes)
            
            for rrr in range(qqq+1, n_frames):
            
                TC_frames_rrr = values_TC[n_nodes*rrr : n_nodes*(rrr+1)]
                
                TC_other_nodes = [TC_frames_rrr[i] for i in all_nodes]
                matrix_other_nodes[qqq,rrr] = np.mean(TC_other_nodes)
        else:

            TC_patches = [TC_frames_qqq[i] for i in patches]
            matrix_patches[qqq, qqq] = np.mean(TC_patches)
        
            # computes nodes in the first shell
            A = get_graph_A(U, qqq, list_molecules, list_nodes, boundary, dist) 
            ls = []
            for bb in patches:                            
                B = A[bb,:]; 
                a = csr_matrix.nonzero(B)
                ls = ls + list(a[1])                
            first_shell = list(set(ls) - set(patches))
            
            TC_first_shell = [TC_frames_qqq[i] for i in first_shell]
            matrix_first_shell[qqq,qqq] = np.mean(TC_first_shell)
            
            # computes the remaining nodes
            all_nodes = list(np.arange(n_nodes))
            first_diff = list(set(all_nodes) - set(first_shell))
            other_nodes = list(set(first_diff) - set(patches))
            
            TC_other_nodes = [TC_frames_qqq[i] for i in other_nodes]
            matrix_other_nodes[qqq,qqq] = np.mean(TC_other_nodes)
            
            for rrr in range(qqq+1, n_frames):
            
                TC_frames_rrr = values_TC[n_nodes*rrr : n_nodes*(rrr+1)]
                
                TC_patches = [TC_frames_rrr[i] for i in patches]
                matrix_patches[qqq, rrr] = np.mean(TC_patches)
            
                # computes nodes in the first shell
                A = get_graph_A(U, rrr, list_molecules, list_nodes, boundary, dist) 
                ls = []
                for bb in patches:                            
                    B = A[bb,:]; 
                    a = csr_matrix.nonzero(B)
                    ls = ls + list(a[1])                
                first_shell = list(set(ls) - set(patches))
                
                TC_first_shell = [TC_frames_rrr[i] for i in first_shell]
                matrix_first_shell[qqq,rrr] = np.mean(TC_first_shell)
                
                # computes the remaining nodes
                all_nodes = list(np.arange(n_nodes))
                first_diff = list(set(all_nodes) - set(first_shell))
                other_nodes = list(set(first_diff) - set(patches))
                
                TC_other_nodes = [TC_frames_rrr[i] for i in other_nodes]
                matrix_other_nodes[qqq,rrr] = np.mean(TC_other_nodes)
            
    mean_patches =  matrix_patches.sum(0)/(matrix_patches!=0).sum(0).astype(float)
    mean_first_shell = matrix_first_shell.sum(0)/(matrix_first_shell!=0).sum(0).astype(float)
    mean_other_nodes = matrix_other_nodes.sum(0)/(matrix_other_nodes!=0).sum(0).astype(float)
        
    pp = np.nan_to_num(np.array(mean_patches).squeeze())     
    ff = np.nan_to_num(np.array(mean_first_shell).squeeze())    
    oo = np.nan_to_num(np.array(mean_other_nodes).squeeze())    
    
    
    with open(os.path.join(path, 'TC_nodes_patches_NO_agg_n_frames_100.txt'), 'a+') as fobj:
        for i in range(n_frames):
            fobj.write('{}\n'.format(pp[i]))
        
    with open(os.path.join(path, 'TC_nodes_first_shell_NO_agg_n_frames_100.txt'), 'a+') as fobj:
        for i in range(n_frames):
            fobj.write('{}\n'.format(ff[i]))
        
    with open(os.path.join(path, 'TC_remaining_nodes_NO_agg_n_frames_100.txt'), 'a+') as fobj:
        for i in range(n_frames):
            fobj.write('{}\n'.format(oo[i]))
            
    fig1, ax = plt.subplots(); 
    x = np.arange(n_frames)
    plt.plot(x, pp)   #np.array(matrix_patches[i,:]).squeeze()
    plt.xlabel('frames')
    plt.ylabel('TC nodes in patches')
    plt.axhline(93, xmin=0,  color='black', linestyle='dashed')
    plt.ylim(0, 500)

    fig2, ax = plt.subplots(); 
    x = np.arange(n_frames)
    plt.plot(x, ff)
    plt.xlabel('frames')
    plt.ylabel('TC nodes in first shell')
    plt.axhline(93, xmin=0,  color='black', linestyle='dashed')
    plt.ylim(0, 500)
    
    fig3, ax = plt.subplots(); 
    x = np.arange(n_frames)
    plt.plot(x, oo)
    plt.xlabel('frames')
    plt.ylabel('TC nodes in first shell')
    plt.axhline(93, xmin=0,  color='black', linestyle='dashed')
    plt.ylim(0, 500)
    plt.show()      

    
    fig4, ax = plt.subplots(); 
    x = np.arange(n_frames)
    plt.plot(x, pp, label = 'nodes in patches')
    plt.plot(x, ff, label = 'nodes in first shell')
    plt.plot(x, oo, label = 'nodes in first shell')
    plt.xlabel('frames')
    plt.ylabel('TC')
    plt.axhline(93, xmin=0,  color='black', linestyle='dashed')
    plt.legend()
    plt.ylim(0, 500) 
    plt.show() 
    
    return
    
    
def res_8_info_patches(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold):
    
    with open(os.path.join(path, 'nodes_patches_first_shell.txt'), 'a+') as fobj:
        fobj.write('The id of the molecules has already been increased by 1, so the values range from 1 to 710. \n')
        fobj.write('\n')
    
    
    for qqq in range(n_frames):
        print(qqq)
        
        with open(os.path.join(path, 'nodes_patches_first_shell.txt'), 'a+') as fobj:
            fobj.write('FRAME {}: \n'.format(qqq))
    
        _, patches_m_larger_5, other_patches, _ = nodes_patches(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
        
        n_patches_larger_5 = len(patches_m_larger_5)
        
        for i, patch in enumerate(patches_m_larger_5):
        
            with open(os.path.join(path, 'nodes_patches_first_shell.txt'), 'a+') as fobj:
                fobj.write('number of molecules in the patch {} = {} \n'.format(i, len(patch)))
                fobj.write('id molecules in the patch {} = {} \n'.format(i, [hh+1 for hh in patch]))
                
            A = get_graph_A(U, qqq, list_molecules, list_nodes, boundary, dist) 
            ls = []
            for bb in patch:                            
                B = A[bb,:]; 
                a = csr_matrix.nonzero(B)
                ls = ls + list(a[1])                
            first_shell = list(set(ls) - set(patch))
             
            with open(os.path.join(path, 'nodes_patches_first_shell.txt'), 'a+') as fobj:
                fobj.write('number of molecules in the first shell of patch {} = {} \n'.format(i, len(first_shell)))
                fobj.write('id molecules in the first shell of patch {} = {} \n'.format(i, [hh+1 for hh in first_shell]))
                
        for ii, patch2 in enumerate(other_patches):
        
            with open(os.path.join(path, 'nodes_patches_first_shell.txt'), 'a+') as fobj:
                fobj.write('number of molecules in the patch {} = {} \n'.format(ii+n_patches_larger_5, len(patch2)))
                fobj.write('id molecules in the patch {} = {} \n'.format(ii+n_patches_larger_5, [hh+1 for hh in patch2]))
        
        with open(os.path.join(path, 'nodes_patches_first_shell.txt'), 'a+') as fobj:
             fobj.write('\n')    
    
        
            

    return
    
    
def res_8_lifetime_patches(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold):


    ans_initial_frame = 'From which frame do you want to start the analysis? '
    initial_frame = input(ans_initial_frame)
    n_frame = int(initial_frame)
    
    ans_final_frame = 'At what frame do you want to stop the analysis? (Type enter if you want to consider the end of the trajectory.) '
    resp_final_frame = input(ans_final_frame)
    if not resp_final_frame:
        final_frame = n_frames
    else:
        final_frame = int(resp_final_frame)
    
    fmt_jump = 'Insert the jump(ps): '
    jump = int(input(fmt_jump))

        
    y=[]
    cl_vect = []
    
    ans_threshold_dim_biggest_path = 'Do you also want to insert a threshold on the size of the biggest patch? (yes/no) '
    resp_dim = input(ans_threshold_dim_biggest_path)
    
    if resp_dim == 'yes':
        threshold_dim_biggest_path = int(input('Insert threshold on the dimension of the biggest patch = '))
        
        
    
    while n_frame <  final_frame - 1:
        print(n_frame)
    
        y_int = []
        cl_int = []
        
        biggest_patches, _, _, _  = nodes_patches(U, n_frame, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
        
        if resp_dim == 'yes':
            if len(biggest_patches) >= threshold_dim_biggest_path:
                
                with open(os.path.join(path, 'info_lifetime_jump_'+str(jump)+'ps_from_frame_'+str(initial_frame)+'to frame'+str(final_frame)+'.txt'), 'a+') as fobj:
                        fobj.write('frame = {}\n'.format(n_frame))
                        fobj.write(' {}\n'.format([ii+1 for ii in biggest_patches]))
        
                while (len(biggest_patches)) > 0:
                    y_int.append(len(biggest_patches))
                    n_frame += 1
                    if n_frame <= final_frame - 1:
                        _, patches_greater_5, other_patches, _  = nodes_patches(U, n_frame, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
                        
                        all_clusters = patches_greater_5 + other_patches  # creo una lista contenente tutti i patches
                        
                        inters = []
                        for i in range(len(all_clusters)):
                            int_new = list(set(biggest_patches) & set(all_clusters[i]))
                            inters.append(int_new)
                            
                        tuple_inters = max(enumerate(inters), key = lambda x: len(x[1]))
                        biggest_patches = tuple_inters[1]
                        cl = tuple_inters[0]
                        if len(biggest_patches) > 0:
                            cl_int.append(cl)
                            
                            with open(os.path.join(path, 'info_lifetime_jump_'+str(jump)+'ps_from_frame_'+str(initial_frame)+'to frame'+str(final_frame)+'.txt'), 'a+') as fobj:
                                fobj.write('frame = {}\n'.format(n_frame))
                                fobj.write(' {}\n'.format([ii+1 for ii in biggest_patches]))
                            
                    else:
                        break
                
                y.append(y_int)
                cl_vect.append(cl_int)
                n_frame += jump    # jump
                
            else:
                n_frame += 1
                
        else:
        
            while (len(biggest_patches)) > 0:
                y_int.append(len(biggest_patches))
                n_frame += 1
                if n_frame <= final_frame - 1:
                    _, patches_greater_5, other_patches, _  = nodes_patches(U, n_frame, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
                    
                    all_clusters = patches_greater_5 + other_patches  # creo una lista contenente tutti i patches
                    
                    inters = []
                    for i in range(len(all_clusters)):
                        int_new = list(set(biggest_patches) & set(all_clusters[i]))
                        inters.append(int_new)
                        
                    tuple_inters = max(enumerate(inters), key = lambda x: len(x[1]))
                    biggest_patches = tuple_inters[1]
                    cl = tuple_inters[0]
                    if len(biggest_patches) > 0:
                        cl_int.append(cl)
                        
                else:
                    break
                
            y.append(y_int)
            cl_vect.append(cl_int)
            n_frame += jump    # jump
        
        with open(os.path.join(path, 'info_lifetime_jump_'+str(jump)+'ps_from_frame_'+str(initial_frame)+'to frame'+str(final_frame)+'.txt'), 'a+') as fobj:
            fobj.write('\n')

    lifetime = [len(i) for i in y]
    print('max lifetime = ', max(lifetime))
    
    print('average lifetime of clusters = ', np.mean(lifetime))
    
    with open(os.path.join(path, 'lifetime_jump_'+str(jump)+'ps_from_frame_'+str(initial_frame)+'to frame'+str(final_frame)+'.txt'), 'a+') as fobj:
        fobj.write('Lifetime = {}\n'.format(lifetime))
        for i in range(len(y)):
            fobj.write('Number of nodes in the max overlap (first value is dim. of biggest patch) = {}\n'.format(y[i]))
            fobj.write('Patch with max overlap = {}\n'.format(cl_vect[i]))
            fobj.write('\n')


    return
    
def res_8_follow_node_TC (U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, res_threshold):

    print('Insert file with TC values')
    root = Tk()
    root.withdraw()
    nfile_TC = filedialog.askopenfilename()
    root.destroy()
    file_TC=open(nfile_TC,'r').read()
    print('Insert file with Patches')
    root = Tk()
    root.withdraw()
    nfile_patches = filedialog.askopenfilename()
    root.destroy()        
    file_patches=open(nfile_patches,'r').read()
    print('Which FRAME are you interested in ?')
    FRAME=input()
    FRAMEplus1=str(int(FRAME)+1)
    print('Which patch are you intersted in  ?')
    Patch=input()
    ##################################################    	
    TC_all=file_TC.splitlines()
    lines=2
    TC_values=[]
    nmol=710
    i=0
    TC_line=''
    while lines<len(TC_all):
        i=0
        while i<nmol:
           TC_line=TC_all[lines].split()
           TC_values.append(TC_line[14])
           i=i+1
           lines=lines+1
        lines=lines+1
    ##################################################Questo blocco ha generato una lista 1D di lunghezza nmol*nframe
    ##################################################Ciascun elemento è il valore della TC dei nodi in ordine
    lines_patches=file_patches.splitlines()
    lines=2
    arr_patch=[]
    while lines<len(lines_patches):
       if str(lines_patches[lines])=='FRAME '+FRAME+': ':
          lines=lines+1
          while str(lines_patches[lines+1])!='FRAME '+FRAMEplus1+': ':
             arr_patch.append(lines_patches[lines])
             lines=lines+1
          break
       lines=lines+1
    i=0
    sep_equal=[]
    sep_equal_nodes=[]
    dim_patch=0
    node_list=[]
    while i<len(arr_patch):
        sep_equal=arr_patch[i].split('=')
        if sep_equal[0]=='number of molecules in the patch '+Patch+' ' :
           dim_patch=int(sep_equal[1].strip(' '))
           sep_equal_nodes=arr_patch[i+1].split('=')
           node_list=sep_equal_nodes[1].strip(' ][,')
           break
        i=i+1
    nodes=node_list.split()
    i=0
    arr_nodes=[]
    while i<len(nodes):
       node=str(nodes[i].replace(",",""))
       arr_nodes.append(node)
       i=i+1
    print('Nodes in the patch   '+ str(arr_nodes))
    ##################################################Questo blocco ha generato la lista di nodi del patch di interesse
    TC_values=TC_values
    arr_nodes=arr_nodes
    i=0
    with open(os.path.join(path, 'TC_nodes_from_patch_'+Patch+'_at_FRAME_'+FRAME+'_.txt'), 'a+') as fobj:
        fobj.write('#Each column is different NODE of the patch '+Patch+ 'at_FRAME'+FRAME+'\n')
        fobj.write('\n')
        while i<len(arr_nodes):
           fobj.write('#NODE '+str(arr_nodes[i])+'\n')
           j=int(arr_nodes[i])-1
           while j<len(TC_values):
              fobj.write(TC_values[j]+'\n')
              j=j+nmol
           fobj.write('\n')
           i=i+1

def res_8_patches_under_TC_constraint (U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, res_threshold):
    
    kk=0
    n_frames=400
    print('Insert file with TC values')
    root = Tk()
    root.withdraw()
    nfile_TC = filedialog.askopenfilename()
    root.destroy()
    file_TC=open(nfile_TC,'r').read()	
    TC_all=file_TC.splitlines()
    lines=2
    TC_values=[]
    nmol=710
    i=0
    TC_line=''
    while lines<len(TC_all):
        i=0
        while i<nmol:
           TC_line=TC_all[lines].split()
           TC_values.append(TC_line[14])
           i=i+1
           lines=lines+1
        lines=lines+1
    ##################################################Questo blocco ha generato una lista 1D di lunghezza nmol*nframe
    ##################################################Ciascun elemento è il valore della TC dei nodi in ordine
    TC_frame=[]
    lmn=0
    with open(os.path.join(path, 'nodes_patches_first_shell_TC_threshold.txt'), 'a+') as fobj:
       fobj.write('The id of the molecules has already been increased by 1, so the values range from 1 to 710. \n')
       fobj.write('\n')    
    while kk<n_frames:
        print(kk)
        TC_frame=[]
        j=0
        while j<nmol:
            TC_frame.append(TC_values[lmn])
            lmn=lmn+1
            j=j+1
        with open(os.path.join(path, 'nodes_patches_first_shell_TC_threshold.txt'), 'a+') as fobj:
            fobj.write('FRAME {}: \n'.format(kk))
        _, patches_m_larger_25, other_patches=nodes_patches_TC_constraint(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, res_threshold,TC_frame)       
        
        n_patches_larger_25 = len(patches_m_larger_25)
        
        for i, patch in enumerate(patches_m_larger_25):
        
            with open(os.path.join(path, 'nodes_patches_first_shell_TC_threshold.txt'), 'a+') as fobj:
                fobj.write('number of molecules in the patch {} = {} \n'.format(i, len(patch)))
                fobj.write('id molecules in the patch {} = {} \n'.format(i, [hh+1 for hh in patch]))
                
            A = get_graph_A(U, kk, list_molecules, list_nodes, boundary, dist) 
            ls = []
            for bb in patch:                            
                B = A[bb,:]; 
                a = csr_matrix.nonzero(B)
                ls = ls + list(a[1])                
            first_shell = list(set(ls) - set(patch))
             
            with open(os.path.join(path, 'nodes_patches_first_shell_TC_threshold.txt'), 'a+') as fobj:
                fobj.write('number of molecules in the first shell of patch {} = {} \n'.format(i, len(first_shell)))
                fobj.write('id molecules in the first shell of patch {} = {} \n'.format(i, [hh+1 for hh in first_shell]))
                
        for ii, patch2 in enumerate(other_patches):
        
            with open(os.path.join(path, 'nodes_patches_first_shell_TC_threshold.txt'), 'a+') as fobj:
                fobj.write('number of molecules in the patch {} = {} \n'.format(ii+n_patches_larger_25, len(patch2)))
                fobj.write('id molecules in the patch {} = {} \n'.format(ii+n_patches_larger_25, [hh+1 for hh in patch2]))
        
        with open(os.path.join(path, 'nodes_patches_first_shell.txt'), 'a+') as fobj:
             fobj.write('\n')    
        kk=kk+1
        
            

    return       
       
def res_8_info_patches_second_shell(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold):
    
    with open(os.path.join(path, 'nodes_patches_first_second_shell.txt'), 'a+') as fobj:
        fobj.write('The id of the molecules has already been increased by 1, so the values range from 1 to 710. \n')
        fobj.write('\n')
    
    
    for qqq in range(n_frames):
        print(qqq)
        
        with open(os.path.join(path, 'nodes_patches_first_and_second_shell.txt'), 'a+') as fobj:
            fobj.write('FRAME {}: \n'.format(qqq))
    
        _, patches_m_larger_5, other_patches, _ = nodes_patches(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
        
        n_patches_larger_5 = len(patches_m_larger_5)
        
        for i, patch in enumerate(patches_m_larger_5):

        
            with open(os.path.join(path, 'nodes_patches_first_and_second_shell.txt'), 'a+') as fobj:
                fobj.write('number of molecules in the patch {} = {} \n'.format(i, len(patch)))
                fobj.write('id molecules in the patch {} = {} \n'.format(i, [hh+1 for hh in patch]))
                
            A = get_graph_A(U, qqq, list_molecules, list_nodes, boundary, dist) 
            ls = []
            for bb in patch:                            
                B = A[bb,:]; 
                a = csr_matrix.nonzero(B)
                ls = ls + list(a[1])                
            first_shell = list(set(ls) - set(patch))
            ls2 = []
            for bbb in first_shell:
                B = A[bbb,:]
                a = csr_matrix.nonzero(B)
                ls2 = ls2 + list(a[1])
            second_shell =list(set(ls2) - set (patch) -set (ls))
             
            with open(os.path.join(path, 'nodes_patches_first_and_second_shell.txt'), 'a+') as fobj:
                fobj.write('number of molecules in the first shell of patch {} = {} \n'.format(i, len(first_shell)))
                fobj.write('id molecules in the first shell of patch {} = {} \n'.format(i, [hh+1 for hh in first_shell]))
                fobj.write('number of molecules in the second shell of patch {} = {} \n'.format(i, len(second_shell)))
                fobj.write('id molecules in the second shell of patch {} = {} \n'.format(i, [hhh+1 for hhh in second_shell]))
                
        for ii, patch2 in enumerate(other_patches):
        
            with open(os.path.join(path, 'nodes_patches_first_and_second_shell.txt'), 'a+') as fobj:
                fobj.write('number of molecules in the patch {} = {} \n'.format(ii+n_patches_larger_5, len(patch2)))
                fobj.write('id molecules in the patch {} = {} \n'.format(ii+n_patches_larger_5, [hh+1 for hh in patch2]))
        
        with open(os.path.join(path, 'nodes_patches_first_and_second_shell.txt'), 'a+') as fobj:
             fobj.write('\n')    
    
        
            

    return
    
          
       
       
def res_8_follow_node_TC_from_prompt (U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, res_threshold):

    print('Insert file with TC values')
    root = Tk()
    root.withdraw()
    nfile_TC = filedialog.askopenfilename()
    root.destroy()
    file_TC=open(nfile_TC,'r').read()
    ##################################################    	
    TC_all=file_TC.splitlines()
    lines=2
    TC_values=[]
    nmol=710
    i=0
    TC_line=''
    while lines<len(TC_all):
        i=0
        while i<nmol:
           TC_line=TC_all[lines].split()
           TC_values.append(TC_line[14])
           i=i+1
           lines=lines+1
        lines=lines+1
    ##################################################Questo blocco ha generato una lista 1D di lunghezza nmol*nframe
    ##################################################Ciascun elemento è il valore della TC dei nodi in ordine
    print ('Insert nodes indexes separated by a blank space')
    str_nodes=input()
    split_str_nodes=str_nodes.split()
    node_list=[]
    i=0
    while i <len(split_str_nodes):
       node_list.append(int(split_str_nodes[i]))
       i=i+1
    fobj=open('TC_nodes_from_prompt_.txt','w')
    fobj.write('#Each NODE is separated by a space')
    fobj.write('\n')
    i=0
    while i <len(node_list):
       j=(node_list[i]-1)
       while j<len(TC_values):
          fobj.write(str(TC_values[j])+'\n')
          j=j+nmol
       i=i+1
       fobj.write('\n')     
    return

def res_8_gen_vmd_input(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold):
    thresh_dim_cluster=25
    nmol=710
    output=open('user.dat','w')
    line=""
    qqq=0
    print("do you wish to write TC values HDL? y/n")
    TCyn=input()
    while qqq <(n_frames):
        print(str(qqq+1))
    
        _, _, _, lst = nodes_patches(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
        
        _, G, _, _, _, _ = get_graph(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save)     
        dict_nodes = range_nodes_to_save(U, qqq, list_molecules, list_nodes, nodes_to_save)
        range_nodes = list(dict_nodes.keys())
        n_nodes = len(range_nodes)
        degree = G.degree
        
        list_set_nodes_WITHOUT_1 = [list(set(x)) for x in lst if len(x)>1]
        NC = len(list_set_nodes_WITHOUT_1)
        TC = []
        mol_clust_list=[]
        TC_mol_clust_list=[]
        with open(file_txt_centrality, 'r') as data:
            data.seek(0)
            for ff in range(((n_nodes)+1)*qqq +2):
                data.readline()


            for i in range(n_nodes):
                cent = data.readline()
                 
                nn_for_read = 0
                while cent[nn_for_read] != '>':
                    nn_for_read +=1

                cent = float(cent[nn_for_read+1 : nn_for_read+25])
                TC.append(cent)    
        for c in lst:
            values_tc = [TC[i] for i in c]
            index = [i for i in c]
            if len(index)>=25:
               j=0
               while j<len(index):
                  mol_clust_list.append(index[j])
                  TC_mol_clust_list.append(values_tc[j])
                  j=j+1
        arr_ord=[]
        arr_ord=mol_clust_list
        arr_ord.sort()
        ###############################################Fino a qua funziona
        if len(mol_clust_list)==0:
           j=0
           while j<nmol:
              line=line+" 0"
              j=j+1
           output.write(line+'\n')
           line=""
           qqq=qqq+1
           continue
        qwe=0
        jj=0
        mol_arr_ind=0
        while qwe<nmol:
           if jj<len(arr_ord) and qwe==int(arr_ord[jj]):
              if TCyn=='y':
                 mol_arr_ind=mol_clust_list.index(arr_ord[jj])
                 line=line+" "+str(TC_mol_clust_list[mol_arr_ind])
              else:
                 line=line+" 1"
              jj=jj+1
           else:
              line=line+" 0"      
           qwe=qwe+1
        output.write(line+'\n')
        line=""
        qqq=qqq+1
    return



def res_8_plot_cardinality(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold):
    nmol=710
    output=open('plot_card_biggest_second_biggest.xvg','w')
    line=""
    qqq=0
    arr_bigg=[]
    arr_sbigg=[]
    card=[]
    while qqq <(n_frames):
        print(str(qqq+1))
    
        _, _, _, lst = nodes_patches(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
        
        _, G, _, _, _, _ = get_graph(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save)     
        dict_nodes = range_nodes_to_save(U, qqq, list_molecules, list_nodes, nodes_to_save)
        range_nodes = list(dict_nodes.keys())
        n_nodes = len(range_nodes)
        degree = G.degree
        
        list_set_nodes_WITHOUT_1 = [list(set(x)) for x in lst if len(x)>1]
        NC = len(list_set_nodes_WITHOUT_1)
        TC = []
        mol_clust_list=[]
        TC_mol_clust_list=[]
        with open(file_txt_centrality, 'r') as data:
            data.seek(0)
            for ff in range(((n_nodes)+1)*qqq +2):
                data.readline()


            for i in range(n_nodes):
                cent = data.readline()
                 
                nn_for_read = 0
                while cent[nn_for_read] != '>':
                    nn_for_read +=1

                cent = float(cent[nn_for_read+1 : nn_for_read+25])
                TC.append(cent)    
        for c in lst:
            values_tc = [TC[i] for i in c]
            index = [i for i in c]
            card.append(len(index))    
        arr_ord=[]
        arr_ord=card
        card=[]
        arr_ord.sort()
        if len(arr_ord)==0:
            arr_bigg.append(0)
            arr_sbigg.append(0)
            qqq=qqq+1
            continue
        elif len(arr_ord)==1:
            arr_bigg.append(arr_ord[0])
            arr_sbigg.append(0)
            qqq=qqq+1
            continue
        else:
            arr_bigg.append(arr_ord[-1])
            arr_sbigg.append(arr_ord[-2])
            qqq=qqq+1
            continue
    i=0
    while i<len(arr_bigg):
        output.write(str(arr_bigg[i])+"\n")
        i=i+1
    output.write("\n")
    i=0
    while i<len(arr_sbigg):
        output.write(str(arr_sbigg[i])+"\n")
        i=i+1
    return
    
def res_8_gen_vmd_input_bigg_2bigg(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold):
       nmol=710
       output= open('userb2b.dat','w')
       line=""
       qqq=0
       while qqq <(n_frames):
          print(str(qqq+1))
          _, _, _, lst = nodes_patches(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
          _, G, _, _, _, _ = get_graph(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save)     
          dict_nodes = range_nodes_to_save(U, qqq, list_molecules, list_nodes, nodes_to_save)
          range_nodes = list(dict_nodes.keys())
          n_nodes = len(range_nodes)
          degree = G.degree
        
          list_set_nodes_WITHOUT_1 = [list(set(x)) for x in lst if len(x)>1]
          NC = len(list_set_nodes_WITHOUT_1)
          TC = []
          mol_clustb_list=[]
          TC_mol_clustb_list=[]
          mol_clust2b_list=[]
          TC_mol_clust2b_list=[]
          with open(file_txt_centrality, 'r') as data:
             data.seek(0)
             for ff in range(((n_nodes)+1)*qqq +2):
                data.readline()


             for i in range(n_nodes):
                cent = data.readline()
                 
                nn_for_read = 0
                while cent[nn_for_read] != '>':
                   nn_for_read +=1

                cent = float(cent[nn_for_read+1 : nn_for_read+25])
                TC.append(cent)    
          c=0
          clu=[]
          while c<2:
             if len(lst)==0:
                break
             if len(lst)==1:
                clu=lst[0]
                index= [i for i in clu]
                values_tc=[TC[i] for i in clu]
                j=0
                while j<len(index):
                   mol_clustb_list.append(index[j])
                   TC_mol_clustb_list.append(values_tc[j])
                   j=j+1
                break
             clu=lst[c]
             values_tc = [TC[i] for i in clu]
             index = [i for i in clu]
             if c<1:
                j=0
                while j<len(index):
                   mol_clustb_list.append(index[j])
                   TC_mol_clustb_list.append(values_tc[j])
                   j=j+1       
             else:
                j=0
                while j<len(index):
                   mol_clust2b_list.append(index[j])
                   TC_mol_clust2b_list.append(values_tc[j])
                   j=j+1
             c=c+1    
          arr_ordb=[]
          arr_ord2b=[]
          arr_ordb=mol_clustb_list
          arr_ord2b=mol_clust2b_list
          arr_ordb.sort()
          arr_ord2b.sort()
          if len(mol_clustb_list)==0:
             j=0
             while j<nmol:
                line=line+" 0"
                j=j+1
          elif len(mol_clustb_list)>0 and len(mol_clust2b_list)==0:
             qwe=0
             jj=0
             mol_arr_ind=0
             while qwe<nmol:
                if jj<len(arr_ordb) and qwe==int(arr_ordb[jj]):
                   line=line+" 1"
                   jj=jj+1
                else:
                   line=line+" 0"      
                qwe=qwe+1             
          else:
             qwe=0
             jb=0
             j2b=0
             mol_arr_indb=0
             mol_arr_ind2b=0
             while qwe<nmol:
                if jb<len(arr_ordb) and qwe==int(arr_ordb[jb]):
                   line=line+" 1"
                   jb=jb+1
                elif j2b<len(arr_ord2b) and qwe==int(arr_ord2b[j2b]):
                   line=line+" 2"
                   j2b=j2b+1
                else:
                   line=line+" 0"      
                qwe=qwe+1                   
          output.write(line+'\n')
          line=""
          qqq=qqq+1
       return

def res_8_HDLclu_vs_tot(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, file_txt_centrality, res_threshold):

    thresh_dim_cluster=25
    nmol=710
    output=open('ratioHDLclu_HDLtot.txt','w')
    weights=open('weights_clusters.txt','w')
    hdl_tot=open('HDL_tot.txt','w')
    hdl_clu=open('HDL_clu.txt','w')
    fraz=open('HDLfraction.txt','w')
    line=""
    qqq=0
    print('Insert file with TC values')
    root = Tk()
    root.withdraw()
    nfile_TC = filedialog.askopenfilename()
    root.destroy()
    file_TC=open(nfile_TC,'r').read()
    ##################################################    	
    TC_all=file_TC.splitlines()
    lines=2
    TC_values=[]
    nmol=710
    ratio=0.0
    i=0
    TC_line=''
    while lines<len(TC_all):
        i=0
        while i<nmol:
           TC_line=TC_all[lines].split()
           TC_values.append(TC_line[14])
           i=i+1
           lines=lines+1
        lines=lines+1
    ##################################################Questo blocco ha generato una lista 1D di lunghezza nmol*nframe
    sc=0
    arr_wgt=[]
    while qqq <(n_frames):
        print(str(qqq+1))
        HDL_tot=0
        i=0
        while i<nmol:
            if float(TC_values[sc])>=93:
                HDL_tot=HDL_tot+1
            sc=sc+1
            i=i+1
        _, _, _, lst = nodes_patches(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold)
        _, G, _, _, _, _ = get_graph(U, qqq, list_molecules, list_nodes, boundary, dist, nodes_to_save)     
        dict_nodes = range_nodes_to_save(U, qqq, list_molecules, list_nodes, nodes_to_save)
        range_nodes = list(dict_nodes.keys())
        n_nodes = len(range_nodes)
        degree = G.degree
        #list_set_nodes_WITHOUT_1 = [list(set(x)) for x in lst if len(x)>1]
        #NC = len(list_set_nodes_WITHOUT_1)
        dim_clu=[]   
        HDL_clu=0 
        wgt=[]
        hdl_tot.write(str(HDL_tot)+"\n")
        for c in lst:
            index = [i for i in c]
            if len(index)>=thresh_dim_cluster:
                dim_clu.append(len(index))
        nclu=0
        while nclu<len(dim_clu):
        	HDL_clu=HDL_clu+dim_clu[nclu]    
        	nclu=nclu+1    
        if HDL_tot==0:
           output.write(str(0)+"\n")
           fraz.write(str(0)+"\n")
           hdl_clu.write(str(0)+"\n")
           qqq=qqq+1
           wgt.append(0)
           arr_wgt.append(wgt)
           continue
        elif HDL_clu==0:
           output.write(str(0)+"\n")
           hdl_clu.write(str(0)+"\n")
           wgt.append(1)
           arr_wgt.append(wgt)
           qqq=qqq+1
           continue   
        else:
           HDL_clu=sum(dim_clu)
           hdl_clu.write(str(HDL_clu)+"\n")
           HDL_noclu=HDL_tot-HDL_clu
           wgt.append(HDL_noclu/HDL_tot)
           nclu=0
           dim_clu.sort()
           while nclu<len(dim_clu):
                wgt.append(dim_clu[nclu]/HDL_tot)
                nclu=nclu+1
           arr_wgt.append(wgt)     
           ratio=HDL_clu/HDL_tot
           output.write(str(ratio)+"\n")
           fraz.write(str(HDL_tot/nmol)+"\n")
           qqq=qqq+1
    weights.write(str(arr_wgt))
    return


def res_10_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):

    fmt = 'Do you want normalized TC? \n \
    1) no \n \
    2) yes, w.r.t. number of edges \n'
    
    normalization = int(input(fmt))
    
    fmt = 'Value of beta (beta > 0) = '
    beta = float(input(fmt))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt = 'Value of beta (beta > 0) = '
        beta = float(input(fmt))
    
    
    with open(os.path.join(path, 'directed_TC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(normalization) + '_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        fobj.write('Type of normalization = {},  molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(normalization, list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for gggg in range(n_frames):
    
        compute_TC_directed(U, gggg, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, normalization, dist)


    return
    
    
def res_10_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
                
    with open(os.path.join(path, 'directed_DEGREE'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for hhhh in range(n_frames):
    
        compute_directed_degree(U, hhhh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        
    return
    
    
def res_10_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):
                
    with open(os.path.join(path, 'directed_EIG'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for hhhh in range(n_frames):
    
        compute_directed_eigenvector(U, hhhh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        
    return

    
def res_10_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):
    
    with open(os.path.join(path, 'directed_HITS'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for hh in range(n_frames):
    
        compute_directed_HITS(U, hh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    

    return

    
def res_10_5(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):
    
    with open(os.path.join(path, 'directed_gTC'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for hh in range(n_frames):
    
        compute_directed_gTC(U, hh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    

    return
    
    
def res_10_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    fmt = 'Do you want normalized NTC? \n \
    1) no \n \
    2) yes, w.r.t. number of edges \n'
    
    normalization = int(input(fmt))
    
    fmt = 'Value of beta (beta > 0) = '
    beta = float(input(fmt))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt = 'Value of beta (beta > 0) = '
        beta = float(input(fmt))
    
    
    with open(os.path.join(path, 'directed_bTC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(normalization) + '_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        fobj.write('Type of normalization = {},  molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(normalization, list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for gggg in range(n_frames):
    
        compute_bTC_directed(U, gggg, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, normalization, dist)
        
    return
    
    
