# Chiara Faccio, Scuola Normale Superiore, Pisa
# chiara.faccio@sns.it
# March 2023

import numpy as np
import numpy.linalg
from numpy.linalg import norm
from numpy import ones, zeros, arange, inf

import random
from collections import deque

import networkx as nx
import dynetx as dnx

import sys, os, math
from math import exp, log

from tkinter import Tk, filedialog
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import scipy
import scipy.sparse.linalg
import scipy.io
from scipy.sparse import csr_matrix, csc_matrix, find
from scipy.spatial import distance_matrix
from scipy.sparse.linalg import norm as scipy_norm
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply, cg
from scipy.sparse import identity, diags, bmat, kron
from scipy.sparse.linalg import spsolve
from scipy.spatial.distance import pdist

import itertools
import MDAnalysis as mda



def graph_slice(G, t):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
    
    
    Author: Alberto Bucci
    
    extract a slice/snapshot of the Dynamic graph, that is a snapshot of the graph :math:`G` at time :math:`t` in \
    NetworkX format.
    Parameters
    __________
    G: DynGraph or DynDiGraph object.
    t: float; snapshot time.
    :return s: (NetworkX Graph object)
        Snapshot at time :math:`t` of :math:`G`.
    """

    """
    Examples
    --------
    >>> import dynetx as dn
    >>> import networkx as nx
    >>> G = dn.DynGraph()
    >>> G.add_interaction(1, 2, 2)
    >>> G.add_interaction(1, 2, 2, e=6)
    >>> G.add_interaction(1, 2, 7, e=11)
    >>> h = graph_slice(G, 3)
    >>> print(nx.adjacency_matrix(h))
    >>>    (0, 1)	1
    >>>    (1, 0)	1
    """

    node_list = G.nodes()
    slice_t = list(dnx.interactions(G, t=t))
    edge_list = ([e[0], e[1]] for e in slice_t)
    sliced_graph = nx.Graph()
    sliced_graph.add_nodes_from(node_list)
    sliced_graph.add_edges_from(edge_list)
    return sliced_graph


def broadcast_centrality(G, alpha=None, conj_grad_maxiter=100, conj_grad_tol=1e-7):
    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
    
    
    Author: Alberto Bucci
    
    Computes the broadcast centrality of the dynamic graph :math:`G`.
    Denoting with :math:`A_t` the adjacency matrix at time :math:`t` and with :math:`\\mathbf{1}` the vector of all ones
    broadcast centrality is :math:`bc = (I-\\alpha A_1)^{-1}(I-\\alpha A_2)^{-1}...(I-\\alpha A_k)^{-1} \\mathbf{1}`\
     [1]_.
    Recall that to ensure that each resolvent :math:`(I-\\alpha A_t)` can be expressed as a power series in the matrix,
    the parameter :math:`\\alpha` must satisfy\
     :math:`0<\\alpha< \\frac{1}{\\rho^*}`, where :math:`\\rho^* = \\max_t \\rho(A_t)`.
    Parameters
    __________
    G: DynGraph object
        a dynamic graph.
    alpha: float, optional
        parameter, if None computed to ensure that each resolvent can be expressed as a power series in the matrix, default: None.
    conj_grad_maxiter: integer
        Maximum number of iterations for solving :math:`Ax=b` with conjugate gradient method. \
        Iterations will stop after maxiter steps even if the specified tolerance has not been achieved, default: 100.
    conj_grad_tol: float, optional
        Tolerance for solving :math:`Ax=b` with conjugate gradient method.
         ``norm(residual) <= conj_grad_tol*norm(b)`` , default: 1e-7.
    Returns
    _______
    bc: dict
        broadcast centrality.
    alpha: float
        alpha parameter.
    Examples
    ________
    .. code:: python
     >>>  from networksns import centrality_measures as cm
     >>>  import dynetx as dn
    Create dynamic graph :math:`G`
    .. code:: python
     >>>    G = dn.DynGraph()
     >>>    G.add_interaction(1, 2, 2, 5)
     >>>    G.add_interaction(1, 3, 2, 5)
     >>>    G.add_interaction(2, 3, 4)
             EdgeView([(1, 2), (1, 3), (2, 3)])
    Compute broadcast centrality
     >>>   bc, alpha = cm.broadcast_centrality(G)
    References
    ----------
    .. [1] Michele Benzi, Isabel Chen, Howard H. Chang, Vicki S. Hertzberg (2017),
           Dynamic communicability and epidemic spread: a case study on an empirical dynamic contact network,
           Journal of Complex Networks, Volume 5, Issue 2,
           https://doi.org/10.1093/comnet/cnw017
    """

    time_snapshots = G.temporal_snapshots_ids()
    n = G.number_of_nodes()
    e = ones(n)  # vector of all ones
    bc = e
    if alpha is None:
        spectral_bound = 0
        for t in time_snapshots[::-1]:
            G_t = graph_slice(G, t)
            Adj_t = nx.adjacency_matrix(G_t)
            row_sum = scipy_norm(Adj_t, ord=inf)
            spectral_bound = max(spectral_bound, row_sum)
        if spectral_bound != 0:
            alpha = 0.9/spectral_bound  # 0 < alpha < 0.9/max(rho(Adj_t))
    for t in time_snapshots[::-1]:
        G_t = graph_slice(G, t)
        Adj_t = nx.adjacency_matrix(G_t)
        bc = cg(identity(n) - alpha * Adj_t, bc, maxiter=conj_grad_maxiter, tol=conj_grad_tol, atol=0)
        if bc[1] != 0:
            raise ValueError('convergence not achieved')
        bc = bc[0]
        if any(bc != 0):
            bc = bc / norm(bc)
            
    G_0 = graph_slice(G, time_snapshots[0])
    node_list = G.nodes()
    bc = dict(zip(node_list, bc))
    return bc, alpha


def receive_centrality(G, alpha=None, conj_grad_maxiter=100, conj_grad_tol=1e-7):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
    
    
    Author: Alberto Bucci
    
    Computes the receive centrality of the dynamic graph :math:`G`.
    Denoting with :math:`A_t` the adjacency matrix at time :math:`t` and with :math:`\\mathbf{1}` the vector of all ones
    receive centrality is :math:`rc = (I-\\alpha A_k)^{-1}(I-\\alpha A_{k-1})^{-1}...(I-\\alpha A_1)^{-1} \\mathbf{1}`\
     [1]_.
    Recall that to ensure that each resolvent :math:`(I-\\alpha A_t)` can be expressed as a power series in the matrix,
    the parameter :math:`\\alpha` must satisfy \
    :math:`0<\\alpha< \\frac{1}{\\rho^*}`, where :math:`\\rho^* = \\max_t \\rho(A_t)`.
    Parameters
    __________
    G: DynGraph object
        a dynamic graph.
    alpha: float, optional
        parameter, if None computed to ensure that each resolvent can be expressed as a power series in the matrix. Default: None.
    conj_grad_maxiter: integer
        Maximum number of iterations for solving :math:`Ax=b` with conjugate gradient method.\
         Iterations will stop after maxiter steps even if the specified tolerance has not been achieved, default: 100.
    conj_grad_tol: float, optional
        Tolerance for solving :math:`Ax=b` with conjugate gradient method.
         ``norm(residual) <= conj_grad_tol*norm(b)`` , default: 1e-7.
    Returns
    _______
    rc: dict
        receive centrality.
    alpha: float
        alpha parameter.
    Examples
    ________
    .. code:: python
     >>>  from networksns import centrality_measures as cm
     >>>  import dynetx as dn
    Create dynamic graph :math:`G`.
    .. code:: python
     >>>    G = dn.DynGraph()
     >>>    G.add_interaction(1, 2, 2, 5)
     >>>    G.add_interaction(1, 3, 2, 5)
     >>>    G.add_interaction(2, 3, 4)
             EdgeView([(1, 2), (1, 3), (2, 3)])
    Compute receive centrality
     >>>   rc, alpha = cm.receive_centrality(G)
    References
    ----------
    .. [1] Michele Benzi, Isabel Chen, Howard H. Chang, Vicki S. Hertzberg (2017),
           Dynamic communicability and epidemic spread: a case study on an empirical dynamic contact network,
           Journal of Complex Networks, Volume 5, Issue 2,
           https://doi.org/10.1093/comnet/cnw017
    """

    time_snapshots = G.temporal_snapshots_ids()
    n = G.number_of_nodes()
    e = ones(n)  # vector of all ones
    rc = e
    if alpha is None:
        spectral_bound = 0
        for t in time_snapshots:
            G_t = graph_slice(G, t)
            Adj_t = nx.adjacency_matrix(G_t)
            row_sum = scipy_norm(Adj_t, ord=inf)
            spectral_bound = max(spectral_bound, row_sum)
        if spectral_bound != 0:
            alpha = 0.9/spectral_bound  # 0 < alpha < 0.9/max(rho(Adj_t))
    for t in time_snapshots:
        G_t = graph_slice(G, t)
        Adj_t = nx.adjacency_matrix(G_t)
        rc = cg(identity(n) - alpha * Adj_t, rc, maxiter=conj_grad_maxiter, tol=conj_grad_tol, atol=0)
        if rc[1] != 0:
            raise ValueError('convergence not achieved')
        rc = rc[0]
        if any(rc != 0):
            rc = rc / norm(rc)
    G_0 = graph_slice(G, time_snapshots[0])
    node_list = list(G_0.nodes)
    rc = dict(zip(node_list, rc))
    return rc, alpha


def exponential_symmetric_quadrature(A, u, tol=1e-7, maxit=50):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    
    Author : Alberto Bucci
    
    Computes :math:`q=u^Te^Au`.
    The computation is done by means of Lanczos method according to [1]_.
    Parameters
    __________
    A: array_like
        sparse/dense symmetric matrix.
    u: array
        vector.
    tol: float,optional
        tolerance for convergence, relative accuracy, default: 1e-7.
    maxit: integer, optional
     maximum number of Lanczos iterations, default: 50.
    :return: **q**: (float)
     value of the quadratic form :math:`u^Te^Au`.
    Examples
    ________
    .. code:: python
     >>>  from networksns import centrality_measures as cm
     >>>  import numpy as np
    Create symmetric matrix :math:`A`
    .. code:: python
     >>>    A = np.arange(0, 9, 1)
     >>>    A = A.reshape(3, 3)
     >>>    A = A + A.transpose()
            array([[ 0,  4,  8],
                   [ 4,  8, 12],
                   [ 8, 12, 16]])
    Create vector :math:`u`
        .. code:: python
    >>>     u = np.arange(0, 3)
            array([0, 1, 2])
    Compute :math:`q=u^T e^A u`.
     >>>    q = cm.exponential_symmetric_quadrature(A, u)
    References
    ----------
    .. [1] G. H. Golub, and G. Meurant (2010)
           "Matrices, Moments and Quadrature with Applications",
           Princeton University Press, Princeton, NJ.
    """
    quadrature = 1
    quadrature_old = 2
    old_err = 1
    err = 1
    omega = []
    gamma = []
    u_norm = norm(u)
    if u_norm == 0:
        return 0
    else:
        x_0 = u/u_norm
        #  computing Lanczos matrix J
        omega.append(x_0.dot(A.dot(x_0)))
        r = A.dot(x_0) - omega[0] * x_0
        r_norm = norm(r)
        if r_norm == 0:
            return exp(omega[0]) * u_norm ** 2
        gamma.append(r_norm)
        x_1 = r / r_norm
        it = 1  # iterations
        while err > tol and old_err > tol and it < maxit:
            z = A.dot(x_1) - gamma[it - 1] * x_0
            omega.append(x_1.dot(z))
            x_2 = z - omega[it] * x_1
            if norm(x_2) == 0:
                gamma.append(0)  # variable used only to avoid deletion in line: eta = gamma[:-1]
                eta = gamma[:-1]
                J = diags(omega)
                J = J + diags(eta, 1)
                J = J + diags(eta, -1)
                e_1 = np.zeros(len(omega))
                e_1[0] = 1
                quadrature = e_1.dot(expm_multiply(J, e_1)) * (norm(u)) ** 2
                break
            gamma.append(norm(x_2))
            x_0 = x_1
            x_1 = x_2 / gamma[it]
            eta = gamma[:-1]
            J = diags(omega)
            J = J + diags(eta, 1)
            J = J + diags(eta, -1)
            e_1 = np.zeros(len(omega))
            e_1[0] = 1
            quadrature_very_old = quadrature_old
            quadrature_old = quadrature
            quadrature = e_1.dot(expm_multiply(J, e_1))*(norm(u))**2
            old_err = err
            err = max(abs((quadrature_old-quadrature))/quadrature_old,
                      abs((quadrature_very_old-quadrature_old))/quadrature_very_old)
            it = it+1
    return quadrature
    
    
def subgraph_centrality_Bucci(G, t=1, tol=1e-7, maxit=50):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    
    Author : Alberto Bucci
    
    Computes the subgraph centrality of all the nodes in graph :math:`G`.
    The subgraph centrality of the :math:`i^{th}` node is given by :math:`[e^{tA}]_{ii}=e_i^T (e^{tA})e_i`,
    where :math:`e_i` and :math:`A` denote respectively the :math:`i^{th}` vector of the canonical basis and the adjacency matrix of the graph [1]_.
    Parameters
    __________
    G: Graph or DiGraph object
        a graph.
    t: scalar, optional
     when exponentiating multiply the adjacency matrix by t, default: 1.
    tol: float,optional
     tolerance for convergence, relative accuracy, default: 1e-7.
    maxit: integer, optional
     maximum number of Lanczos iterations, default: 50.
    :return: **sc** (dict)  subgraph centrality of all nodes in :math:`G`.
    Examples
    ________
    .. code:: python
     >>>  from networksns import centrality_measures as cm
     >>>  import networkx as nx
    Create graph :math:`G`.
    .. code:: python
     >>>    G = nx.Graph()
     >>>    G.add_edge(1, 2)
     >>>    G.add_edge(2, 3)
            EdgeView([(1, 2), (2, 3)])
    Compute the subgraph centrality.
     >>>    sc = cm.subgraph_centrality(G)
    References
    ----------
    .. [1] Ernesto Estrada and Juan A. Rodríguez-Velázquez (2005)
           Subgraph centrality in complex networks,
           Physical Review, Volume 71, Issue 5,
           https://doi.org/10.1103/PhysRevE.71.056103
    """
    
    n = G.number_of_nodes()
    node_list = list(G.nodes)
    subgraph_centrality = np.zeros(n)
    for i in range(n):
        subgraph_centrality[i] = node_subgraph_centrality(G, node_list[i], t, tol, maxit)
    centrality = dict(zip(node_list, subgraph_centrality))
    return centrality
    

def node_subgraph_centrality(G, u, t=1, tol=1e-7, maxit=50):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    
    Author : Alberto Bucci
    
    Computes the subgraph centrality of node :math:`u`.
    If node :math:`u` is the :math:`i^{th}` node of the graph, the subgraph centrality of node :math:`u` is given by :math:`[e^{tA}]_{ii}=e_i^T (e^{tA})e_i`,
    where :math:`e_i` and :math:`A` denote respectively the :math:`i^{th}` vector of the canonical basis and the adjacency matrix of the graph [1]_.
    Parameters
    __________
    G: Graph object
        a graph.
    u: node_id
        node in G.
    t: scalar, optional
     when exponentiating multiply the adjacency matrix by t, default: 1.
    tol: float,optional
     tolerance for convergence, relative accuracy, default: 1e-7.
    maxit: integer, optional
     maximum number of Lanczos iterations, default: 50.
     :return: **sc_u** (float) subgraph centrality of node :math:`u`.
    Examples
    ________
    .. code:: python
     >>>  from networksns import centrality_measures as cm
     >>>  import networkx as nx
    Create graph :math:`G`.
    .. code:: python
     >>>    G = nx.Graph()
     >>>    G.add_edge(1, 'u')
     >>>    G.add_edge('u', 2)
            EdgeView([(1, 'u'), ('u', 2)])
    Compute the node total communicability
     >>>    sc_u = cm.node_subgraph_centrality(G, 'u')
    References
    ----------
    .. [1] Ernesto Estrada and Juan A. Rodríguez-Velázquez (2005)
           Subgraph centrality in complex networks,
           Physical Review, Volume 71, Issue 5,
           https://doi.org/10.1103/PhysRevE.71.056103
    """
    
    n = G.number_of_nodes()
    node_list = G.nodes
    enumerated_nodes = dict(zip(node_list, np.arange(n)))
    node_position = enumerated_nodes[u]
    e_node = np.zeros(n)
    e_node[node_position] = 1
    Adj = nx.adjacency_matrix(G)
    if t != 1:
        Adj = Adj * t
    subgraph_centrality = exponential_symmetric_quadrature(Adj, e_node, tol, maxit)
    return subgraph_centrality


def total_communicability(G, t=1):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    
    Author : Alberto Bucci
    
    Computes the total communicability of :math:`G`.

    Total communicability is defined as the row sum of the exponential of the adjacency matrix [1]_, so denoting with
    :math:`A` the adjacency matrix of graph :math:`G` and with :math:`\\mathbf{1}` the vector of all ones,
    we have :math:`tc= e^A \\mathbf{1}`.

    Parameters
    __________
    G: Graph or DiGraph object
        a graph.
    t: scalar, optional
        exponentiate multiply the adjacency matrix by :math:`t`, default None.


    :return: **tc** (dict)  total communicability of graph :math:`G`.

    Examples
    ________
    .. code:: python

     >>>  from sobigdatainit.sobigdata_pkg import sobigdata as sbd
     >>>  import networkx as nx

    Create graph :math:`G`

    .. code:: python

     >>>    G = nx.Graph()
     >>>    G.add_edge(1, 2)
     >>>    G.add_edge(2, 3)
            EdgeView([(1, 2), (2, 3)])

    Compute :math:`tc= e^A \\mathbf{1}`

     >>>    tc = sbd.total_communicability(G)


    References
    ----------
    .. [1] Michele Benzi, Christine Klymko (2013)
           Total communicability as a centrality measure,
           Journal of Complex Networks, Volume 1, Issue 2, Pages 124–149,
           https://doi.org/10.1093/comnet/cnt007
    """

    n = G.number_of_nodes()
    node_list = G.nodes
    e = ones(n)  # vector of all ones
    Adj = nx.adjacency_matrix(G)
    if t != 1:
        Adj = Adj*t
    tot_communicability = expm_multiply(Adj, e)
    centrality = dict(zip(node_list, tot_communicability))
    return centrality


def total_directed_communicability(G, t=None):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
    
    
    Author: Alberto Bucci
    
    Computes broadcast and receive communicability of a directed graph :math:`G`.

    Let :math:`A` be the adjacency matrix of :math:`G`, then\
     :math:`\\mathcal{A}=\\begin{pmatrix} 0 & A \\\\ A^T & 0 \\end{pmatrix}`
    is the adjacency matrix of the associated undirected bipartite graph [1]_.

    Broadcast communicability  (hub centrality) is given by :math:`bc = (e^{\\mathcal{A}}\\mathbf{1})_{1:n}`.


    Receive communicability (authority centrality) is given by :math:`rc = (e^{\\mathcal{A}}\\mathbf{1})_{n+1:2n}`.

    Parameters
    __________

    G: DiGraph object
        a directed graph.
    t: scalar, optional
     when exponentiate multiply the adjacency matrix by t, default None.

    Returns
    ________

    bc: dict
     broadcast communicability
    rc: dict
     receive communicability

    Examples
    ________

    .. code:: python

     >>>  from sobigdatainit.sobigdata_pkg import sobigdata as sbd
     >>>  import networkx as nx

    Create graph :math:`G`

    .. code:: python

     >>>    G = nx.DiGraph()
     >>>    G.add_edge(1, 2)
     >>>    G.add_edge(1, 3)
     >>>    G.add_edge(2, 3)
     >>>    G.add_edge(3, 1)
            OutEdgeView([(1, 2), (1, 3), (2, 3), (3, 1)])

    Compute broadcast and receive communicability

     >>>    bc, rc = sbd.total_directed_communicability(G)

    References
    ----------
    .. [1] Michele Benzi, Paola Boito (2020),
           Matrix functions in network analysis,
           GAMM‐Mitteilungen, Volume 43, Issue 3,
           https://onlinelibrary.wiley.com/doi/abs/10.1002/gamm.202000012
    """

    n = G.number_of_nodes()
    node_list = G.nodes
    e = ones(2*n)  # vector of all ones
    Adj = nx.adjacency_matrix(G, nodelist = range(0,n))
    Bip_Adj = bmat([[None, Adj], [Adj.transpose(), None]])
    if t is not None:
        Bip_Adj = Bip_Adj*t
    tot_communicability = expm_multiply(Bip_Adj, e)
    bc = dict(list(zip(node_list, tot_communicability[0:n])))  # hub centrality
    rc = dict(list(zip(node_list, tot_communicability[n:2*n])))  # authority centrality
    return bc, rc


def total_network_communicability(G, t=None, tol=1e-7, maxit=50):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
    
    Author : Alberto Bucci
    
    Computes the total network communicability of :math:`G`.

    Total network communicability is defined as the sum over all elements of the exponential of the adjacency matrix\
     [1]_, so denoting with
    :math:`A` the adjacency matrix of graph :math:`G` and with :math:`\\mathbf{1}` the vector of all ones,
    we have :math:`tnc = \\mathbf{1}^T e^A \\mathbf{1}`.

    Parameters
    __________
    G: Graph object
        an undirected graph.
    t: scalar, optional
        exponentiate multiply the adjacency matrix by t, default None.
    tol: float,optional
     tolerance for convergence, relative accuracy, default: 1e^7.
    maxit: integer, optional
     maximum number of Lanczos iterations, default: 50


    :return: **tnc** (float)  total network communicability of graph :math:`G`.

    Examples
    ________
    .. code:: python

     >>>  from sobigdatainit.sobigdata_pkg import sobigdata as sbd
     >>>  import networkx as nx

    Create graph :math:`G`

    .. code:: python

     >>>    G = nx.Graph()
     >>>    G.add_edge(1, 'u')
     >>>    G.add_edge('u', 2)
            EdgeView([(1, 'u'), ('u', 2)])

    Compute the total network communicability

     >>>    tnc = sbd.total_network_communicability(G)

    References
    ----------
    .. [1] Michele Benzi, Christine Klymko (2013)
           Total communicability as a centrality measure,
           Journal of Complex Networks, Volume 1, Issue 2, Pages 124–149,
           https://doi.org/10.1093/comnet/cnt007
    """

    n = G.number_of_nodes()
    Adj = nx.adjacency_matrix(G, nodelist = range(0,n))
    if t is not None:
        Adj = Adj*t
    e = ones(n)
    tot_net_communicability = exponential_symmetric_quadrature(Adj, e, tol, maxit)
    return tot_net_communicability


def select_topology_file():

    '''
    Selects the tolopogy file
    '''
    
    print("Choose the topology file.")
    root = Tk()
    root.withdraw()
    fname = filedialog.askopenfilename()
    root.destroy()

    name_file = os.path.splitext(os.path.basename(fname))[0]
    extension = os.path.splitext(os.path.basename(fname))[-1]
    print("You have chosen the file {}{} ".format(name_file, extension))
    
    path = os.path.split(fname)[0]
        
    return path, fname

    
def select_trajectory_file():

    '''
    Selects the trajectory file
    '''
    
    print("Choose the trajectory file.")
    root = Tk()
    root.withdraw()
    fname = filedialog.askopenfilename()
    root.destroy()

    name_file = os.path.splitext(os.path.basename(fname))[0]
    extension = os.path.splitext(os.path.basename(fname))[-1]
    print("You have chosen the file {}{} ".format(name_file, extension))
        
    return fname


def parameters_graph_undirected(topol):

    '''
    Returns the parameters that allow you to build the undirected graph: 
    1) dist = the threshold on the distance between two nodes (nm),
    2) boundary = indicates whether or not there are periodic boundary conditions,
    3) list_molecules = which molecules consider,
    4) list_nodes = for each molecule which atom to consider as a node,
    5) nodes_to_save = the name of the molecules for which you would like to save the values of centrality measures of its nodes
    '''

    # We choose the distance
    fmt_dist = 'Threshold to define the edge (nm): '
    #dist = input(fmt_dist)
    dist = 0.35 #float(dist)
    
    while dist < 0.0:
        print('ERROR: the distance must be positive')
        dist = input(fmt_dist)
        dist = float(dist)
        
    # We impose the periodic boundary conditions
    fmt_bound = 'Do you want to impose periodic boundary conditions: (1) no, (2) yes: '
    boundary = '2' #input(fmt_bound)
    
    while (boundary != str(1)) and (boundary != str(2)) :
        print('ERROR: not valid value') 
        boundary = input(fmt_bound)
    
    # We choose the nodes
    ML = 'tmp'
    list_molecules =['water']#[]
    #fmt_ml = 'Enter the names of the MOLECULE to be considered (press ENTER to end the insertion): '
    #while ML:
    #    ML = input(fmt_ml)
    #    if ML:
    #        list_molecules.append(ML)
    
    U = mda.Universe(topol)
    list_nodes = [['OW1']]
#    for ii in range(len(list_molecules)):
#        molecules = U.select_atoms('resname '+str(list_molecules[ii]))
#        name_atoms = list(set(molecules.names))
#        fmt_nodes = 'For the {} molecule, whom do you want to consider as a node? Possibility: {}, enter the name of the ATOM without quotation marks (press ENTER to end the insertion) #'.format(list_molecules[ii], name_atoms)
#        AT = 'tmp'
    list_atoms =['OW1']#[]
    #    while AT:
    #        AT = input(fmt_nodes)
    #        if AT:
    #            list_atoms.append(AT)
    #    list_nodes.append(list_atoms)
    #nd = 'tmp'
    nodes_to_save = ['water']#[]
     #fmt_cm = 'Please enter the name of the MOLECULE for which you would like to save the values of centrality measures of its nodes (press ENTER to end the insertion): '
    #while nd:
    #    nd = input(fmt_cm)
    #    if nd:
    #        nodes_to_save.append(nd)

    return {'dist':dist, 'boundary':boundary, 'list_molecules':list_molecules, 'list_nodes':list_nodes, 'nodes_to_save':nodes_to_save}


def parameters_graph_directed(topol):

    '''
    Returns the parameters that allow you to build the undirected graph: 
    1) dist = the threshold on the distance between two nodes (nm),
    2) boundary = indicates whether or not there are periodic boundary conditions,
    3) list_molecules = which molecules consider,
    4) list_nodes = for each molecule which atom to consider as a node,
    5) nodes_to_save = the name of the molecule for which you would like to save the values of centrality measures of its nodes
    '''
    
    # We choose the distance
    fmt_dist = 'Threshold to define the edge (nm): '
    dist = input(fmt_dist)
    dist = float(dist)
    
    while dist < 0.0:
        print('ERROR: the distance must be positive')
        dist = input(fmt_dist)
        dist = float(dist)

    print('Two nodes are connected if their distance is less than {} and if the angle is between + - 30 degrees\n'.format(dist))
    
    # We impose the periodic boundary conditions
    fmt_bound = 'Do you want to impose periodic boundary conditions: (1) no, (2) yes: '
    boundary = input(fmt_bound)
    
    while (boundary != str(1)) and (boundary != str(2)) :
        print('ERROR: invalid value, choose 1 or 2') 
        boundary = input(fmt_bound)
    
    # We choose the nodes
    ML = 'tmp'
    list_molecules = []
    fmt_ml = 'Enter the name of the MOLECULE to be considered (press ENTER to end the insertion): '
    while ML:
        ML = input(fmt_ml)
        if ML:
            list_molecules.append(ML)
    
    U = mda.Universe(topol)
    list_nodes = []
    for ii in range(len(list_molecules)):
        molecules = U.select_atoms('resname '+str(list_molecules[ii]))
        name_atoms = list(set(molecules.names))
        fmt_nodes = 'For the {} molecule, whom do you want to consider as a node? Possibility: {}, enter the name of the ATOM without quotation marks (press ENTER to end the insertion) '.format(list_molecules[ii], name_atoms)
        AT = 'tmp'
        list_atoms = []
        while AT:
            AT = input(fmt_nodes)
            if AT:
                list_atoms.append(AT)
        list_nodes.append(list_atoms)
                
        
    nd = 'tmp'
    nodes_to_save = []
    fmt_cm = 'Please enter the name of the MOLECULE for which you would like to save the values of centrality measures of its nodes (press ENTER to end the insertion): '
    while nd:
        nd = input(fmt_cm)
        if nd:
            nodes_to_save.append(nd)

    return {'dist':dist, 'boundary':boundary, 'list_molecules':list_molecules, 'list_nodes':list_nodes, 'nodes_to_save':nodes_to_save}


def extract_data(U, kk, list_molecules, list_nodes):
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
 
 
def range_nodes_to_save(U, kk, list_molecules, list_nodes, nodes_to_save):

    dict_nodes = {}
    
    for tr in U.trajectory[kk:kk+1]:

        hh = 0
        for ii in range(len(list_molecules)):
            ml = U.select_atoms('resname '+str(list_molecules[ii]))
            for jj in range(len(list_nodes[ii])):
                at = ml.select_atoms('name '+str(list_nodes[ii][jj]))
                len_at = len(at)
                for ll in range(len_at):
                    if str(list_molecules[ii]) in nodes_to_save:
                        dict_nodes[hh] = at[ll]
                    hh += 1
            
    return dict_nodes


def get_adjacency_matrix(dist, coord):

    '''
    Returns the adjacency matrix without periodic boundary conditions  
    '''

    coord = np.asmatrix(coord)    
    D = distance_matrix(coord, coord)
    
    D1 = D <= dist #boolean matrix whose entry {i,j} is True if the distance between node v_i and node v_j is less or equal than dist
    
    A = D1*1
    A = A - scipy.sparse.diags(A.diagonal())
    
    return csr_matrix(A)
    
    
def get_adjacency_matrix_pbc(dist, coord, box):

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


def get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist):

    '''
    Construts the adjacency matrix, but  return only the graph and the number of edges
    '''
    
    data = extract_data(U, kk, list_molecules, list_nodes)
    coord = data['coord']
    box = data['box']
    
    if int(boundary) == 1:
        A = get_adjacency_matrix(dist, coord)
    elif int(boundary) == 2:
        A = get_adjacency_matrix_pbc(dist, coord, box)
    else:
        print('Value not valid')
        sys.exit(1)
    
    G = nx.from_scipy_sparse_matrix(A)
    
    n_edges = G.number_of_edges()
        
    return G, n_edges
    

def angle_between(H, O1, O2):

    '''
    Returns the angle in radians between the points H O1 O2
    '''
    
    O1H =  H - O1
    O1O2 = O2 - O1
    V1 = O1H/np.linalg.norm(O1H)
    V2 = O1O2/np.linalg.norm(O1O2)
    res = np.clip(np.dot(V1, V2), -1.0, 1.0)
    angle = np.arccos(res)
    
    return angle
 

def oriented_adjacency_matrix(dist, coord_O, coord_H1, coord_H2, box, boundary):

    '''
    Returns the sparse directed adjacency matrix
    '''
    
    if int(boundary) == 1:
        A_dist = get_adjacency_matrix(dist, coord_O)
    elif int(boundary) == 2:
        A_dist = get_adjacency_matrix_pbc(dist, coord_O, box)

    n = np.shape(A_dist)[0]
    dim = 3
    
    lower = ((-30)*np.pi)/180
    upper = ((30)*np.pi)/180
    
 
    for i in range(n):
        oxyg_1 = np.squeeze(np.array(coord_O[i,:]))
        
        if int(boundary) == 1:
            hyd_1 = np.squeeze(np.array(coord_H1[i,:]))
            hyd_2 = np.squeeze(np.array(coord_H2[i,:]))
            
        elif int(boundary) == 2:
            hyd_1 = np.squeeze(np.array(coord_H1[i,:]))
            hyd_2 = np.squeeze(np.array(coord_H2[i,:]))
            
            # I select the coordinates of H1 among all the periodic coordinates. I chose the one with the shortest distance from O.
            
            for d in range(dim):
                dist = abs(oxyg_1[d] - hyd_1[d])
                if dist > box[d]/2:
                    if oxyg_1[d] >= hyd_1[d]:
                        hyd_1[d] = box[d] + hyd_1[d]
                    else:
                        hyd_1[d] = hyd_1[d] - box[d]
            
            # I select the coordinates of H2 among all the periodic coordinates. I chose the one with the shortest distance from O.
            for d in range(dim):
                dist = abs(oxyg_1[d] - hyd_2[d])
                if dist > box[d]/2:
                    if oxyg_1[d] >= hyd_2[d]:
                        hyd_2[d] = box[d] + hyd_2[d]
                    else:
                        hyd_2[d] = hyd_2[d] - box[d]
        
        vect = A_dist[i,:]
        _,J,_ = find(vect)
        
        for j in J:
            # I select the coordinates of O2 among all the periodic coordinates. I chose the one with the shortest distance from O.
            oxyg_2 = np.squeeze(np.array(coord_O[j,:]))
            for d in range(dim):
                dist = abs(oxyg_1[d] - oxyg_2[d])
                if dist > box[d]/2:
                    if oxyg_1[d] >= oxyg_2[d]:
                        oxyg_2[d] = box[d] + oxyg_2[d]
                    else:
                        oxyg_2[d] = oxyg_2[d] - box[d]
                        
            deg_1 = angle_between(hyd_1, oxyg_1, oxyg_2)
            deg_2 = angle_between(hyd_2, oxyg_1, oxyg_2) 
            
            if lower <= deg_1 <= upper or lower <= deg_2 <= upper:
                pass
            else:
                A_dist[i,j] = 0
            #endif

    A_dist.eliminate_zeros()

    return A_dist
  

def get_graph_directed(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save):

    '''
    Construts the adjacency matrix and the graph
    '''
        
    data_O = extract_data(U, kk, list_molecules, list_nodes)
    coord_O = data_O['coord']
    box = data_O['box']
    
    data_H1 = extract_data(U, kk, ['water'], ['HW2'], ['water'])
    coord_H1 = data_H1['coord']
    
    data_H2 = extract_data(U, kk, ['water'], ['HW3'], ['water'])
    coord_H2 = data_H2['coord']


    A = oriented_adjacency_matrix(dist, coord_O, coord_H1, coord_H2, box, boundary)
    
    G = nx.from_scipy_sparse_matrix(A, create_using=nx.DiGraph)
    
    n_edges = G.number_of_edges()
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    pos={}
    for i in range_nodes:
        position = np.squeeze(np.array(coord_O[i, 0:2]))
        pos[i] = position
        
    
    return A, G, pos, coord_O, n_edges, box, dict_nodes
 
 
def get_graph(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save):

    '''
    Construts the adjacency matrix and the graph (undirected case). Returns:
    - A = adjacency matrix, 
    - G = the graph, 
    - pos = phisical position of the nodes in the central box,
    - coord = coordinates of the nodes with PBC,
    - n_edges = number of edges,
    - box = values of the box.
    '''
        
    data = extract_data(U, kk, list_molecules, list_nodes)
    coord = data['coord']
    box = data['box']
    diag_box = data['diag_box']
    
    if int(boundary) == 1:
        A = get_adjacency_matrix(dist, coord)
    elif int(boundary) == 2:
        A = get_adjacency_matrix_pbc(dist, coord, box)
    else:
        print('Value not valid')
        sys.exit(1) 
    
    G = nx.from_scipy_sparse_matrix(A)
    
    n_edges = G.number_of_edges()
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    pos={}
    for i in range_nodes:
        position = np.squeeze(np.array(coord[i, 0:2]))
        pos[i] = position
        
    return A, G, pos, coord, n_edges, diag_box
    
    
def compute_NTC(U, kk, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, normalization, dist):

    '''
    Returns NTC
    '''    
    
    if normalization == 1:     
        G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist)    
        
        # TOTAL COMMUNICABILITY    
        sub_tot_comm = total_communicability(G,  t = beta)  
        
            
    elif normalization == 2:
        G, n_edges = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
        nodelist = list(G.nodes())
        
        # TOTAL COMMUNICABILITY  
        tot_comm = total_communicability(G, t = beta)          
        sub_tot_comm = {}
        for i in nodelist:
            sub_tot_comm[i] = tot_comm[i]/n_edges
            
            
    elif normalization == 3:
        G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
        nodelist = list(G.nodes())
        
        # TOTAL COMMUNICABILITY 
        tot_comm = total_communicability(G, t = beta)  
        tnc = total_network_communicability(G)
        sub_tot_comm = {}
        for i in nodelist:
            sub_tot_comm[i] = tot_comm[i]/tnc
            
    elif normalization == 4:
        A, G, _, _, _, _,  = get_graph(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save)
        nodelist = list(G.nodes())
        
        # TOTAL COMMUNICABILITY    
        tot_comm = total_communicability(G, t = beta) 
        
        w = list(scipy.linalg.eigh(A.todense(), eigvals_only=True))
        w.sort(reverse=True)
        spectral_radius = abs(w[0])
        print(spectral_radius)
        print(w[0])
        print('')
        sub_tot_comm = {}
        for i in nodelist:
            sub_tot_comm[i] = tot_comm[i]/math.exp(spectral_radius * beta)
        with open(os.path.join(path, 'spectral_radius.txt'), 'a+') as fobj:
            fobj.write('\n')
            fobj.write('{:=18.6f} \n'.format(spectral_radius))
        
    elif normalization == 5:
        G, n_edges  = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
        nodelist = list(G.nodes())
        n_nodes = len(nodelist)
        
        # TOTAL COMMUNICABILITY  
        tot_comm = total_communicability(G, t = beta)          
        sub_tot_comm = {}
        for i in nodelist:
            sub_tot_comm[i] = tot_comm[i]/n_nodes
        
    elif normalization == 6:
        G, _  = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
        nodelist = list(G.nodes())
        n_nodes = len(nodelist)
        
        # TNC of a complete graph
        TNC_complete = n_nodes * np.exp((n_nodes - 1) * beta)
        
        # TOTAL COMMUNICABILITY  
        tot_comm = total_communicability(G, t = beta)          
        sub_tot_comm = {}
        for i in nodelist:
            sub_tot_comm[i] = tot_comm[i]/TNC_complete
    
    elif normalization == 7:
        test = 100
        G, m  = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
        n = len(G)
        nodelist = list(G.nodes())
        TNC_random = []
        
        for ii in range(test):
            G_random_nm = nx.gnm_random_graph(n, m)
            tmp = total_network_communicability(G_random_nm, t=beta)
            TNC_random.append(tmp)
            
        mean_TNC_random = np.mean(TNC_random)
            
        # TOTAL COMMUNICABILITY  
        tot_comm = total_communicability(G, t = beta)          
        sub_tot_comm = {}
        for i in nodelist:
            sub_tot_comm[i] = tot_comm[i]/mean_TNC_random
        
    elif normalization == 8:
    
        G, _  = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
        n = len(G)
        nodelist = list(G.nodes())
        
        # DEGREE 
        deg = G.degree
        deg_vect = [np.exp(deg[i]) for i in nodelist]
            
        mean_degree = np.mean(deg_vect)            
        
        # TOTAL COMMUNICABILITY  
        tot_comm = total_communicability(G, t = beta)          
        sub_tot_comm = {}
        for i in nodelist:
            sub_tot_comm[i] = tot_comm[i]/deg_vect[i]  #mean_degree
    
    else:
        print('Value not valid')
        sys.exit(1)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    with open(os.path.join(path, 'NTC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(normalization) + '_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        for ii in range_nodes:
            fobj.write('{}      {:=18.10e}\n'.format(dict_nodes[ii], sub_tot_comm[ii]))
        fobj.write('\n')
        
    tc_list = [sub_tot_comm[i] for i in range_nodes]
    
    return tc_list


def compute_NBT_total_communicability(U, kk, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Computes the non-backtracking total communicability
    '''

    A, G, _, _, _, _ = get_graph(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save)
    nodelist = list(G.nodes())
    n = len(G)
     
    # DEGREE      
    degree = G.degree
    deg = [degree[i] for i in nodelist ]
    min_deg = min(deg)
    
    if min_deg >= 2:   # if the graph does not contain leaves (nodes with degree 1)
    
        deg_matrix = np.asmatrix(deg).T
        
        one_plus_deg_minus_one = [1 + 1/(deg[i]-1) for i in range(n)]
        one_plus_deg_minus_one_matrix = np.asmatrix(one_plus_deg_minus_one).T
        
        deg_minus_one = [1/(deg[i]-1) for i in range(n)]
        deg_minus_one_matrix = np.asmatrix(deg_minus_one).T
        
        D = scipy.sparse.diags(deg)
        
        I = scipy.sparse.identity(n)
        
        O = csr_matrix((n,n))
        
        I_D = I - D

        Y = scipy.sparse.bmat([[None, I], [I_D, A]])
        
        Y = beta*Y
        
        I_O = scipy.sparse.hstack((I,O))
        
        b = np.concatenate((one_plus_deg_minus_one_matrix, deg_matrix), axis=0)
        
        vect0 = expm_multiply(Y, b)
        vect = I_O @ vect0
        
        NBT_TC = vect - deg_minus_one_matrix
        list_NBT_TC = lista = np.array(NBT_TC)[:,0]
        
        centrality = dict(zip(nodelist, list_NBT_TC))
        
    else:
    
        deg_matrix = np.asmatrix(deg).T
        D = scipy.sparse.diags(deg)
        
        I = scipy.sparse.identity(n)
        
        O = csr_matrix((n,n))
        
        I_D = I - D

        Y = scipy.sparse.bmat([[None, I], [I_D, A]])
        
        one = np.asmatrix(np.ones(n)).T
        A_2_D = (A@A - D)@ one
        
        v = np.concatenate((deg_matrix, A_2_D), axis=0)
        
        approx_Y = scipy.sparse.bmat([[beta*Y, v], [None, 0]])
        
        e_2n_plus_1 = np.asmatrix(np.zeros(2*n+1)).T
        e_2n_plus_1[-1] = 1
        
        vect_2n_plus_1 = expm_multiply(approx_Y, e_2n_plus_1)
        
        NBT_TC = beta*vect_2n_plus_1[0:n] + np.ones((n,1))
        
        list_NBT_TC = np.array(NBT_TC)[:,0]
        
        centrality = dict(zip(nodelist, list_NBT_TC))
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    with open(os.path.join(path, 'NBT_TC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        for ii in range_nodes:
            fobj.write('{}      {:=18.10e}\n'.format(dict_nodes[ii], centrality[ii]))
        fobj.write('\n')
        
    cent = [centrality[ii] for ii in range_nodes]
    
    return cent


def compute_subgraph(U, kk, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Computes subgraph centrality
    '''

    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
    nodelist = list(G.nodes())    
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    # SUBGRAPH     
    subgraph = subgraph_centrality_Bucci(G, t=beta)
       
    
    sub_subgraph = {}
    for i in range_nodes:
        sub_subgraph[i] = subgraph[i]
    
    with open(os.path.join(path, 'SUBGRAPH_beta_'+ str(beta) + '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range_nodes:
            fobj.write('{}      {:=18.6f}\n'.format(dict_nodes[ii], sub_subgraph[ii]))
    
    return sub_subgraph
         
    
def compute_degree(U, kk, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns degree centrality
    '''

    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
    n = G.number_of_nodes()
    nodelist = list(G.nodes())
    
    # DEGREE      
    degree = G.degree 
    
    sub_deg = {}
    for i in nodelist:
        sub_deg[i] = degree[i]
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    with open(os.path.join(path, 'DEGREE'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range_nodes:
            fobj.write('{}      {:=18.1f}\n'.format(str(dict_nodes[ii]), sub_deg[ii]))
    
    deg = [sub_deg[ii] for ii in range_nodes]
    
    return deg
   

def compute_CLOSENESS(U, kk, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns the values of closeness centrality
    '''

    G, _  = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
    cl = nx.closeness_centrality(G, u=None, distance=None,wf_improved=True)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    with open(os.path.join(path, 'CLOSENESS_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range_nodes:
            fobj.write('{}      {:=18.6f}\n'.format(dict_nodes[ii], cl[ii]))
    
    closeness = [cl[ii] for ii in range_nodes]
    
    return closeness

 
def compute_BETWEENNESS(U, kk, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns the values of the betweenness centrality 
    '''

    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
    nodelist = G.nodes()
    bw = nx.betweenness_centrality(G, normalized=False)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
        
    with open(os.path.join(path, 'BETWEENNESS_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range_nodes:
            fobj.write('{}      {:=18.6e}\n'.format(dict_nodes[ii], bw[ii]))
            
    betweenness = [bw[ii] for ii in range_nodes]
    
    return betweenness
 

def katz_direct_solver_sparse_matrix(G, alpha):

    '''
    Returns Katz centrality using a direct solver for sparse matrices
    '''

    n = G.number_of_nodes()
    nodelist = list(G.nodes())
    A = nx.adjacency_matrix(G, nodelist = nodelist, dtype=float) 
    
    b = np.ones((n, 1))    
    
    centrality = spsolve(scipy.sparse.identity(n) - (alpha * A), b)
    kz = dict(zip(nodelist, map(float, centrality)))
    
    return kz
 

def compute_KATZ(U, kk, alpha, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns the values of the Katz centrality using a direct solver for sparse matrices
    '''
    
    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
    nodelist = G.nodes()
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())


    # direct method with dense matrix
    #kz = nx.katz_centrality_numpy(G, alpha = alpha, beta = 1.0, normalized=False) 
    
    # powers method
    #kz = nx.katz_centrality(G, alpha = alpha, beta = 1.0, normalized=False) 
    
    # direct method with sparse matrix    
    kz = katz_direct_solver_sparse_matrix(G, alpha)

    with open(os.path.join(path, 'KATZ_'+str(alpha)+ '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range_nodes:
            fobj.write('{}      {:=18.6f}\n'.format(dict_nodes[ii], kz[ii]))
            
    katz = [kz[ii] for ii in range_nodes]
    
    return katz


def compute_eigenvector(U, kk, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns eigenvector centrality
    '''
    
    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist)  
    comp = nx.connected_components(G)
    sub_eig= {}
    nodes = len(G.nodes())

    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    for c in comp:
        G2 = G.subgraph(c)
        nodelist = list(G2.nodes())
        
        if len(nodelist) == 1:
            ii = nodelist[0]
            sub_eig[ii] = 1/nodes
        else:
            eig = nx.eigenvector_centrality(G2, max_iter=1000, tol=1e-06)
            #eig = nx.eigenvector_centrality_numpy(G2)
            sum_eig = sum([eig[i] for i in nodelist])
            for i in nodelist:
                sub_eig[i] = eig[i]/sum_eig

            for i in nodelist:
                sub_eig[i] = sub_eig[i]*(len(nodelist)/nodes)
    
    with open(os.path.join(path, 'EIGENVECTOR'+ '_boundary_' + str(boundary) + '.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range_nodes:
            fobj.write('{}      {:=18.6e}\n'.format(dict_nodes[ii], sub_eig[ii]))
            
    return
                

def get_graph_A(U, kk, list_molecules, list_nodes, boundary, dist):

    '''
    Returns only the adjacency matrix
    '''
        
    data = extract_data(U, kk, list_molecules, list_nodes)
    coord = data['coord']
    box = data['box']
    
    if int(boundary) == 1:
        A = get_adjacency_matrix(dist, coord)
    elif int(boundary) == 2:
        A = get_adjacency_matrix_pbc(dist, coord, box)
    else:
        print('Value not valid')
        sys.exit(1)
        
    return A
    
    
def count_cycles(U, kk, list_molecules, list_nodes, path, boundary, dist):

    '''
    Counts number of edges and cycles of length 3,4 and 5
    '''

    A = get_graph_A(U, kk, list_molecules, list_nodes, boundary, dist)   
    G = nx.from_scipy_sparse_matrix(A)
    
    n_edges = G.number_of_edges()    
    n,_  = np.shape(A)
    
    with open(os.path.join(path, 'number_edges.txt'), 'a+') as fobj:
        fobj.write('Number edges = {:d} \n'.format(n_edges))
        
    A2 = A@A
    
    n_nodes, _ = np.shape(A)
    range_nodes = np.arange(0, n_nodes)
    deg = G.degree
    degree = [deg[node] for node in range_nodes]
    f1 = [degree[i]*(degree[i]-1) for i in range(n_nodes)]
    F1 = sum(f1)/2   
    
    A3 = A2@A
    a3 = A3.diagonal().sum()
    F2 = a3/6
    
    A4 = A3@A
    m = n_edges
    a4 = A4.diagonal().sum()
    F5 = (a4 - 4*F1 - 2*m)/8
    
    f6 = [A3[i,i]*(degree[i]-2)/2 for i in range(n_nodes) if degree[i] > 2]
    F6 = sum(f6)
    
    A5 = A4@A
    a5 = A5.diagonal().sum()
    F8 = (1/10)*(a5 -30*F2 -10*F6)
    
    with open(os.path.join(path, 'cycles_3.txt'), 'a+') as fobj:
        fobj.write('Number cycles 3 = {:f} \n'.format( F2))
        
    with open(os.path.join(path, 'cycles_4.txt'), 'a+') as fobj:
        fobj.write('Number cycles 4 = {:f} \n'.format(F5))
        
    with open(os.path.join(path, 'cycles_5.txt'), 'a+') as fobj:
        fobj.write('Number cycles 5 = {:f} \n'.format(F8))
        
    return n_edges, F2, F5, F8


def estrada_index(U, kk, list_molecules, list_nodes, path, boundary, dist):

    '''
    Returns the Estrada index
    '''

    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist)   
    estrada_index_values = nx.estrada_index(G)
    
    with open(os.path.join(path, 'estrada_index_dist_'+str(dist)+'unweighted.txt'), 'a+') as fobj:
        fobj.write('{:=18.8e} \n'.format(estrada_index_values))
    
    return estrada_index_values


def watts_strogatz(U, kk, list_molecules, list_nodes, path, boundary, dist):

    '''
    Returns the Wattz-Strogatz coefficient clustering
    '''

    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist)   
    watts_strogatz_values = nx.average_clustering(G)
    
    with open(os.path.join(path, 'watts_strogatz.txt'), 'a+') as fobj:
        fobj.write('Watts-Strogatz coeff. clustering = {:f} \n'.format(watts_strogatz_values))
        
    return watts_strogatz_values
  
  
def transitivity_index(U, kk, list_molecules, list_nodes, path, boundary, dist):

    '''
    Returns the trasitivity index
    '''

    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
    n = len(G)

    number_of_triangles = sum(nx.triangles(G).values())/3
    
    range_nodes = np.arange(0, n)
    
    degree_norm = G.degree
    degree = [degree_norm[node] for node in range_nodes]
    P2 = sum([x*(x-1)/2 for x in degree])
    
    transitivity_index_values = 3*number_of_triangles/P2
    
    with open(os.path.join(path, 'transitivity_index.txt'), 'a+') as fobj:
        fobj.write('Transitivity index = {:f} \n'.format(transitivity_index_values))
        
    return transitivity_index_values
            

def assortativity(U, kk, list_molecules, list_nodes, nodes_to_save, path, boundary, dist): 

    '''
    Returns the assortativity coefficient of the graph G
    '''

    A, G, _, _, n_edges, _ = get_graph(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save)
    
    assortativity_values = nx.degree_assortativity_coefficient(G)

    A2 = A@A
    
    n_nodes, _ = np.shape(A)
    range_nodes = np.arange(0, n_nodes)
    deg = G.degree
    degree = [deg[node] for node in range_nodes]
    degree_minusone = np.asmatrix([i-1 for i in degree])
    f1 = [degree[i]*(degree[i]-1) for i in range(n_nodes)]
    F1 = sum(f1)/2   
    
    A3 = A2@A
    a3 = A3.diagonal().sum()
    F2 = a3/6
    
    m = n_edges    
    
    path2_path1 = F1/m       
    
    S2 = (degree_minusone@A@(degree_minusone.T))/2 - 3*F2
    
    path3_path2 = S2/F1
    
    with open(os.path.join(path, 'assortativity.txt'), 'a+') as fobj:
        fobj.write('Assortativity = {:f} \n'.format(assortativity_values))
        fobj.write('P2/P1 = {:f} \n'.format(path2_path1))
        fobj.write('P3/P2 = {:f} \n'.format(path3_path2[0,0]))
        fobj.write('\n')
        
    return assortativity_values, path2_path1, path3_path2[0,0]    


def bipartivity(U, kk, list_molecules, list_nodes, path, path3, boundary, dist): 

    '''
    Returns the bipartivity value of the graph G. Firstly it computes the eigenvalues and it saves them in a folder, then it computes the bipartivity measure
    '''
   
    eig = eigenvalues(U, kk, list_molecules, list_nodes, path3, boundary, dist)
    eig.sort(reverse=True)
    num = sum(np.cosh(eig))
    den = sum(np.exp(eig))
    bipartivity_values = num/den    
    
    with open(os.path.join(path, 'bipartivity.txt'), 'a+') as fobj:
        fobj.write('Bipartivity = {:f} \n'.format(bipartivity_values))
        
    return bipartivity_values


def eigenvalues(U, kk, list_molecules, list_nodes, path3, boundary, dist):

    '''
    Computes the eigenvalues of the adjacency matrix. It can take a long time
    '''

    A = get_graph_A(U, kk, list_molecules, list_nodes, boundary, dist) 
    n = np.shape(A)[0]
    
    w = scipy.linalg.eigh(A.todense(), eigvals_only=True)
     
    with open(os.path.join(path3, 'eigenvalues_'+str(kk)+'.txt'), 'a+') as fobj:
        for ii in range(n):
            fobj.write('{:f} \n'.format(w[ii]))
        fobj.write('\n')
        
    return list(w)


def bipartivity_for_files(path, n): 

    '''
    Returns the bipartivity value given the eigenvalues. The eigenvalues are in a folder. In the folder, there is a file.txt, with the eigenvalues, for each frame of the trajectory
    '''
    
    bipartivity_vect = []
    
    print('Select the folder that contains the eigenvalues \n')
    
    root = Tk()
    root.withdraw()
    path2 = filedialog.askdirectory()
    root.destroy()
    
    sorted_list_files = [os.path.join(path2, str(file)) for file in os.listdir(path2)]
    sorted_list_files.sort(key=os.path.getctime) #sorted by creation time

    for file in sorted_list_files:    
        eig = []
        
        with open(file, 'r') as data:
            data.seek(0)
            for line in range(n):
                d = data.readline()
                eig.append(float(d))
        eig.sort(reverse=True)
        num = sum(np.cosh(eig))
        den = sum(np.exp(eig))
        bipartivity_values = num/den
        
        with open(os.path.join(path, 'bipartivity.txt'), 'a+') as fobj:
            fobj.write('Bipartivity = {:f} \n'.format(bipartivity_values))
            
        bipartivity_vect.append(bipartivity_values)
                        
    print('Average number of bipartivity value = {} \n'.format(np.mean(bipartivity_vect)))
        
    return 


def sparsity_A(U, frame, list_molecules, list_nodes, boundary, dist, path):

    '''
    Sparsity plot of th adjacency matrix 
    '''
 
    A = get_graph_A(U, frame, list_molecules, list_nodes, boundary, dist) 
    n = np.shape(A)[0]

    fig1 = plt.figure(frame)
    plt.spy(A, markersize=0.3, color = 'k')
    fig1.savefig(os.path.join(path, 'sparsity A_frame' + str(frame) + '.png'))
    plt.show()
    plt.close()
    return


def entropy_subgraph(U, kk, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns the entropy of a graph based on the subgraph centrality: S(G, \beta) = - \sum_{i=1}^N p_i *ln(p_i), where p_i = \frac{exp(\beta A)_{i i} }{Tr(exp(\beta A))}
    '''
    
    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
    
    largest_cc = max(nx.connected_components(G), key=len)
    G1 = G.subgraph(largest_cc)    
    node_list = list(G1.nodes())
    n = len(node_list)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())

    
    if beta == 1.0:
        sub_values = nx.subgraph_centrality(G1)

    else:
        sub_values = subgraph_centrality_Bucci(G1, t=beta)
    
    with open(os.path.join(path, 'SUBGRAPH_beta_'+ str(beta) + '_boundary_' + str(boundary) +'max_connected_component.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range_nodes:
            if ii in node_list:
                fobj.write('{}      {:=18.6f}\n'.format(dict_nodes[ii], sub_values[ii]))

    EE = sum(sub_values.values())
    
    p_vect = [sub_values[ii]/EE for ii in node_list]
    
    H = [p_vect[ii]*log(p_vect[ii]) for ii in range(n) ]
    
    entropy_values = -sum(H)
    
    with open(os.path.join(path, 'entropy_subgraph_dist'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Entropy = {:f} \n'.format(entropy_values))

    return entropy_values
    
    
def entropy_TC(U, kk, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns the entropy of a graph based on the total communicability : TC(G, \beta) = - \sum_{i=1}^N p_i *ln(p_i), where p_i = \frac{(exp(\beta A) 1)_{i} }{1^T (exp(\beta A)) 1}
    '''

    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
    
    largest_cc = max(nx.connected_components(G), key=len)
    G1 = G.subgraph(largest_cc)    
    node_list = list(G1.nodes())
    n = len(node_list)

    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    TC_values = total_communicability(G1, t = beta)  

    NTC = sum(TC_values.values())
    
    p_vect = [TC_values[ii]/NTC for ii in node_list]
    
    H = [p_vect[ii]*log(p_vect[ii]) for ii in range(n) ]
    
    entropy_values = -sum(H)
    
    with open(os.path.join(path, 'entropy_TC_dist'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Entropy = {:f} \n'.format(entropy_values))
        
    with open(os.path.join(path, 'NTC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(1) + '_dist_'+ str(dist) +'max_connected_component.txt'), 'a+') as fobj:
        for ii in node_list:
            if ii in node_list:
                fobj.write('{}      {:=18.10e}\n'.format(dict_nodes[ii], TC_values[ii]))
        fobj.write('\n')

    return entropy_values
    
    
def entropy_von_neumann(U, kk, list_molecules, list_nodes, path, boundary, dist):

    '''
    Returns the Von-Neumann entropy: E(G) = - Tr( p * ln(p)), where p is the density matrix, p = L/Tr(L), and L is the Laplacian graph
    '''
    
    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
    
    largest_cc = max(nx.connected_components(G), key=len)
    G1 = G.subgraph(largest_cc)    
    A = nx.adjacency_matrix(G1, nodelist=list(G1), dtype=float)
     
    n = A.shape[0]
    e1 = np.asmatrix(np.ones((n,1)))
    degree = A@e1
    deg = np.squeeze(np.array(degree))
    
    L = scipy.sparse.diags(deg) - A
    
    trace_L = sum(deg)
    
    Den = L/trace_L
    
    w,_ = scipy.sparse.linalg.eigsh(Den, k=n-1, which='LM')
    
    H = [lambda_i * log(lambda_i) for lambda_i in w ]
    
    entropy_values = -sum(H)
    
    with open(os.path.join(path, 'entropy_Von_Neumann_dist'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Entropy = {:f} \n'.format(entropy_values))

    return entropy_values


def max_min_eigenvalues(n): 

    '''
    Returns the max and min eigenvalue
    '''
    
    print('Select the folder that contains the eigenvalues \n')
    
    root = Tk()
    root.withdraw()
    path2 = filedialog.askdirectory()
    root.destroy()
        
    min_eig = []
    max_eig = []
    diff = []
    sorted_list_files = [os.path.join(path2, str(file)) for file in os.listdir(path2)]
    sorted_list_files.sort(key=os.path.getctime) #sorted by creation time

    for file in sorted_list_files:    # the order in which the files are opened does not matter to me
        eig = []
        with open(file, 'r') as data:
            data.seek(0)
            for line in range(n):
                d = data.readline()
                eig.append(float(d))
        
        mi_e = min(eig)
        ma_e = max(eig)
        
        min_eig.append(mi_e)
        max_eig.append(ma_e)
        diff.append(ma_e + mi_e)
        
    print('Mean value of the maximum eigenvalue = {} \n'.format(np.mean(max_eig)))
    print('Mean value of the minimum eigenvalue = {} \n'.format(np.mean(min_eig)))
    print('Mean value of the differnce = {} \n'.format(np.mean(diff)))
    print('Maximum eigenvalue = {} \n'.format(max(max_eig)))
    
    return 
    

def density(U, kk, list_molecules, list_nodes, n):

    '''
    Returns the density of the central box, defining as n_nodes/volume(nm)
    '''

    for tr in U.trajectory[kk, kk+1]:
        volume = (tr.volume)/1000
        dens = n/volume
        
    return dens
    
    
def density_graph(U, kk, list_molecules, list_nodes, boundary, dist, path):

    '''
    Returns the density of the graph G
    '''

    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 

    dens = nx.density(G)
    
    with open(os.path.join(path, 'graph_density_dist_'+str(dist)+'nm.txt'), 'a+') as fobj:
        fobj.write('graph density = {:f} \n'.format(dens))
    
    return dens


def energy(path, n):

    '''
    Returns the energy of the graph E(G) = \sum_{i=1}^N |\lambda_i|
    '''

    print('Select the folder that contains the eigenvalues \n')
    
    root = Tk()
    root.withdraw()
    path2 = filedialog.askdirectory()
    root.destroy()
        
    energy_vect = []
    sorted_list_files = [os.path.join(path2, str(file)) for file in os.listdir(path2)]
    sorted_list_files.sort(key=os.path.getctime) #sorted by creation time
    
    for file in sorted_list_files:    # the order in which the files are opened does not matter to me
        eig = []
        with open(file, 'r') as data:
            data.seek(0)
            for line in range(n):
                d = data.readline()
                eig.append(float(d))
                
        energy = np.sum(np.abs(eig))
        
        with open(os.path.join(path, 'energy_graph.txt'), 'a+') as fobj:
            fobj.write('{:f} \n'.format(energy))
            fobj.write('\n')
        
        energy_vect.append(energy)
        
        
    print('Mean value of the energy = {} \n'.format(np.mean(energy_vect)))
    
    return


def ASPL_diameter_isolated(U, kk, list_molecules, list_nodes, boundary, dist):

    '''
    Returns the average shortest path length, the diameter of the graph and the number of isolated points. We consider the largest connected component
    '''
    
    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist)
    
    n_nodes_G = len(G.nodes())
    largest_cc = max(nx.connected_components(G), key=len)
    G_s = G.subgraph(largest_cc)
    
    ASPL = nx.average_shortest_path_length(G_s)
    
    diam = nx.diameter(G_s)
    
    comp = nx.connected_components(G)
    n_G2 = 0
    for c in comp:
        G2 = G.subgraph(c)
        if len(G2) > 1:
            n_G2 += G2.number_of_nodes()

    n_isolated = n_nodes_G - n_G2
    
    
    return ASPL, diam, n_isolated


def algebraic_connectivity(U, kk, list_molecules, list_nodes, boundary, dist):

    '''
    Returns the algebraic connectivity of the graph
    '''
    
    G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
    AG = nx.algebraic_connectivity(G)
    
    return AG
    
    
def plot_network_3d(coord, dist, path, name_img, color = None, vmin = None, vmax = None): 

    '''
    Plot 3d of a frame
    '''

    n = np.shape(coord)[0]
    
    if color == None:
        color = np.ones(n)
  
    if vmin == None: 
        vmin = min(color)
        vmax = max(color)
                    

    fig = plt.figure(0, figsize=(10,7))
    ax1 = fig.add_subplot(1, 1, 1, projection='3d')
    plot1 = ax1.scatter(coord[:,0], coord[:,1], coord[:,2], c=color, cmap=plt.cm.jet, vmin=vmin, vmax=vmax, s=30  )
    

    for i in range(n):
        for j in range(i+1, n):
            matrix_i_j = np.linalg.norm(coord[i]-coord[j])
            if matrix_i_j <= dist:
                ax1.plot([coord[i,0], coord[j,0]], [coord[i,1], coord[j,1]],[coord[i,2], coord[j,2]], 'k-', linewidth=0.3)
    
    if any(color != np.ones(n)):
        fig.colorbar(plot1, ax = ax1)
        
    ax1.set_rasterized(True)
    fig.savefig(os.path.join(path, 'fig_' + str(name_img) + '.png'))
    #fig.savefig(os.path.join(path, 'fig_' + str(name_img) + '.eps'))
    plt.show()
    plt.close()
    return
    
      
def plot_network_2d(pos, G, range_nodes, path, name_img, color = None, vmin = None, vmax = None): 

    '''
    Plot 2d of a frame
    '''
    
    n = len(range_nodes)
    
    if color == None:
        color = np.ones(n)
  
    if vmin == None: 
        vmin = min(color)
        vmax = max(color)

    fig, axes = plt.subplots(nrows=1, ncols=1,  figsize=(10, 5)) 
    
    ec = nx.draw_networkx_edges(G, pos, alpha=0.2, ax = axes)
    nc = nx.draw_networkx_nodes(G, pos, nodelist=range_nodes, node_color=color, cmap=plt.cm.jet, vmin=vmin, vmax=vmax, ax = axes)
    lc = nx.draw_networkx_labels(G, pos, font_size=6, ax = axes)

    if any(color != np.ones(n)):
        fig.colorbar(nc, ax = axes)
    
    axes.set_rasterized(True)
    fig.savefig(os.path.join(path, 'fig_' + str(name_img) + '.png'))
    #fig.savefig(os.path.join(path, 'fig_' + str(name_img) + '.eps'))
    plt.show()
    plt.close()  

    return
    

def save_matrix_matlab_format(U, kk, list_molecules, list_nodes, boundary, dist, new_path):

    '''
    Saves a matrix in matlab format A.mat
    '''

    A = get_graph_A(U, kk, list_molecules, list_nodes, boundary, dist) 
    
    scipy.io.savemat(os.path.join(new_path,'A_' + str(kk) +'.mat'), {'A':A})   
    
    return
    
    
def katz_dynamic_graph(U, list_molecules, list_nodes, boundary, dist, n_frames, path, alpha, nodes_to_save):

    '''
    Computes broadcast and receive centrality of the dynamic graph
    '''
    
    dynamic_graph = dnx.DynGraph()
    
    for kk in range(n_frames): 
        G, _ = get_graph_G(U, kk, list_molecules, list_nodes, boundary, dist) 
        dynamic_graph.add_nodes_from(G.nodes())
        dynamic_graph.add_interactions_from(G.edges(data=True), t=kk)
        
    bc, _ = broadcast_centrality(dynamic_graph, alpha=alpha, conj_grad_maxiter=100, conj_grad_tol=1e-7)
    rc, _ = receive_centrality(dynamic_graph, alpha=alpha, conj_grad_maxiter=100, conj_grad_tol=1e-7)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    with open(os.path.join(path, 'dynamic_communicability_'+ str(alpha)+'.txt'), 'a+') as fobj:
        for ii in range_nodes:
            fobj.write('{}      {:=18.6e}      {:=18.6e}\n'.format(dict_nodes[ii], bc[ii], rc[ii]))
        fobj.write('\n')
    
    return
    
    
def aggregated_degree(U, list_molecules, list_nodes, boundary, dist, n_frames, path, nodes_to_save):

    '''
    Computes aggregated degree  of the dynamic graph
    '''
    
    A, _, _, _, _, _= get_graph(U, 0, list_molecules, list_nodes, boundary, dist, nodes_to_save)
    
    for kk in range(1, n_frames):
        A = A + get_graph_A(U, kk, list_molecules, list_nodes, boundary, dist) 
        
    n = np.shape(A)[0]
        
    e = np.ones(n)
    e = np.asmatrix(e).T
    ad = A@e
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    with open(os.path.join(path, 'aggregated_degree.txt'), 'a+') as fobj:
        for ii in range_nodes:
            fobj.write('{}      {:=18.6e}      \n'.format(dict_nodes[ii], ad[ii,0]))
        fobj.write('\n')
    
    return


def get_weighted_adjacency_matrix_pbc(dist, coord, box, weight_edges = 1):

    '''
    Returns the weighted adjacency matrix with periodic boundary conditions
    '''
    
    coord = np.squeeze(np.array(coord))
    
    N = np.shape(coord)[0]
    
    M = np.zeros((N,N))
    
    triu = np.triu_indices_from(M, k=1)
    
    self_distance_pbc = mda.lib.distances.self_distance_array(coord, box = box)
    
    M[triu] = self_distance_pbc
    M.T[triu] = self_distance_pbc
    
    if weight_edges == 1:
        M[np.where(M > dist)] = M[np.where(M > dist)] - M[np.where(M > dist)]
        M[np.where((M <= dist)&(M > 0 ))] = np.power(M[np.where((M <= dist)&(M > 0 ))] , -1)
    elif weight_edges == 2:
        M[np.where(M > dist)] = M[np.where(M > dist)] - M[np.where(M > dist)]
        M[np.where((M <= dist)&(M > 0 ))] = np.power(M[np.where((M <= dist)&(M > 0 ))] , -2)
    elif weight_edges == 3:
        M[np.where(M > dist)] = M[np.where(M > dist)] - M[np.where(M > dist)]
        M[np.where((M <= dist)&(M > 0 ))] = np.power(M[np.where((M <= dist)&(M > 0 ))] , -3)
    elif weight_edges ==4:
        M[np.where(M > dist)] = M[np.where(M > dist)] - M[np.where(M > dist)]
        M[np.where((M <= dist)&(M > 0 ))] = np.power(M[np.where((M <= dist)&(M > 0 ))] , -0.5)
    
    A = M - scipy.sparse.diags(M.diagonal())
    
    return csr_matrix(A)

   
def get_weighted_adjacency_matrix(dist, coord, weight_edges=1):

    '''
    Returns the weighted adjacency matrix without periodic boundary conditions
    '''
    
    coord = np.asmatrix(coord)
    
    M = distance_matrix(coord, coord)
    
    if weight_edges == 1:
        M[np.where(M > dist)] = M[np.where(M > dist)] - M[np.where(M > dist)]
        M[np.where((M <= dist)&(M > 0 ))] = np.power(M[np.where((M <= dist)&(M > 0 ))] , -1)
    elif weight_edges == 2:
        M[np.where(M > dist)] = M[np.where(M > dist)] - M[np.where(M > dist)]
        M[np.where((M <= dist)&(M > 0 ))] = np.power(M[np.where((M <= dist)&(M > 0 ))] , -2)
    elif weight_edges == 3:
        M[np.where(M > dist)] = M[np.where(M > dist)] - M[np.where(M > dist)]
        M[np.where((M <= dist)&(M > 0 ))] = np.power(M[np.where((M <= dist)&(M > 0 ))] , -3)
    elif weight_edges ==4:
        M[np.where(M > dist)] = M[np.where(M > dist)] - M[np.where(M > dist)]
        M[np.where((M <= dist)&(M > 0 ))] = np.power(M[np.where((M <= dist)&(M > 0 ))] , -0.5)
    
    A = M - scipy.sparse.diags(M.diagonal())
    
    return csr_matrix(A)


def compute_weighted_NTC_degree(U, kk, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist, weight_edges= 1):

    '''
    Returns  NTC and degree for a weighted graph
    '''
    
    data = extract_data(U, kk, list_molecules, list_nodes)
    coord = data['coord']
    box = data['box']
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    if int(boundary) == 1:
        A = get_weighted_adjacency_matrix(dist, coord, weight_edges = weight_edges)
    elif int(boundary) == 2:
        A = get_weighted_adjacency_matrix_pbc(dist, coord, box, weight_edges = weight_edges)
        
    G = nx.from_scipy_sparse_matrix(A)
    nodelist = G.nodes()
    
    if weight_edges == 1:
        str_weight = '1fracd'
    elif weight_edges ==2:
        str_weight = '1fracd2'
    elif weight_edges ==3:
        str_weight = '1fracd3'
    elif weight_edges ==4:
        str_weight = '1fracsqrt(d)'
    
    
    # TOTAL COMMUNICABILITY        
    tot_comm = total_communicability(G, t = beta)                  
    
    with open(os.path.join(path, 'weighted_NTC_beta_'+ str(beta)+'_undirected_boundary_'+str(boundary)+'_dist_'+str(dist)+'_weight_'+str(str_weight)+'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range_nodes:
            fobj.write('{}      {:=18.6e}\n'.format(dict_nodes[ii], tot_comm[ii]))

    # DEGREE      
    with open(os.path.join(path, 'weighted_DEGREE_undirected_boundary_'+str(boundary)+'_dist_'+str(dist)+'_weight_'+str(str_weight)+'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range_nodes:
            fobj.write('{}      {:=18.6e}\n'.format(dict_nodes[ii], np.sum(A[ii,:])))
    return


def compute_weighted_Estrada_index(U, kk, beta, list_molecules, list_nodes, path, boundary, dist, method, weight_edges = 1):

    '''
    Returns Estrada index for a weighted graph
    '''
    
    data = extract_data(U, kk, list_molecules, list_nodes)
    coord = data['coord']
    box = data['box']
    
    if int(boundary) == 1:
        A = get_weighted_adjacency_matrix(dist, coord, weight_edges = weight_edges)
    elif int(boundary) == 2:
        A = get_weighted_adjacency_matrix_pbc(dist, coord, box, weight_edges = weight_edges)
        
        
    if weight_edges == 1:
        str_weight = '1fracd'
    elif weight_edges ==2:
        str_weight = '1fracd2'
    elif weight_edges ==3:
        str_weight = '1fracd3'
    elif weight_edges ==4:
        str_weight = '1fracsqrt(d)'
    
    if method == 1:
        w = scipy.linalg.eigh(A.todense(), eigvals_only=True)
        exp_w = [np.exp(beta*x) for x in w]
        EE = sum(exp_w)
        
        with open(os.path.join(path, 'estrada_index_weighted_dist_'+str(dist)+'_using_eigenvalues_weight_'+str(str_weight)+'.txt'), 'a+') as fobj:
            fobj.write('{:=18.6e}\n'.format(EE))
    
    elif method == 2:
        G = nx.from_scipy_sparse_matrix(A)
        nodelist = G.nodes()
        
        # Compute the subgraph centrality of each node      
        sub_values = subgraph_centrality_Bucci(G, t=beta)    
        EE = 0
        for ii in nodelist:
            EE += sub_values[ii]
            
        with open(os.path.join(path, 'estrada_index_weighted_dist_'+str(dist)+'_using_subgraph_weight_'+str(str_weight)+'.txt'), 'a+') as fobj:
            fobj.write('{:=18.6e}\n'.format(EE))
            
    return EE
    
    
def weighted_entropy_TC(U, kk, list_molecules, list_nodes, dist, path, boundary, weight_edges = 1):

    '''
    Returns the entropy of a graph based on the total communicability : TC(G, \beta) = - \sum_{i=1}^N p_i *ln(p_i), where p_i = \frac{(exp(\beta A) 1)_{i} }{1^T (exp(\beta A)) 1}
    '''
    
    data = extract_data(U, kk, list_molecules, list_nodes)
    coord = data['coord']
    box = data['box']
    
    if weight_edges == 1:
        str_weight = '1fracd'
    elif weight_edges ==2:
        str_weight = '1fracd2'
    elif weight_edges ==3:
        str_weight = '1fracd3'
    elif weight_edges ==4:
        str_weight = '1fracsqrt(d)'
    
    if int(boundary) == 1:
        A = get_weighted_adjacency_matrix(dist, coord, weight_edges = weight_edges)
    elif int(boundary) == 2:
        A = get_weighted_adjacency_matrix_pbc(dist, coord, box, weight_edges = weight_edges)
        
    G = nx.from_scipy_sparse_matrix(A)
    
    largest_cc = max(nx.connected_components(G), key=len)
    
    G1 = G.subgraph(largest_cc)    
    node_list = list(G1.nodes())
    n = len(node_list)
    
    TC_values = total_communicability(G1, t = 1.0)  

    NTC = sum(TC_values.values())
    
    p_vect = [TC_values[ii]/NTC for ii in node_list]
    
    H = [p_vect[ii]*log(p_vect[ii]) for ii in range(n) ]
    
    entropy_values = -sum(H)
    
    with open(os.path.join(path, 'weighted_entropy_TC_dist'+str(dist)+'_weight_'+str(str_weight)+'_max_connected_components.txt'), 'a+') as fobj:
        fobj.write('Entropy = {:f} \n'.format(entropy_values))

    return entropy_values
    
    
def weighted_entropy_von_neumann(U, kk, list_molecules, list_nodes, dist, path, boundary, weight_edges = 1):

    '''
    Returns the Von-Neumann entropy: E(G) = - Tr( p * ln(p)), where p is the density matrix, p = L/Tr(L), and L is the Laplacian graph
    '''
    
    data = extract_data(U, kk, list_molecules, list_nodes)
    coord = data['coord']
    box = data['box']
    
    if weight_edges == 1:
        str_weight = '1fracd'
    elif weight_edges ==2:
        str_weight = '1fracd2'
    elif weight_edges ==3:
        str_weight = '1fracd3'
    elif weight_edges ==4:
        str_weight = '1fracsqrt(d)'
    
    if int(boundary) == 1:
        A = get_weighted_adjacency_matrix(dist, coord, weight_edges = weight_edges)
    elif int(boundary) == 2:
        A = get_weighted_adjacency_matrix_pbc(dist, coord, box, weight_edges = weight_edges)
        
    G = nx.from_scipy_sparse_matrix(A)
    largest_cc = max(nx.connected_components(G), key=len)
    G1 = G.subgraph(largest_cc)    
    A = nx.adjacency_matrix(G1,  nodelist=list(G1), dtype=float)
     
    n = A.shape[0]
    e1 = np.asmatrix(np.ones((n,1)))
    degree = A@e1
    deg = np.squeeze(np.array(degree))
    
    L = scipy.sparse.diags(deg) - A
    
    trace_L = sum(deg)
    
    Den = L/trace_L
    
    w,_ = scipy.sparse.linalg.eigsh(Den, k=n-1, which='LM')
    
    H = [lambda_i * log(lambda_i) for lambda_i in w ]
    
    entropy_values = -sum(H)
    
    with open(os.path.join(path, 'weighted_entropy_Von_Neumann_dist'+str(dist)+'_weight_'+str(str_weight)+'.txt'), 'a+') as fobj:
        fobj.write('Entropy = {:f} \n'.format(entropy_values))

    return entropy_values
    
    
def compute_TC_directed(U, gggg, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, normalization, dist):

    '''
    Returns bc and rc for TC. We use the total communicability e^{\bata A} \mathbf{1}
    '''
    
    _, G, _, _, n_edges, _= get_graph_directed(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    # TOTAL DIRECTED COMMUNICABILITY        
    bc = total_communicability(G,  t = beta)  
    rc = total_communicability(G.reverse(),  t = beta)  
    
    if normalization == 1:      
        sub_bc = {}
        sub_rc = {}
        for i in range_nodes:
            sub_bc[i] = bc[i]
            sub_rc[i] = rc[i]
            
    elif normalization == 2:
        sub_bc = {}
        sub_rc = {}
        for i in range_nodes:
            sub_bc[i] = bc[i]/n_edges
            sub_rc[i] = rc[i]/n_edges
    
    else:
        print('Value not valid')
        sys.exit(1)
            
    with open(os.path.join(path, 'directed_TC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(normalization) + '_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        for ii in range_nodes:
            fobj.write('{}      {:=18.6e}      {:=18.6e}\n'.format(dict_nodes[ii], sub_bc[ii],sub_rc[ii]))
        fobj.write('\n')
    
    return


def compute_bTC_directed(U, kk, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, normalization, dist):

    '''
    Returns bc and rc for TC in the case of a directed graph. We use bipartization form of A
    '''
    
    _, G, _, _, n_edges, _= get_graph_directed(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    # TOTAL DIRECTED COMMUNICABILITY        
    bc, rc = total_directed_communicability(G, t=beta)
    
    if normalization == 1:      
        sub_bc = {}
        sub_rc = {}
        for i in range_nodes:
            sub_bc[i] = bc[i]
            sub_rc[i] = rc[i]
            
    elif normalization == 2:
        sub_bc = {}
        sub_rc = {}
        for i in range_nodes:
            sub_bc[i] = bc[i]/n_edges
            sub_rc[i] = rc[i]/n_edges
    
    else:
        print('Value not valid')
        sys.exit(1)
            
    with open(os.path.join(path, 'directed_bTC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(normalization) + '_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        for ii in range_nodes:
            fobj.write('{}      {:=18.6e}      {:=18.6e}\n'.format(dict_nodes[ii], sub_bc[ii],sub_rc[ii]))
        fobj.write('\n')
    
    return
    
    
def compute_directed_gTC(U, hh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns bc and rc for gTC in the case of a directed graph. 
    '''
    
    A, G, _, _, n_edges, _= get_graph_directed(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save)
    
    U, sigma, Vh = scipy.linalg.svd(A.todense(), full_matrices=False)
    
    sinh_sigma = [math.sinh(x) for x in sigma]
    
    func_gen = U @ scipy.sparse.diags(sinh_sigma) @ Vh
    
    one_vector = np.asmatrix(ones(n))
    
    C_h = np.flatten(func_gen @ one_vector)
    C_a = np.flatten((one_vector.T) @ func_gen)
    
    bc = dict(list(zip(nodelist, C_h)))  # hub centrality
    rc = dict(list(zip(nodelist, C_a)))  # authority centrality
 
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())

    with open(os.path.join(path, 'directed_gTC'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        for ii in range_nodes:
            fobj.write('{}      {:=18.6e}      {:=18.6e}\n'.format(dict_nodes[ii], bc[ii], rc[ii]))
        fobj.write('\n')
        
    return
    
    
def compute_directed_eigenvector(U, hhhh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns bc and rc for eigenvector centrality in the case of a directed graph
    '''
    _, G, _, _, n_edges, _ = get_graph_directed(U, hhhh, list_molecules, list_nodes, boundary, dist)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    nodelist = list(G.nodes())
    
    sub_bc = {}
    sub_rc = {}
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())

    # EIGENVECTOR  

    # We check if the graph is strongly connected
    st_conn = nx.is_strongly_connected(G)
    if st_conn == True:
    
        bc = nx.eigenvector_centrality_numpy(G.reverse())
        sum_bc = sum([bc[i] for i in nodelist])
        for i in nodelist:
            sub_bc[i] = bc[i]/sum_bc
            
        rc = nx.eigenvector_centrality_numpy(G)
        sum_rc = sum([rc[i] for i in nodelist])
        for i in nodelist:
            sub_rc[i] = rc[i]/sum_rc
    
        
        with open(os.path.join(path, 'directed_EIG'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
            for ii in range_nodes:
                fobj.write('{}      {:=18.6e}      {:=18.6e}\n'.format(dict_nodes[ii], sub_bc[ii],sub_rc[ii]))
            fobj.write('\n')
    
    else:
        print('The directed graph in frame {} is not strongly conneted. Eig. centrality is not computed'.format(hhhh))
    
    return
    
    
def compute_directed_HITS(U, hh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns bc and rc for HITS algorithm in the case of a directed graph
    '''
    _, G, _, _, _, _ = get_graph_directed(U, hh, list_molecules, list_nodes, boundary, dist)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())

    conn = nx.is_strongly_connected(G)
    
    if conn == True:   
    
        hub, auth = nx.hits(G)

        with open(os.path.join(path, 'directed_HITS'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
            for ii in range_nodes:
                    fobj.write('{}      {:=18.6e}      {:=18.6e}\n'.format(dict_nodes[ii], hub[ii], auth[ii]))
            fobj.write('\n')
            
    else:
        print('The directed graph in frame {} is not strongly connected. HITS is not computed'.format(hh))

    return
  
    
def compute_directed_degree(U, kk, list_molecules, list_nodes, nodes_to_save, path, boundary, dist):

    '''
    Returns bc and rc for degree centrality in the case of a directed graph
    '''
    _, G, _, _, n_edges, _ = get_graph_directed(U, kk, list_molecules, list_nodes, boundary, dist)
    n = len(G)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    
    # DEGREE      
    bc = G.out_degree #np.squeeze(np.array(A@np.ones((n,1))))
    rc = G.in_degree #np.squeeze(np.array(np.ones((1,n))@A))
    
    sub_deg_bc = {}
    for i in range_nodes:
        sub_deg_bc[i] = bc[i]
        
    sub_deg_rc = {}
    for i in range_nodes:
        sub_deg_rc[i] = rc[i]
    
    with open(os.path.join(path, 'directed_DEGREE'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        for ii in range_nodes:
            fobj.write('{}      {:=18.1f}      {:=18.1f}\n'.format(dict_nodes[ii], sub_deg_bc[ii], sub_deg_rc[ii]))
        fobj.write('\n')
    
    return
    
    
def nodes_patches(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, file_txt_centrality, res_threshold):

    '''
    Computes patches of nodes with CM >= threshold using periodic boundary conditions.
    Returns:
    - biggest patches
    - patches with dimensions >=5
    - the remaining patches
    If such patches do not exist, the algorithm returns empty sets
    '''

    values_cluster = []
    index_clusters = []
    
    _, G, _, _, _, _ = get_graph(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save)     
    
    dict_nodes = range_nodes_to_save(U, kk, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    n_nodes = len(range_nodes)
    
    CMs = []
                
    # import the values of CMs
    with open(file_txt_centrality, 'r') as data:
        data.seek(0)
        for ff in range(((n_nodes)+1)*kk +2):
            data.readline()

        for i in range(n_nodes):
            cent = data.readline()
             
            nn_for_read = 0
            while cent[nn_for_read] != '>':
                nn_for_read +=1

            cent = float(cent[nn_for_read+1 : nn_for_read+25])
            CMs.append(cent)    
            
    for i in range_nodes:
        cent = CMs[i]
        if cent >= res_threshold:
            values_cluster.append(cent)
            index_clusters.append(i)
            
    clusters_G = G.subgraph(index_clusters)
    connected_components_clusters = [c for c in sorted(nx.connected_components(clusters_G), key=len, reverse=True)]
    lst = [list(x) for x in connected_components_clusters]
            
    lst.sort(key=len, reverse = True)
    
    

    if len(lst) == 0:
        patches_m_larger_5 = [[]]
        other_patches = [[]]
        biggest = []
    else:
        biggest = lst[0]
        patches_m_larger_5 = [list(set(x)) for x in lst if len(x) > 4]
        other_patches = [list(set(x)) for x in lst if len(x) <= 4]
        
    
    return biggest, patches_m_larger_5, other_patches, lst
    
def nodes_patches_TC_constraint(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save, path, res_threshold,TC_frame):

    '''
    Computes patches of nodes with degree >= threshold using periodic boundary conditions.
    Returns only the patches with dimensions >=5. If such patches do not exist, the algorithm returns 0
    '''

    values_cluster = []
    index_clusters = []
        
    A, G, _, _, _, _ = get_graph(U, kk, list_molecules, list_nodes, boundary, dist, nodes_to_save) 
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    n_nodes = len(range_nodes)
                
    # TC      
    n = len(G.nodes())
    i=0
    while i<len(TC_frame):
       if float(TC_frame[i])>= int(res_threshold):
          values_cluster.append(TC_frame[i])
          index_clusters.append(i)
       i=i+1
            
    clusters_G = G.subgraph(index_clusters)

    connected_components_clusters = [c for c in sorted(nx.connected_components(clusters_G), key=len, reverse=True)]
    lst = [list(x) for x in connected_components_clusters]
            
    lst.sort(key=len, reverse = True)
    
    

    if len(lst) == 0:
        patches_m_larger_25 = [[]]
        other_patches = [[]]
        biggest = []
    else:
        biggest = lst[0]
        patches_m_larger_25 = [list(set(x)) for x in lst if len(x) > 24]
        other_patches = [list(set(x)) for x in lst if len(x) <= 24]
        


    return biggest, patches_m_larger_25, other_patches
 
