import igraph as ig
import random
import numpy as np
#import matplotlib.pyplot as plt
from collections import Counter
try:
    from ctns.steps import step
except ImportError as e:
    from steps import step

def fix_distribution_node_number(distribution, n_nodes):
    """
    Make the sum of the elements in the distribution be equal to the number of nodes
    
    Parameters
    ----------
    distribution: list of int
        Distribution of int

    n_nodes: int
        Number of nodes in The contact network
        
    Return
    ------
    distribution: list of int
        Distribution where the sum of elements is equal to n_nodes

    """
    refined_distribution = list()
    while True:
        refined_distribution.append(distribution.pop())
        if (np.sum(refined_distribution) > n_nodes):
            refined_distribution.pop()
            break
    refined_distribution.append(n_nodes - np.sum(refined_distribution))    
    return refined_distribution  

def reset_network(G):
    """
    Reset network status, removing all edges and setting all the infection 
    related attributes of each node to the default value
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network
        
    Return
    ------
    None

    """
    G.delete_edges(list(G.es))
    for node in G.vs:
        node["agent_status"] = 'S'
        node["infected"] = False
        node["days_from_infection"] = 0
        node["probability_of_being_infected"] = 0.0
        node["quarantine"] = 0
        node["test_validity"] = 0
        node["test_result"] = -1
        node["symptoms"] = list()  

def dump_network(G, path):
    """
    Dump the contact network to a file
    N.B.: if you need to reload the network after dump, use only pickle, not json 
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    path: string
        The file path
    Return
    ------
    None

    """
    G.write_picklez("{}.pickle".format(path))

def compute_IR(G, R_0, infection_duration, incubation_days):
    """
    Compute the infection rate of the disease in the network.
    The factor is computed as R_0 / (average_weighted_degree * (infection_duration - incubation_days))
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    R_0: float
        R_0 of the disease

    infection_duration: int
        Average total duration of the disease

    incubation_days: int
        Average number of days where the patient is not infective

    Return
    ------
    infection_rate: float
        The infection rate for the network

    """
    avr_deg = list()
    # compute average weighted degree on 20 steps
    for i in range (20):
        step(G, i, 0, 0, 0, 0, 0, 0, False, list(), 0, "Random", 0)

        degrees = G.strength(list(range(len(G.vs))), weights = "weight")
        avr_deg.append(sum(degrees) / len(degrees))
    #reset network status
    reset_network(G)
    return R_0 /((infection_duration - incubation_days) * (sum(avr_deg) / len(avr_deg)))

'''

def plot_degree_dist(G, node_sociality = None, edge_category = None, title = None):
    """
    Print degree probability distribution of the network.
    node_sociality and edge_category can be used to filter only certain types of
    nodes and edges
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    node_sociality: string
        Consider only nodes that have sociality == node_sociality

    edge_category: string
        Consider only edges that have category == edge_category

    title: string
        The title of the plot

    Return
    ------
    None

    """
    #plt.figure(figsize = (8, 6), dpi = 300)
    if edge_category != None:
        G = G.copy()
        toRemove = []
        for edge in G.es:
            if edge["category"] != edge_category:
                toRemove.append(edge.index)
        G.delete_edges(toRemove)


    if node_sociality == None:
        degs = G.degree()
    else:
        selected_nodes = list()
        for node in G.vs:
            if node['sociality'] == node_sociality:
                selected_nodes.append(node)
        degs = G.degree(selected_nodes) 

    conuted_degs = Counter(degs)
    for key in conuted_degs.keys():
        conuted_degs[key] = conuted_degs[key] / len(G.vs())
        plt.scatter(key, conuted_degs[key], marker = '.', color = "red", s = 10)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.xlabel("Degree", fontsize = 15)
    plt.ylabel("Probability", fontsize = 15)
    if title == None:
        plt.title("Degree distribution", fontsize = 20)
    else:
        plt.title(title, fontsize = 20)
    plt.tight_layout()
    plt.show()  

# Print degree distribution
def print_degree_summary(G):
    """
    Print degree probability distribution dummary of the network.
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    Return
    ------
    None

    """
    plot_degree_dist(G, title = "Degree distribution")
    plot_degree_dist(G, node_sociality = "low", title = "Degree distribution low sociality")
    plot_degree_dist(G, node_sociality = "medium", title = "Degree distribution medium sociality")
    plot_degree_dist(G, node_sociality = "high", title = "Degree distribution high sociality")

    plot_degree_dist(G, edge_category = "family_contacts", title = "Degree distribution family contacts")
    plot_degree_dist(G, edge_category = "frequent_contacts", title = "Degree distribution frequent contacts")
    plot_degree_dist(G, edge_category = "occasional_contacts", title = "Degree distribution occasional contacts")
    plot_degree_dist(G, edge_category = "random_contacts", title = "Degree distribution random contacts")

'''