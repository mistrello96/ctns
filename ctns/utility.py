import igraph as ig
import random, pickle
from pathlib import Path
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
    refined_distribution.append(int(n_nodes - np.sum(refined_distribution)))
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

def dump_simulation(nets, path, dump_type):
    """
    Dump the simulation to a pickle file
    
    Parameters
    ----------
    nets: list
        List of ig.Graph()

    path: string
        The file path

    dump_type: string
        Can be either ["full", "light"]. In the first case, the full simulation is dumped.
        Otherwise, only a report about node status is saved

    Return
    ------
    None. The function saves a picke file containing:
        - a list of ig.Graph() if dump_type is full
        - a dict[i][j][attribute] where i is the step index, j is the node index and attribute can either be [agent_status, test_result, quarantine]

    """

    if dump_type == "full":
        try:
            with open(Path(path + ".pickle"), "wb") as f:
                pickle.dump(nets, f, protocol = pickle.HIGHEST_PROTOCOL)
        except:
            print("Simulation is too big to be dumped, try dumping in separated files. Note that this operation will be slow")
            try:
                for i in range(len(nets)):
                    net.write_picklez(fname = Path(path + "_network{}.pickle".format(i)), version = pickle.HIGHEST_PROTOCOL)
            except:
                print("Cannot dump network, the single network is too big. Please deactivate the dump option")
                sys.exit(-1)
    else:
        to_dump = dict()
        for i in range(len(nets)):
            to_dump[i] = dict()
            for j in range(len(list(nets[i].vs))):
                to_dump[i][j] = dict()
                to_dump[i][j]["agent_status"] = nets[i].vs[j]["agent_status"]
                to_dump[i][j]["test_result"] = nets[i].vs[j]["test_result"]
                to_dump[i][j]["quarantine"] = nets[i].vs[j]["quarantine"]
        with open(Path(path + ".pickle"), "wb") as f:
            pickle.dump(to_dump, f, protocol = pickle.HIGHEST_PROTOCOL)

def compute_TR(G, R_0, infection_duration, incubation_days):
    """
    Compute the transmission rate of the disease in the network.
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
    transmission_rate: float
        The transmission rate for the network

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

def plot_degree_dist(G, node_sociability = None, edge_category = None, title = None, path = None):
    """
    Print degree probability distribution of the network.
    node_sociability and edge_category can be used to filter only certain types of
    nodes and edges
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    node_sociability: string
        Consider only nodes that have sociability == node_sociability

    edge_category: string
        Consider only edges that have category == edge_category

    title: string
        The title of the plot

    path: string
        Path to the folder where to save images

    Return
    ------
    None

    """

    plt.figure(figsize = (8, 6), dpi = 300)
    if edge_category != None:
        G = G.copy()
        toRemove = []
        for edge in G.es:
            if edge["category"] != edge_category:
                toRemove.append(edge.index)
        G.delete_edges(toRemove)


    if node_sociability == None:
        degs = G.degree()
    else:
        selected_nodes = list()
        for node in G.vs:
            if node['sociability'] == node_sociability:
                selected_nodes.append(node)
        degs = G.degree(selected_nodes) 

    counted_degs = Counter(degs)
    for key in counted_degs.keys():
        counted_degs[key] = counted_degs[key] / len(G.vs())
        plt.scatter(key, counted_degs[key], marker = '.', color = "red", s = 10)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.xlabel("Degree", fontsize = 15)
    plt.ylabel("Probability", fontsize = 15)
    if title == None:
        plt.title("Degree distribution", fontsize = 20)
    else:
        plt.title(title, fontsize = 20)
    plt.tight_layout()
    plt.savefig(Path(path + "/" + title + ".png"))
    plt.savefig(Path(path + "/" +  title + ".pdf"))

# Print degree distribution
def print_degree_summary(G, path):
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
    
    plot_degree_dist(G, title = "Degree distribution", path = path)
    plot_degree_dist(G, node_sociability = "low", title = "Degree distribution low sociability", path = path)
    plot_degree_dist(G, node_sociability = "medium", title = "Degree distribution medium sociability", path = path)
    plot_degree_dist(G, node_sociability = "high", title = "Degree distribution high sociability", path = path)

    plot_degree_dist(G, edge_category = "family_contacts", title = "Degree distribution family contacts", path = path)
    plot_degree_dist(G, edge_category = "frequent_contacts", title = "Degree distribution frequent contacts", path = path)
    plot_degree_dist(G, edge_category = "occasional_contacts", title = "Degree distribution occasional contacts", path = path)
    plot_degree_dist(G, edge_category = "random_contacts", title = "Degree distribution random contacts", path = path)

'''