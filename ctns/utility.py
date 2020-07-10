import igraph as ig
import random
import numpy as np
from collections import Counter
#import matplotlib.pyplot as plt
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
        node["prob_inf"] = 0.0
        node["quarantine"] = 0
        node["test_validity"] = 0
        node["test_result"] = -1
        node["symptoms"] = list()  

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
        step(G, i, 0, 0, 0, 0, 0, 0, False, list(), 0, "Random", 0, False, 0.5, 0.5)

        degrees = G.strength(list(range(len(G.vs))), weights = "weight")
        avr_deg.append(sum(degrees) / len(degrees))
    #reset network status
    reset_network(G)
    return R_0 /((infection_duration - incubation_days) * (sum(avr_deg) / len(avr_deg)))

def update_dump_report(to_dump, net):
    """
    Update the simulation dump in case light dump is selected
    
    Parameters
    ----------
    to_dump: dict
        Old dump doctionary

    net: ig.Graph()
        The contact network

    Return
    ------
    to_dump: dict
        Updated dump doctionary

    """

    network_report = Counter(net.vs["agent_status"])
    tested = 0
    positive = 0
    quarantined = 0
    for node in net.vs:
        if node["test_result"] != -1:
            tested += 1
        if node["test_result"] == 1:
            positive += 1
        if node["quarantine"] != 0:
            quarantined += 1

    to_dump['S'].append(network_report['S'])
    to_dump['E'].append(network_report['E'])
    to_dump['I'].append(network_report['I'])
    to_dump['R'].append(network_report['R'])
    to_dump['D'].append(network_report['D'])
    to_dump['quarantined'].append(quarantined)
    to_dump['positive'].append(positive)
    to_dump['tested'].append(tested)
    to_dump['total'].append(sum(network_report.values()))

    return to_dump