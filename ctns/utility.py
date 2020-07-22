import numpy as np
from collections import Counter

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
        node["quarantine"] = 0
        node["test_validity"] = 0
        node["test_result"] = -1
        node["symptoms"] = list()
        if "prob_inf" in G.vs.attributes():
            node["prob_inf"] = 0.0

def update_dump_report(to_dump, net, new_positive_counter, use_probabilities):
    """
    Update the simulation dump in case light dump is selected
    
    Parameters
    ----------
    to_dump: dict
        Old dump doctionary

    net: ig.Graph()
        The contact network

    new_positive_counter: int
        Number of new positive nodes found on this iteration

    use_probabilities: bool
        Enables probabilities of being infected estimation

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
    to_dump['new_positive_counter'].append(new_positive_counter)
    if use_probabilities:
        to_dump['avg_prob_inf'].append(np.mean(net.vs['prob_inf']))
    to_dump['total'].append(sum(network_report.values()))

    return to_dump

def compute_sd_reduction(step_index, initial_day_restriction, restriction_duration, social_distance_strictness):
    """
    Calculate the decreased social distance strictness index
    
    Parameters
    ----------
        
    step_index : int 
        index of the step
    
    initial_day_restriction: int
        Day index from when social distancing measures are applied

    restriction_duration: int
        How many days the social distancing last. Use -1 to make the restriction last till the end of the simulation

    social_distance_strictness: int
        How strict from 0 to 4 the social distancing measures are. 
        Represent the portion of contact that are dropped in the network (0, 25%, 50%, 75%, 100%)
        Note that family contacts are not included in this reduction

    Return
    ------
    social_distance_strictness: int
        Updated value for social_distance_strictness
        
    """
    
    default_days = restriction_duration // social_distance_strictness
    spare_days = restriction_duration - default_days * social_distance_strictness

    strictness_sequence = list()

    for i in range(social_distance_strictness):
        updated_days = default_days
        if spare_days > 0:
            updated_days += 1
            spare_days -= 1
        strictness_sequence.append(updated_days)


    for i in range(1, len(strictness_sequence)):
      strictness_sequence[i] += strictness_sequence[i - 1]

    social_distance_reduction = None

    for i in range(len(strictness_sequence)):
      if step_index - initial_day_restriction < strictness_sequence[i]:
        social_distance_reduction = i
        break

    return social_distance_strictness - social_distance_reduction

def perform_test(node, incubation_days, use_probabilities):
    """
    Perform tampon and syerological test on the node
    
    Parameters
    ----------
        
    node: igraph.Vertex
        Node to test

    incubation_days: int
        Average number of days where the patient is not infective

    use_probabilities: bool
        Enables probabilities of being infected estimation

    Return
    ------
    positive: bool
        If the node tested is found positive
        
    """

    if node["infected"]:
        node["test_result"] = 1
        node["quarantine"] = 14
        node["test_validity"] = 14
        if use_probabilities:
            node["prob_inf"] = 1
        return True
    else:
        node["test_result"] = 0
        node["test_validity"] = incubation_days
        if use_probabilities:
            node["prob_inf"] = 0
        return False

    return None

def retrive_to_test_quarantine(nodes, values, n_new_test, reverse = True):
    """
    Test some nodes of the network and put the in quarantine if needed
    
    Parameters
    ----------

    nodes: list
        List of all nodes that needs a test

    values: list
        List of values used to sort the previous list

    n_new_test: int
        Number of new avaiable tests

    reverse: bool
        Sort using descending order

    Return
    ------
    selected: list
        List of nodes that will be tested

    """

    selected = list()

    zipped_lists = zip(values, nodes)
    sorted_pairs = sorted(zipped_lists, reverse = reverse)
    tuples = zip(*sorted_pairs)
    sorted_values, sorted_nodes = [list(tuple) for tuple in  tuples]
    selected = sorted_nodes[:min(n_new_test, len(nodes))]

    return selected