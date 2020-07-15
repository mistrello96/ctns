import igraph as ig
import random, math
import numpy as np
try:
    from ctns.utility import fix_distribution_node_number, reset_network
    from ctns.steps import step
except ImportError as e:
    from utility import fix_distribution_node_number, reset_network
    from steps import step

def generate_node_list_attribute(G, attribute_name, distribution):
    """
    Add to each node a list named attribute_name containing edge representing the
    possible social contanct to each node.
    This list is used to draw the new edges for the node in the step function
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    attribute_name: string
        Name of the node attribute

    distribution: list of int
        Distribution of the "community" dimension
        
    Return
    ------
    None

    """

    # extract node list and randomize
    node_list = list(G.vs)
    random.shuffle(node_list)
    # iterate over community
    for element in distribution:
        community_nodes = list()
        # extract community members
        for i in range(0, element):
            community_nodes.append(node_list.pop())
        # for each node, pick a random number (proportional to sociability) of 
        # community member and create a possible link
        for node in community_nodes:
            tmp = int(math.ceil(len(community_nodes) / 3))
            tmp2 = int(math.ceil(2 * len(community_nodes) / 3))
            if node["sociability"] == "low":
                n_contact = 1 + int(random.random() * (tmp - 1))
            if node["sociability"] == "medium":
                n_contact = tmp + int(random.random() * (tmp2 - tmp))
            if node["sociability"] == "high":
                n_contact = tmp2 + int(random.random() * (len(community_nodes) + 1 - tmp2))

            targets = random.sample(community_nodes, n_contact)
            for target in targets:
                if target != node and not (node.index, target.index) in node[attribute_name]:
                    node[attribute_name].append((node.index, target.index))
                    target[attribute_name].append((target.index, node.index))  

def generate_node_family_attribute(G, distribution):
    """
    Add to each node a list named family_contacts containing edge representing the
    interaction with family and a family ID.
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    distribution: list of int
        Distribution of the "community" dimension
        
    Return
    ------
    None

    """

    node_counter = 0;
    for i in range(len(distribution)):
        family_nodes = list(range(node_counter, node_counter + distribution[i]))
        node_counter = node_counter + distribution[i]
        # for ech node in the family, generate a complete graph with other family members
        for source in family_nodes:
            G.vs[source]["family_id"] = i
            for target in family_nodes:
                if source != target:
                    G.vs[source]["family_contacts"].append((source, target))  

def generate_network(n_of_families, use_probabilities):
    """
    Generte contact network nodes
    
    Parameters
    ----------
    n_of_families: int
        Number of families in the network

    use_probabilities: bool
        Enables probabilities of being infected estimation
        
    Return
    ------
    G: ig.Graph()
        The contact network

    """

    #DISTRIBUTIONS
    # family distribution
    family_distribution = np.random.normal(2, 2, n_of_families).round().astype(int)
    family_distribution = [x if x > 1 else 1 for x in family_distribution]

    number_of_nodes = np.sum(family_distribution)

    # frequent contact distribution
    frequent_contact_distribution = np.random.normal(12, 3, n_of_families).round().astype(int)
    frequent_contact_distribution = list(frequent_contact_distribution[frequent_contact_distribution >= 0])

    # refine to have sum = node number
    frequent_contact_distribution = fix_distribution_node_number(frequent_contact_distribution, number_of_nodes)

    # occasional contact distribution
    occasional_contact_distribution = np.random.normal(24, 6, n_of_families).round().astype(int)
    occasional_contact_distribution = list(occasional_contact_distribution[occasional_contact_distribution >= 0])

    # refine to have sum = node number
    occasional_contact_distribution = fix_distribution_node_number(occasional_contact_distribution, number_of_nodes)

    # random edge distribution -> see step function

    # GENERATE NODES
    G = ig.Graph()
    G.add_vertices(number_of_nodes)
    for node in G.vs:
        node["sex"] = np.random.choice(["man", "woman"], p = [0.487, 0.513])
        # age is the represented by the lower bound, so if age is 20, the person has age [20-29]
        node["age"] = np.random.choice([0, 10, 20, 30, 40, 50, 60, 70, 80, 90], 
            p = [0.084, 0.096, 0.102, 0.117, 0.153, 0.155, 0.122, 0.099, 0.059, 0.013]) 

        node["family_id"] = -1
        node["family_contacts"] = list()
        node["frequent_contacts"] = list()
        node["occasional_contacts"] = list()

        # setting sociability
        if node["age"] < 6:
            node["sociability"] = "low"
            node["pre_existing_conditions"] = 0
        elif node["age"] > 70:
            node["sociability"] = np.random.choice(
                ['low', 'medium', 'high'], p=[0.75, 0.23, 0.02])
            node["pre_existing_conditions"] = np.random.choice(
                [0, 1, 2, 3], p=[0.1, 0.4, 0.3, 0.2])
        else:
            node["sociability"] = np.random.choice(
                ['low', 'medium', 'high'], p=[0.6, 0.3, 0.1])
            node["pre_existing_conditions"] = np.random.choice(
                [0, 1, 2, 3], p=[0.6, 0.2, 0.1, 0.1])
        
        if node["age"] < 10:
            node["death_rate"] = 0.002
        elif node["age"] < 20:
            node["death_rate"] = 0.001
        elif node["age"] < 30:
            node["death_rate"] = 0.001
        elif node["age"] < 40:
            node["death_rate"] = 0.004
        elif node["age"] < 50:
            node["death_rate"] = 0.009
        elif node["age"] < 60:
            node["death_rate"] = 0.026
        elif node["age"] < 70:
            node["death_rate"] = 0.10
        elif node["age"] < 80:
            node["death_rate"] = 0.249
        elif node["age"] < 90:
            node["death_rate"] = 0.308
        elif node["age"] < 100:
            node["death_rate"] = 0.261

        node["agent_status"] = 'S'
        node["infected"] = False
        node["days_from_infection"] = 0
        node["quarantine"] = 0
        node["test_validity"] = 0
        node["test_result"] = -1
        node["symptoms"] = list()
        node["needs_IC"] = False
        if use_probabilities:
            node["prob_inf"] = 0.0

    # GENERATE NODE LIST ATTRIBUTES
    # generate families contact attribute in nodes
    generate_node_family_attribute(G, family_distribution)
    # generate frequent contact attribute in nodes
    generate_node_list_attribute(G, "frequent_contacts", frequent_contact_distribution)
    # generate occasional contact attribute in nodes
    generate_node_list_attribute(G, "occasional_contacts", occasional_contact_distribution)

    return G  

def init_infection(G, n_initial_infected_nodes):
    """
    Make random nodes infected in the network
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    n_initial_infected_nodes: int
        Number of infected nodes to put on the network

    Return
    ------
    None

    """
    
    for node in random.sample(list(G.vs), n_initial_infected_nodes):
        node["agent_status"] = "E"
        node["infected"] = True
        node["days_from_infection"] = 1

    if "prob_inf" in G.vs.attributes():
        number_of_nodes = len(list(G.vs))
        for node in G.vs:
            node["prob_inf"] = n_initial_infected_nodes / number_of_nodes

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
        step(G, i, 0, 0, 0, 0, 0, 0, False, list(), 0, "Random", 0, 0, False, 0, 0, 0)

        degrees = G.strength(list(range(len(G.vs))), weights = "weight")
        avr_deg.append(sum(degrees) / len(degrees))
    #reset network status
    reset_network(G)
    return R_0 /((infection_duration - incubation_days) * (sum(avr_deg) / len(avr_deg)))