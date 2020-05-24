import igraph as ig
import random
import numpy as np
from collections import Counter

def generate_family_edges(G):
    """
    Generate family edges. All edges between nodes of the same family are created
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network
        
    Return
    ------
    None

    """

    toAdd = []

    for node in G.vs:
        for edge in node["family_contacts"]:
            if G.vs[edge[0]]["agent_status"] != "D" \
            and G.vs[edge[1]]["agent_status"] != "D" \
            and not G.vs[edge[0]]["quarantine"] \
            and not G.vs[edge[1]]["quarantine"]:
            #and not (edge[0], edge[1]) in toAdd
            #and not (edge[1], edge[0]) in toAdd:
                toAdd.append(edge)

    weights = np.random.randint(3, 8, len(toAdd))
    G.add_edges(toAdd)
    for (edge_index, i) in zip(G.get_eids(toAdd, directed = False), range(len(weights))):
      G.es[edge_index]["weight"] = weights[i]
      G.es[edge_index]["category"] = "family_contacts"  

def generate_occfreq_edges(G, edge_category, restriction_value):
    """
    Create edges from the node of type edge_category
    The number of edges is chosen according to the sociability of the node
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    edge_category: string
        Category of the edge attribute
    
    restriction_value: float
        How many edges are dropped in proportion to normal condition?
        
    Return
    ------
    None

    """

    toAdd = []

    for node in G.vs:
        possible_edges = node[edge_category].copy()
        random.shuffle(possible_edges)
        tmp = int(len(possible_edges) / 3)
        tmp2 = int(2 * len(possible_edges) / 3)
        if node["sociability"] == "low":
            n_edges = int(random.random() * tmp)
        if node["sociability"] == "medium":
            n_edges = tmp + int(random.random() * (tmp2 - tmp))
        if node["sociability"] == "high":
            n_edges = tmp2 + int(random.random() * (len(possible_edges) + 1 - tmp2))
        n_edges = int(n_edges * restriction_value)

        for i in range (0, n_edges):
            edge = possible_edges[i]
            if G.vs[edge[0]]["agent_status"] != "D" \
            and G.vs[edge[1]]["agent_status"] != "D" \
            and not G.vs[edge[0]]["quarantine"] \
            and not G.vs[edge[1]]["quarantine"] \
            and not G[edge[0], edge[1]]:
            #and not (edge[0], edge[1]) in toAdd \
            #and not (edge[1], edge[0]) in toAdd:
                toAdd.append(edge)

    weights = np.random.randint(1, 6, len(toAdd))
    G.add_edges(toAdd)
    for (edge_index, i) in zip(G.get_eids(toAdd, directed = False), range(len(weights))):
      G.es[edge_index]["weight"] = weights[i]
      G.es[edge_index]["category"]= edge_category  

def generate_random_edges(G, number_of_random_edges, restriction_value):
    """
    Create number_of_random_edges random edges in the contact network
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network
        
    number_of_random_edges : int
        Number of edges to create
    
    restriction_value: float
        How many edges are dropped in proportion to normal condition?
        
    Return
    ------
    None

    """

    toAdd = []
    number_of_random_edges = int(number_of_random_edges * restriction_value)
    edge_list = np.random.randint(0, len(list(G.vs)) - 1, 2 * number_of_random_edges)
    for i in range (0, 2 * number_of_random_edges, 2):
        source = edge_list[i]
        target = edge_list[i + 1]
        if not G[source, target] \
        and source != target \
        and G.vs[source]["agent_status"] != "D" \
        and G.vs[target]["agent_status"] != "D" \
        and not G.vs[source]["quarantine"] \
        and not G.vs[target]["quarantine"] \
        and not G[source, target]:
        #and not (source, target) in toAdd \
        #and not (target, source) in toAdd:
            toAdd.append((source, target))

    G.add_edges(toAdd)
    for edge_index in G.get_eids(toAdd, directed = False):
      G.es[edge_index]["weight"] = 1
      G.es[edge_index]["category"]= "random_contacts"  

def step_edges(G, restriction_value):
    """
    Removes old edges and creates new edges
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    restriction_value: float
        How many edges are dropped in proportion to normal condition?
        
    Return
    ------
    None

    """
    
    G.delete_edges(list(G.es))

    generate_family_edges(G)

    generate_occfreq_edges(G, "frequent_contacts", restriction_value)
    generate_occfreq_edges(G, "occasional_contacts", restriction_value)

    random_contact_total = len(list(G.vs)) + random.random() * (7 * len(list(G.vs)) - len(list(G.vs)))
    generate_random_edges(G, random_contact_total, restriction_value)

    # since edge generation produce a multigraph but a single edge has attributes as wanted,
    # all other edges are removed. This produce a simple graph
    # this operation is way much faster than checking befor adding new edge
    
    toRemove = list()
    for edge in G.es:
        if edge["weight"] == None:
            toRemove.append(edge)
    G.delete_edges(toRemove)  

def step_spread(G, incubation_days, infection_duration, transmission_rate):
    """
    Make the infection spread across the network
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    incubation_days: int
        Average number of days where the patient is not infective

    infection_duration: int
        Average total duration of the disease
    
    transmission_rate: float
        Value of the transmission rate for the disease in the network

    Return
    ------
    None
    
    """

    for node in G.vs:
        # update parameters if node is infected
        if node["infected"] == True:
            node["days_from_infection"] += 1
            if node["days_from_infection"] == infection_duration:
                # if infection is over, it will be dead of recovered
                node["agent_status"] = np.random.choice(["R", "D"], p = (1 - node["death_rate"], node["death_rate"]))
                node["infected"] = False
                node["symptoms"] = list()
                node["days_from_infection"] = 0
                node["needs_IC"] = False
                if node["agent_status"] == "D":
                    node["quarantine"] = 0
                    node["test_result"] = - 1
                
            # if it is still infective, spread the infection with his contacts
            if node["agent_status"] == "I":
                for contact in G.neighborhood(node)[1:]:
                    if G.vs[contact]["agent_status"] == "S":
                        prob = transmission_rate * G[node, contact] # access to the weight of the edge
                        # has the new node been infected?
                        G.vs[contact]["agent_status"] = np.random.choice(["S", "E"], p = (1 - prob, prob))
                        if G.vs[contact]["agent_status"] == "E":
                            G.vs[contact]["infected"] = True
                            G.vs[contact]["days_from_infection"] = 1

            # if the node become I, pick some symptoms
            if node["agent_status"] == "E" and node["days_from_infection"] == incubation_days:
                node["agent_status"] = "I"
                #if mild case
                case = random.uniform(0, 1)
                if case < 0.8:
                    if case < 0.05:
                        node["symptoms"].append("Loss of taste or smell")
                    if case < 0.2:
                        node["symptoms"].append("Fever")
                    if case < 0.2:
                        node["symptoms"].append("Cough")
                    if case < 0.2:
                        node["symptoms"].append("Tiredness")
                # if severe
                else:
                    if case < 0.99:
                        node["symptoms"].append("Fever")
                    if case < 0.7:
                        node["symptoms"].append("Tiredness")
                    if case < 0.6:
                        node["symptoms"].append("Cough")
                    if case < 0.3:
                        node["symptoms"].append("Dyspnea")
                
                if random.uniform(0, 1) < 0.02:
                    node["needs_IC"] = True
   
def step_test(G, nets, incubation_days, n_new_test, policy_test, contact_tracking_efficiency):
    """
    Test some nodes of the network and put the in quarantine if needed
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network

    nets: list of ig.Graph()
        History of the network
    
    incubation_days: int
        Average number of days where the patient is not infective

    n_new_test: int
        Number of new avaiable tests

    policy_test: string
        Test strategy
        Can be ["Random, Degree Centrality, Betweenness Centrality"]
    
    contact_tracking_efficiency: float
        The percentage of contacts successfully traced

    Return
    ------
    None

    """

    # create pool of nodes to test
    high_priority_test_pool = set()
    low_priority_test_pool = set()
    for node in G.vs:
        # update quarantine
        if node["quarantine"] > 0:
            node["quarantine"] -= 1
            # if node has been found positive and quarantine is over, re-test the node
            if node["quarantine"] == 0 and node["test_result"] == 1:
                high_priority_test_pool.add(node.index)
        # update test validity
        if node["test_validity"] > 0:
            node["test_validity"] -= 1
        # test also detect if the node is recovered (sierological test)
        # if node is not dead, if test validity is expired, if node is not a known recovered, add to low priority test pool
        if node["agent_status"] != "D" and node["test_validity"] <= 0 \
         and not (node["test_result"] == 0 and node["agent_status"] == "R"):
            low_priority_test_pool.add(node.index)

    low_priority_test_pool = low_priority_test_pool - high_priority_test_pool
    found_positive = set()
    #cannot directly put nodes in a set cause node object are not supported
    high_priority_test_pool = [G.vs[i] for i in high_priority_test_pool]
    low_priority_test_pool = [G.vs[i] for i in low_priority_test_pool]
    for node in high_priority_test_pool:
        if node["infected"]:
            node["test_result"] = 1
            node["quarantine"] = 14
            node["test_validity"] = 14
            found_positive.add(node.index)
        else:
            node["test_result"] = 0
            node["test_validity"] = incubation_days

    if policy_test == "Random" and n_new_test > 0:
        to_test = random.sample(low_priority_test_pool, min(len(low_priority_test_pool), n_new_test))
        for node in to_test:
            if node["infected"]:
                node["test_result"] = 1
                node["quarantine"] = 14
                node["test_validity"] = 14
                found_positive.add(node.index)
            else:
                node["test_result"] = 0
                node["test_validity"] = incubation_days

    if policy_test == "Degree Centrality" and n_new_test > 0:
        low_priority_test_pool_index = [x.index for x in low_priority_test_pool]
        degree_results = G.strength(low_priority_test_pool_index, weights = "weight")
        zipped_lists = zip(degree_results, low_priority_test_pool)
        sorted_pairs = sorted(zipped_lists, reverse = True)
        tuples = zip(*sorted_pairs)
        _, sorted_nodes = [ list(tuple) for tuple in  tuples]
        for i in range(min(n_new_test, len(low_priority_test_pool))):
            node = sorted_nodes[i]
            if node["infected"]:
                node["test_result"] = 1
                node["quarantine"] = 14
                node["test_validity"] = 14
                found_positive.add(node.index)
            else:
                node["test_result"] = 0
                node["test_validity"] = incubation_days

    if policy_test == "Betweenness Centrality" and n_new_test > 0:
        low_priority_test_pool_index = [x.index for x in low_priority_test_pool]
        betweenness_results = G.betweenness(low_priority_test_pool_index, 
                                            directed = False, weights = "weight",
                                            cutoff = None)
        zipped_lists = zip(betweenness_results, low_priority_test_pool)
        sorted_pairs = sorted(zipped_lists, reverse = True)
        tuples = zip(*sorted_pairs)
        _, sorted_nodes = [ list(tuple) for tuple in  tuples]
        for i in range(min(n_new_test, len(low_priority_test_pool))):
            node = sorted_nodes[i]
            if node["infected"]:
                node["test_result"] = 1
                node["quarantine"] = 14
                node["test_validity"] = 14
                found_positive.add(node.index)
            else:
                node["test_result"] = 0
                node["test_validity"] = incubation_days
    
    # to_quarantine will contain family contacts (quarantine 100%), 
    # possibly_quarantine will contain other contacts, quarantine influenced by contact tracing efficiency

    if len(found_positive) > 0:
        to_quarantine = set()
        possibly_quarantine = set()
        # track to max 14 days before
        for i in range(1, min(len(nets) + 1, 15)):
            net = nets[-i]
            for edge in net.es:
                if edge["category"] == "family_contacts" \
                and (edge.source in found_positive or edge.target in found_positive):
                    if edge.source in found_positive:
                        to_quarantine.add(G.vs[edge.target].index)
                    else:
                        to_quarantine.add(G.vs[edge.source].index)
                else:
                    if edge.source in found_positive:
                        possibly_quarantine.add(G.vs[edge.target].index)
                    if edge.target in found_positive:
                        possibly_quarantine.add(G.vs[edge.source].index)
        
        # set diff to remove double contacts
        possibly_quarantine = possibly_quarantine - to_quarantine
        if len(possibly_quarantine) > int(len(possibly_quarantine) * contact_tracking_efficiency):
            possibly_quarantine = random.sample(possibly_quarantine, int(len(possibly_quarantine) * contact_tracking_efficiency))
        else:
            possibly_quarantine = list(possibly_quarantine)
        to_quarantine = list(to_quarantine) + possibly_quarantine

        to_quarantine = [G.vs[i] for i in to_quarantine]

        # put them in quarantine
        # Need test?
        if to_quarantine != list() and to_quarantine != None:
            for node in to_quarantine:
                node["quarantine"] = 14  


'''
def step_vaccine(G, n_vacc, policy_vacc, vacc_pool, agent_status_report):
    """
    Vaccine some individuals in the network
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network
    
    n_vacc: int
        Number of vaccines avaiable that day
    
    policy_vacc: string
        How to distribute the daily avaiable vaccines
        Can be ["Random, Degree Centrality, Betweenness Centrality, Improved Degree Centrality, Improved Betweenness Centrality"]
    
    agent_status_report: list
        Report of nodes status

    Return
    ------
    agent_status_report: list
        Updated report of nodes status

    """

    if len(vacc_pool) > 0:
        
        # random policy
        if policy_vacc == "Random":
            if n_vacc >= len(vacc_pool):
                for index in vacc_pool:
                    agent_status_report[index] = "V"
                    G.nodes[index]["agent_status"] = "V"
            else:
                random_selected_to_vaccinate = random.sample(vacc_pool, n_vacc)
                for index in random_selected_to_vaccinate:
                    agent_status_report[index] = "V"
                    G.nodes[index]["agent_status"] = "V"
        
        # degree policy
        elif policy_vacc == "Degree Centrality":
            degree_results = G.degree(vacc_pool)
            order_by_degree = sorted(degree_results, key=lambda x: x[1], reverse = True)

            for node, value in order_by_degree[:n_vacc]:
                agent_status_report[node] = "V"
                G.nodes[node]["agent_status"] = "V"

        elif policy_vacc == "Improved Degree Centrality":
            Z = G.subgraph(vacc_pool)
            degree_results = Z.degree()
            order_by_degree = sorted(degree_results, key=lambda x: x[1], reverse = True)

            for node, value in order_by_degree[:n_vacc]:
                agent_status_report[node] = "V"
                G.nodes[node]["agent_status"] = "V"

        # Betweenness policy
        elif policy_vacc == "Betweenness Centrality":
            betweenness_results = {node: val for node, val in ig.betweenness_centrality(G).items() if node in vacc_pool}

            order_by_betweenness = sorted(betweenness_results, key=betweenness_results.get, reverse = True)
            
            for node in order_by_betweenness[:n_vacc]:
                agent_status_report[node] = "V"
                G.nodes[node]["agent_status"] = "V"

        elif policy_vacc == "Improved Betweenness Centrality":
            betweenness_results = {node: val for node, val in ig.betweenness_centrality_subset(G, sources = vacc_pool, targets = vacc_pool).items() if node in vacc_pool}

            order_by_betweenness = sorted(betweenness_results, key=betweenness_results.get, reverse = True)
            
            for node in order_by_betweenness[:n_vacc]:
                agent_status_report[node] = "V"
                G.nodes[node]["agent_status"] = "V"
    return agent_status_report  

'''

def step(G, step_index, incubation_days, infection_duration, transmission_rate,
         initial_day_restriction, restriction_duration, social_distance_strictness,
         restriction_decreasing, nets, n_test, policy_test, contact_tracking_efficiency):
    """
    Advance the simulation of one step
    
    Parameters
    ----------
    G: ig.Graph()
        The contact network
        
    step_index : int 
        index of the step

    incubation_days: int
        Average number of days where the patient is not infective

    infection_duration: int
        Average total duration of the disease
    
    transmission_rate: float
        Value of the transmission rate for the disease in the network
    
    initial_day_restriction: int
        Day index from when social distancing measures are applied

    restriction_duration: int
        How many days the social distancing last. Use -1 to make the restriction last till the end of the simulation

    social_distance_strictness: int
        How strict from 0 to 4 the social distancing measures are. 
        Represent the portion of contact that are dropped in the network (0, 25%, 50%, 75%, 100%)
        Note that family contacts are not included in this reduction

    restriction_decreasing: bool
        If the social distancing will decrease the strictness during the restriction_duration period

    nets: list of ig.Graph()
        History of the network

    n_test: int
        Number of avaiable tests

    policy_test: string
        Test strategy
    
    contact_tracking_efficiency: float
        The percentage of contacts successfully traced

    Return
    ------
    report: dict
        Counter of nodes status

    G: ig.Graph()
        The contact network
        
    """

    # generate new edges
    if not restriction_duration:
        if step_index >= initial_day_restriction:
            step_edges(G, 1 - (25 * social_distance_strictness / 100))
        else:
            step_edges(G, 1)
    else:
        if step_index >= initial_day_restriction and step_index < initial_day_restriction + restriction_duration:
            if restriction_decreasing:
                social_distance_strictness = compute_sd_reduction(step_index, initial_day_restriction, restriction_duration, social_distance_strictness)
                step_edges(G, 1 - (25 * social_distance_strictness / 100))
            else:
                step_edges(G, 1 - (25 * social_distance_strictness / 100))
        else:
            step_edges(G, 1)
    
    # spread infection
    step_spread(G, incubation_days, infection_duration, transmission_rate)
          
    # make some test on nodes     
    step_test(G, nets, incubation_days, n_test, policy_test, contact_tracking_efficiency)

    '''
    # if avaiable, vaccinate some nodes
    if step_index >= initial_day_vaccination:
        vacc_pool = [i for i in range(len(agent_status_report)) if agent_status_report[i] == "S"]
        agent_status_report = step_vaccine(G, n_vacc, policy_vacc, vacc_pool, agent_status_report)
    '''
    agent_status_report = list()
    for node in G.vs:
        agent_status_report.append(node["agent_status"])

    return Counter(agent_status_report), G  

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