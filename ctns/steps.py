import random
import numpy as np
try:
    from ctns.utility import retrive_to_test_quarantine, perform_test, compute_sd_reduction
except ImportError as e:
    from utility import retrive_to_test_quarantine, perform_test, compute_sd_reduction

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

def step_spread(G, incubation_days, infection_duration, transmission_rate, use_probabilities, alpha, gamma):
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

    use_probabilities: bool
        Enables probabilities of being infected estimation

    alpha: float
        Parameter to regulate probability of being infected contact decay. Domain = (0, 1). Lower values corresponds to higher probability decay

    gamma: float
        Parameter to regulate probability of being infected contact diffusion. Domain = (0, +inf). Higher values corresponds to stronger probability diffusion

    Return
    ------
    None
    
    """

    if use_probabilities:
        old_prob =  G.vs["prob_inf"]

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
                    if use_probabilities:
                        node["prob_inf"] = 0
                
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

    # update prob of being infected
    if use_probabilities:
        nodes_contact_probs = [1] * len(list(G.vs))   
        
        for edge in G.es:
            weight = edge["weight"]
            source = edge.source
            target = edge.target
            nodes_contact_probs[source] *= 1 - (1 - 1 / (1 + alpha * 21)) * old_prob[target] * (1 - np.e**(- gamma * weight))
            nodes_contact_probs[target] *= 1 - (1 - 1 / (1 + alpha * 21)) * old_prob[source] * (1 - np.e**(- gamma * weight))

        for node in G.vs:
            if node["agent_status"] != "D" and not (node["test_result"] == 0 and node["agent_status"] == "R"):
                index = node.index
                node["prob_inf"] = 1 - (1 / (1 + alpha * 21)) * old_prob[index] - (1 - old_prob[index]) * nodes_contact_probs[index]
   
def step_test(G, nets, incubation_days, n_new_test, policy_test, contact_tracing_efficiency, quarantine_efficiency, use_probabilities, lambdaa):
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
        Can be ["Random, Degree Centrality, Betweenness Centrality, PBI"]
    
    contact_tracing_efficiency: float
        The percentage of contacts successfully traced

    quarantine_efficiency: float
        The percentage of contacts put under quarantine

    use_probabilities: bool
        Enables probabilities of being infected estimation

    lambdaa: float
        Parameter to regulate influence of contacts with a positive. Domain = (0, 1). Higher values corresponds to stronger probability diffusion

    Return
    ------
    new_positive_counter: int
        Number of new positive nodes found on this iteration

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

    low_priority_test_pool -= high_priority_test_pool
    found_positive = set()
    #cannot directly put nodes in a set cause node object are not supported
    high_priority_test_pool = [G.vs[i] for i in high_priority_test_pool]
    low_priority_test_pool = [G.vs[i] for i in low_priority_test_pool]
    for node in high_priority_test_pool:
        result = perform_test(node, incubation_days, use_probabilities)
        if result:
            found_positive.add(node.index)

    to_test = list()

    new_positive_counter = 0

    if n_new_test:
        if policy_test == "Random":
            to_test = random.sample(low_priority_test_pool, min(len(low_priority_test_pool), n_new_test))

        if policy_test == "Degree Centrality":
            low_priority_test_pool_index = [x.index for x in low_priority_test_pool]
            degree_results = G.strength(low_priority_test_pool_index, weights = "weight")
            to_test = retrive_to_test_quarantine(low_priority_test_pool, degree_results, n_new_test)

        if policy_test == "Betweenness Centrality":
            low_priority_test_pool_index = [x.index for x in low_priority_test_pool]
            betweenness_results = G.betweenness(low_priority_test_pool_index, 
                                                directed = False, weights = "weight",
                                                cutoff = None)
            to_test = retrive_to_test_quarantine(low_priority_test_pool, betweenness_results, n_new_test)

        if policy_test == "PBI":
            probs_infected = [i["prob_inf"] for i in low_priority_test_pool]
            to_test = retrive_to_test_quarantine(low_priority_test_pool, probs_infected, n_new_test)

    for node in to_test:
        result = perform_test(node, incubation_days, use_probabilities)
        if result:
            found_positive.add(node.index)
            new_positive_counter += 1
            
    # tracked will contain family contacts (quarantine 100%), 
    # possibly_tracked will contain other contacts, quarantine influenced by contact tracing efficiency
    ct_nets = list(nets) + [G.copy()]

    if len(found_positive) > 0:
        tracked = set()
        possibly_tracked = set()
        # trace contacts
        for net in ct_nets:
            for edge in net.es:
                if edge["category"] == "family_contacts" \
                and (edge.source in found_positive or edge.target in found_positive):
                    if edge.source in found_positive:
                        tracked.add(G.vs[edge.target].index)
                    else:
                        tracked.add(G.vs[edge.source].index)
                else:
                    if edge.source in found_positive:
                        possibly_tracked.add(G.vs[edge.target].index)
                    if edge.target in found_positive:
                        possibly_tracked.add(G.vs[edge.source].index)
        
        # set diff to remove double contacts
        possibly_tracked -= tracked
        if int(len(possibly_tracked) * contact_tracing_efficiency) > 0:
            possibly_tracked = random.sample(possibly_tracked, int(len(possibly_tracked) * contact_tracing_efficiency))
        else:
            possibly_tracked = list(possibly_tracked)
        tracked = list(tracked)
        if use_probabilities:
            # update prob of being infected of current contact 
            for node in tracked + possibly_tracked:
                for contact in G.neighborhood(node)[1:]:
                    contact_node = G.vs[contact]
                    if contact_node["agent_status"] != "D" and not (contact_node["test_result"] == 0 and contact_node["agent_status"] == "R"):
                        current_contact_weight = G[node, contact]
                        contact_node["prob_inf"] = contact_node["prob_inf"] \
                                                    + lambdaa * np.e**(-(1 / current_contact_weight)) * (1 - contact_node["prob_inf"])

            
            # update prob of being infected of past tracked contact 
                for net_index in range(1, len(nets) + 1):
                    net = nets[- net_index]
                    for contact in net.neighborhood(node)[1:]:
                        contact_node = G.vs[contact]
                        if contact_node["agent_status"] != "D" and not (contact_node["test_result"] == 0 and contact_node["agent_status"] == "R"):
                            if contact in tracked + possibly_tracked:
                                current_contact_weight = net[node, contact]
                                # use net_index instead of net_index + 1 since net index is already +1 from line above 
                                contact_node["prob_inf"] = contact_node["prob_inf"] \
                                                            + lambdaa * np.e**(- (net_index) * (1 / current_contact_weight)) * (1 - contact_node["prob_inf"])
        
        possibly_quarantine = list()

        if quarantine_efficiency == 1:
            possibly_quarantine = list(possibly_tracked)
        elif quarantine_efficiency == 0 or possibly_tracked == list():
            possibly_quarantine = list()
        else:
            if policy_test == "Random":
                possibly_quarantine = random.sample(possibly_tracked, int(len(possibly_tracked) * quarantine_efficiency))
            if policy_test == "Degree Centrality":
                values = G.strength(possibly_tracked, weights = "weight")
                possibly_quarantine = retrive_to_test_quarantine(possibly_tracked, values, int(len(possibly_tracked) * quarantine_efficiency))
            if policy_test == "Betweenness Centrality":
                values = G.betweenness(possibly_tracked, directed = False, weights = "weight", cutoff = None)
                possibly_quarantine = retrive_to_test_quarantine(possibly_tracked, values, int(len(possibly_tracked) * quarantine_efficiency))
            if policy_test == "PBI":
                values = [G.vs[x]["prob_inf"] for x in possibly_tracked]
                possibly_quarantine = retrive_to_test_quarantine(possibly_tracked, values, int(len(possibly_tracked) * quarantine_efficiency))

        to_quarantine = [G.vs[i] for i in tracked + possibly_quarantine]
        
        # put them in quarantine
        for node in to_quarantine:
            node["quarantine"] = 14

    return new_positive_counter

def step(G, step_index, incubation_days, infection_duration, transmission_rate,
         initial_day_restriction, restriction_duration, social_distance_strictness,
         restriction_decreasing, nets, n_test, policy_test, contact_tracing_efficiency,
         quarantine_efficiency, use_probabilities, alpha, gamma, lambdaa):
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
    
    contact_tracing_efficiency: float
        The percentage of contacts successfully traced

    quarantine_efficiency: float
        The percentage of contacts put under quarantine

    use_probabilities: bool
        Enables probabilities of being infected estimation

    alpha: float
        Parameter to regulate probability of being infected contact decay. Domain = (0, 1). Lower values corresponds to higher probability decay

    gamma: float
        Parameter to regulate probability of being infected contact diffusion. Domain = (0, +inf). Higher values corresponds to stronger probability diffusion

    lambdaa: float
        Parameter to regulate influence of contacts with a positive. Domain = (0, 1). Higher values corresponds to stronger probability diffusion

    Return
    ------
    G: ig.Graph()
        The contact network

    new_positive_counter: int
        Number of new positive nodes found on this iteration

        
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
    step_spread(G, incubation_days, infection_duration, transmission_rate, use_probabilities, alpha, gamma)
          
    # make some test on nodes
    new_positive_counter = step_test(G, nets, incubation_days, n_test, policy_test, contact_tracing_efficiency, quarantine_efficiency, use_probabilities, lambdaa)

    return G, new_positive_counter