import numpy as np
from pathlib import Path
from collections import deque, Counter
import sys, random, time, pickle
try:
    from ctns.generator import generate_network, init_infection, compute_TR
    from ctns.steps import step
    from ctns.utility import update_dump_report
except ImportError as e:
    from generator import generate_network, init_infection, compute_TR
    from steps import step
    from utility import update_dump_report

def run_simulation(n_of_families = 250,
    use_steps = True,
    number_of_steps = 150,
    incubation_days = 5,
    infection_duration = 21,
    initial_day_restriction = 50,
    restriction_duration = 21,
    social_distance_strictness = 2,
    restriction_decreasing = True,
    n_initial_infected_nodes = 3,
    R_0 = 2.9,
    n_test = 3,
    policy_test = "Random",
    contact_tracing_efficiency = 0.8,
    contact_tracing_duration = 14,
    quarantine_efficiency = 0.4,
    use_random_seed = True,
    seed = 42,
    use_probabilities = False,
    alpha = 0.5,
    gamma = 0.003,
    lambdaa = 0.02,
    dump_type = None,
    path = None):
    """
    Execute the simulation and dump/return resulting networks

    Parameters
    ----------
    n_of_families: int
        Number of families in the network
        
    use_steps : bool 
        Use a fixed number of steps or keep the simulation going untill the spreading is not over

    number_of_steps : int
        Number of simulation step to perform

    incubation_days: int
        Number of days where the patient is not infective

    infection_duration: int
        Total duration of the disease per patient

    initial_day_restriction: int
        Day index from when social distancing measures are applied

    restriction_duration: int
        How many days the social distancing last. Use 0 to make the restriction last till the end of the simulation

    social_distance_strictness: int
        How strict from 0 to 4 the social distancing measures are. 
        Represent the portion of contact that are dropped in the network (0, 25%, 50%, 75%, 100%)
        Note that family contacts are not included in this reduction

    restriction_decreasing: bool
        If the social distancing will decrease the strictness during the restriction_duration period

    n_initial_infected_nodes: int
        Number of nodes that are initially infected

    R_0: float
        The R0 facotr of the disease

    n_test: int
        Number of avaiable tests

    policy_test: string
        Strategy with which test are made. Can be Random, Degree Centrality, Betweenness Centrality or PBI

    contact_tracing_efficiency: float
        The percentage of contacts successfully traced back in the past contact_tracing_duration days

    contact_tracing_duration: int
        How many days the contacts of a node are tracked back

    quarantine_efficiency: float
        The percentage of contacts put under quarantine

    use_random_seed: bool
        Id use or not a fixed random seed

    seed: int
        The random seed value

    use_probabilities: bool
        Enables probabilities of being infected estimation

    alpha: float
        Parameter to regulate probability of being infected contact decay. Domain = (0, 1). Lower values corresponds to higher probability decay

    gamma: float
        Parameter to regulate probability of being infected contact diffusion. Domain = (0, +inf). Higher values corresponds to stronger probability diffusion

    lambdaa: float
        Parameter to regulate influence of contacts with a positive. Domain = (0, 1). Higher values corresponds to stronger probability diffusion

    dump_type: string
        Can be either ["full", "light"]. In the first case, the full simulation (all nets structure) is dumped.
        Otherwise, only a report about node status is saved.
        NB, full method will use significally more RAM than light; also the dump will have much bigger size.
        The dumped file will have the following structure:
        - a dict containig simulation parameters and a list of ig.Graph() if dump_type is full
        - a dict[class] where class can be [S, E, I, R, D, quarantined, positive, tested, total] and value is a list of the corresponding attribute value on day i

    path: string
        The path to the file/folder where the networks will be saved

    Return
    ------
    None

    """

    # generate new edges
    if use_random_seed:
        np.random.seed(seed = seed)
        random.seed(seed)
    else:
        np.random.seed(int(time.time()))
        random.seed(time.time())
    
    # check values
    if n_of_families < 10:
        print("Invalid number of families. Use at least 10 families")
        sys.exit()
    if use_steps:
        if number_of_steps < 0:
            print("Invalid number of steps")
            sys.exit()
    if infection_duration < 0:
        print("Invalid infection duration")
        sys.exit()
    if incubation_days < 0 or incubation_days >= infection_duration:
        print("Invalid incubation duration")
        sys.exit()
    if initial_day_restriction < 0 :
        print("Invalid initial day social distancing")
        sys.exit()
    if social_distance_strictness < 0 or social_distance_strictness > 4:
        print("Invalid social distancing value")
        sys.exit()
    if n_initial_infected_nodes < 0 or n_initial_infected_nodes > n_of_families:
        print("Invalid number of initial infected nodes")
        sys.exit()
    if R_0 < 0:
        print("Invalid value of R0")
        sys.exit()
    if n_test < 0:
        print("Invalid number of test per day")
        sys.exit()
    if not (policy_test == "Random" or policy_test == "Degree Centrality" or policy_test == "Betweenness Centrality" or "PBI"):
        print("Invalid test strategy")
        sys.exit()
    if contact_tracing_efficiency < 0 or contact_tracing_efficiency > 1:
        print("Invalid contact tracing efficiency")
        sys.exit()
    if quarantine_efficiency < 0 or quarantine_efficiency > 1:
        print("Invalid quarantine efficiency")
        sys.exit()
    if contact_tracing_duration < 0:
        print("Invalid contatc tracing duration")
        sys.exit()
    if restriction_duration < 0:
        print("Invalid restriction restriction_duration")
        sys.exit()
    if restriction_duration == 0:
        restriction_decreasing = False
        social_distance_strictness = 0
    if social_distance_strictness == 0:
        restriction_decreasing = False
        restriction_duration = 0
    if gamma < 0 :
        print("Value of gamma must greater than 0")
        sys.exit()
    if lambdaa < 0 or lambdaa > 1:
        print("Value of lambda must be 0 < lambda < 1")
        sys.exit()
    if dump_type != "full" and dump_type != "light":
        print("Invalid dump type")
        sys.exit()
    if path == None:
        print("Invalid path")
        sys.exit()
    if not use_probabilities and policy_test == "PBI":
        print("Cannot use PBI if probability of being infected is not enabled")
        sys.exit()

    # making parameters consistent
    if restriction_duration == 0 or social_distance_strictness == 0:
        restriction_decreasing = False
        social_distance_strictness = 0
        initial_day_restriction = 0
        restriction_duration = 0

    if n_test == 0:
        policy_test = None
        contact_tracing_efficiency = 0
        contact_tracing_duration = 0

    if contact_tracing_duration == 0 or contact_tracing_efficiency == 0:
        contact_tracing_duration = 0
        contact_tracing_efficiency = 0

    if not use_probabilities:
        gamma = 0
        lambdaa = 0

    config = locals()
    a = time.perf_counter()
    # init network
    G = generate_network(n_of_families, use_probabilities)
    transmission_rate = compute_TR(G, R_0, infection_duration, incubation_days)
    init_infection(G, n_initial_infected_nodes)

    nets = deque(maxlen = contact_tracing_duration)
    if dump_type == "full":
        to_dump = dict()
        to_dump["nets"] = list()
        to_dump["parameters"] = config
    if dump_type == "light":
        to_dump = dict()
        to_dump['S'] = list()
        to_dump['E'] = list()
        to_dump['I'] = list()
        to_dump['R'] = list()
        to_dump['D'] = list()
        to_dump['quarantined'] = list()
        to_dump['positive'] = list()
        to_dump['tested'] = list()
        to_dump['total'] = list()
        to_dump['new_positive_counter'] = list()
        to_dump['avg_prob_inf'] = list()
        to_dump['parameters'] = config

    if use_steps:
        for sim_index in range (0, number_of_steps):
            net, new_positive_counter = step(G, sim_index, incubation_days, infection_duration, transmission_rate,
                             initial_day_restriction, restriction_duration, social_distance_strictness, 
                             restriction_decreasing, nets, n_test, policy_test, contact_tracing_efficiency,
                             quarantine_efficiency, use_probabilities, alpha, gamma, lambdaa)
            nets.append(net.copy())
            if dump_type == "full":
                to_dump["nets"].append(net.copy())
            if dump_type == "light":
                to_dump = update_dump_report(to_dump, net, new_positive_counter, use_probabilities)
    else:
        exposed = n_initial_infected_nodes
        infected = 0
        sim_index = 0
        while((infected + exposed) != 0):
            net, new_positive_counter = step(G, sim_index, incubation_days, infection_duration, transmission_rate,
                             initial_day_restriction, restriction_duration, social_distance_strictness, 
                             restriction_decreasing, nets, n_test, policy_test, contact_tracing_efficiency,
                             quarantine_efficiency, use_probabilities, alpha, gamma, lambdaa)
            nets.append(net.copy())
            sim_index += 1

            if dump_type == "full":
                to_dump["nets"].append(net.copy())
            if dump_type == "light":
                to_dump = update_dump_report(to_dump, net, new_positive_counter, use_probabilities)
            
            network_report = Counter(net.vs["agent_status"])
            infected = network_report["I"]
            exposed = network_report["E"]
            if infected + exposed == 0:
                break

    with open(Path(path + ".pickle"), "wb") as f:
        pickle.dump(to_dump, f, protocol = pickle.DEFAULT_PROTOCOL)
    b = time.perf_counter()
    print("\n Simulation ended successfully \n")
    print("time elapsed " + str(b-a))

def main():

    dump_type = None
    path = None

    user_interaction = int(input("Press 0 to load the default values or 1 to manually input the configuration for the simulation: "))
    # get values from user
    if user_interaction:
        n_of_families = int(input("Please insert the number of families in the simulation: "))
        use_steps = int(input("Press 1 use a fixed number of steps or 0 to run the simulation untill the infection is over: "))
        if use_steps:
            number_of_steps = int(input("Please insert the number of the simulation step: "))
        incubation_days = int(input("Please insert the disease incubation duration: "))
        infection_duration = int(input("Please insert the disease duration: "))
        initial_day_restriction = int(input("Please insert the step index from which the social distance is applied: "))
        restriction_duration = int(input("Please insert the number of days which the social distance last. Insert 0 to make the restriction last for all the simulation: "))
        social_distance_strictness = int(input("Please insert a value between 0 and 4 to set the social distance strictness: "))
        restriction_decreasing = int(input("Press 1 to make the strictness of the social distance decrease during the simulation or 0 to keep it fixed: "))
        n_initial_infected_nodes = int(input("Please insert the number of initial infected individuals: "))
        R_0 = float(input("Please insert the value of R0: "))
        n_test = int(input("Please insert the number of available test per day: "))
        policy_test = input("Please insert strategy with which test are made. Can be Random, Degree Centrality, Betweenness Centrality or PBI (probability of being infected): ")
        contact_tracing_efficiency = float(input("Please insert a value between 0 and 1 to set the contact tracing efficiency: "))
        contact_tracing_duration = int(input("Please insert for how many days the contact tracing is computed: "))
        quarantine_efficiency = float(input("Please insert a value between 0 and 1 to set the quarantine efficiency: "))
        use_random_seed = int(input("Press 1 use a fixed a random seed or 0 to pick a random seed: "))
        if use_random_seed:
            seed = int(input("Please insert the random seed: "))
        use_probabilities = int(input("Press 1 to enable estimation of probabilities of being infected, 0 otherwise: "))
        if use_probabilities:
            alpha = float(input("Please insert the value of alpha: "))
            gamma = float(input("Please insert the value of gamma: "))
            lambdaa = float(input("Please insert the value of lambdaa: "))
        dump_type = input("Please insert the dump type. Can be either full of light: ")
        path = input("Please insert the path with the file to dump. Please omit file type, that will be set automatically: ")

        run_simulation(n_of_families, use_steps, number_of_steps, incubation_days, infection_duration,
            initial_day_restriction, restriction_duration, social_distance_strictness, restriction_decreasing,
            n_initial_infected_nodes, R_0, n_test, policy_test, contact_tracing_efficiency, contact_tracing_duration,
            quarantine_efficiency, use_random_seed, seed, use_probabilities, alpha, gamma, lambdaa, dump_type, path)
    else:
        dump_type = input("Please insert the dump type. Can be either full of light: ")
        path = input("Please insert the path with the file to dump. Please omit file type, that will be set automatically: ")
        run_simulation(path = path, dump_type = dump_type)

if __name__ == "__main__":
    main()