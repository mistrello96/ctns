import igraph as ig
import numpy as np
from pathlib import Path
import sys, random, time
try:
    from ctns.generator import generate_network, init_infection
    from ctns.steps import step
    from ctns.utility import dump_simulation, compute_TR
except ImportError as e:
    from generator import generate_network, init_infection
    from steps import step
    from utility import dump_simulation, compute_TR

def run_simulation(n_of_families = 500,
    use_steps = True,
    number_of_steps = 150,
    incubation_days = 5,
    infection_duration = 21,
    initial_day_restriction = 50,
    restriction_duration = 21,
    social_distance_strictness = 2,
    restriction_decreasing = True,
    n_initial_infected_nodes = 10,
    R_0 = 2.9,
    n_test = 5,
    policy_test = "Random",
    contact_tracking_efficiency = 0.8,
    use_random_seed = None,
    seed = None,
    dump = False,
    dump_type = "full",
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
        Strategy with which test are made. Can be Random, Degree Centrality, Betweenness Centrality

    contact_tracking_efficiency: float
        The percentage of contacts successfully traced back in the past 14 days

    use_random_seed: bool
        Id use or not a fixed random seed

    seed: int
        The random seed value

    dump: bool
        If true dump networks to file, otherwise return the networks list
        NB, if network simulation is too heavy, the dump can crash

    dump_type: string
        Can be either ["full", "light"]. In the first case, the full simulation is dumped.
        Otherwise, only a report about node status is saved
        The dumped file will have the following structure:
        - a list of ig.Graph() if dump_type is full
        - a dict[class] where class can be [S, E, I, R, D, quarantined, positive, tested, total] and value is a list of the corresponding attribute value on day i

    path: string
        The path to the file/folder where the networks will be saved

    Return
    ------
    nets:list
        List of nets of the simulation

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
    if not (policy_test == "Random" or policy_test == "Degree Centrality" or policy_test == "Betweenness Centrality"):
        print("Invalid test strategy")
        sys.exit()
    if contact_tracking_efficiency < 0 or contact_tracking_efficiency > 1:
        print("Invalid contact tracing efficiency")
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
    if dump_type != "full" and dump_type != "light":
        print("Invalid dump type")
        sys.exit()

    config = locals()

    # init network
    G = generate_network(n_of_families)
    transmission_rate = compute_TR(G, R_0, infection_duration, incubation_days)
    init_infection(G, n_initial_infected_nodes)
    nets = list()

    if use_steps:
        for sim_index in range (0, number_of_steps):
            report, net= step(G, sim_index, incubation_days, infection_duration, transmission_rate,
                             initial_day_restriction, restriction_duration, social_distance_strictness, 
                             restriction_decreasing, nets, n_test, policy_test, contact_tracking_efficiency)
            nets.append(net.copy())
    else:
        exposed = n_initial_infected_nodes
        infected = 0
        sim_index = 0
        while((infected + exposed) != 0):
            report, net= step(G, sim_index, incubation_days, infection_duration, transmission_rate,
                             initial_day_restriction, restriction_duration, social_distance_strictness, 
                             restriction_decreasing, nets, n_test, policy_test, contact_tracking_efficiency)
            nets.append(net.copy())
            sim_index += 1
            infected = report["I"]
            exposed = report["E"]

            if infected + exposed == 0:
                break
    if dump:
        dump_simulation(nets, path, dump_type, config)
  
    print("\n Simulation ended successfully \n")

    return nets

def main():

    n_of_families = None
    use_steps = None
    number_of_steps = None
    incubation_days = None
    infection_duration = None
    initial_day_restriction = None
    restriction_duration = None
    social_distance_strictness = None
    restriction_decreasing = None
    n_initial_infected_nodes = None
    R_0 = None
    n_test = None
    policy_test = None
    contact_tracking_efficiency = None
    use_random_seed = None
    seed = None
    dump = True
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
        restriction_decreasing = int(input("Press 0 to make the strictness of the social distance decrease during the simulation or 1 to keep it fixed: "))
        n_initial_infected_nodes = int(input("Please insert the number of initial infected individuals: "))
        R_0 = float(input("Please insert the value of R0: "))
        n_test = int(input("Please insert the number of available test per day: "))
        policy_test = input("Please insert strategy with which test are made. Can be Random, Degree Centrality, Betweenness Centrality: ")
        contact_tracking_efficiency = float(input("Please insert a value between 0 and 1 to set the contact tracing efficiency: "))
        use_random_seed = int(input("Press 1 use a fixed a random seed or 0 to pick a random seed: "))
        if use_random_seed:
            seed = int(input("Please insert the random seed: "))
        dump_type = input("Please insert the dump type. Can be either full of light: ")
        path = input("Please insert the path with the file to dump. Please omit file type, that will be set automatically: ")

        run_simulation(n_of_families, use_steps, number_of_steps, incubation_days, infection_duration,
            initial_day_restriction, restriction_duration, social_distance_strictness, restriction_decreasing,
            n_initial_infected_nodes, R_0, n_test, policy_test, contact_tracking_efficiency, use_random_seed,
            seed, dump, dump_type, path)
    else:
        dump_type = input("Please insert the dump type. Can be either full of light: ")
        path = input("Please insert the path with the file to dump. Please omit file type, that will be set automatically: ")
        run_simulation(path = path, dump = True, dump_type = dump_type)

if __name__ == "__main__":
    main()