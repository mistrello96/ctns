# CTNS

Contact Tracing Network Simulator, a tool to simulate digital contact networks beetween people in a population where a disease is spreading.
The simulation is highly customizable and will return (or dump) a time series of networks. We suggest to use (not jet implemented) to analize the results.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Installing

To install this module you simply need to run:

```
pip install ctns
```

on you terminal.

You can also clone this repo and run the contact_network_simulator file in the ctns folder to get an interactive window to customize and run your simulation.

### Usage

To run a simulation, you first need to import this library in your python code:

```
from ctns.contact_network_simulator import run_simulation
```

and then simply call the function run_simulation.

If you would like to launch the tool directly from Terminal/CMD, you can just type

```
ctns
```

from Linux terminal or

```
ctns.bat
```

from Windows CMD/Power Shell.

Remember to specify a path to a file for the network dump if needed.

You can alternatively clone the repo, navigate to the ctns/ctns folder and run

```
python3 contact_network_simulator.py
```

## Network structure

The networks are igraph networks that will have nodes (representing people) and weighed edges (representing contacts).
Each node will have the following relevant attributes:
- sex -> either man or woman
- age -> a number representing the age slice of the node (e.g., 40 represent a node with age between [40,50))
- family_id -> id of the family of the node
- sociability -> a value that can be either low, medium or high and represent how social (average number of edge) a node is
- pre_existing_conditions -> number of existing pathologies of the node
- agent_status -> the status of the node in the simulation. It can be either S(susceptible), E(Exposed), I(Infective), R(Recovered) or D(Dead)
- infected -> if the node is currently infected or not
- quarantine -> 0 if the node is not in quarantine, an int representing how many days of quarantine remains otherwise
- test_result -> -1 if no test has been made on the node, 0 if the test result is negative, 1 if the test result is positive
- symptoms -> a list of symptoms shown by the node
- need_IC -> if the node needs or not treatment in intensive care

Please note that you will find other attributes used in the simulation process.

Each edge will have a weight, representing the duration of the contact between two people and a category, that can be either family, frequent, occasional or random.