import bluepysnap
import pandas as pd

circuit_path = "/home/dhuruva/projects/ctb-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker/circuit_config.json"
circuit = bluepysnap.Circuit(circuit_path)

edge = circuit.edges['S1nonbarrel_neurons__S1nonbarrel_neurons__chemical']
print("Available properties:", edge.property_names)

# Fetch one synapse to see values
synapses = edge.get(limit=1, properties=list(edge.property_names))
print(synapses.T)
