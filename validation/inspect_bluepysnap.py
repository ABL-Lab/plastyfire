import bluepysnap

circuit_path = "/home/dhuruva/projects/ctb-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker/circuit_config.json"
circuit = bluepysnap.Circuit(circuit_path)
edge = circuit.edges['S1nonbarrel_neurons__S1nonbarrel_neurons__chemical']

print("Methods of edge object:")
print(dir(edge))

print("\nHelp on pair_edges:")
print(help(edge.pair_edges))

print("\nHelp on afferent_edges:")
print(help(edge.afferent_edges))

