from collections import defaultdict
import joblib
import pandas as pd
import networkx as nx


def get_pairs(data):
    data.reset_index(inplace=True)
    # Separate out the variants from the base VCF and add new columns of the base/comp ids
    base = data[data['state'].isin(['tpbase'])].copy()
    base['base_id'] = base['MatchId'].apply(lambda x: x[0])
    base['comp_id'] = base['MatchId'].apply(lambda x: x[1])

    # Separate out the variants from the comparison VCF and add new columns of the base/comp ids
    comp = data[data['state'].isin(['tp'])].copy()
    comp['base_id'] = comp['MatchId'].apply(lambda x: x[0])
    comp['comp_id'] = comp['MatchId'].apply(lambda x: x[1])

    # Merge the base/comparison variants
    combined = pd.merge(base, comp, left_on='base_id', right_on='comp_id', suffixes=('_base', '_comp'))
    ret = combined[["hash_base", "hash_comp"]]
    return ret.values

if True: # Turn off/on making data
    graph = nx.Graph()

    data = joblib.load("kanpig_1kgp_multi/data.jl")
    base = list(data[data['state'].isin(['tpbase', 'fn'])].index)
    print("kanpig total:", len(base))
    for i in base:
        graph.add_node('k_' + i, prog='kanpig')

    comp = list(data[data['state'].isin(['tp', 'fp'])].index)
    print("1kgp total:", len(comp))
    for i in comp:
        graph.add_node('o_' + i, prog='1kgp')

    for i,j in get_pairs(data):
        graph.add_edge('k_' + i, 'o_' + j)

    data = joblib.load("snif_1kgp_multi/data.jl")
    base = list(data[data['state'].isin(['tpbase', 'fn'])].index)
    print("snif total:", len(base))
    for i in base:
        graph.add_node('s_' + i, prog='snif')

    for i,j in get_pairs(data):
        graph.add_edge('s_' + i, 'o_' + j)

    data = joblib.load("snif_kanpig_multi/data.jl")
    for i,j in get_pairs(data):
        graph.add_edge('s_' + i, 'k_' + j)

    joblib.dump(graph, 'graph.jl')

import joblib
import networkx as nx
from collections import defaultdict, Counter

G = joblib.load("graph.jl")

# Dictionaries to track membership in each set
membership = defaultdict(set)

# Populate membership based on edges and node properties
for i, subg in enumerate(nx.connected_components(G)):
    m_members = set()
    for node in subg:
        m_members.add(G.nodes[node]['prog'])
    membership[i] = m_members

print("Total members:", len(membership))
# Count nodes in each category for the Venn diagram
exclusive_A = exclusive_B = exclusive_C = 0
shared_AB = shared_AC = shared_BC = shared_ABC = 0

membership_counts = Counter()

for node, sets in membership.items():
    if sets == {'kanpig'}:
        exclusive_A += 1
    elif sets == {'snif'}:
        exclusive_B += 1
    elif sets == {'1kgp'}:
        exclusive_C += 1
    elif sets == {'kanpig', 'snif'}:
        shared_AB += 1
    elif sets == {'kanpig', '1kgp'}:
        shared_AC += 1
    elif sets == {'snif', '1kgp'}:
        shared_BC += 1
    elif sets == {'kanpig', 'snif', '1kgp'}:
        shared_ABC += 1
    
    for name in sets:
        membership_counts[name] += 1

# Display results
for name in ['kanpig', 'snif', '1kgp']:
    print(f"Total members {name}:", membership_counts[name])

print("Exclusive to kanpig:", exclusive_A)
print("Exclusive to snif:", exclusive_B)
print("Exclusive to 1kgp:", exclusive_C)
print("Shared between kanpig and snif (not 1kgp):", shared_AB)
print("Shared between kanpig and 1kgp (not snif):", shared_AC)
print("Shared between snif and 1kgp (not kanpig):", shared_BC)
print("Shared among kanpig, snif, and 1kgp:", shared_ABC)
