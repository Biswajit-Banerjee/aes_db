import matplotlib.pyplot as plt
import networkx as nx
import math
import numpy as np
import json 

def visualize_rna_graph(dot_bracket, sequence, ax, title=""):
    G = nx.Graph()
    base_colors = {
        'A': '#FFA07A', 'C': '#98FB98', 'G': '#87CEFA', 'U': '#DDA0DD', '-': '#757575'
    }
    
    for i, (dot, base) in enumerate(zip(dot_bracket, sequence)):
        label = f"{base}{i+1}"
        G.add_node(i, base=base, label=label)
    
    for i in range(len(sequence) - 1):
        G.add_edge(i, i+1, color='gray', weight=1, style='solid')
    
    stack = []
    for i, dot in enumerate(dot_bracket):
        if dot == '(':
            stack.append(i)
        elif dot == ')' and stack:
            j = stack.pop()
            G.add_edge(i, j, color='#8BDDBB', weight=2, style='solid')
    
    pos = nx.kamada_kawai_layout(G)
    
    node_colors = [base_colors[data['base']] for _, data in G.nodes(data=True)]
    nx.draw(G, pos, ax=ax, with_labels=False, node_color=node_colors, 
            node_size=300, font_size=6, font_weight='bold', linewidths=1, edgecolors='black')
    
    labels = nx.get_node_attributes(G, 'label')
    nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=6, font_weight='bold')
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=[G[u][v]['color'] for u, v in G.edges()], 
                           width=[G[u][v]['weight'] for u, v in G.edges()])
    
    ax.axis('off')
    ax.set_title(title, fontsize=16)

def get_base_pairs(dot_bracket):
    stack = []
    base_pairs = set()
    for i, dot in enumerate(dot_bracket):
        if dot == '(':
            stack.append(i)
        elif dot == ')' and stack:
            j = stack.pop()
            base_pairs.add((j, i))
    return base_pairs

def combine_structures(structure1, structure2, mode='union'):
    bp1 = get_base_pairs(structure1)
    bp2 = get_base_pairs(structure2)
    
    if mode == 'union':
        combined_bp = bp1.union(bp2)
    elif mode == 'intersection':
        combined_bp = bp1.intersection(bp2)
    elif mode == 'fr3d_backbone':
        combined_bp = bp1.union(bp2 - bp1)
    elif mode == 'ipknot_backbone':
        combined_bp = bp2.union(bp1 - bp2)
    
    result = ['.' for _ in range(len(structure1))]
    for i, j in combined_bp:
        result[i] = '('
        result[j] = ')'
    return ''.join(result)

def visualize_rna_structures(data, storage_path):
    fig, axes = plt.subplots(2, 3, figsize=(30, 20))
    axes = axes.flatten()
    
    fr3d = data['FR3D']
    ipknot = data['IPKNOT']
    sequence = data['SEQUENCE']
    
    structures = [
        (fr3d, "FR3D Structure"),
        (ipknot, "IPKNOT Structure"),
        (combine_structures(fr3d, ipknot, 'union'), "FR3D ∪ IPKNOT"),
        (combine_structures(fr3d, ipknot, 'intersection'), "FR3D ∩ IPKNOT"),
        (combine_structures(fr3d, ipknot, 'fr3d_backbone'), "FR3D backbone + IPKNOT"),
        (combine_structures(fr3d, ipknot, 'ipknot_backbone'), "IPKNOT backbone + FR3D")
    ]
    
    for (structure, title), ax in zip(structures, axes):
        visualize_rna_graph(structure, sequence, ax, title)
    
    fig.suptitle(data['TITLE'], fontsize=24)
    plt.tight_layout()
    
    plt.savefig(storage_path, dpi=300)
    plt.close()


if __name__ == "__main__":
    from tqdm.auto import tqdm
    import os
    
    with open("test_1.json") as f:
        all_data = json.load(f)
    
    for data in tqdm(all_data):
        aes = data["TITLE"].replace("AES ", "")
        storage_path = os.path.join("cm_builder_storage", "results", f"aes_{aes}", "plot.png")
        os.makedirs(os.path.dirname(storage_path), exist_ok=True)
        visualize_rna_structures(data, storage_path)