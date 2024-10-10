import matplotlib.pyplot as plt
from Bio import AlignIO
import networkx as nx
import math
import numpy as np

def visualize_rna_graph(dot_bracket, sequence, ax, title=""):
    G = nx.Graph()
    base_colors = {
        'A': '#FFA07A', 'C': '#98FB98', 'G': '#87CEFA', 'U': '#DDA0DD', '-': '#757575'
    }
    
    for i, (dot, base) in enumerate(zip(dot_bracket, sequence)):
        label = f"{base}{i+1}" if base != "-" else str(i+1)
        G.add_node(i, base=base, label=label)
    
    for i in range(len(sequence) - 1):
        G.add_edge(i, i+1, color='gray', weight=1, style='solid')
    
    stack = []
    for i, dot in enumerate(dot_bracket):
        if dot == '(':
            stack.append(i)
        elif dot == ')' and stack:
            j = stack.pop()
            if sequence[i] == '-' and sequence[j] == '-':
                G.add_edge(i, j, color='#FFFFFF', weight=2, style='dashed')
            else:
                G.add_edge(i, j, color='#8BDDBB', weight=2, style='solid')
    
    pos = nx.kamada_kawai_layout(G)
    
    node_colors = [base_colors[data['base']] for _, data in G.nodes(data=True)]
    nx.draw(G, pos, ax=ax, with_labels=False, node_color=node_colors, 
            node_size=500, font_size=6, font_weight='bold', linewidths=1.5, edgecolors='black')
    
    labels = nx.get_node_attributes(G, 'label')
    nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=8, font_weight='bold')
    
    # Draw edges with different styles
    solid_edges = [(u, v) for (u, v, d) in G.edges(data=True) if d['style'] == 'solid']
    dashed_edges = [(u, v) for (u, v, d) in G.edges(data=True) if d['style'] == 'dashed']
    
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=solid_edges, 
                           edge_color=[G[u][v]['color'] for u, v in solid_edges], 
                           width=[G[u][v]['weight'] for u, v in solid_edges])
    
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=dashed_edges, 
                           edge_color=[G[u][v]['color'] for u, v in dashed_edges], 
                           width=[G[u][v]['weight'] for u, v in dashed_edges],
                           style='dashed')
    
    ax.axis('off')
    ax.text(0.5, -0.1, title, transform=ax.transAxes, ha='center', va='center', fontsize=20, wrap=True)

def visualize_multiple_rna_graphs(dot_bracket, mapping, image_title=""):
    n = len(mapping)
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)

    fig, axs = plt.subplots(rows, cols, figsize=(10*cols, 6*rows), dpi=300)
    fig.suptitle(image_title, fontsize=40, y=0.98)

    axs = axs.flatten() if isinstance(axs, np.ndarray) else [axs]

    for i, (title, sequence) in enumerate(mapping.items()):
        visualize_rna_graph(dot_bracket, sequence, axs[i], title)

    for i in range(n, len(axs)):
        fig.delaxes(axs[i])

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    # plt.figtext(0.5, 0.01, f"Dot-bracket: {dot_bracket}", ha="center", fontsize=12, 
    #             bbox={"facecolor":"white", "alpha":0.5, "pad":5})
    
    if not title:
        title = 'rna_structures'
    
    plt.savefig(f'{image_title}.png', dpi=300, bbox_inches='tight')
    # plt.show()

def parse_stockholm(file_path):
    alignment = AlignIO.read(file_path, "stockholm")

    # Dictionary to store sequences with IDs as keys
    sequences = {}

    # Iterate over each sequence record in the alignment
    for record in alignment:
        sequences[record.id] = str(record.seq)
    
    ss_cons = alignment.column_annotations["secondary_structure"]
    
    return sequences, ss_cons

def show_stockholm(file_path, image_title=""):
    sequences, dot_bracket = parse_stockholm(file_path)
    if not image_title:
        image_title = file_path.replace("stockholm", "png")
        
    visualize_multiple_rna_graphs(dot_bracket, sequences, image_title=image_title)
    

if __name__ == "__main__":

    # Call the function to visualize all RNA graphs
    show_stockholm("/home/sumon/repos/aes_db/cm_builder_storage/results/aes_30/msa.stockholm")