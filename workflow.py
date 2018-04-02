

import networkx as nx
import matplotlib.pyplot as plt


G = nx.MultiDiGraph()

G.add_edges_from(['P1','P2'])
G.add_edges_from(['P1','P3'])

nx.draw(G)
plt.show()



