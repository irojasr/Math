from sage.all import *

def test():
    G = DiGraph(multiedges=True, loops=True)
    weights = [3, 1, -1, -3]
    G.add_vertices(weights)
    for w in weights:
        G.add_edge(w, w, 'H')
        if w + 2 in weights:
            G.add_edge(w, w+2, 'E')
        if w - 2 in weights:
            G.add_edge(w, w-2, 'F')

    pos = {w: (w, 0) for w in weights}
    try:
        # Looking at SageMath DiGraph plot parameters
        p = G.plot(pos=pos, edge_labels=True, color_by_label={'E':'green', 'F':'red', 'H':'blue'})
        print("Plotted successfully. Dir elements:", dir(p))
    except Exception as e:
        print("Error:", e)

test()
