import pickle

from plplpl.NoisyOr import NoisyOrBayesianNetwork


"""
Returns a NoisyOrBayesianNetwork

modelFolder: path to directory
dataFolder: path to directory
modelname: unique name for this model/computation
colour_min: minimum number of timesteps for colour to appear
colour_max: maximum number of timesteps for colour to appear
maturation_min: minimum number of timesteps for maturation 
maturation_max: maximum number of timesteps for maturation
save: if we should save the model to a file (pickle)
"""

#c_min = 10
#c_max = 36
#m_min = 4
#m_max = 29
def setupGraph(modelFolder, dataFolder, modelname, colour_min, colour_max, maturation_min, maturation_max, save=True, debug=False):
    # Initialize a NoisyOrBayesianNetwork.
    model = NoisyOrBayesianNetwork()

    # Get list of human friendly names to use as variable names.
    with open(dataFolder + modelname + "_humanFriendlyName.pickle", "rb") as f:
        name = pickle.load(f)

    # Save the list of unique cell ids for later.
    cells = list(name.keys())

    # Add the nodes to the graph.
    for cell in cells:
        # Add three nodes per cell per timestep.
        # g[cellId]_[step] is 1 if cellId has the gene at timestep step, else 0.
        # c[cellId]_[step] is 1 if cellId is coloured (not green) at timestep step, else 0.
        # m[cellId]_[step] is 1 if cellId is matured (can conjugate) at timestep step, else 0.
        model.add_nodes_from(["g"+name[cell], "c"+name[cell], "m"+name[cell]])
        if debug:
            print("Added cell " + name[cell])

    # Add edges to the graph.
    # First: get the children of a cell, either true children or itself at the next timestep.
    # We have to add edges one timestep forward for the three properties, gene, colour, and maturation.
    with open(dataFolder + modelname + "_forwardLinks.pickle", "rb") as f:
        forwardLinks = pickle.load(f)
    # Second: get the neighbours of a cell.
    # We have to add edges one timestep forward for each possible mate.
    with open(dataFolder + modelname + "_neighbours.pickle", "rb") as f:
        neighbours = pickle.load(f)

    # Add the edges, one cell at a time.
    for cell in cells:
        # First: add the simple downward edges. (g -> g, c -> c, m -> m)
        for child in forwardLinks[cell]:
            model.add_edge("g"+name[cell], "g"+name[child])
            model.add_edge("c"+name[cell], "c"+name[child])
            model.add_edge("m"+name[cell], "m"+name[child])
        # Second: add inter-cell edges to neighbours. (m -> g)
        for neighbour in neighbours[cell]:
            model.add_edge("m"+name[cell], "g"+name[neighbour])
        # Third: add the intra-cell edges to children. (g -> c, g -> m)
        children = forwardLinks[cell].copy()
        depth = 1
        while depth <= max(colour_max, maturation_max):
            # If we're in the correct range for colour, add edges.
            if depth >= colour_min and depth <= colour_max:
                for child in children:
                    model.add_edge("g"+name[cell], "c"+name[child])
            # If we're in the correct range for maturation, add edges.
            if depth >= maturation_min and depth <= maturation_max:
                for child in children:
                    model.add_edge("g"+name[cell], "m"+name[child])
            # In any case, get all of the children's children and repeat until we hit the required depth.
            newChildren = []
            for child in children:
                newChildren.extend(forwardLinks[child])
            children = newChildren
            depth += 1

        if debug:
            print("Cell " + name[cell] + " complete.")

    if save:
        with open(modelFolder + modelName + "_graph.pickle", "wb") as f:
            pickle.dump(model, f)

    return model

