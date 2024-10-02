import math

from collections import defaultdict

"""
We are no longer using pgmpy and networkx.
The benefits of those libraries are not worth the additional complexity.
"""
class BayesianNetwork(object):
    """
    Low overhead Bayesian network object for use in plplpl v2.

    Optional Arguments:
        vertices: converted to set.
        edges: set of tuples of directed edges (v1, v2). 
        cpd: dict of vertices to CPDs
        constants: dict (arbitrary)
    """
    def __init__(self, vertices=set(), edges=set(), cpds=dict(), constants=dict()):
        # Initialize vertices.
        try:
            self._vertices = set(vertices)
        except Exception as e:
            raise TypeError("Tried to interpret vertices as a set and failed.") from e

        # Initialize edges.
        self._children = defaultdict(set)
        self._parents = defaultdict(set)
        self._edges = set()
        for edge in edges:
            if (not isinstance(edge, tuple)) or (len(edge) != 2):
                raise ValueError("Expected edges to be tuples of size 2.")
            if (edge[0] not in self._vertices) or (edge[1] not in self._vertices):
                raise ValueError("Edge " + str(edge) + " references a vertex not in the graph.")
            self._edges.add(edge)
            self._children[edge[0]].add(edge[1])
            self._parents[edge[1]].add(edge[0])

        # Initialize CPDs.
        if not isinstance(cpds, dict):
            raise TypeError("Expected cpds to be dict but got " + str(type(cpds)) + ".")
        self._cpds = cpds

        # Initialize constants
        if not isinstance(constants, dict):
            raise TypeError("Expected constants to be dict but got " + str(type(constants)) + ".")
        self.constants = constants

    def __iter__(self):
        return self._vertices.__iter__()

    def __len__(self):
        return len(self._vertices)

    #
    # Variable Access
    #

    # Get vertices.
    def vertices(self):
        return self._vertices

    # Add a new vertex.
    def add_vertex(self, vertex):
        self._vertices.add(vertex)

    # Remove a vertex.
    # If soft=True, then check that there are no edges or cpds to lose.
    def remove_vertex(self, vertex, soft=False):
        # Check if very safe to remove if soft=True.
        if soft:
            if len(self._children[vertex]) > 0:
                raise ValueError("Vertex " + str(vertex) + " still had outgoing edges.")
            if len(self._parents[vertex]) > 0:
                raise ValueError("Vertex " + str(vertex) + " still has incoming edges.")
            if vertex in self._cpds:
                raise ValueError("Vertex " + str(vertex) + " still has an associated CPD.")
        
        # Remove outgoing edges.
        for child in list(self._children[vertex]):
            self.remove_edge((vertex, child))
        del self._children[vertex]

        # Remove incoming edges.
        for parent in list(self._parents[vertex]):
            self.remove_edge((parent, vertex))
        del self._parents[vertex]

        # Remove an associated CPD if it exists.
        if vertex in self._cpds:
            del self._cpds[vertex]

    # Get edges.
    def edges(self):
        return self._edges

    # Add a new edge.
    def add_edge(self, edge):
        if edge not in self._edges:
            if (not isinstance(edge, tuple)) or (len(edge) != 2):
                raise ValueError("Expected edges to be tuples of size 2.")
            if (edge[0] not in vertices) or (edge[1] not in vertices):
                raise ValueError("Edge " + str(edge) + " references a vertex not in the graph.")
            self._edges.add(edge)
            self._children[edge[0]].add(edge[1])
            self._parents[edge[1]].add(edge[0])

    # Remove an existing edge.
    def remove_edge(self, edge, fromCPD=False):
        if edge not in self._edges:
            raise ValueError("Edge " + str(edge) + " not in graph.")
        self._edges.remove(edge)
        self._children[edge[0]].remove(edge[1])
        self._parents[edge[1]].remove(edge[0])
        if fromCPD:
            if edge[1] in self._cpds:
                self._cpds[edge[1]].remove_evidence(edge[0])

    # Get all vertices pointed to by a vertex.
    def children(self, vertex):
        return self._children[vertex]

    # Get all vertices pointing to a vertex.
    def parents(self, vertex):
        return self._parents[vertex]

    # Get the direct parent, i.e. the same object one timestep previous.
    def directParent(self, vertex):
        candidates = []
        for parent in self._parents[vertex]:
            if parent[0] == vertex[0]:
                candidates.append(parent)
        if len(candidates) == 1: # If it's unambiguous, return it.
            return candidates[0]
        else: # Otherwise, return None.
            return None

    # Return a simple dictionary of all direct parents.
    def getDirectParents(self):
        directParents = dict()
        for vertex in self._vertices:
            directParents[vertex] = self.directParent(vertex)
        return directParents

    # Get the direct children, i.e. the same object one timestep later. (There could be more than one.)
    def directChildren(self, vertex):
        directChildren = []
        for child in self._children[vertex]:
            if child[0] == vertex[0]:
                directChildren.append(child)
        return directChildren
    
    # Return a simple dictionary of all direct children.
    def getDirectChildren(self):
        directChildren = dict()
        for vertex in self._vertices:
            directChildren[vertex] = self.directChildren(vertex)
        return directChildren

    # Add a CPD for a vertex.
    # If safe, checks the fails if a cpd already exists.
    def add_cpd(self, vertex, cpd, safe=False):
        if safe and (vertex in self._cpds):
            raise ValueError("Vertex " + vertex + " already has a CPD.")
        self._cpds[vertex] = cpd

    # Remove a CPD for a vertex.
    def remove_cpd(self, vertex):
        if vertex in self._cpds:
            del self.cpds[vertex]

    # Return a CPD for a vertex.
    def cpd(self, vertex):
        return self._cpds[vertex]
    def get_cpd(self, vertex):
        return self._cpds[vertex]

    # Return all CPDs.
    def get_cpds(self):
        return self._cpds 

    #
    # Self Checks
    #

    # Ensure that self._edges, self._parents, and self._children all agree and have no invalid vertices.
    # Returns True if no errors, False otherwise.
    # Uses tqdm for progressBar.
    def checkStructure(self, progressBar=False):
        if progressBar:
            import tqdm
        # Check that self._children and self._parents are correct for all vertices.
        if progressBar:
            iterator = tqdm.tqdm(self._vertices)
        else:
            iterator = self._vertices
        edges = set()
        for vertex in iterator:
            if progressBar:
                iterator.set_description("Checking " + str(vertex))
            if vertex in self._children:
                for child in self._children[vertex]:
                    if child not in self._vertices:
                        return False
                    if (vertex, child) not in self._edges:
                        return False
                    if vertex not in self._parents[child]:
                        return False
                    edges.add((vertex, child))
            if vertex in self._parents:
                for parent in self._parents[vertex]:
                    if parent not in self._vertices:
                        return False
                    if (parent, vertex) not in self._edges:
                        return False
                    if vertex not in self._children[parent]:
                        return False

        # If we got here, then every parent/child edge is in edges.
        # So if the length of the ones we found is the same as the original, then the edges look right.
        if len(edges) != len(self._edges):
            return False

        #TODO
        # Check if acyclic

        # Everything looks correct.
        return True

    # Ensure that each vertex has a correct CPD associated.
    def checkCPDs(self, progressBar=False):
        if progressBar:
            import tqdm
        if progressBar:
            iterator = tqdm.tqdm(self._vertices)
        else:
            iterator = self._vertices
        for vertex in iterator:
            if progressBar:
                iterator.set_description("Checking " + str(vertex))
            # Check that all nodes have CPDs
            if vertex not in self._cpds:
                return False
            # Check that the CPD is defined over the vertex associated.
            if self._cpds[vertex]._variable != vertex:
                return False
            # Check that all parents are in evidence.
            if set(self._cpds[vertex]._evidence) != self._parents[vertex]:
                return False
            # Check that the CPD implements get_value.
            if not (hasattr(self._cpds[vertex], "get_value") and callable(self._cpds[vertex].get_value)):
                return False
        return True
    
    def checkCPDWeights(self, progressBar=False):
        if progressBar:
            import tqdm
        if progressBar:
            iterator = tqdm.tqdm(self._vertices)
        else:
            iterator = self._vertices
        for vertex in iterator:
            if progressBar:
                iterator.set_description("Checking " + str(vertex))
            # Check if the weights are appropriate.
            for weight in self._cpds[vertex]._evidence_noise:
                if (weight < 0) or (weight > 1) or (weight in [math.inf, -math.inf, math.nan]):
                    return False
        return True
