import pathlib

from plplpl.base_functions import BaseConjugationFunction

# determining the quality of contact between cells

# input: list of all cells, dictionary of all neighbours, later will include other info
# output: dictionary mapping cell pairs to the weight of the edge from cell1 to cell2 (cell1,cell2)
def getContactWeights(all_cells, dictNeighbours):

    contact_weights = dict()

    # iterate through each cell
    for cell1 in all_cells:
        # and check the weight for all of its neighbors
        for cell2 in dictNeighbours[cell1]:
            contact_weights[(cell1,cell2)] = calcWeight(cell1,cell2)

    return contact_weights

# input: two cells, later will include other info
# output: weight of edge from first cell to second cell
def calcWeight(cell1,cell2):
    return 1/1000

class ContactWeightsFixedNaive(BaseConjugationFunction):
    def __init__(self):
        super().__init__("contactWeightsFixedNaive", 0, ["g"], [])

    def weight(self, cell1, cell2, debug=0, **kwargs):
        return 1/1000
