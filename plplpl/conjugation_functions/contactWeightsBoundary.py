import pathlib
import shapely.geometry as sg

from plplpl.base_functions import BaseConjugationFunction
from plplpl.conjugation_functions.boundingBoxes import calcBox, extendDonorEnds

# for the edge cell1 --> cell2
# input: endpoints + widths of two cells to compare, epsilon to extend cell1 by, tolerance for min area
# output: length of cell2 boundary in range of cell1
def getContactBoundary(endpoints1, width1, endpoints2, width2, epsilon=0.4,tol=1e-3):

    # cell 1 is acting like donor, so bound with hexagon
    # get the corners of the rectangular portion, then add the ends
    corners1 = calcBox(endpoints1,width1,epsilon)
    ends = extendDonorEnds(endpoints1,epsilon)
    
    # add the ends to the corners in the correct location to get a hexagon
    corners1.insert(1,ends[0])
    corners1.insert(4,ends[1])

    # cell 2 is acting like recipient, so bound with rectangle
    corners2 = calcBox(endpoints2,width2,0)
    
    # make the shapely polygons for donor, recipient and check if they intersect
    p1 = sg.Polygon(corners1)
    p2 = sg.Polygon(corners2)

    boundary = p2.boundary
    #print(p1,p2)
    #print(boundary)

    result = p1.intersection(boundary).length
    
    if result < 1e-3:
        return 0
    else:
        return result

class ContactWeightsBoundary(BaseConjugationFunction):
    def __init__(self):
        super().__init__("contactWeightsBoundary", 0, ["g"], ["ends", "width"])

    def weight(self, cell1, cell2, debug=0, **kwargs):
        endpoints1 = self.data["ends"][cell1]
        endpoints2 = self.data["ends"][cell2]
        width1 = self.data["width"][cell1]
        width2 = self.data["width"][cell2]

        return getContactBoundary(endpoints1,width1,endpoints2,width2)
