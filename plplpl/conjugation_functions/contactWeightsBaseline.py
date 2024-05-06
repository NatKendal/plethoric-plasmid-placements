import pathlib
import shapely.geometry as sg

from plplpl.base_functions import BaseConjugationFunction
from plplpl.conjugation_functions.boundingBoxes import calcBox, extendDonorEnds

# for the edge cell1 --> cell2
# input: endpoints of long axes as [[x1,y1],[x2,y2]], widths of short axes, epsilon to extend cell1 by
# output: 0 if cells do not overlap; 1 if they do
def checkContact(endpoints1, width1, endpoints2, width2, epsilon=0.3):

    # cell 1 is acting like donor, so bound with hexagon
    # get the corners of the rectangular portion, then add the ends
    corners1 = calcBox(endpoints1,width1,epsilon)
    ends = extendDonorEnds(endpoints1,epsilon)
    
    # since order of box corners is top left, bottom left, bottom right, top right
    # and order of ends is left, right
    # first end goes in second position (index 1) and second end goes in fifth position (index 4)
    corners1.insert(1,ends[0])
    corners1.insert(4,ends[1])

    # cell 2 is acting like recipient, so bound with rectangle
    corners2 = calcBox(endpoints2,width2)
    
    # make the shapely polygons for donor, recipient and check if they intersect
    p1 = sg.Polygon(corners1)
    p2 = sg.Polygon(corners2)

    result = p1.intersects(p2)
    
    if result:
        return 1
    else:
        return 0

class ContactWeightsBaseline(BaseConjugationFunction):
    def __init__(self):
        super().__init__("contactWeightsBaseline", 0, ["g"], ["ends", "width"])

    def weight(self, cell1, cell2, debug=0, **kwargs):
        endpoints1 = self.data["ends"][cell1]
        endpoints2 = self.data["ends"][cell2]
        width1 = self.data["width"][cell1]
        width2 = self.data["width"][cell2]

        return checkContact(endpoints1,width1,endpoints2,width2)
