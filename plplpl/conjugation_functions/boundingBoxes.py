# get the bounding boxes for two cells

import numpy as np
import shapely.geometry as sg

# input: endpoints of long axis as [[x1,y1],[x2,y2]], width of short axis
# optionally epsilon to extend the box by (defaults to zero)
# output: list with four corners of a bounding box, ordered for shapely
# (top left, bottom left, bottom right, top right)
def calcBox(endpoints, width, epsilon=0):

    x1 = endpoints[0][0]
    y1 = endpoints[0][1]
    x2 = endpoints[1][0]
    y2 = endpoints[1][1]

    # make sure end1 has the smaller x coordinate
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1

    # calculate the angle relative to x-axis
    # must first check it is not exactly vertical
    if x1 == x2:
        angle = np.pi/2
    else:
        angle = abs(np.arctan((y2-y1)/(x2-x1)))

    # calculate the distance in x,y to each corner from the endpoint
    dx = (width/2 + epsilon) * np.sin(angle)
    dy = (width/2 + epsilon) * np.cos(angle)

    # find the corners - depends on relative position of endpoints
    # order is top left, bottom left, bottom right, top right
    if y1 <= y2:
        c1 = [x1-dx, y1+dy]
        c2 = [x1+dx, y1-dy]
        c3 = [x2+dx, y2-dy]
        c4 = [x2-dx, y2+dy]

    else:
        c1 = [x1+dx, y1+dy]
        c2 = [x1-dx, y1-dy]
        c3 = [x2-dx, y2-dy]
        c4 = [x2+dx, y2+dy]

    return [c1,c2,c3,c4]

# input: endpoints of long axis as [[x1,y1],[x2,y2]], epsilon to extend by
# output: endpoints extended by epsilon along the long axis in same form, with x1 <= x2
def extendDonorEnds(endpoints, epsilon):

    x1 = endpoints[0][0]
    y1 = endpoints[0][1]
    x2 = endpoints[1][0]
    y2 = endpoints[1][1]

    # make sure end1 has the smaller x coordinate
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1

    # calculate the angle relative to x-axis
    # must first check it is not exactly vertical
    if x1 == x2:
        angle = np.pi/2
    else:
        angle = abs(np.arctan((y2-y1)/(x2-x1)))

    # calculate the distance in x,y to move each endpoint
    dx = epsilon * np.cos(angle)
    dy = epsilon * np.sin(angle)

    # find the new endpoints - depends on relative position of endpoints
    if y1 <= y2:
        e1 = [x1-dx,y1-dy]
        e2 = [x2+dx,y2+dy]

    else:
        e1 = [x1-dx,y1+dy]
        e2 = [x2+dx,y2-dy]

    return [e1,e2]

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

# for the edge cell1 --> cell2
# input: endpoints + widths of two cells to compare, epsilon to extend cell1 by, tolerance for min area
# output: area of cell2 in range of cell1
def getContactArea(endpoints1, width1, endpoints2, width2, epsilon=0.3,tol=1e-3):

    # cell 1 is acting like donor, so bound with hexagon
    # get the corners of the rectangular portion, then add the ends
    corners1 = calcBox(endpoints1,width1,epsilon)
    ends = extendDonorEnds(endpoints1,epsilon)
    
    # add the ends to the corners in the correct location to get a hexagon
    corners1.insert(1,ends[0])
    corners1.insert(4,ends[1])

    # cell 2 is acting like recipient, so bound with rectangle
    corners2 = calcBox(endpoints2,width2)
    
    # make the shapely polygons for donor, recipient and check intersecting area
    p1 = sg.Polygon(corners1)
    p2 = sg.Polygon(corners2)

    result = p1.intersection(p2).area
    
    # snap anything below tolerance to 0
    if result < tol:
        return 0
    else:
        return result
    
# for the edge cell1 --> cell2
# input: endpoints + widths of two cells to compare, epsilon to extend cell1 by, tolerance for min area
# output: length of cell2 boundary in range of cell1
def getContactBoundary(endpoints1, width1, endpoints2, width2, epsilon=0.3,tol=1e-3):

    # cell 1 is acting like donor, so bound with hexagon
    # get the corners of the rectangular portion, then add the ends
    corners1 = calcBox(endpoints1,width1,epsilon)
    ends = extendDonorEnds(endpoints1,epsilon)
    
    # add the ends to the corners in the correct location to get a hexagon
    corners1.insert(1,ends[0])
    corners1.insert(4,ends[1])

    # cell 2 is acting like recipient, so bound with rectangle
    corners2 = calcBox(endpoints2,width2)
    
    # make the shapely polygons for donor, recipient and check if they intersect
    p1 = sg.Polygon(corners1)
    p2 = sg.Polygon(corners2)

    boundary = p2.boundary
    print(p1,p2)
    print(boundary)

    result = p1.intersection(boundary).length
    
    if result < 1e-3:
        return 0
    else:
        return result
 
    
