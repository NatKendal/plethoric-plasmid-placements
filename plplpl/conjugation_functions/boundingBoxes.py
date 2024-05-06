# get the bounding boxes for two cells

import numpy as np

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

 
    
