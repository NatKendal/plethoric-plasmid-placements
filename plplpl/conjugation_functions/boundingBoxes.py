# get the bounding boxes for two cells

import numpy as np

# input: endpoints as [[x1,y1],[x2,y2]], width
# output: four corners of a tight bounding box
def calcBox(endpoints,width):

    x1 = endpoints[0][0]
    y1 = endpoints[0][1]
    x2 = endpoints[1][0]
    y2 = endpoints[1][1]

    # make sure end1 has the smaller x coordinate
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1

    # calculate the angle between the cells
    angle = abs(np.arctan((y2-y1)/(x2-x1)))

    # calculate the distance in x,y to each corner from the endpoint
    dx = width/2 * np.sin(angle)
    dy = width/2 * np.cos(angle)

    # find the corners - depends on relative position of endpoints
    if y1 <= y2:
        c1 = [x1 - dx, y1 + dy]
        c2 = [x1 + dx, y1 - dy]
        c3 = [x2 - dx, y2 + dy]
        c4 = [x2 + dx, y2 - dy]

    else:
        c1 = [x1 - dx, y1 - dy]
        c2 = [x1 + dx, y1 - dy]
        c3 = [x2 - dx, y2 - dy]
        c4 = [x2 + dx, y2 - dy]

    return [c1,c2,c3,c4]

# need to write a function to modify endpoints of donor box

# also need a function to check if donor/sender?? then calculate boxes for both

# and then point in polygon can go here or separate