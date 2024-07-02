from collections import deque
def getFlaggedChildren(raw, forwardLinks, cellUID):
    queue = deque()
    queue.append(cellUID)
    firstFlagged = []
    while queue:
        cell = queue.popleft()
        for child in forwardLinks[cell]:
            if raw[child]["flag"] == 1:
                firstFlagged.append(child)
            else:
                queue.append(child)
    return firstFlagged
