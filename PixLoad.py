import numpy as np

def loadFrame(frame, frameData):
    result = ""
    with open(frameData, "r") as stream:
        result = stream.read()

    dim = [int(a) for a in result.split()]

    nb = []

    with open(frame, 'r') as stream:
        nb = np.fromfile(stream)
    
    #return np.reshape(nb, (dim[2], dim[1], dim[0]))
    """print("p")
    for x in np.nditer(nb):
        if x < 0 or x>1:
            print(x)"""
    nb = nb.astype(np.float32)
    #print(dim[2], dim[1], dim[0])
    return nb, dim

