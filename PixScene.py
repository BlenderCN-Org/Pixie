from math import radians

class Scene():
    def __init__(self):
        self.camera = None
        self.light = None
        self.volume = None


class Camera():
    def __init__(self, pos, rot, resolution, focalLength, sensorWidth):
        self.pos = pos
        self.rot = (radians(rot[0]), radians(rot[1]), radians(rot[2]))
        self.resolution = resolution
        self.angle = 0
        self.sensorWidth = sensorWidth
        self.focalLength = focalLength


class Volume():
    def __init__(self, pos, scale, dim, grid, densityMul):
        self.pos = pos
        self.scale = scale
        self.dim = dim
        self.grid = grid
        self.bounds = [[pos[a] - scale[a]/2 for a in range(3)],
                       [pos[a] + scale[a]/2 for a in range(3)]]
        self.densityMul = densityMul


class Light():
    def __init__(self, pos, intensity):
        self.pos = pos
        self.intensity = intensity

