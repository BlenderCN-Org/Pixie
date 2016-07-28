from PIL import Image
import numpy as np
from PixLoad import loadFrame
from PixScene import Scene, Camera, Volume, Light
from PixRender import renderScene

import time
import math

l = Light((.5, .5, 1), 1)
s = Scene()

s.light = l

x = 0
y = -1
z = -.1

downScale = 1

startFrame = 1
lastFrame = 30

t = time.time()
for f in range(startFrame, lastFrame+1, 1):
    print(f)
    dat, dim = loadFrame("/media/sf_Test_smoke_cache/frame{}".format(f), "/media/sf_Test_smoke_cache/frameData.txt")
    # dim in format z, y, x
    v = Volume((0, 0, 0), (.5, .5, .5), dim, dat, .5)
    s.volume = v
    
    #x = math.sin(math.pi*f/30)
    #y = -math.cos(math.pi*f/30)
    c = Camera((x, y, z), (0, 180, 0), (1920//downScale, 1080//downScale), 30/1000, 36/1000) # (... focal length, sensor width)
    s.camera = c
    im = renderScene(s)
    #im.transpose(Image.FLIP_TOP_BOTTOM)
    #big = im.resize((500, int(500*im.size[1]/im.size[0])), Image.NEAREST)
    im.save("testRender/out{}.png".format(f))

totTime = time.time() - t
print("Completed " + str(lastFrame-startFrame+1) + " in " + str(round(totTime, 1)) + " seconds. " 
      + str(round(totTime/(lastFrame-startFrame+1), 1)) + " seconds per frame.")

#im = renderScene(s)
#print(time.time() - t)
#im = Image.new("RGB", (dim[2], dim[1]))

#big = im.resize((500, int(500*im.size[1]/im.size[0])), Image.NEAREST)
# Image.ANTIALIAS
#big.show()

