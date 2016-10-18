import bpy
import numpy as np


with open("C:/Users/Peter/Documents/blender/Test smoke cache/frameData.txt", 'w') as stream:
    a = bpy.context.active_object
    ds = a.modifiers["Smoke"].domain_settings
    x, y, z = ds.domain_resolution
    if ds.use_high_resolution:
        x *= 1 + ds.amplify
        y *= 1 + ds.amplify
        z *= 1 + ds.amplify
    stream.write(str(x) + " " + str(y) + " " + str(z))

def frameChange(scene):
    
    a = scene.objects["Smoke Domain"]
    ds = a.modifiers["Smoke"].domain_settings
    n = np.array(ds.density_grid)
    with open("C:/Users/Peter/Documents/blender/Test smoke cache/frame{}".format(scene.frame_current), 'wb') as stream:
        n.tofile(stream)
    print("frame {} exported".format(scene.frame_current))
    
    if scene.frame_current == ds.point_cache.frame_end:
        bpy.ops.screen.animation_play()
        bpy.app.handlers.frame_change_post.clear()
        
bpy.app.handlers.frame_change_post.append(frameChange)
bpy.context.scene.frame_set(ds.point_cache.frame_start)
bpy.ops.screen.animation_play()
    
"""nb = []

with open("C:/Users/Peter/Documents/blender/Test smoke cache/frame", 'r') as stream:
    nb = np.fromfile(stream)

print(np.array_equal(n, nb))"""