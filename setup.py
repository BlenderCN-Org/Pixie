from distutils.core import setup, Extension

RayCastAccel = Extension('RayCastAccel',
                    sources = ['PixRenderAccel.c'])

setup (name = 'RayCastAccel',
       version = '0.0.1',
       description = 'Functions for accelerating raycasting',
       ext_modules = [RayCastAccel])
