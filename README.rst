
===========
`pyemblite`
===========

.. start long description.

Python wrapper for Embree-3. Source code adapted from
`pyembree <https://github.com/scopatz/pyembree>`_ Embree-2 wrapper.

.. end long description.

Quick Start
===========

Example::

   import numpy as np
   import trimesh
   from trimesh.primitives import Sphere
   from pyemblite.mesh_construction import TriangleMesh
   from pyemblite.rtcore_scene import EmbreeScene

   # Create Embree scene which holds meshes.
   scene = EmbreeScene()

   # Create a mesh using trimesh (https://github.com/mikedh/trimesh).
   tmesh = Sphere(radius=5.0, subdivisions=1)

   # Create Embree triangle mesh geometry
   emesh = TriangleMesh(scene, tmesh.vertices, tmesh.faces)

   # Commit the scene (builds spatial acceleration structures).
   scene.commit()

   # Generate ray origins and ray directions
   ray_orgs = (
       np.zeros((tmesh.vertices.shape[0], 3), dtype=np.float32)
       +
       tmesh.centroid
   ).astype(np.float32)
   ray_dirs = (tmesh.vertices - tmesh.centroid).astype(np.float32)
   ray_dirs /= np.linalg.norm(ray_dirs, axis=1)[np.newaxis, 1]

   # Query the index of the first face which gets hit by ray
   # (index of -1 indicates ray did not hit a face)
   primID = scene.run(ray_orgs, ray_dirs, query='INTERSECT')

   # Query the distance from the ray origin where face which gets hit by ray
   # Intersection points are ray_orgs + tfar * ray_dirs
   tfar = scene.run(ray_orgs, ray_dirs, query='DISTANCE')
   print(tfar)

   # Query all info, intersect_info is a dict with keys:
   # ['u', 'v', 'Ng', 'tfar', 'primID', 'geomID']
   intersect_info = scene.run(ray_orgs, ray_dirs, output=True)


Installation
============

Install from latest github source:

   ``python -m pip install --no-deps --no-build-isolation --user git+git://github.com/AppliedMathematicsANU/pyemblite.git#egg=pyemblite``

If you're on windows, you need to have embree3 installed and in the default location (C:\Program Files\Intel\Embree3\) before you install pyemblite. To complicate matters, the most recent versions no longer come with an installer, the last version to do so is this one: https://github.com/embree/embree/releases/download/v3.13.2/embree-3.13.2.x64.vc14.msi .
If you want to use the most recent version instead, make sure you unzip the contents of the zip file into the aforementioned folder.
You also still need to have build tools installed (some kind of C/C++ compiler). One way to achieve this is to install Visual Studio Build tools. Visual studio build tools likely require the installation of visual studio community edition first. This link should (hopefully) get you started: https://visualstudio.microsoft.com/downloads/

Requirements
============

Requires:

- python-3 version `>= 3.4`,
- `numpy <http://www.numpy.org/>`_ version `>= 1.7`,
- `embree <https://embree.github.io>`_ `>= 3.0` (`Latest release <https://github.com/embree/embree/releases/latest>`_)


Testing
=======

Run tests (unit-tests and doctest module docstring tests) using::

   python -m pyemblite.test


Latest source code
==================

Source at github:

   https://github.com/AppliedMathematicsANU/pyemblite


License information
===================

See the file `LICENSE.txt <https://github.com/AppliedMathematicsANU/pyemblite/blob/dev/LICENSE.txt>`_
for terms & conditions, for usage and a DISCLAIMER OF ALL WARRANTIES.

