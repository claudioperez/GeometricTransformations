import opensees.openseespy as ops

# Set up the model
ops.wipe()
ops.model("basic", "-ndm", 3)

ops.node(1, 0, 0, 0)
ops.node(2, 0, 5, 0)
ops.node(3, 0, 5, 5)

fixed = 6 * [1]
ops.fix(1, *fixed)

# Define geometry transforms
ops.geomTransf("Linear", 1, (-1, 0, 0))
ops.geomTransf("Linear", 2, ( 0, 0, 1))             # <---- note subtly wrong vecxz


ops.section("Elastic", 1, E=29e9, A=0.18, Iz=0.005, Iy=0.005, G=12e9, J=0.01)
# define elements

element = "ForceFrame"
ops.element(element, 1, (1, 2), section=1, transform=1)
ops.element(element, 2, (2, 3), section=1, transform=2)   # <-- Raises exception

