
import veux
import xara
import numpy as np
import matplotlib.pyplot as plt
from veux.motion import Motion


def create_prism(length:    float,
                 element:   str,
                 section:   dict,
                 boundary:  tuple,
                 orient:    tuple  = (0, 1, 1),
                 transform: str = None,
                 divisions: int = 1,
                 rotation = None):

    L  = length

    # Number of elements discretizing the column
    ne = divisions

    nn = ne + 1

    model = xara.Model(ndm=3, ndf=6)

    for i in range(1, nn+1):
        x = (i-1)/float(ne)*L
        location = (x, x/15, -x/15)

        if rotation is not None:
            location = tuple(rotation@location)

        model.node(i, location)


    # Define boundary conditions
    model.fix( 1, boundary[0])
    model.fix(nn, boundary[1])

    #
    # Define cross-section 
    #
    sec_tag = 1
    properties = []
    for k,v in section.items():
        properties.append("-" + k)
        properties.append(v)

    model.section("FrameElastic", sec_tag, *properties)

    # Define geometric transformation
    geo_tag = 1

    if rotation is not None:
        orient = tuple(map(float, rotation@orient))

    model.geomTransf(transform, geo_tag, *orient)


    # Define elements
    for i in range(1, ne+1):
        model.element(element, i, (i, i+1),
                    section=sec_tag,
                    shear=1,
                    transform=geo_tag)


    model.pattern("Plain", 1, "Linear", load={
        ne+1: [0, 0, 0,   0, 0, 0]
    })
    return model


def create_prism_openseespy(length:    float,
                 element:   str,
                 section:   dict,
                 boundary:  tuple,
                 orient:    tuple  = (0, 1, 1),
                 divisions: int = 1,
                 rotation = None):

    L  = length

    # Number of elements discretizing the column
    ne = divisions

    nn = ne + 1

    import openseespy.opensees as ops
    model = ops 
    model.model("-ndm", 3, "-ndf", 6)

    for i in range(1, nn+1):
        x = (i-1)/float(ne)*L
        location = (x, x/15, -x/15)

        if rotation is not None:
            location = tuple(rotation@location)

        model.node(i, *location)


    # Define boundary conditions
    model.fix( 1, *boundary[0])
    model.fix(nn, *boundary[1])

    #
    # Define cross-section 
    #
    sec_tag = 1
    properties = []
    for k,v in section.items():
        properties.append("-" + k)
        properties.append(v)

    model.section("Elastic", sec_tag, *properties)

    # Define geometric transformation

    if rotation is not None:
        orient = tuple(map(float, rotation@orient))

    model.geomTransf("Corotational", 1, *orient)


    # Define elements
    for i in range(1, ne+1):
        model.element(element, i, (i, i+1),
                      section=sec_tag,
                      transform=1,
                      shear=1)


    model.timeSeries("Linear", 1)
    model.pattern("Plain", 1,1)
    model.load( ne+1, 0, 0, 0, 0, 0, 0)
    return model


def analyze_moment(model, steps=1, ne=5):

    model.system("BandGen")
    model.test("EnergyIncr", 1e-10, 10, 1)
    model.integrator("DisplacementControl", ne+1, 2, 0.1)
    model.analysis("Static")


    artist = veux.create_artist(model)
    artist.draw_axes()
    artist.draw_outlines()
    motion = Motion(artist)

    u = []
    for i in range(steps):

        if model.analyze(1) != 0:
            break

        motion.draw_sections(position=model.nodeDisp,
                             rotation=getattr(model, "nodeRotation", None))
        
        motion.advance(time=i/10)
        u.append(np.linalg.norm(model.nodeDisp(ne+1,4))/length)

    motion.add_to(artist.canvas)
    return model,u,artist


if __name__ == "__main__":

    E  = 10000
    I  = 200.0
    length = 100
    model = create_prism(
        length = length,
        section = dict(
            E   = E,
            G   = E/2.0,
            A   = 20.0,
            J   = 20.0,
            Iy  = I,
            Iz  = I,
            Ay  = 20.0,
            Az  = 20.0),
        boundary = ((1,1,1,  1,1,1),
                    (0,0,0,  0,0,0)),
        divisions=3,
        transform="Corotational",
        element="PrismFrame"
    )

    m,u,artist = analyze_moment(model, steps=200)
    
    for node in m.getNodeTags():
        print(f"Node {node}: {np.linalg.norm(m.nodeDisp(node))}")

    
    veux.serve(artist)
    plt.plot(u,'.')
    plt.show()

