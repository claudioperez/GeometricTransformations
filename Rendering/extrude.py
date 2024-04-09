import sys
import numpy as np
from scipy.spatial.transform import Rotation

import shps.curve
from shps.model.frame import FrameState

class Solution: # FieldState
    location:    iter
    rotation:    iter
    def __init__(self, model, conv_hist=None, iter_hist=None, transform=None):
        self._model = model
        if transform is None:
            transform = np.eye(6)

        self._elements  = model["assembly"]
        self._conv_hist = [
            {
                "U":   {k: transform@s["U"][k]   for k in s["U"]},
                "DU":  {k: transform@s["DU"][k]  for k in s["DU"]},
                "DDU": {k: transform@s["DDU"][k] for k in s["DDU"]},
                "Time": s["Time"]
            } for s in conv_hist
        ]

        self._iter_hist = [
            {
                "U":   {k: transform@s["U"][k]   for k in s["U"]},
                "DU":  {k: transform@s["DU"][k]  for k in s["DU"]},
                "DDU": {k: transform@s["DDU"][k] for k in s["DDU"]},
                "Time": s["Time"]
            } for s in iter_hist
        ]
        self._transform = transform
        self._locations = None
        self._rotations = None
        self._reference = {
              "initial": [FrameState(elem) for elem in model["assembly"].values()]
        }

    def iter_conv(self):
        n = len(self._iter_hist)
        for i,state in enumerate(self._iter_hist):
            if (i+1) < n and self._iter_hist[i+1]["Time"] == state["Time"]:
                continue
            else:
                yield (i,state)

    def displacements(self, time):
        for state in self._conv_hist:
            if state["Time"] == time:
                return state["U"]

#   @abstractmethod
    def timeline(self):
        pass

#   @abstractmethod
    def element_solutions(self,  time, elems=None):
        pass

    def locations(self, time, points=None, elems=None):
        for elem in self.element_solutions(time, elems):
            yield [elem.location(x) for x in X]

    def rotations(self, time, points=None, elems=()):
        for elem in self.element_solutions(time):
            if not elems or elem.name in elems:
#               yield [elem.rotation(x).as_matrix() for x in elem.linspace(points)]
                yield [elem._rotation[i].as_matrix() for i in range(points)]

class SpinSolution(Solution):
    def timeline(self):
        states = iter(self._iter_hist)
        for state in states:
            time  = state["Time"]
            yield time

            try:
                while time == state["Time"]:
                    state = next(states)
            except StopIteration:
                return

    def element_solutions(self, time, elems=None):
        for elem in self._elements.values():

            elem_state = FrameState(elem)

            for i,solution in enumerate(self._iter_hist):
                if solution["Time"] > time: break

                spin       = [solution["DDU"][n][3:] for n in elem["nodes"]]
                elem_state = elem_state.update(axis=spin)

            location = [solution[  "U"][n][:3] for n in elem["nodes"]]

            yield elem_state.update(location=location)

class QuatSolution(Solution):
    total : bool

class AxisSolution(Solution):
    total : bool
    def timeline(self):
        for _,state in self.iter_conv():
            yield state["Time"]

    def element_solutions(self, time, elems=None):
        for elem in self._elements.values():

            elem_state = FrameState(elem)

            for i,state in self.iter_conv():

                axis       = [state["DU"][n][3:] for n in elem["nodes"]]
                elem_state = elem_state.update(axis=axis)

                if state["Time"] == time: break

            location = [state[  "U"][n][:3] for n in elem["nodes"]]

            yield elem_state.update(location=location)


def plot_extruded_frames(renderer, solution=None, time=None, options=None):
    if solution is not None:
        displ = solution.displacements(time)

    sections = sees.get_section_geometries(renderer.model, renderer.config)

    nodes = renderer.model["nodes"]

    coords = []
    triang = []

    I = 0
    for i,el in enumerate(renderer.model["assembly"].values()):
        # TODO: probably better to loop over sections dict rather
        #       than elements
        try:
            section = sections[el["name"]]
        except KeyError:
            if "default_section" in renderer.config:
                section = renderer.config["default_section"]
            else:
                continue

        section = section*options["objects"]["sections"]["scale"]

        N  = len(el["nodes"])

        ne = len(section)
        if solution is not None:
            glob_displ = np.array([
                u for n in el["nodes"] #[el["nodes"][0], el["nodes"][-1]]
                    for u in displ[nodes[n]["name"]]
            ])
            nn = N # number of nodes
            X = shps.curve.displace(el["crd"], glob_displ.reshape(nn,6)[:,:3], N)
            X = X.T
            R = next(solution.rotations(time, N, (el["name"],)))
        else:
            section = section*0.98
            X = np.array(el["crd"])
            R = [sees.orientation(el["crd"], el["trsfm"]["yvec"]).T]*N

        # loop over sample points along element length to assemble
        # `coord` and `triang` arrays
        for j in range(N):
            # loop over section edges
            for k,edge in enumerate(section):
                # append rotated section coordinates to list of coordinates
                coords.append(X[j, :] + R[j]@[0, *edge])

                if j == 0:
                    # skip the first section
                    continue

                elif k < ne-1:
                    triang.extend([
                        # tie two triangles to this edge
                        [I+    ne*j + k,   I+    ne*j + k + 1,    I+ne*(j-1) + k],
                        [I+ne*j + k + 1,   I+ne*(j-1) + k + 1,    I+ne*(j-1) + k]
                    ])
                else:
                    # elif j < N-1:
                    triang.extend([
                        [I+    ne*j + k,    I + ne*j , I+ne*(j-1) + k],
                        [      I + ne*j, I + ne*(j-1), I+ne*(j-1) + k]
                    ])

        I += N*ne

    renderer.canvas.plot_mesh(coords, triang,
                              color   = "gray" if solution is not None else "white",
                              opacity = None   if solution is not None else 0.0
                             )

    show_edges = True

    if not show_edges:
        return

    IDX = np.array((
        (0, 2),
        (0, 1)
    ))

    nan = np.zeros(renderer.ndm)*np.nan
    coords = np.array(coords)
    if "extrude.sections" in options["show_objects"]:
        tri_points = np.array([
            coords[idx]  if (j+1)%3 else nan for j,idx in enumerate(np.array(triang).reshape(-1))
        ])
    else:
        tri_points = np.array([
            coords[i]  if j%2 else nan for j,idx in enumerate(np.array(triang)) for i in idx[IDX[j%2]]
        ])

    renderer.canvas.plot_lines(tri_points,
                               color="black" if solution is not None else "#808080",
                               width=4)

class so3:
    @classmethod
    def exp(cls, vect):
        return Rotation.from_rotvec(vect).as_matrix()

def _add_moment(renderer, loc=None, axis=None):
    import meshio
    loc = [1.0, 0.0, 0.0]
    axis = [0, np.pi/2, 0]
    mesh_data = meshio.read('chrystals_moment.stl')
    coords = mesh_data.points

    coords = np.einsum('ik, kj -> ij',  coords,
                       so3.exp([0, 0, -np.pi/4])@so3.exp(axis))
    coords = 1e-3*coords + loc
#   for node in coords:
#       node = so3.exp(axis)@node
    for i in mesh_data.cells:
        if i.type == "triangle":
            triangles =  i.data #mesh_data.cells['triangle']

    renderer.canvas.plot_mesh(coords, triangles)

def _add_solutions(renderer, res_file=None, rparam=None, **opts):
    soln = sees.read_displacements(res_file)


    if "IterationHistory" in soln:
        if "iter" in opts["show_objects"]:
            solution = SpinSolution(renderer.model,
                                    transform=renderer.dofs2plot,
                                    conv_hist=soln["ConvergedHistory"],
                                    iter_hist=soln["IterationHistory"])
        if "incr" in opts["show_objects"]:
            solution = AxisSolution(renderer.model,
                                    transform=renderer.dofs2plot,
                                    conv_hist=soln["ConvergedHistory"],
                                    iter_hist=soln["IterationHistory"])


    time = list(solution.timeline())[-1] if "time" not in opts else int(opts["time"])
#   renderer.plot_chords(renderer.model["assembly"], displ = solution.displacements(time))
    plot_extruded_frames(renderer, solution, time, opts)
#   renderer.canvas.fig.update_traces(contours_z=dict(show=True, usecolormap=True,
#                                 highlightcolor="limegreen", project_z=True))



import sees
from sees import RenderError, read_model
def _render(sam_file, res_file=None, noshow=False, **opts):
    # Configuration is determined by successively layering
    # from sources with the following priorities:
    #      defaults < file configs < kwds 

    config = sees.config.Config()


    if sam_file is None:
        raise RenderError("ERROR -- expected positional argument <sam-file>")

    # Read and clean model
    if not isinstance(sam_file, dict):
        model = sees.read_model(sam_file)
    else:
        model = sam_file

    if "RendererConfiguration" in model:
        sees.apply_config(model["RendererConfiguration"], config)

    sees.apply_config(opts, config)

    renderer = sees.SkeletalRenderer(model, **config)

    plot_extruded_frames(renderer, options=opts)
    # -----------------------------------------------------------
    _add_solutions(renderer, sam_file, **opts)
    # -----------------------------------------------------------
    _add_moment(renderer)
    # -----------------------------------------------------------

    camera = dict(
      up=dict(x=0, y=0, z=1),
      center=dict(x=0, y=0, z=0),
      eye=dict(x=1.25, y=1.25, z=1.25)
    )

    #fig.update_layout(scene_camera=camera, title=name)

    # write plot to file if file name provided
    if config["write_file"]:
        renderer.plot()
        renderer.write(config["write_file"])

    else:
        renderer.plot()
        if not noshow:
            renderer.canvas.show()
        # renderer.repl()

    return renderer


if __name__ == "__main__":
    import sees.__main__
    config = sees.__main__.parse_args(sys.argv)

    try:
        _render(**config)

    except (FileNotFoundError, RenderError) as e:
        # Catch expected errors to avoid printing an ugly/unnecessary stack trace.
        print(e, file=sys.stderr)
        print("         Run '{NAME} --help' for more information".format(NAME=sys.argv[0]), file=sys.stderr)
        sys.exit()

