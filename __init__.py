from __future__ import annotations
import bmesh
import bpy
from bpy_extras.io_utils import ImportHelper
from bpy.props import BoolProperty, CollectionProperty, FloatProperty, IntProperty
import math
from numbers import Number
import numpy as np
import os
from osgeo import gdal, osr
from xml.etree import ElementTree as ET
import time
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

bl_info = {
    "name": "Import CityGML",
    "author": "Dealga McArdle, ppaawweeuu, Simon Nieswand",
    "version": (0, 5, 0),
    "blender": (2, 80, 0),
    "category": "Import-Export",
    "description": "Import geometry from CityGML file(s)",
    "wiki_url": "https://github.com/ppaawweeuu/Import_CityGML",
}


def get_namespaces(
    filename: str, allow_duplicates: bool = False, exit_early: bool = True
) -> Dict[str, str]:
    """Get namespaces used in a GML/XML file.

    Args:
        filename: Path to the GML/XML file to get the namespaces from.
        allow_duplicates: Whether a namespace is allowed to appear multiple
            times with different URIs. If True, the returned map will only
            contain the URI of the last occurrence of the duplicate namespace.
        exit_early: Whether the look-up of namespaces should be stopped after
            the root element of the GML/XML tree was handled. If False, the
            input file will be parsed fully which can result in long run times
            depending on the file size.

    Returns:
        Dictionary mapping namespaces to URIs.

    Raises:
        KeyError: If @allow_duplicates is False and a namespace is encounterd
            at least twice with different URIs.
    """
    namespaces = {}
    for event, element in ET.iterparse(filename, events=("start-ns", "start")):
        if event == "start-ns":
            namespace, uri = element
            saved_uri = namespaces.get(namespace)
            if saved_uri is not None and saved_uri != uri and not allow_duplicates:
                raise KeyError(
                    f"Duplicate namespace '{namespace}' with different URIs found. "
                    f"Saved URI: {saved_uri}, new URI: {uri}.")
            namespaces[namespace] = uri
        elif event == "start" and exit_early:
            break
    return namespaces


def substitute_uri(
    tag: str, namespaces: Optional[Dict[str, str]] = None
) -> str:
    """Subtitute the namespace prefix in an GML/XML element tag with the
    corresponding URI or prepend the URI of the default namespace if no prefix
    is given.

    Args:
        tag: GML/XML element tag as a string. Usually in the form of
            'namespace_prefix:some_tag'.
        namespaces: A dictionary mapping namespace prefixes to URIs. If None or
            empty, @tag will be returned unchanged. The URI of the default
            namespace should be included with an empty string ('') as the key.

    Returns:
        Element tag with namespace prefix replaced with the corresponding URI.
        Or simply @tag if @namespaces is None or empty.

    Raises:
        KeyError: If the namespace prefix found in @tag is not present as a key
            in @namespaces or if no namespace prefix was detected and @namespaces
            has no entry for the default namespace.
        AttributeError: If multiple colons are found in @tag.
    """
    if namespaces is None or len(namespaces) == 0:
        return tag
    splits = tag.split(":")
    if len(splits) == 1:
        uri = namespaces.get("")
        if uri is None:
            raise KeyError(
                f"Tag '{tag}' doesn't include a namespace prefix and "
                "no default namespace found in provided prefix map.")
    elif len(splits) == 2:
        prefix, tag = splits
        uri = namespaces.get(prefix)
        if uri is None:
            raise KeyError(
                f"Prefix '{prefix}' not found in provided prefix map.")
    else:
        raise AttributeError(f"Multiple colons in tag '{tag}' are invalid.")
    return f"{{{uri}}}{tag}"


def string_to_list(
    string: str, type: Optional[Callable] = None, delimiter: str = " "
) -> List[Any]:
    """Convert a string containing multiple elements into a list of those elements.
    The elements can be converted to a custom type and will be filtered for
    empty / None elements (before the conversion).

    Args:
        string: String containing elements separated by @delimiter.
        type: Callable to convert elements from string by calling type(element: str).
            If None, the elements will not be converted.
        delimiter: String that is used as delimiter between list elements.

    Returns:
        List of elements filtered for empty / None elements and possibly converted.
    """
    splits = string.split(delimiter)
    splits = filter(None, splits)
    if type is None:
        return list(splits)
    return list(map(type, filter(None, splits)))


def unflatten_sequence(
    sequence: Sequence[Any], n_dims: int = 3
) -> List[Tuple[Any]]:
    """Convert a flat (one-dimensional) sequence into a list of tuples.

    Args:
        sequence: One-dimensional sequence. The length of @sequence has to be
            an integer multiple of @n_dims.
        n_dims: Length of the tuples that @sequence is to be divided into.

    Returns:
        Elements of @sequence as a list of tuples of length @n_dims.

    Raises:
        IndexError: If the length of @sequence is not an integer multiple of
            @n_dims.
    """
    return [
        tuple(
            sequence[i + k] for k in range(n_dims)
        ) for i in range(0, len(sequence), n_dims)
    ]


def offset_coords(
    coords: Sequence[Sequence[Number]], offset: Sequence[Number]
) -> List[Tuple[Number]]:
    """Translate each coordinate in a sequence by a given offset.

    Args:
        coords: Sequence of coordinates.
        offset: Offset to translate each element in @coords by. All elements
            in @coords must be of the same length which has to match the length
            of @offset.

    Returns:
        List of tuples containing the coordinates in @coords translated by
        @offset.

    Raises:
        IndexError: If not all elements in @coords have the same length and / or
            if this length does not match the length of @offset.
    """
    if len(coords) == 0:
        return coords
    length = len(coords[0])
    if not all(len(x) == length for x in coords):
        raise IndexError(
            "Not all elements in the sequence of coordinates have the same length")
    if len(offset) != length:
        raise IndexError(
            "The length of the elements in the sequence of coordinates does "
            "not match the length of the given offset")
    return [
        tuple(
            coord[i] + offset[i] for i in range(len(offset))
        ) for coord in coords
    ]


def scale_coords(
    coords: Sequence[Sequence[Number]], scale: Number
) -> List[Tuple[Number]]:
    """Scale each coordinate in a sequence by a given scale.

    Args:
        coords: Sequence of coordinates.
        scale: Factor to scale each element in @coords by.

    Returns:
        List of tuples containing the coordinates in @coords scaled by @scale.
    """
    return [
        tuple(
            dim * scale for dim in coord
        ) for coord in coords
    ]


def transform_coords(
    coords: Sequence[Sequence[Number]], offset: Sequence[Number], scale: Number
) -> List[Tuple[Number]]:
    """Translate and then scale each coordinate in a sequence by a given offset
    and scale.

    Args:
        coords: Sequence of coordinates.
        offset: Offset to translate each element in @coords by. All elements
            in @coords must be of the same length which has to match the length
            of @offset.
        scale: Factor to scale each element in @coords by. The scaling is
            applied AFTER the translation and is applied to all dimensions of
            the coordinates equally.

    Returns:
        List of tuples containing the coordinates in @coords translated by
        @offset and scaled by @scale.

    Raises:
        IndexError: If not all elements in @coords have the same length and / or
            if this length does not match the length of @offset.
    """
    return scale_coords(offset_coords(coords, offset), scale)


def reduce_coords(
    coords: Sequence[Sequence[Number]], n: Optional[int]
) -> List[Tuple[Number]]:
    if n is None or n >= len(coords):
        return coords
    if n <= 0:
        return []
    step = int(len(coords)/n)
    return coords[:n*step:step]


def reproject_coords(
    coords: Sequence[Sequence[Number]],
    transformation: osr.CoordinateTransformation
) -> List[Tuple[Number]]:
    return [transformation.TransformPoint(*coord) for coord in coords]


def get_boundaries(
    element: ET.Element, namespaces
) -> Tuple[np.ndarray, np.ndarray]:
    list_elements = element.findall(".//gml:posList", namespaces)
    single_elements = element.findall(".//gml:pos", namespaces)
    pos_list = ""
    for pos_element in (list_elements + single_elements):
        pos_list += f"{pos_element.text} "
    pos_list = " ".join(pos_list.split()) # remove tabs etc.
    coords = unflatten_sequence(string_to_list(pos_list, type=float))
    return np.amin(coords, axis=0), np.amax(coords, axis=0)


# TODO:
#   Docstrings
#   Triangulation
#   Textures
#   Members could have multiple or no buildings
#   Members could have other sub-elements than buildings
#   Exception handling for ID checks and everything namespace related
#   How to reliably get all surfaces?
#   Surfaces could have multiple polygons / posLists
#   'Inner' polygons cut-outs
#   n_corners <= 2 doesn't make sense / raise exception
class Surface():
    def __init__(
        self, id: Optional[str] = None, verts = List[Tuple[Number]],
    ) -> None:
        self.id = id
        self.verts = verts

    def __str__(self) -> str:
        return f"Surface(id: {self.id}, verts: {self.verts})"

    def __repr__(self) -> str:
        return str(self)

    @staticmethod
    def from_element(
        element: ET.Element, namespaces: Dict[str, str] = {},
        n_corners: int = 0
    ) -> Optional[Surface]:
        face_id = element.attrib.get(substitute_uri("gml:id", namespaces))
        pos_list = element.find(".//gml:posList", namespaces)
        if pos_list is None:
            pos_list = ""
            positions = element.findall(".//gml:pos", namespaces)
            for pos in positions:
                pos_list += f"{pos.text} "
        else:
            pos_list = pos_list.text
        pos_list = " ".join(pos_list.split()) # remove tabs etc.
        coords = unflatten_sequence(string_to_list(pos_list, type=float))
        if coords[0] == coords[-1]:
            coords = coords[:-1]
        if n_corners:
            coords = reduce_coords(coords, n_corners)
        if len(coords) <= 2:
            return
        return Surface(id=face_id, verts=coords)

    def to_mesh(
        self, offset: Sequence[Number] = (0, 0, 0), scale: Number = 1.,
        triangulate: bool = False
    ) -> bmesh.types.BMesh:
        bm = bmesh.new()
        name = "Surface" if self.id is None else self.id
        #coords = reproject_coords(self.verts, trans)
        coords = transform_coords(self.verts, offset=offset, scale=scale)
        for vert in coords:
            bm.verts.new(vert)
        bm.faces.new(bm.verts)
        if triangulate:
            #bm.normal_update()
            bmesh.ops.triangulate(bm, faces=bm.faces)
        mesh = bpy.data.meshes.new(name)
        bm.to_mesh(mesh)
        del bm
        return mesh

    @property
    def lower(self):
        return np.amin(self.verts, axis=0)

    @property
    def upper(self):
        return np.amax(self.verts, axis=0)


class BuildingPart():
    def __init__(
        self, id: Optional[str] = None, faces = List[Surface],
        height: Optional[float] = None
    ) -> None:
        self.id = id
        self.faces = faces
        #if height is None:
        #    lowers = []
        #    uppers = []
        #    for face in faces:
        #        uppers.append(face.lower)
        #        lowers.append(face.upper)
        #    if not self.faces:
        #        height = 0.
        #    else:
        #        height = float(
        #            (np.amax(uppers, axis=0) - np.amin(lowers, axis=0))[-1])
        self.height = height

    def __str__(self) -> str:
        return (
            f"BuildingPart(id: {self.id}, faces: {self.faces}, "
            f"height: {self.height})")

    def __repr__(self) -> str:
        return str(self)

    @staticmethod
    def from_element(
        element: ET.Element, namespaces: Dict[str, str] = {},
        only_footprint: bool = False, n_corners: int = 0
    ) -> Optional[BuildingPart]:
        part_id = element.attrib.get(substitute_uri("gml:id", namespaces))
        faces: List[Surface] = []
        bounds = element.findall('.//bldg:boundedBy', namespaces)
        #lowers = []
        #uppers = []
        for bound in bounds:
            #face = Surface.from_element(bound, namespaces, n_corners)
            #if face is None:
            #    continue
            #lowers.append(face.lower)
            #uppers.append(face.upper)
            #ground = bound.findall(".//bldg:GroundSurface", namespaces)
            #if not ground and only_footprint:
            #    continue
            #faces.append(face)
            ground = bound.findall(".//bldg:GroundSurface", namespaces)
            if not ground and only_footprint:
                continue
            face = Surface.from_element(bound, namespaces, n_corners)
            if face is None:
                continue
        #lower = np.amin(lowers, axis=0)
        #upper = np.amax(uppers, axis=0)
        #height = float((upper - lower)[-1])
        height = None
        if not faces:
            return
        return BuildingPart(id=part_id, faces=faces, height=height)

    def to_mesh(
        self, offset: Sequence[Number] = (0, 0, 0), scale: Number = 1.,
        merge_vertices: bool = False, merge_distance: float = 0.0001,
        extrude_footprint: bool = False, triangulate: bool = False
    ) -> bmesh.types.BMesh:
        bm = bmesh.new()
        name = "BuildingPart" if self.id is None else self.id
        for face in self.faces:
            face_mesh = face.to_mesh(
                offset=offset, scale=scale, triangulate=triangulate)
            bm.from_mesh(face_mesh)
            del face_mesh
        if merge_vertices:
            bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=merge_distance)
        if extrude_footprint:
            # FIXME: Handle extrusion with multiple faces and triangulation
            bm.normal_update()
            ext = bmesh.ops.extrude_face_region(bm, geom=bm.faces)
            ext_verts = [e for e in ext["geom"] if isinstance(e, bmesh.types.BMVert)]
            bmesh.ops.translate(bm, vec=[0, 0, self.height], verts=ext_verts)
        mesh = bpy.data.meshes.new(name)
        bm.to_mesh(mesh)
        del bm
        return mesh


class Building():
    def __init__(
        self, id: Optional[str] = None, parts = List[BuildingPart],
    ) -> None:
        self.id = id
        self.parts = parts

    def __str__(self) -> str:
        return f"Building(id: {self.id}, parts: {self.parts})"

    def __repr__(self) -> str:
        return str(self)

    @staticmethod
    def from_element(
        element: ET.Element, namespaces: Dict[str, str] = {},
        only_footprint: bool = False, n_corners: int = 0
    ) -> Optional[Building]:
        bldg_id = element.attrib.get(substitute_uri("gml:id", namespaces))
        parts: List[BuildingPart] = []
        part_elements = element.findall(".//bldg:BuildingPart", namespaces)
        if not part_elements:
            part_elements = [element]
        for part_element in part_elements:
            part = BuildingPart.from_element(
                part_element, namespaces, only_footprint, n_corners)
            if part is None:
                continue
            parts.append(part)
        if not parts:
            return
        return Building(id=bldg_id, parts=parts)

    def to_mesh(
        self, offset: Sequence[Number] = (0, 0, 0), scale: Number = 1.,
        merge_vertices: bool = False, merge_distance: float = 0.0001,
        extrude_footprints: bool = False, triangulate: bool = False
    ) -> bmesh.types.BMesh:
        bm = bmesh.new()
        name = "Building" if self.id is None else self.id
        for part in self.parts:
            part_mesh = part.to_mesh(
                offset=offset, scale=scale,
                merge_vertices=merge_vertices,
                merge_distance=merge_distance,
                extrude_footprint=extrude_footprints,
                triangulate=triangulate)
            bm.from_mesh(part_mesh)
            del part_mesh
        mesh = bpy.data.meshes.new(name)
        bm.to_mesh(mesh)
        del bm
        return mesh


def main(
    filename: str, origin: Sequence[Number] = (0, 0, 0), scale: Number = 1,
    only_footprints: bool = False, n_corners: int = 0,
    extrude_footprints: bool = False, triangulate_faces: bool = False,
    merge_vertices: bool = False, merge_distance: float = 0.001,
    recalculate_view: bool = False
) -> None:
    namespaces = get_namespaces(filename)
    offset = tuple([-1 * dim for dim in origin])

    #gdal.UseExceptions()
    #crs_src = osr.SpatialReference()
    #crs_src.ImportFromEPSG(7415)
    #crs_tar = osr.SpatialReference()
    #crs_tar.ImportFromEPSG(25832)
    #trans = osr.CoordinateTransformation(crs_src, crs_tar)

    try:
        member_tag = substitute_uri("core:cityObjectMember", namespaces)
    except KeyError:
        member_tag = substitute_uri("cityObjectMember", namespaces)
    #member_tag = substitute_uri("ogr:featureMember", namespaces)

    buildings: List[Building] = []
    for _, element in ET.iterparse(filename, events=("end",)):
        if element.tag != member_tag:
            continue
        building_elements = element.findall(".//bldg:Building", namespaces)
        #building_elements = element.findall(".//ogr:entities", namespaces)
        for building_element in building_elements:
            building = Building.from_element(
                building_element, namespaces, only_footprints, n_corners)
            if building is not None:
                buildings.append(building)
        element.clear()

    if not buildings:
        raise ValueError("No valid geometry found in GML file.")

    #bldg_ids = []
    #poly_ids = []
    bm = bmesh.new()
    for building in buildings:
        building_mesh = building.to_mesh(
            offset=offset, scale=scale,
            merge_vertices=merge_vertices,
            merge_distance=merge_distance,
            extrude_footprints=extrude_footprints,
            triangulate=triangulate_faces)
        bm.from_mesh(building_mesh)
        del building_mesh

#    for _, element in ET.iterparse(filename, events=("end",)):
#        if element.tag != member_tag:
#            continue
#
#        #print(f"element: {element}")
#        bldgs = element.findall(".//bldg:Building", namespaces)
#        #bldgs = element.findall(".//ogr:entities", namespaces)
#        for bldg in bldgs:
#            #print(f"\tbuilding: {bldg}")
#            ##bldg_id = bldg.attrib.get(substitute_uri("gml:id", namespaces))
#            ##if bldg_id is not None and bldg_id in bldg_ids:
#            ##    continue
#            ##bldg_ids.append(bldg_id)
#            #print(f"\tid: {bldg_id}")
#            bldg_verts = []
#            bldg_faces = []
#            #parts = bldg.findall(".//bldg:BuildingPart", namespaces)
#            parts = []
#            if len(parts) == 0:
#                parts = [bldg]
#            for part in parts:
#                polygons = part.findall('.//gml:Polygon', namespaces)
#                patches = part.findall('.//gml:PolygonPatch', namespaces)
#                triangles = part.findall('.//gml:Triangle', namespaces)
#                faces = polygons + triangles + patches
#                #print(f"\t{len(faces)}")
#                for face in faces:
#                    ##poly_id = face.attrib.get(substitute_uri("gml:id", namespaces))
#                    ##if poly_id is not None and poly_id in poly_ids:
#                    ##    continue
#                    ##poly_ids.append(poly_id)
#                    face_verts = []
#                    pos_list = face.find(".//gml:posList", namespaces)
#                    if pos_list is None:
#                        pos_list = ""
#                        positions = face.findall(".//gml:pos", namespaces)
#                        for pos in positions:
#                            pos_list += f"{pos.text} "
#                    else:
#                        pos_list = pos_list.text
#                    pos_list = " ".join(pos_list.split()) # remove tabs etc.
#                    coords = unflatten_sequence(string_to_list(pos_list, type=float))
#                    if coords[0] == coords[-1]:
#                        coords = coords[:-1]
#                    #coords = reproject_coords(coords, trans)
#                    coords = transform_coords(coords, offset=offset, scale=scale)
#                    for vert in coords:
#                        face_verts.append(bm.verts.new(vert))
#                    if len(face_verts) < 3:
#                        continue
#                    bldg_faces.append(bm.faces.new(face_verts))
#                    bldg_verts.extend(face_verts)
#            if triangulate_faces:
#                bm.normal_update()
#                bmesh.ops.triangulate(bm, faces=bldg_faces)
#            if merge_vertices:
#                bmesh.ops.remove_doubles(
#                    bm, verts=bldg_verts, dist=merge_distance)
#        element.clear()

    obj_name = os.path.basename(filename)
    mesh = bpy.data.meshes.new(obj_name)
    bm.to_mesh(mesh)
    obj = bpy.data.objects.new(obj_name, mesh)
    bpy.context.collection.objects.link(obj)
    del bm
    del mesh
    del buildings

    if recalculate_view:
        viewport = None
        screen = bpy.context.screen
        for area in screen.areas:
            if area.type == "VIEW_3D":
                viewport = area
        if viewport is not None:
            pow10 = 1 + math.floor(math.log10(obj.dimensions.length))
            viewport.spaces[0].clip_end = max(10**(pow10 + 3), 100)


class CityGMLSelector(bpy.types.Operator, ImportHelper):
    @classmethod
    def description(cls, context, properties):
        return "Import geometry from CityGML file(s)"

    bl_idname = "wm.citygml_selector"
    bl_label = "Pick GML file(s)"

    filename_ext = ".gml"
    use_filter_folder = True

    files: CollectionProperty(type=bpy.types.PropertyGroup)
    import_scale: FloatProperty(
        name="Import Scale",
        description="Scale imported geometry by this factor.",
        soft_min=0.0, soft_max=1.0,
        precision=3,
        default=1.0,
    )
    origin_x: FloatProperty(
        name="Origin Point X",
        description="X value from imported file to be 0",
        precision=1,
        default=0.0,
    )
    origin_y: FloatProperty(
        name="Origin Point Y",
        description="Y value from imported file to be 0",
        precision=1,
        default=0.0,
    )
    origin_z: FloatProperty(
        name="Origin Point Z",
        description="Z value from imported file to be 0",
        precision=1,
        default=0.0,
    )
    recalculate_view: BoolProperty(
        name="Recalculate View Clip End",
        description=(
            "Recalculate view clip end in current viewport "
            "based on size of imported geometry"),
        default=False,
    )
    only_footprints: BoolProperty(
        name="Import only footprints",
        description="Import only 'GroundSurface' faces",
        default=False,
    )
    reduce_footprints: IntProperty(
        name="Reduce footprints",
        description="Reduce footprint polygons",
        default=0,
        min=0,
        subtype="UNSIGNED",
    )
    extrude_footprints: BoolProperty(
        name="Extrude footprints",
        description="Extrude footprints to height of buildings",
        default=False,
    )
    triangulate_faces: BoolProperty(
        name="Triangulate Faces",
        description="Triangulate faces (slow)",
        default=False,
    )
    merge_vertices: BoolProperty(
        name="Merge Vertices",
        description="Merge duplicate vertices of buildings (slow)",
        default=False,
    )
    merge_distance: FloatProperty(
        name="Merge Distance",
        description="Vertex merge distance",
        precision=4,
        default=0.001,
    )

    def draw(self, context):
        layout = self.layout

        box = layout.box()
        row = box.row(align=True)
        row.label(text="Origin Point:")
        row = box.row(align=True)
        row.prop(self, "origin_x", text="X:")
        row = box.row(align=True)
        row.prop(self, "origin_y", text="Y:")
        row = box.row(align=True)
        row.prop(self, "origin_z", text="Z:")

        box = layout.box()
        row = box.row(align=True)
        row.label(text="Import Scale:")
        row = box.row(align=True)
        row.prop(self, "import_scale", text="")

        box = layout.box()
        row = box.row(align=True)
        row.label(text="Footprint Options:")
        row = box.row(align=True)
        row.prop(self, "only_footprints", text="Import Footprints Only")
        row = box.row(align=True)
        row.prop(self, "reduce_footprints", text="Reduce to N-Gons (0 = off)")
        row = box.row(align=True)
        row.prop(self, "extrude_footprints", text="Extrude Footprints")

        box = layout.box()
        row = box.row(align=True)
        row.label(text="Mesh Operations:")
        row = box.row(align=True)
        row.prop(self, "triangulate_faces", text="Triangulate Faces (slow)")
        row = box.row(align=True)
        row.prop(self, "merge_vertices", text="Merge Vertices (slow)")
        row = box.row(align=True)
        row.prop(self, "merge_distance", text="Merge Distance")

        box = layout.box()
        row = box.row(align=True)
        row.prop(self, "recalculate_view", text="Set View Clip End")

    def execute(self, context):
        folder = (os.path.dirname(self.filepath))
        n = len(self.files)
        pow10 = 1 + math.floor(math.log10(n))
        for i, file in enumerate(self.files):
            path_to_file = os.path.join(folder, file.name)
            start_time = time.monotonic()
            print(f"[{i+1:0{pow10}}/{n}] Processing '{file.name}'")
            try:
                main(
                    filename=path_to_file,
                    origin=(
                        self.origin_x,
                        self.origin_y,
                        self.origin_z,
                    ),
                    scale=self.import_scale,
                    only_footprints=self.only_footprints,
                    n_corners=self.reduce_footprints,
                    extrude_footprints=self.extrude_footprints,
                    triangulate_faces=self.triangulate_faces,
                    merge_vertices=self.merge_vertices,
                    merge_distance=self.merge_distance,
                    recalculate_view=self.recalculate_view,
                )
                duration = time.monotonic() - start_time
                print(f"\tSuccesfully imported in {duration:.2f} seconds.")
            except Exception as exc:
                print(f"\tError while importing: [{type(exc).__name__}] {exc}")
        return{"FINISHED"}


def menu_import(self, context):
    self.layout.operator(CityGMLSelector.bl_idname, text="CityGML (.gml)")


def register():
    bpy.utils.register_class(CityGMLSelector)
    bpy.types.TOPBAR_MT_file_import.append(menu_import)


def unregister():
    bpy.types.TOPBAR_MT_file_import.remove(menu_import)
    bpy.utils.unregister_class(CityGMLSelector)
