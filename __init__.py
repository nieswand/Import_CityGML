import bmesh
import bpy
from bpy_extras.io_utils import ImportHelper
from bpy.props import BoolProperty, FloatProperty, CollectionProperty
import math
from numbers import Number
import os
from xml.etree import ElementTree as ET
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


def main(
    filename: str, origin: Sequence[Number] = (0, 0, 0), scale: Number = 1,
    merge_vertices: bool = False, merge_distance: float = 0.001,
    recalculate_view: bool = False
) -> None:
    namespaces = get_namespaces(filename)
    offset = tuple([-1 * dim for dim in origin])

    try:
        member_tag = substitute_uri("core:cityObjectMember", namespaces)
    except KeyError:
        member_tag = substitute_uri("cityObjectMember", namespaces)

    bm = bmesh.new()
    for _, element in ET.iterparse(filename, events=("end",)):
        if element.tag != member_tag:
            continue

        bldgs = element.findall(".//bldg:Building", namespaces)
        for bldg in bldgs:
            bldg_verts = []
            parts = bldg.findall(".//bldg:BuildingPart", namespaces)
            if len(parts) == 0:
                parts = [bldg]
            for part in parts:
                polygons = part.findall('.//gml:Polygon', namespaces)
                triangles = part.findall('.//gml:Triangle', namespaces)
                faces = polygons + triangles
                for face in faces:
                    face_verts = []
                    pos_list = face.find(".//gml:posList", namespaces)
                    if pos_list is None:
                        pos_list = ""
                        positions = face.findall(".//gml:pos", namespaces)
                        for pos in positions:
                            pos_list += f"{pos.text} "
                    else:
                        pos_list = pos_list.text
                    pos_list = " ".join(pos_list.split()) # remove tabs etc.
                    coords = unflatten_sequence(string_to_list(pos_list, type=float))
                    if coords[0] == coords[-1]:
                        coords = coords[:-1]
                    coords = transform_coords(coords, offset=offset, scale=scale)
                    for vert in coords:
                        face_verts.append(bm.verts.new(vert))
                    bm.faces.new(face_verts)
                    bldg_verts.extend(face_verts)
            if merge_vertices:
                bmesh.ops.remove_doubles(
                    bm, verts=bldg_verts, dist=merge_distance)
        element.clear()

    if len(bm.verts) == 0:
        del bm
        raise ValueError("No geometry found in GML file.")

    obj_name = os.path.basename(filename)
    mesh = bpy.data.meshes.new(obj_name)
    bm.to_mesh(mesh)
    obj = bpy.data.objects.new(obj_name, mesh)
    bpy.context.collection.objects.link(obj)

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
        row.prop(self, "merge_vertices", text="Merge Vertices (slow)")
        row = box.row(align=True)
        row.prop(self, "merge_distance", text="Merge Distance")

        box = layout.box()
        row = box.row(align=True)
        row.prop(self, "recalculate_view", text="Set View Clip End")

    def execute(self, context):
        folder = (os.path.dirname(self.filepath))
        for file in self.files:
            path_to_file = os.path.join(folder, file.name)
            try:
                main(
                    filename=path_to_file,
                    origin=(
                        self.origin_x,
                        self.origin_y,
                        self.origin_z,
                    ),
                    scale=self.import_scale,
                    merge_vertices=self.merge_vertices,
                    merge_distance=self.merge_distance,
                    recalculate_view=self.recalculate_view,
                )
                print(f"'{file.name}' imported")
            except Exception as exc:
                print(
                    f"Error while importing '{file.name}': "
                    f"[{type(exc).__name__}] {exc}")
        return{"FINISHED"}


def menu_import(self, context):
    self.layout.operator(CityGMLSelector.bl_idname, text="CityGML (.gml)")


def register():
    bpy.utils.register_class(CityGMLSelector)
    bpy.types.TOPBAR_MT_file_import.append(menu_import)


def unregister():
    bpy.types.TOPBAR_MT_file_import.remove(menu_import)
    bpy.utils.unregister_class(CityGMLSelector)
