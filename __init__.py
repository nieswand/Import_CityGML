bl_info = {
    "name": "Import CityGML",
    "author": "Dealga McArdle, ppaawweeuu, Simon Nieswand",
    "version": (0, 5, 0),
    "blender": (2, 80, 0),
    "category": "Import-Export",
    "description": "Import geometry from CityGML file(s)",
    "wiki_url": "https://github.com/ppaawweeuu/Import_CityGML",
}


import bmesh
import bpy
from bpy_extras.io_utils import ImportHelper
from bpy.props import BoolProperty, FloatProperty, CollectionProperty
import math
from numbers import Number
import os
from xml.etree import ElementTree as ET
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple


def get_prefix_map(
    filename: str, allow_duplicates: bool = False, exit_early: bool = True
) -> Dict[str, str]:
    prefix_map = {}
    for event, element in ET.iterparse(filename, events=("start-ns", "start")):
        if event == "start-ns":
            prefix, uri = element
            saved_uri = prefix_map.get(prefix)
            if saved_uri is not None and saved_uri != uri and not allow_duplicates:
                raise KeyError(
                    f"Duplicate prefix '{prefix}' with different URIs found. "
                    f"Saved URI: {saved_uri}, new URI: {uri}.")
            prefix_map[prefix] = uri
        elif event == "start" and exit_early:
            break
    return prefix_map


def substitute_uri(
    tag: str, prefix_map: Optional[Dict[str, str]] = None
) -> str:
    if prefix_map is None or len(prefix_map) == 0:
        return tag
    splits = tag.split(":")
    if len(splits) == 1:
        return tag
    elif len(splits) == 2:
        prefix, tag = splits
        uri = prefix_map.get(prefix)
        if uri is None:
            raise KeyError(
                f"Prefix '{prefix}' not found in provided prefix map.")
        return f"{{{uri}}}{tag}"
    else:
        raise AttributeError(f"Multiple colons in tag '{tag}' are invalid.")


def string_to_list(
    string: str, typ: Optional[Callable] = None, delimiter: str = " "
) -> List[Any]:
    splits = string.split(delimiter)
    splits = filter(None, splits)
    if typ is None:
        return list(splits)
    return list(map(typ, filter(None, splits)))


def unflatten_poslist(
    poslist: Sequence[Number], n_dims: int = 3
) -> List[Tuple[Number]]:
    return [
        tuple(
            poslist[i + k] for k in range(n_dims)
        ) for i in range(0, len(poslist), n_dims)
    ]


def offset_coords(
    coords: Sequence[Sequence[Number]], offset: Sequence[Number]
) -> List[Tuple[Number]]:
    return [
        tuple(
            coord[i] + offset[i] for i in range(len(offset))
        ) for coord in coords
    ]


def scale_coords(
    coords: Sequence[Sequence[Number]], scale: Number
) -> List[Tuple[Number]]:
    return [
        tuple(
            dim * scale for dim in coord
        ) for coord in coords
    ]


def transform_coords(
    coords: Sequence[Sequence[Number]], offset: Sequence[Number], scale: Number
) -> List[Tuple[Number]]:
    return scale_coords(offset_coords(coords, offset), scale)


def main(
    filename: str, origin: Sequence[Number] = (0, 0, 0), scale: Number = 1,
    merge_vertices: bool = False, merge_distance: float = 0.001,
    recalculate_view: bool = False
) -> None:
    prefix_map = get_prefix_map(filename)
    offset = tuple([-1 * dim for dim in origin])

    bm = bmesh.new()
    for _, element in ET.iterparse(filename, events=("end",)):
        if element.tag != substitute_uri("core:cityObjectMember", prefix_map):
            continue

        bldgs = element.findall(".//bldg:Building", prefix_map)
        for bldg in bldgs:
            bldg_verts = []
            parts = bldg.findall(".//bldg:BuildingPart", prefix_map)
            if len(parts) == 0:
                parts = [bldg]
            for part in parts:
                polygons = part.findall('.//gml:Polygon', prefix_map)
                triangles = part.findall('.//gml:Triangle', prefix_map)
                faces = polygons + triangles
                for face in faces:
                    face_verts = []
                    pos_list = face.find(".//gml:posList", prefix_map)
                    if pos_list is None:
                        pos_list = ""
                        positions = face.findall(".//gml:pos", prefix_map)
                        for pos in positions:
                            pos_list += f"{pos.text} "
                    else:
                        pos_list = pos_list.text
                    pos_list = " ".join(pos_list.split()) # remove tabs etc.
                    coords = unflatten_poslist(string_to_list(pos_list, typ=float))
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
