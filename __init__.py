bl_info = {
    "name": "Import CityGML",
    "author": "",
    "version": (0, 4, 2),
    "blender": (2, 90, 0),
    "category": "Import-Export",
    "description": "Import geometry from CityGML file(s)",
    "wiki_url": "https://github.com/ppaawweeuu/Import_CityGML",    
#    "tracker_url": "http://  "
}

import os
from xml.etree import ElementTree as et

import bpy
from bpy_extras.io_utils import ImportHelper

import re

from bpy.props import (BoolProperty,
                       FloatProperty,
                       StringProperty,
                       EnumProperty,
                       CollectionProperty,
                       FloatVectorProperty)

def main(filename, scale, origin):

    def unflatten(coords, s=scale):
        return [((coords[i]-origin[0]) * s, (coords[i + 1]-origin[1]) * s, (coords[i + 2]-origin[2]) * s) for i in range(0, len(coords), 3)]

    tree = et.parse(filename)
    polygons = []
    polygons = tree.findall('.//{http://www.opengis.net/gml}Polygon')
    triangles = []
    triangles = tree.findall('.//{http://www.opengis.net/gml}Triangle')
    faces = []
    faces = polygons + triangles
        
    master_verts = []
    master_faces = []
    extend_verts = master_verts.extend
    append_faces = master_faces.append
    
    test = len(faces[0].findall('.//{http://www.opengis.net/gml}posList'))

    if test == 0:
        for i, p in enumerate(faces):
            texts = []
            for poslist in p.findall('.//{http://www.opengis.net/gml}pos'):
                positionfull = poslist.text
                position = positionfull.strip() #remove first and last whitespace
                texts.append(position)
            text = ' '.join(texts)

            coords = [float(i) for i in text.split(' ')]
            verts = unflatten(coords)
            
            start_idx = len(master_verts)
            end_idx = start_idx + len(verts)
            extend_verts(verts)
            append_faces([i for i in range(start_idx, end_idx)])
            
    else:
        for i, p in enumerate(faces):
            poslist = p.find('.//{http://www.opengis.net/gml}posList')
            textfull = poslist.text
            textreduce1 = re.sub('\n', ' ', textfull) #remove line breaks
            textreduce2 = re.sub('\t', ' ', textreduce1) #remove tabs
            textreduce3 = re.sub(' +', ' ', textreduce2) # remove multiple whitespaces
            text = textreduce3.strip() #remove first and last whitespace

            coords = [float(i) for i in text.split(' ')]
            verts = unflatten(coords)

            start_idx = len(master_verts)
            end_idx = start_idx + len(verts)
            extend_verts(verts)
            append_faces([i for i in range(start_idx, end_idx)])

    ob_name = os.path.basename(filename)
    
    mesh = bpy.data.meshes.new(ob_name)
    mesh.from_pydata(master_verts, [], master_faces)
    mesh.update()

    obj = bpy.data.objects.new(ob_name, mesh)

    scene = bpy.context.scene
    current_collection = bpy.context.collection.name
    scene.collection.children[current_collection].objects.link(obj)


# pick folder and import from... maybe this is a bad name.
class CityGMLDirectorySelector(bpy.types.Operator, ImportHelper):
    bl_idname = "wm.citygml_folder_selector"
    bl_label = "pick an xml file(s)"

    filename_ext = ".xml"
    use_filter_folder = True
    
    files: CollectionProperty(type=bpy.types.PropertyGroup)
    scale_setting: FloatProperty(
        name="Import Scale",
        description="1 for meters, 0.001 for kilometers",
        min=0.0, max=1.0,
        soft_min=0.0, soft_max=1.0,
        precision = 3,
        default=1.0
        )
    origin_setting: FloatVectorProperty(
        name="Origin Point",
        description="can be helpful for reduce large coordinates values",
        precision = 1,
        size = 3,
        default=(0.0, 0.0, 0.0)
        )    
    
    def draw(self, context):
        layout = self.layout
        
        box = layout.box()
        row = box.row(align=True)
        row.label(text='Import Scale:')
        row = box.row(align=True)
        row.prop(self, "scale_setting", text='')
        
        box = layout.box()
        row = box.row(align=True)
        row.label(text='Origin Offset (X,Y,Z):')
        row = box.row(align=True)
        row.prop(self, "origin_setting", text='')

    def execute(self, context):
        folder = (os.path.dirname(self.filepath))
        for i in self.files:
            path_to_file = (os.path.join(folder, i.name))
            try:
                main(
                    filename = path_to_file,
                    scale = self.scale_setting,
                    origin = self.origin_setting)
                print(str(i.name) + " imported")
            except:
                print(str(i.name) + " error, no valid geometry")
        return{'FINISHED'}

def menu_import(self, context):
    self.layout.operator(CityGMLDirectorySelector.bl_idname, text="cityGML (.xml)")


def register():
    bpy.utils.register_class(CityGMLDirectorySelector)
    bpy.types.TOPBAR_MT_file_import.append(menu_import)


def unregister():
    bpy.types.TOPBAR_MT_file_import.remove(menu_import)
    bpy.utils.unregister_class(CityGMLDirectorySelector)
