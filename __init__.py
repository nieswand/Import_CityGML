bl_info = {
    "name": "CityGML Import Basic",
    "author": "",
    "version": (0, 4),
    "blender": (2, 80, 0),
    "category": "Import-Export"
}

import os
from xml.etree import ElementTree as et

import bpy
from bpy_extras.io_utils import ImportHelper

import re

def main(filename):

    def unflatten(coords, s=0.001):
        return [(coords[i] * s, coords[i + 1] * s, coords[i + 2] * s) for i in range(0, len(coords), 3)]

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
    bl_label = "pick an xml file"

    filename_ext = ".xml"
    use_filter_folder = True

    def execute(self, context):
        fdir = self.properties.filepath
        main(fdir)
        return{'FINISHED'}


def menu_import(self, context):
    self.layout.operator(CityGMLDirectorySelector.bl_idname, text="cityGML (.xml)")


def register():
    bpy.utils.register_class(CityGMLDirectorySelector)
    bpy.types.TOPBAR_MT_file_import.append(menu_import)


def unregister():
    bpy.types.TOPBAR_MT_file_import.remove(menu_import)
    bpy.utils.unregister_class(CityGMLDirectorySelector)
