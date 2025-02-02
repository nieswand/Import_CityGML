# Import CityGML

This is a [Blender](https://www.blender.org/) add-on to import geometry (only geometry) from 3D city models encoded in CityGML files.
It was started by Dealga McArdle ([zeffii](https://github.com/zeffii/)) in 2015 as [io_cityGML_basic](https://github.com/zeffii/io_cityGML_basic) add-on for Blender 2.7. In 2020, the add-on was updated for Blender 2.8 and 2.9 by [ppaawweeuu](https://github.com/ppaawweeuu).
In 2023, the add-on was updated again by Simon Nieswand ([nieswand](https://github.com/nieswand)).

## Basic Functionality
This add-on imports geometry from CityGML file (or files) as one mesh/object (per file).
> basic means basic :)

There can be much non-geometrical information in CityGML files, for example materials, textures, building part classifications (wall, roof, floor, furniture), number of floors, roof type, adresses, ocupancy and much more. Properties can be assigned to objects (building, road, tunnel) or parts (wall, roof), or geometry (solids, polygons, triangles). There is no limitation in amount of information stored in CityGML file and there was no intention to import more than geometry.

## Installation

1. Download repository as zip (through GitHub this can be done through the `Clone or download` > `Download ZIP` green button).

2. Run `Blender` and go to `Edit > Preferences > Add-ons` and press the `Install...` button.

3. Select the downloaded ZIP and press `Install Add-on from File...`.

4. Enable the add-on from the list, by ticking the empty box next to the add-on's name. (Optional: If this is an update, un-tick and tick the add-on to reload it!)

## Features:

### Menu

Import menu item (import CityGML .gml)

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/01_menu.png)


### As One

Imports all geometry as one object

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/04_single.png)


### Multi-file

Single or multiply file import (multiple files = multiple objecs)

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/03_multiply.png)


### Names

Imported file name = object & mesh name

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/02_obj_mesh_names.png)


### Options

Imported object are located in (0, 0, 0) (their origin) but mesh points can be for example x=12342222, y=23343322, z=23. City models can be huge. That is the reason why there are some problems in Blender with displaying geometry. It means problems with zooming, rotating, shadows, view clipping etc.

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/06_problems.gif)

There are 3 workarounds:

1. Origin Point

We can overwrite imported geometry origin point to be closer to Blender world origin (0, 0, 0).

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/06_options_origin.png)

CityGML files are readable in Notepad ([Notepad++](https://notepad-plus-plus.org/) is much better, especially for very big text files). We can find series of coordinates between `<gml:posList>` or `<gml:pos>` tags.

Below example when we could set import Origin Point x=750000, y=380000, z=0

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/07_size.png)



2. Import Scale

We can scale all imported geometry to be closer to World Origin.

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/06_options_scale.png)



3. Recalculate View Clip End

When we work on big objects we need to change the view clipping range. This option can recalulate the view clip end based on the size of the imported geometry.

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/06_options_clips.png)



### Other

Other features:
* Imports into current collection
* Imports `<gml:posList>` and lists of `<gml:pos>` tags
* Imports ground triangles (if any)
* Exception handling if GML file has no geometry & console progress info

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/05_console.png)



## Testing Datasets

Add-on was tested on files from polish [geoportal.gov.pl](https://mapy.geoportal.gov.pl/)

Here are some direct links to CityGML files for 3 cities:

[Warszawa - map link](https://mapy.geoportal.gov.pl/imap/Imgp_2.html?locale=pl&gui=new&sessionID=5214135) / [Warszawa CityGML LoD1](https://integracja.gugik.gov.pl/Budynki3D/pobierz.php?d=2&plik=powiaty/lod1/1465_gml.zip) / [Warszawa CityGML LoD2](https://integracja.gugik.gov.pl/Budynki3D/pobierz.php?d=2&plik=powiaty/1465_gml.zip)

[Lublin - map link](https://mapy.geoportal.gov.pl/imap/Imgp_2.html?locale=pl&gui=new&sessionID=5214138) / [Lublin CityGML LoD1](https://integracja.gugik.gov.pl/Budynki3D/pobierz.php?d=2&plik=powiaty/lod1/0663_gml.zip) / [Lublin CityGML LoD2](https://integracja.gugik.gov.pl/Budynki3D/pobierz.php?d=2&plik=powiaty/0663_gml.zip)

[Krakrów - map link](https://mapy.geoportal.gov.pl/imap/Imgp_2.html?locale=pl&gui=new&sessionID=5214142) /
[Kraków CityGML LoD1](https://integracja.gugik.gov.pl/Budynki3D/pobierz.php?d=2&plik=powiaty/lod1/1261_gml.zip) /
[Kraków CityGML LoD2](https://integracja.gugik.gov.pl/Budynki3D/pobierz.php?d=2&plik=powiaty/1261_gml.zip)

Add-on was also tested on other [CityGML datasets](https://www.cityjson.org/datasets/#datasets-converted-from-citygml) from [CityJSON.org](https://www.cityjson.org/) site (for example [New York](https://www1.nyc.gov/site/doitt/initiatives/3d-building.page)).



## Limitations
There is general issue with polygons with holes. Polygons in CityGML can be created by two (or more) rings. In that case inner ring is interpreted as void and is substracted from outer ring. [Rules_for_Validating_GML_Geometries_in_CityGML](http://en.wiki.quality.sig3d.org/index.php/Modeling_Guide_for_3D_Objects_-_Part_1:_Basics_(Rules_for_Validating_GML_Geometries_in_CityGML)#gml:Polygon).

Blender sees two coplanar faces with connected vericles.

![](https://github.com/ppaawweeuu/Import_CityGML/blob/master/img/11_problems_holes.png)
