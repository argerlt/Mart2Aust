# **Mart2Aust**
---
### ***WARNING: this is still a work in progress. Use with caution, as things might break and are often inadequately documented**
---

Mart2Aus reconstructs Parent grain structures from Electron BackScatter Diffraction (EBSD) scans of child structures. The code is provided in Matlab, specifically to take advantage of the crystallographic functions available through MTEX.

Users should note, while the crystallographic calculations and distribution functions are performed using MTEX, this algorithm is fundamentally differnent than the default parent-child reconstruction tool provided in MTEX. Whereas that tool is based on a clustering algorithm, this code uses a variation on image graph cut theory, which was originally developed by (cite). A Mincut-Maxflow algorithm is used to seperate out pixels suspected of originating from the same parent grain. 

While this tool has been set up with Martensitic steel in mind, generalizing it to any other system (such as $\alpha/\beta$ titanium should be semi-trivial

## Installation

Mart2Aust Requires Matlab 2018a or higher, and MTEX 5.7.0 (NOTE: MTEX is not always backwards compatable. Users hoping to use newer features in the future may need to downgrade)

To install, copy this repo and unzip the mtex.zip file into the "Functions" folder. 
To run a reconstruction, open the "Examples" folder and chose any of the files to run. Most focus on a small subset of scans taken from AF96, details of which can be found here :(insert link) along with several other EBSD scans of the same material.


## Background information

Users seeking to understand exactly **WHAT** this code is doing are encouraged to explor the Literature folder, which contains the publications relevant to this work. However broadly speaking, the most important concepts are:

- **Texture and orientation analysis**
- **Graph cut Theory**
- **Parent/Child Oreintation relationships**
  
TODO: add actual links for places to learn about these things

## Example Reconstructions

**Initial:**
---
![Martensitic AF96 scan](/Resources/EBSD/images/AF96_Large.png)


**Reconstructed:**
---
![Reconstructed Austenite](/Resources/EBSD/images/AF96_Large_out.jpeg)



## Credits

**Implementation of the Graph Cut Algorithm**
Alexander Brust
Austin Gerlt
Steve Niezgoda

**Orientation relationship determination**
Eric Payton
Victoria Yardley

**Collection of the example EBSD data**
Vikas Sinha
Jared Shank

**Code cleanup, modernization, bugfixing, etc**
Austin Gerlt


## TODO
- Cleanup this Readme
- Add proper references to graph cut, MTEX, MATLAB, AF96 data, etc.
- fix example files
- generalize to allow for Titanium and/or any Parent/child system
- document AusRecon precursor
- link publications
- add variant plotting functions
- add ASTM grain size