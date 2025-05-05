# A Surface-based Appearance Model for Pennaceous Feathers

## Authors

[Juan Raul Padron Griffe](https://juanraul8.github.io/),
[Dario Lanza](https://dariolanza95.github.io/),
[Adrian Jarabo](http://giga.cps.unizar.es/~ajarabo/),
[Adolfo Muñoz](http://webdiis.unizar.es/~amunoz/es/)<br>

## Introduction

The appearance of a real-world feather is the result of a complex light interaction with its multi-scale biological structure including the central shaft, branching barbs and interlocking barbules on those barbs. In this work, we propose a practical surface-based appearance model for feathers. We represent the far-field appearance of feathers using a BSDF that implicitly represents the light scattering from the main biological structures of a feather such as the shaft, barb and barbules. Our model accounts for the particular characteristics of feather barbs such as the non-cylindrical cross-sections and the scattering media via a numerically-based BCSDF. To model the relative visibility between barbs and barbules, we derive a masking term for the differential projected areas of the different components of the feather’s microgeometry, which allows to analytically compute the masking between barbs and barbules. As opposed to previous works, our model uses a lightweight representation of the geometry based on a 2D texture, and does not require to explicitly represent the barbs as curves. We show the flexibility and potential of our appearance model approach to represent the most important visual features of several pennaceous feathers.\\ 

Implementation of [A Surface-based Appearance Model for Pennaceous Feathers](https://graphics.unizar.es/projects/FeathersAppearance_2024/) in Mitsuba 0.6.

![Teaser render](https://github.com/juanraul8/PennaceousFeathersRendering/blob/main/resources/feathers_appearance_teaser.png)

## Mitsuba Installation

First clone the [Mitsuba](http://mitsuba-renderer.org/download.html) repository and compile it following its instructions. 

Copy the configuration file to build:

``` cp installation/config.py ./mitsuba/ ```

You can set up a conda environment with all dependencies using the following command:

``` conda env create -f mitsuba_environment.yml ```

Activate anaconda environment for scons compilation:

``` conda activate mitsuba ```  

Compile mitsuba by running the following comand inside the mitsuba folder: 

``` scons -j20 ```

-j indicate the number of threads. The compilation usually takes around 2 minutes. For more details, you can check the -mitsuba_utils.sh- script or the Mitsuba documentation.

Once mitsuba is compiled, copy the content inside feather_plugin inside the mitsuba folder and compile again.

## Mitsuba Rendering

Once the compilation is ready, update the paths:

``` source setpath.sh ```

Run the following command to render the teaser scene with our complete scattering model:

``` mitsuba teaser.xml ./results/teaser_ours.exr -p 20 -Dspp=1024 -Dmax_depth=65 -Dsensor_scene=./sensors/high_res_paper_sensor.xml -o  -Duse_default_sampling_ours=$default_sampling -Duse_default_sampling_ours=true ```

This experiment should take around 25 minutes for a 1280x720 image using 20 threads with 1024 spp. 

Use the following command to transform the render from .exr to tonemapped .png:

``` mtsutil tonemap scenes/teaser/teaser.exr ``` 

## Figures Reproduction

The scenes can be download from the project website including the scene configuration file, the 3D assets and spectral data.  

### Figure 01 (Teaser)

Go into the experiment folder:
``` cd scenes/teaser/ ```

Execute the following script to create the renders of the teaser figure:

``` time bash render_teaser.sh ``` 

This experiment takes around 37 minutes for 1280x720 images and 1024 spp.

### Figure 08 (Ablation Studies, Fiber BCSDF)



### Figure 08 (Ablation Studies, Feather BSDF)

Run the following script inside the folder scenes/ablation_BSDF/ to reproduce the renders show in Figure 08 with feather pelts for red, blue, green and black feathers:

``` time bash render_ablation_studies.sh ``` 

Inside the script you can set up the spp. This experiment takes around 1h39m for 1280x70 images and 1024 spp. Using the following command we create renders using the bird and the wing of feathers:

``` time bash render_ablation_studies_wing.sh ``` 

You can also find these renders in the project website.  

### Figure 10 (Parameter Exploration, Supplemental Figures)


Perform the appearance exploration for two parameters using the following script to replicate the Figure 10:

``` time python render_appearance_exploration.py --conf_file=figure12.json ```

This experiment takes around 13 minutes for 512x512 images and 256 spp.

Another option is setting the mitsuba_scene to single_feather_appearance_exploration.xml to perform the parameter exploration in single feathers. All the parameter exploration can be executed with the corresponding configuration file, like the supp_figure9 configuration file for barb eumelanin and barb pheomelanin variation:


``` time python render_appearance_exploration.py --conf_file=supp_figure9.json ```

This experiment takes around 28 minutes for 512x512 images and 256 spp.


### Figure 11 (Appearance Matching: Amazon Parrot, Goniochromatic Effect)

Execute the following script inside the folder scenes/appearance_matching_single_feather/amazon_parrot_matching to create the renders for the appearance matching of the Amazon parrot feather:

``` time bash goniochromatic_appearance.sh ``` 

This experiment takes around 7 minutes for 256x492 images and 256 spp.

A video of the goniochromatic effect can be created with the following script:

``` time bash goniochromatic_appearance_video.sh ``` 

The full video (512x1024, 1024 spp) take around 702 minutes (11 hours). The video from the images can be created using the following command inside the video folder:

```  ffmpeg -framerate 15 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p goniochromatic_example.mp4 ``` 


### Figure 12 (Appearance Matching: Black Goose)

Execute the following script inside the folder scenes/appearance_matching_single_feather/black_goose_matching to create the renders for the appearance matching of the black goose feather:

``` time bash transparent_appearance.sh ``` 

This experiment takes around 4 minutes for 300x700 images and 256 spp.

### Relighting 

Activate the python_enviroment:

 ``` conda activate python_enviromnent``` 

We use the experiment in Figure 08 to create a relighting of the bird and the wing of feathers using the following command:

 ``` time python relight_feathers.py --feather_type=blue ``` 

The full video (2048x2048, 1024 spp) take around X minutes (X hours). The video from the images can be created using the following command inside the video folder:

```  ffmpeg -framerate 15 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p relighting_example.mp4 ``` 

Inside the script the type of feather can be selected (red, green, blue or black) and the type of light trajectory (azimuthal, polar). 

### TO DO

* 2D Fiber Tracer in Python
* Port solution to Mitsuba 3
* Fiber with medulla modeled as a subsurface scattering
* Iridescent feathers using thin film interference 

## Ackowledgements

* Physically-based renderer: [Mitsuba 0.6](https://github.com/mitsuba-renderer/mitsuba)
* Modelling tool: [Blender](https://www.blender.org/)  
* Initial steps: [Alvaro Romeo Arroyo TFG](https://zaguan.unizar.es/record/112309/files/TAZ-TFG-2021-4785.pdf)
* Hair scattering model: [Hair BSDF](https://www.pbrt.org/hair.pdf)
* Fur scattering model: [Fur BSDF](https://sites.cs.ucsb.edu/~lingqi/project_page/fur2/index.html)  
* Simulating 2D Light Transport: [Tantalum](https://benedikt-bitterli.me/tantalum/)  

## Contact

If you have any questions, please feel free to email the authors.   