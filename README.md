# Cell_Kinetics_CellposeSAM_TrackMate

Workflow for converting live-imaging time-lapse data into single-cell trajectories and basic cell-kinetic measurements using **Cellpose-SAM**, **TrackMate in Fiji/ImageJ**, and **MATLAB**.

This repository provides a stepwise example pipeline:

1. Start with a time-lapse image sequence.
2. Segment cells in each frame using Cellpose-SAM.
3. Optionally shrink segmentation masks to improve TrackMate linking.
4. Use TrackMate in Fiji/ImageJ to generate cell tracks.
5. Analyze trajectories in MATLAB to quantify cell movement, persistence, splitting behavior, and mean-squared displacement.

---

## Repository structure

```text
Cell_Kinetics_CellposeSAM_TrackMate/
│
├── Sample_input/
│   └── Example time-lapse image frames in .tif format
│
├── Step1_CellposeSAM/
│   ├── run_Cellpose_SAM_ImageSequence.ipynb
│   └── Example Cellpose-SAM segmentation masks
│
├── Step2_ShrinkMask/
│   ├── shrink_label_masks_in_folder.m
│   └── Example shrunk label masks
│
├── Step3_TrackMate_in_FIJI/
│   ├── 1_Image_spots.csv
│   ├── 1_Image_edges.csv
│   └── Video_Combined Stacks-1.avi
│
├── Step4_Cell_Kinetics/
│   ├── 1_left_epi_roi.csv
│   ├── 1_right_epi_roi.csv
│   └── Trajectory_analysis_plot.m
│
├── LICENSE
└── README.md
```

---

## Workflow overview

### Step 0. Prepare the input image sequence

Place the time-lapse image frames in `Sample_input/`.

The example dataset contains sequential `.tif` frames named with `rot0000`, `rot0001`, etc. Each file represents one frame of the time-lapse movie.

Expected input format:

```text
frame0000.tif
frame0001.tif
frame0002.tif
...
```

The pipeline assumes that the input frames are 2D images or projected 2D images from live-imaging data.

---

### Step 1. Segment cells using Cellpose-SAM

Folder:

```text
Step1_CellposeSAM/
```

Main file:

```text
run_Cellpose_SAM_ImageSequence.ipynb
```

This notebook runs Cellpose-SAM on an image sequence and generates one segmentation label image per frame.

Expected output:

```text
*_masks.tif
```

Each mask image is a label image in which:

* background pixels are labeled as `0`
* each segmented object has a unique integer label
* the label image corresponds to one time point of the movie

The example output files are already included in `Step1_CellposeSAM/`.

---

### Step 2. Shrink segmentation masks before tracking

Folder:

```text
Step2_ShrinkMask/
```

Main file:

```text
shrink_label_masks_in_folder.m
```

This MATLAB function reads Cellpose-SAM label masks and shrinks each labeled object around its centroid. This can help TrackMate detect and link objects more robustly by reducing mask overlap or boundary contact between neighboring cells.

Example MATLAB usage:

```matlab
inDir = 'Step1_CellposeSAM';
outDir = 'Step2_ShrinkMask';
scaleFactor = 0.5;

shrink_label_masks_in_folder(inDir, outDir, scaleFactor);
```

Parameters:

* `inDir`: folder containing the original Cellpose-SAM mask `.tif` files
* `outDir`: folder where shrunk masks will be saved
* `scaleFactor`: fraction by which each mask is scaled around its centroid

  * `0.5` shrinks each mask to 50% of its original size
  * `1.0` keeps the mask unchanged

The script preserves the original label values and writes new `.tif` label images to the output folder.

---

### Step 3. Track cells using TrackMate in Fiji/ImageJ

Folder:

```text
Step3_TrackMate_in_FIJI/
```

Example files:

```text
1_Image_spots.csv
1_Image_edges.csv
Video_Combined Stacks-1.avi
```

Use Fiji/ImageJ and TrackMate to track cells from the segmentation masks or processed movie.

The key TrackMate outputs needed for downstream analysis are:

```text
1_Image_spots.csv
1_Image_edges.csv
```

The spots file contains detected object positions and time points.
The edges file contains the links between spots across time, which define the cell trajectories.

When exporting from TrackMate, make sure to export both:

* **Spots statistics**
* **Edges statistics**

The MATLAB analysis script assumes that the TrackMate CSV files contain standard columns such as spot ID, track ID, position, frame/time, source spot ID, and target spot ID.

---

### Step 4. Analyze cell trajectories and kinetics

Folder:

```text
Step4_Cell_Kinetics/
```

Main file:

```text
Trajectory_analysis_plot.m
```

Example input files:

```text
1_Image_spots.csv
1_Image_edges.csv
1_left_epi_roi.csv
1_right_epi_roi.csv
```

This MATLAB script analyzes TrackMate trajectories within a user-defined region of interest.

The script performs several operations:

1. Loads TrackMate spots and edges CSV files.
2. Loads an ROI polygon from a CSV file.
3. Selects tracks whose starting point lies inside the ROI.
4. Builds a directed graph from TrackMate spot-linking information.
5. Detects track-splitting events.
6. Decomposes branching tracks into root-to-leaf paths.
7. Computes trajectory-level kinetic measurements.
8. Plots tracks, splitting points, velocity, persistence, and MSD.

Key user settings are near the top of the script:

```matlab
spotsFile = "1_Image_spots.csv";
edgesFile = "1_Image_edges.csv";
roiFile   = "1_left_epi_roi.csv";

pixelSize = 0.36;      % microns per pixel
frameInterval = 5;     % minutes per frame
xDirection = 1;        % direction convention for projected velocity
maxTime = 128;         % maximum time/frame cutoff
```

Modify these values to match your own imaging experiment.

---

## ROI file format

The ROI files should be CSV files containing polygon coordinates.

Example:

```text
X,Y
120,80
180,85
200,150
130,160
120,80
```

The MATLAB script uses the ROI polygon to select tracks whose first detected position starts inside the ROI.

The repository includes two example ROI files:

```text
1_left_epi_roi.csv
1_right_epi_roi.csv
```

These can be used to analyze different tissue regions separately.

---

## Main outputs from trajectory analysis

The MATLAB script generates plots and quantitative measurements including:

### Track selection plot

Shows the ROI boundary and the starting positions of selected tracks.

### Time-coded trajectory plot

Shows selected trajectories colored by time and marks splitting points.

### Splitting analysis

For each selected track, the script counts nodes with more than one outgoing edge as splitting events.

### Projected velocity

The script estimates a principal direction of motion and projects track velocity along that direction.

### Persistence

Persistence is calculated as:

```text
net displacement / total path length
```

Values close to `1` indicate straighter motion.
Values close to `0` indicate more wandering or tortuous motion.

### Mean-squared displacement

The script calculates ensemble MSD across time lags and fits the relationship:

```text
MSD = 4D × t^α
```

where:

* `D` is an apparent diffusion/dispersion coefficient
* `α` describes the scaling of motion with time

Interpretation of `α`:

* `α ≈ 1`: random-walk-like motion
* `α > 1`: persistent or directed motion
* `α < 1`: constrained or subdiffusive motion

---

## Software requirements

### Python / notebook environment

Used for Cellpose-SAM segmentation.

Recommended:

* Python
* Jupyter Notebook or Google Colab
* Cellpose / Cellpose-SAM
* NumPy
* tifffile or equivalent image I/O package

GPU acceleration is recommended for Cellpose-SAM.

### Fiji/ImageJ

Used for TrackMate-based cell tracking.

Recommended:

* Fiji
* TrackMate plugin

### MATLAB

Used for mask processing and trajectory analysis.

Required MATLAB functionality includes:

* image reading/writing
* affine image warping
* table import
* graph/digraph analysis
* plotting

The mask-shrinking script uses functions such as `imread`, `imwrite`, `imref2d`, `affine2d`, and `imwarp`.

---

## How to run the full example pipeline

### 1. Segment input images

Open:

```text
Step1_CellposeSAM/run_Cellpose_SAM_ImageSequence.ipynb
```

Run the notebook on the images in:

```text
Sample_input/
```

Save the resulting mask files as:

```text
*_masks.tif
```

---

### 2. Shrink masks in MATLAB

Run:

```matlab
shrink_label_masks_in_folder('Step1_CellposeSAM', 'Step2_ShrinkMask', 0.5);
```

This creates shrunk label masks for TrackMate.

---

### 3. Track cells in Fiji/ImageJ

Open the movie or mask sequence in Fiji/ImageJ.

Use TrackMate to detect and link cells across time.

Export:

```text
1_Image_spots.csv
1_Image_edges.csv
```

Place these files in:

```text
Step3_TrackMate_in_FIJI/
```

or update the MATLAB paths accordingly.

---

### 4. Analyze trajectories in MATLAB

Copy or reference the TrackMate CSV files and ROI file in the analysis folder.

Open:

```text
Step4_Cell_Kinetics/Trajectory_analysis_plot.m
```

Edit the file paths and imaging parameters:

```matlab
spotsFile = "1_Image_spots.csv";
edgesFile = "1_Image_edges.csv";
roiFile   = "1_left_epi_roi.csv";

pixelSize = 0.36;
frameInterval = 5;
maxTime = 128;
```

Run the script.

---

## Notes and limitations

* The example analysis assumes TrackMate CSV files with standard column names.
* TrackMate CSV files may contain metadata rows before the header; the MATLAB script expects the header on line 4.
* The ROI filter selects tracks based on the starting point of each track, not every later position.
* The current analysis is designed for 2D trajectories.
* If the image calibration or frame interval differs from the example dataset, update `pixelSize` and `frameInterval` before interpreting velocity or MSD.
* The mask-shrinking step is optional but can improve tracking when adjacent cell masks touch or overlap.
* The MSD analysis is descriptive and should be interpreted carefully for branching, dividing, or highly constrained cell trajectories.

---

## Suggested citation / acknowledgment

If you use this workflow, please cite or acknowledge the relevant software tools:

* Cellpose / Cellpose-SAM for cell segmentation
  GitHub: https://github.com/MouseLand/cellpose
* Fiji/ImageJ and TrackMate for cell tracking
* MATLAB for downstream trajectory analysis

---

## License

This repository is distributed under the MIT License.
