---
name: MedicalImaging
topic: Medical Image Analysis
maintainer: Brandon Whitcher
email: bwhitcher@gmail.com
version: 2021-12-28
source: https://github.com/cran-task-views/MedicalImaging/
---


### Data Input/Output

*DICOM*

The industry standard format, for data coming off a clinical imaging
device, is [DICOM](http://dicom.nema.org) (Digital Imaging and
Communications in Medicine). The DICOM "standard" is very broad and
very complicated. Roughly speaking each DICOM-compliant file is a
collection of fields organized into two four-byte sequences
(group,element) that are represented as hexadecimal numbers and form a
*tag* . The (group,element) combination announces what type of
information is coming next. There is no fixed number of bytes for a
DICOM header. The final (group,element) tag should be the "data" tag
(7FE0,0010), such that all subsequent information is related to the
image(s).

-   The packages `r pkg("oro.dicom", priority = "core")`,
    `r pkg("divest", priority = "core")`,
    `r pkg("fmri", priority = "core")` and
    `r pkg("tractor.base", priority = "core")` (part of the
    `r gcode("tractor")` project) provide R functions that
    read DICOM files and facilitate their conversion to ANALYZE or NIfTI
    format.

*ANALYZE and NIfTI*

Although the industry standard for medical imaging data is DICOM,
another format has come to be heavily used in the image analysis
community. The ANALYZE format was originally developed in conjunction
with an image processing system (of the same name) at the Mayo
Foundation. An Anlayze (7.5) format image is comprised of two files, the
"hdr" and "img" files, that contain information about the
acquisition and the acquisition itself, respectively. A more recent
adaption of this format is known as
[NIfTI-1](http://nifti.nimh.nih.gov/nifti-1) and is a product of the
Data Format Working Group (DFWG) from the Neuroimaging Informatics
Technology Initiative (NIfTI). The NIfTI-1 data format is almost
identical to the ANALYZE format, but offers a few improvements: merging
of the header and image information into one file (.nii),
re-organization of the 348-byte fixed header into more relevant
categories and the possibility of extending the header information.

-   The packages `r pkg("RNifti", priority = "core")`,
    `r pkg("AnalyzeFMRI", priority = "core")`,
    `r pkg("fmri")`, `r pkg("tractor.base")`
    (part of the `r gcode("tractor")` project),
    `r pkg("oro.nifti", priority = "core")`, and
    `r pkg("neuroim", priority = "core")` all provide
    functions that read/write ANALYZE and NIfTI files.

### Magnetic Resonance Imaging (MRI)

*Diffusion Tensor Imaging (DTI)*

-   The `r pkg("tractor.base")` package (part of the
    [tractor project](http://www.tractor-mri.org.uk/) ) consists of
    functions for reading, writing and visualising MRI images. Images
    may be stored in ANALYZE, NIfTI or DICOM file formats, and can be
    visualised slice-by-slice or in projection. It also provides
    functions for common image manipulation tasks, such as masking and
    thresholding; and for applying arbitrary functions to image data.
    The package is written in pure R.

*Functional MRI*

-   `r pkg("adaptsmoFMRI")` contains R functions for
    estimating the blood oxygenation level dependent (BOLD) effect by
    using functional magnetic resonance imaging (fMRI) data, based on
    adaptive Gauss Markov random fields, for real as well as simulated
    data. Inference of the underlying models is performed by efficient
    Markov Chain Monte Carlo simulation, with the Metropolis Hastings
    algorithm for the non-approximate case and the Gibbs sampler for the
    approximate case. When comparing the results of approximate to the
    non-approximate version the outcome is in favour of the former, as
    the gain of accuracy in estimation, when not approximating, is
    minimal and the computational burden becomes less cumbersome.
-   `r pkg("AnalyzeFMRI")` is a package originally written
    for the processing and analysis of large structural and functional
    MRI data sets under the ANALYZE format. It has been updated to
    include new functionality: complete NIfTI input/output,
    cross-platform visualization based on Tcl/Tk components, and
    spatial/temporal ICA ( [Independent Components
    Analysis](http://en.wikipedia.org/wiki/Independent_component_analysis)
    ) via a graphical user interface (GUI).
-   The R package `r pkg("fmri")` provides tools for the
    analysis of functional MRI data. The core is the implementation of a
    new class of adaptive smoothing methods. These methods allow for a
    significant signal enhancement and reduction of false positive
    detections without, in contrast to traditional non-adaptive
    smoothing methods, reducing the effective spatial resolution. This
    property is especially of interest in the analysis of
    high-resolution functional MRI. The package includes functions for
    input/output of some standard imaging formats (ANALYZE, NIfTI, AFNI,
    DICOM) as well as for linear modelling the data and signal detection
    using [Random Field
    Theory](http://imaging.mrc-cbu.cam.ac.uk/imaging/PrinciplesRandomFields)
    . It also includes ICA and NGCA (non-Gaussian Components Analysis)
    based methods and hence has some overlap with
    `r pkg("AnalyzeFMRI")`.
-   Neuroimage is an R package (currently only available within the
    `r rforge("neuroim")` project on R-Forge) that provides
    data structures and input/output routines for functional brain
    imaging data. It reads and writes NIfTI-1 data and provides S4
    classes for handling multi-dimensional images.

*Structural MRI*

-   The package `r pkg("mritc", priority = "core")` provides
    tools for MRI tissue classification using normal mixture models and
    (partial volume, higher resolution) hidden Markov normal mixture
    models fitted by various methods. Functions to obtain initial values
    and spatial parameters are available. Facilities for visualization
    and evaluation of classification results are provided. To improve
    the speed, table lookup methods are used in various places,
    vectorization is used to take advantage of conditional independence,
    and some computations are performed by embedded C code.

*Visualization*

-   The package `r pkg("brainR")` includes functions for
    creating three-dimensional (3D) and four-dimensional (4D) images
    using WebGL, RGL, and JavaScript commands. This package relies on
    the X ToolKit ( [XTK](https://github.com/xtk/X#readme) ).
-   `r pkg("Morpho", priority = "core")` is a collection of
    tools for statistical shape analysis and visualization of point
    based shape representations (landmarks, meshes). Apart from the core
    functions such as General Procrustes Analysis and sliding of
    semi-landmarks, `r pkg("Morpho")` is sporting a variety
    of statistical procedures to assess group differences and asymmetry,
    most of them based on permutation/bootstrapping methods. For
    registration purposes there are functions to calculate landmark
    transforms (rigid, similarity, affine and thin-plate spline) as well
    as iterative closest point registration and automated alignment
    exploiting the shapes' principal axes. To deal with
    missing/erroneous data there are imputation methods available for
    missing landmarks and interactive outlier detection. For
    visualization there are functions to create interactive 3D plots of
    distance maps as well as visualizing differences between point
    clouds by deforming rectangular grids, both in 2D and 3D.
    Additionally, it includes an algorithm to retrodeform surface meshes
    representing structures that have suffered a series of locally
    affine deformations (e.g. fossils).
-   `r pkg("Rvcg", priority = "core")` interfaces
    [VCGLIB](http://vcg.sourceforge.net) to provide functions for
    manipulating triangular surface meshes; e.g., surfaces generated
    from medical image segmentations. Among those manipulations are
    quadric-edge collapse decimation, smoothing, subsampling, closest
    point search or uniform remeshing. Additionally it allows the
    generation of isosurfaces from 3D arrays. It has capabilities for
    import/export of STL, PLY and OBJ files, both in binary and ASCII
    format.

*Simulation*

-   The package `r pkg("neuRosim", priority = "core")` allows
    users to generate fMRI time series or 4D data. Some high-level
    functions are created for fast data generation with only a few
    arguments and a diversity of functions to define activation and
    noise. For more advanced users it is possible to use the low-level
    functions and manipulate the arguments.

### General Image Processing

-   `r pkg("adimpro", priority = "core")` is a package for 2D
    digital (color and B/W) images, actually not specific to medical
    imaging, but for general image processing.
-   The package `r pkg("bayesImageS")` implements several
    algorithms for segmentation of 2D and 3D images (such as CT and
    MRI). It provides full Bayesian inference for hidden Markov normal
    mixture models, including the posterior distribution for the
    smoothing parameter. The pixel labels can be sampled using
    checkerboard Gibbs or Swendsen-Wang. MCMC algorithms for the
    smoothing parameter include the approximate exchange algorithm
    (AEA), pseudolikelihood (PL), thermodynamic integration (TI), and
    approximate Bayesian computation (ABC-MCMC and ABC-SMC). An external
    field prior can be used when an anatomical atlas or other spatial
    information is available.
-   `r bioc("EBImage")` is an R package which provides
    general purpose functionality for the reading, writing, processing
    and analysis of images. Furthermore, in the context of
    microscopy-based cellular assays, this package offers tools to
    transform the images, segment cells and extract quantitative
    cellular descriptors.
-   The package `r pkg("mmand", priority = "core")`
    (Mathematical Morphology in Any Number of Dimensions) provides
    morphological operations like erode and dilate, opening and closing,
    as well as smoothing and kernel-based image processing. It operates
    on arrays or array-like data of arbitrary dimension.
-   The `r pkg("RNiftyReg", priority = "core")` provides an
    interface to the
    [NiftyReg](http://sourceforge.net/projects/niftyreg/) image
    registration tools. Rigid-body, affine and non-linear registrations
    are available and may be applied in 2D-to-2D, 3D-to-2D and 4D-to-3D
    procedures.
-   The package `r pkg("fslr")` contains wrapper functions
    that interface with the [FMRIB Sofware
    Library](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki) (FSL), a powerful
    and widely-used neuroimaging software library, using system
    commands. The goal with this package is to interface with FSL
    completely in R, where you pass R-based NIfTI objects and the
    function executes an FSL command and returns an R-based NIfTI
    object.

### Positron Emission Tomography (PET)

-   The `r pkg("occ", priority = "core")` package provides a
    generic function for estimating PET neuro-receptor occupancies by a
    drug, from the total volumes of distribution of a set of regions of
    interest (ROI). Fittings methods include the reference region, the
    *ordinary least squares* (OLS, sometimes known as "occupancy
    plot") and the *restricted maximum likelihood estimation* (REML).

### Electroencephalography (EEG)

-   `r pkg("edfReader", priority = "core")` reads some of the
    most popular file formats in EEG recordings.
-   The EEG package (currently only available within the
    `r rforge("eeg")` project on R-Forge) reads in single
    trial EEG (currently only ascii-exported pre-processed and trial
    segmented in Brain Vision Analyzer), computes averages (i.e.,
    event-related potentials or ERP's) and stores ERP's from multiple
    data sets in a `             data.frame           ` like object,
    such that statistical analysis (linear model, (M)ANOVA) can be done
    using the familiar R modeling frame work.
-   `r pkg("eegkit", priority = "core")` includes many useful
    functions for analysing EEG signals (among others, visualizing
    positions of electrodes).
-   `r pkg("PTAk")` is an R package that uses a multiway
    method to decompose a tensor (array) of any order, as a
    generalisation of a singular value decomposition (SVD) also
    supporting non-identity metrics and penalisations. A 2-way SVD with
    these extensions is also available. The package also includes
    additional multiway methods: PCAn (Tucker-n) and PARAFAC/CANDECOMP
    with these extensions. Applications include the analysis of EEG and
    functional MRI data.



### Links
-   [Journal of Statistical Software: Special Volume on Magnetic Resonance Imaging in R](https://www.jstatsoft.org/v44/)
-   [ANTsR is a framework that incorporates ITK and ANTs-based image processing methods into the R programming language.](http://picsl.upenn.edu/antsr)
-   [SimpleITK is a simplified layer built on top of ITK, intended to facilitate its use in rapid prototyping, education, interpreted languages. SimpleITK provides support for 2D and 3D images, and a selected set of pixel types for them. Different image filters may support a different collection of pixel types, in many cases due to computational requirements. The library is wrapped for interpreted languages by using SWIG. In particular, the following wrappings are available: Python, Java, Tcl, Lua, R and Ruby.](http://www.simpleitk.org/)
