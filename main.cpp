/*! \file main.cpp
 * \brief IceGrain's main file.
 */
/*! \mainpage IceGrain Documentation
 * \tableofcontents
 *
 * \section section_licence Licence
 *
 *  Copyright (c) 2013 Tobias Binder.
 *  
 *  This software was developed at the University of Heidelberg by
 *  Tobias Binder, Bjoern Andres, Thorsten Beier and Arthur Kuehlwein.
 *  Enquiries shall be directed to tobias.binder@iwr.uni-heidelberg.de.
 * 
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are met:
 * 
 *  - Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright notice, 
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  - All advertising materials mentioning features or use of this software must 
 *    display the following acknowledgement: ``This product includes the IceGrain
 *    package developed by Tobias Binder and others''.
 *  - The name of the author must not be used to endorse or promote products 
 *    derived from this software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
 *  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 *  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
 *  EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
 *  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
 *  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
 *  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
 *  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * \section synopsis Synopsis
 * The project \b IceGrain builds on the project \b CIS. Based on the segmentation created by \b CIS
 * and the computed boundary features, the microstructure of polar ice is extracted and physical parameters are determined.
 *
 * \section cargs Command line arguments
 * \b Ice Grain's <tt>main.cpp</tt> contains the following command line arguments. Note that there are command line arguments with a suffix \b -gui. These are specifically designed to be used with a GUI and thus not meant to be explicitly called from a terminal.
 * \subsection cargs_predim -predict-image
 * Extracts image structure. The suffix "image" means image batch mode.\n Usage: <tt>./IceGrain -predict-image boundary-features/*.bmp.bin suffix/ </tt>
 * \subsection cargs_predpa -predict-param
 * Extracts image structure in parameter file batch mode.\n Usage: <tt>./IceGrain -predict-param boundary-features/image.bmp.bin parameterfile.txt suffix/</tt>
 * \subsection cargs_correct -corrections
 * Allows correction of the extracted image structure via means of a GUI.\n Usage: <tt>./IceGrain -corrections boundary-features/image.bmp.bin suffix/</tt>
 * \subsection cargs_statistics -statistics
 * Extracts various statistics from an image.\n Usage: <tt>./IceGrain -statistics boundary-features/image.bmp.bin suffix/</tt>
 * \subsection cargs_depthprof -depth-profile
 * Creates a depth profile of the saved parameters.\n Usage: <tt>./IceGrain -depth-profile</tt>
 * \subsection cargs_subgrain -subgrains
 * Try to find subgrains in the given image structure.\n Usage: <tt>./IceGrain -subgrains boundary-features/image.bmp.bin suffix/</tt>
 * 
 * \section includes Libraries
 * This project uses the following (C++) libraries:
 * - <a href="http://www.cimg.sourceforge.net/">CIMG</a> enables easy use, processing and display of images. 
 * - <a href="http://hci.iwr.uni-heidelberg.de/vigra/">VIGRA</a> is a superb collection of algorithms developed in the HCI. It includes algorithms for image analysis and classification with random forests.
 * - <a href="http://www.hdfgroup.org/">HDF5</a> is a data format used for storing scientific data efficiently and flexible. It is used for random forests at the moment, but further data can be integrated.
 * - <a href="http://plplot.sourceforge.net/">plplot</a> is a library to create plots without using additional programs like GNUplot or \b Matlab.
 * 
 * \section dataStructures Data Structures
 * \subsection dataStructures_hdf5 HDF5 Structures
 * To save all obtained structures, the hierarchic data format \b HDF5 was chosen. It enables the storage of a multitude of data structures in a single file. Programs
 * like \b Matlab can read from <B>HDF5</B>-files directly.
 * The tool <a href="http://www.hdfgroup.org/hdf-java-html/hdfview/">HDFView</a> enables viewing and editing of <B>HDF4</B>- and <B>HDF5</B>-files. \b HDFView is available for Windows, Linux and Mac operating systems. In \b HDFView, data sets
 * can be copied into the clipboard or exported as text files.
 * Analogous to the processing of microstructure images, there are three data structures, which build on one another:
 * - \b Segmentation
 * - \b Grain-Boundary-Network
 * - \b Network-Parameters
 * 
 * Note that only the \b Grain-Boundary-Network and \b Network-Parameters data structures will be explained in detail here, as they are part of the \b IceGrain project.
 * The \b Segmentation data structure is part of the \b CIS, and thus further documentation about this structure can be found in the <B>CIS</B>'s code manual.
 * 
 * \subsubsection dataStructures_hdf5_gbn Grain-Boundary-Network
 * Information about the corresponding class (grain boundary, bubble edge, artefact) for each edge as well as which regions represent grains is save din the file <tt>rf-predictions/*.bmp.h5</tt>.
 * The data set "region_labels" (lower left-hand side of \b figure 1) assigns a new region index to each region. In line 14, the region 15 is assigned the index 24. The 
 * data sets "grain_center_of_mass" and "bubble_center_of_mass" state the coordinates of the center of mass of the new regions. Since line 23 has the coordinate (0,0) in both 
 * case, it is neither a grain nor a bubble, but an ignored border region. All other regions can be assigned to a grain or a bubble.
 * \image html hdf5_gbn_1.png "Figure 1: New region indices in the network structure"
 * The data sets "bubble_size" and "grain_size" have zeros in line 23, for all other regions the size here is in pixels.
 * 
 * The group "grain_arcs_index" contains a data set for each grain region. A data set for the entry in "region_labels" has to be opened if one searches a data set corresponding
 * to a grain region. The index 21 is assigned to the former region 54. The data set "21", which can be seen on the upper right-hand corner of <B>figure 2</B>, contains the 
 * seven edge indices of the grain with the center of mass (189, 663). The pixel coordinate of these edges can be found in the "objects" file. For bubble regions, the 
 * edge indices are saved in "bubble_arcs", respectively.
 * \image html hdf5_gbn_2.png "Figure 2: Grain boundary and node indices of the network structure"
 * Furthermore, there are the two data sets "grain_bubble_junctions" and "grain_junctions". The second node point selected in <B>figure 4</B> does neither form a grain
 * boundary triple-point, nor a grain boundary/bubble node in the network structure. Only the nodes mentiond in "grain_bubble_junctions" and "grain_junctions" show up in 
 * the network structure.
 * 
 * \subsubsection dataStructures_hdf5_param Network-Parameters
 * The extracted parameters of the \b Grain-Boundary-Network are compiled in the file <tt>statistics-results/*.bmp.*.h5</tt>. The data set "grain_areas" lists all grain region
 * indices, among them is region 21 (entry 20), which is representedin all other data sets with entry 10. In line 10 of "grain_arc_height" it is stated, that the grain region
 * 21 has a height of 2.08 mm; there is a 0 in line 10 of "grain_arc_number". Of course, each grain has at least one grain boundary, but stating those regions with contact
 * to the border region cannot be faciliated in an exact manner.
 * \image html hdf5_np_1.png "Figure 3: Description of the extracted parameters"
 * 
 * In order to parameterize the \b Grain-Boundary-Network, the segmentation edges are compiled, so that they begin and end on a grain boundary triple-point or a grain/bubble
 * node point. The edges 13 and 14 together, for example, form the network edge 37 (data set "grain_boundary_index/37") on the upper right-hand corner of <B>figure 3</B>. The "-"
 * prefix at the entry 14 states, that the order of pixels stated in the "objects" file has to be reversed in order to create a continuation of edge 13.
 * 
 * The data set "grain_junction_inde" corresponds in principle to the data set "neighborhood-1" of the "objects" file. Here however, only the grain boundary triple-points are listed
 * (stated in "grain_junctions" of the network file). The former assigns the newly compiled edges 53, 54 and 65 to the node 22.
 * The groups "grain_boundary_curvs" and "grain_boundary_angles" indicate a pixel-wise orientation of the normal vector and curvature for all new edges.
 * \image html hdf5_seg_2.png "Figure 4: Objects of the segmentation"
 * 
 */
#include <cgp/cgp_config.hxx>

#include <vigra/edgedetection.hxx>

const int nr_of_classes=3;
const float length_scaling=193.5; //20.0;//one mm is 193.5 Pixels
const float area_scaling=37444.0; //400.0;//one mm is 193.5 Pixels

#include "inc/ParameteredObject.hxx"

#include "path_functions.h"
#include "boundary_data_structure.h"
#include "boundary_probabilities.h"
#include "correct_prediction.h"
#include "statistics.h"
#include "subgrain.cpp"
#include "depth_profile.h"
#include "single_depth_profile.h"
#include "view.h"
#include "analyze.h"

SplitStream sout(std::cout);

int main(int argc, char *argv[])
{
    std::cout<<std::endl<<"IceGrain Calculation based on CIS"<<std::endl<<std::endl;
    std::string command_line_option=argv[1];

    /*
    OPTION -predict-image  (BATCH mode for images)
    */
    if (command_line_option=="-predict-image")
    {
        std::cout<<"Option -predict-image"<<std::endl;

        ParameterFile paramFile;

        if( !paramFile.load("parameters.txt") )//default parameter file
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

       	Parameter<std::string> p_watershed;
       	p_watershed.assign("", "path_watershed", "pixel-classification/watershed-segmentation/");
        p_watershed.load(paramFile,"config");
        std::string path_watershed=p_watershed;

       	Parameter<std::string> path_boundary_rf;
       	path_boundary_rf.assign("", "path_boundary_rf", "boundary-classification/random-forests/test");
        path_boundary_rf.load(paramFile,"config");

       	Parameter<std::string> p_rf_predictions;
       	p_rf_predictions.assign("", "path_rf_predictions", "rf-predictions/");
        p_rf_predictions.load(paramFile,"config");
        std::string path_rf_predictions=p_rf_predictions;

        Parameter<std::string> p_image;
        p_image.assign("","path_image", "images/");
        p_image.load(paramFile,"config");
        std::string path_image=p_image;

        std::string folder=argv[argc-1];
        if (folder=="no") folder="";
        path_watershed.append(folder.c_str());
        path_rf_predictions.append(folder.c_str());
        path_image.append(folder.c_str());

        Parameter<std::string> p_thumbs;
        p_thumbs.assign("", "path_thumbs", "no");
        p_thumbs.load(paramFile,"config");
        std::string path_thumbs=p_thumbs;
        if (path_thumbs!="no") path_thumbs.append(folder.c_str());

        int i=2;
        while(i<argc-1)
        {
            std::cout<<"Image "<<i-1<<"/"<<argc-3<<std::endl;
            extract_boundary_probabilities(argv[i],path_watershed,path_boundary_rf,path_rf_predictions,path_rf_predictions,
                "parameters.txt",paramFile,path_thumbs,path_image);
            i++;
        }
    }

    /*
    OPTION -predict-param  (BATCH mode for parameter files)
    */
    else if (command_line_option=="-predict-param")
    {
        std::cout<<"Option -predict-param"<<std::endl;

        ParameterFile paramFile;

        int i=3;
        while(i<argc-1)
        {
            bool default_rf_trained=false;

            std::cout<<"Parameter file "<<i-2<<" of " <<argc-4<<": "<<argv[i]<<std::endl;

            if( !paramFile.load(argv[i]) )
            {
                std::cout<<"Error: Parameter file could not be found!"<<std::endl;
                return 0;
            }

            Parameter<std::string> p_watershed;
            p_watershed.assign("", "path_watershed", "pixel-classification/watershed-segmentation/");
            p_watershed.load(paramFile,"config");
            std::string path_watershed=p_watershed;

            Parameter<std::string> path_boundary_rf;
            path_boundary_rf.assign("", "path_boundary_rf", "boundary-classification/random-forests/test");
            path_boundary_rf.load(paramFile,"config");

            Parameter<std::string> p_rf_predictions;
            p_rf_predictions.assign("", "path_rf_predictions", "rf-predictions/");
            p_rf_predictions.load(paramFile,"config");
            std::string path_rf_predictions=p_rf_predictions;

            std::string folder=argv[argc-1];
            if (folder=="no") folder="";
            path_watershed.append(folder.c_str());
            path_rf_predictions.append(folder.c_str());

            Parameter<std::string> p_thumbs;
            p_thumbs.assign("", "path_thumbs", "no");
            p_thumbs.load(paramFile,"config");
            std::string path_thumbs=p_thumbs;
            if (path_thumbs!="no") path_thumbs.append(folder.c_str());

            extract_boundary_probabilities(argv[2],path_watershed,path_boundary_rf,path_rf_predictions,path_rf_predictions,
                argv[i],paramFile,path_thumbs);
            i++;
        }
    }

    /*
    OPTION -predict-gui
    */
    else if (command_line_option=="-predict-gui")
    {
        ParameterFile paramFile;

        if( !paramFile.load(argv[argc-3]) )//parameter file defined in gui
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }
        
       	Parameter<std::string> path_boundary_rf;

        /* In order to avoid changing the signature of 'extract_boundary_probabilities', the path to the folders
         * containing the random forests and boundary features is stored in the third to last argument if the option
         * 'Default Classifier' is unchecked in the GUI, in which case 'FALSE' is appended to the end of the argument string.
         * The strings containing the two paths are seperated by a space character.
         */
        std::string arg_path_rf = argv[argc-4];
        std::string filepath_boundary_features;
        if(arg_path_rf.compare(arg_path_rf.size()-5, 5, "FALSE") == 0)
        {            
            arg_path_rf.resize(arg_path_rf.size()-5);            
            size_t found;
            std::string space = " ";
            found = arg_path_rf.find(space);
            if(found != std::string::npos)
            {
                filepath_boundary_features = arg_path_rf;
                filepath_boundary_features = filepath_boundary_features.substr(found+1, arg_path_rf.size()-found);
                std::cout << filepath_boundary_features << std::endl;

                path_boundary_rf = arg_path_rf.substr(0, found);
            }                
        }
        else if(arg_path_rf.compare(arg_path_rf.size()-4, 4, "TRUE") == 0)
        {            
            arg_path_rf.resize(arg_path_rf.size()-4);
            path_boundary_rf.assign("", "path_boundary_rf", "boundary-classification/random-forest/");
            path_boundary_rf.load(paramFile,"config");            
            
            filepath_boundary_features = arg_path_rf;
        }
        else std::cout << "This was not supposed to happen" << std::endl;

        std::string path_watershed=argv[argc-6];
        std::string path_rf_predictions=argv[argc-5];
        std::string path_thumbs=argv[argc-1];

        std::string folder=argv[argc-2];
        if (folder=="no") folder="";
        path_watershed.append(folder.c_str());
        path_rf_predictions.append(folder.c_str());
        if (path_thumbs!= "no") path_thumbs.append(folder.c_str());

        int i=2;
        while(i<argc-6)
        {
            std::string filepath_boundary_features_loop = filepath_boundary_features;
            filepath_boundary_features_loop.append(folder.c_str());
            filepath_boundary_features_loop.append(get_filename(argv[i]));
            filepath_boundary_features_loop.append(".bin");

            extract_boundary_probabilities(filepath_boundary_features_loop,path_watershed,path_boundary_rf(),path_rf_predictions,
                                           path_rf_predictions,argv[argc-3],paramFile,path_thumbs,get_path(argv[i]));
            i++;
        }
    }

    /*
    OPTION -subgrains
    */
    else if (command_line_option=="-subgrains")
    {
        std::cout<<"Option -subgrains"<<std::endl;

        ParameterFile paramFile;

        if( !paramFile.load("parameters.txt") )//default parameter file
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

        Parameter<std::string> p_image;
        p_image.assign("", "path_image", "/home/tobinder/Data/");
        p_image.load(paramFile,"config");
        std::string path_image=p_image;

        Parameter<std::string> p_watershed;
        p_watershed.assign("", "path_watershed", "pixel-classification/watershed-segmentation/");
        p_watershed.load(paramFile,"config");
        std::string path_watershed=p_watershed;

        Parameter<std::string> p_rf_predictions;
        p_rf_predictions.assign("", "path_rf_predictions", "rf-predictions/");
        p_rf_predictions.load(paramFile,"config");
        std::string path_rf_predictions=p_rf_predictions;

        Parameter<std::string> p_subGB_rf_predictions;
        p_subGB_rf_predictions.assign("", "path_subGB_rf_predictions", "rf-predictions/");
        p_subGB_rf_predictions.load(paramFile,"config");
        std::string path_subGB_rf_predictions=p_subGB_rf_predictions;

        Parameter<std::string> p_results;
        p_results.assign("", "path_results", "statistics-results/");
        p_results.load(paramFile,"config");
        std::string path_results=p_results;

        std::string folder=argv[argc-1];
        if (folder=="no") folder="";
        path_image.append(folder.c_str());
        path_watershed.append(folder.c_str());
        path_rf_predictions.append(folder.c_str());
        path_subGB_rf_predictions.append(folder.c_str());
        path_results.append(folder.c_str());

        Parameter<int> low_grain_size;
        low_grain_size.assign("", "low_grain_size", 0);
        low_grain_size.load(paramFile,"config");

        Parameter<int> min_bubble_distance;
        min_bubble_distance.assign("", "min_bubble_distance", 0);
        min_bubble_distance.load(paramFile,"config");

        find_subgrains(argv[argc-2],path_image,path_watershed,path_rf_predictions,path_subGB_rf_predictions,path_results,
                       low_grain_size,min_bubble_distance,"parameters.txt",paramFile);
    }

    /*
    OPTION -subgrains-gui
    */
    else if (command_line_option=="-subgrains-gui")
    {
        std::cout<<"Option -subgrains"<<std::endl;

        ParameterFile paramFile;

        if( !paramFile.load(argv[argc-6]) )//default parameter file
        {            
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

        std::string suffix = argv[argc-7];        
        std::string p_image = argv[argc-5];
        std::string p_watershed = argv[argc-4];
        std::string p_rf_predictions = argv[argc-3];
        std::string p_subGB_rf_predictions = argv[argc-2];
        std::string p_results = argv[argc-1];

        //Append suffix if there is one    
        if(suffix.compare("no") != 0)
        {                                
                p_watershed.append(suffix);
                p_rf_predictions.append(suffix);
                p_subGB_rf_predictions.append(suffix);
                p_results.append(suffix);
        }     
      
        Parameter<int> low_grain_size;
        low_grain_size.assign("", "low_grain_size", 0);
        low_grain_size.load(paramFile,"config");

        Parameter<int> min_bubble_distance;
        min_bubble_distance.assign("", "min_bubble_distance", 0);
        min_bubble_distance.load(paramFile,"config");
        
        int i=2;
        while(i<argc-8)
        {
            std::string name_current = argv[i];
            name_current.append(".bin");

            std::string p_bin = argv[argc-8];
            if(suffix.compare("no") == 0)
            {                
                p_bin.append(name_current);
            }
            else
            {
                p_bin.append(suffix);
                p_bin.append(name_current);
            }                
                        
            find_subgrains(p_bin,p_image,p_watershed,p_rf_predictions,p_subGB_rf_predictions,p_results,
                        low_grain_size,min_bubble_distance,"parameters.txt",paramFile);
            i++;
        }
    }

    /*
    OPTION -corrections
    */
    else if (command_line_option=="-corrections")
    {
        std::cout<<"Option -corrections"<<std::endl;

        ParameterFile paramFile;

        if( !paramFile.load("parameters.txt") )//default parameter file
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

       	Parameter<std::string> p_image;
       	p_image.assign("", "path_image", "/home/tobinder/Data/");
        p_image.load(paramFile,"config");
        std::string path_image=p_image;

       	Parameter<std::string> p_watershed;
       	p_watershed.assign("", "path_watershed", "pixel-classification/watershed-segmentation/");
        p_watershed.load(paramFile,"config");
        std::string path_watershed=p_watershed;

       	Parameter<std::string> p_rf_predictions;
       	p_rf_predictions.assign("", "path_rf_predictions", "rf-predictions/");
        p_rf_predictions.load(paramFile,"config");
        std::string path_rf_predictions=p_rf_predictions;

        std::string folder=argv[argc-1];
        if (folder=="no") folder="";
        path_image.append(folder.c_str());
        path_watershed.append(folder.c_str());
        path_rf_predictions.append(folder.c_str());

        correct_prediction(argv[argc-2],path_image,path_watershed,path_rf_predictions,"parameters.txt",paramFile);
    }

    /*
    OPTION -border-corrections-gui
    */
    else if (command_line_option=="-border-corrections-gui")
    {
        ParameterFile paramFile;

        if( !paramFile.load(argv[argc-2]) )//parameter file defined in gui
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

        std::string path_image=get_path(argv[argc-5]);
        std::string path_watershed=argv[argc-4];

        std::string folder=argv[argc-1];
        if (folder=="no") folder="";
        path_watershed.append(folder.c_str());

        int i=2;
        while(i<argc-4)
        {
            std::string filepath_boundary_features=argv[argc-3];
            filepath_boundary_features.append(folder.c_str());
            filepath_boundary_features.append(get_filename(argv[i]));
            filepath_boundary_features.append(".bin");
            correct_prediction(filepath_boundary_features, path_image, path_watershed, "", argv[argc-2], paramFile, 1);
            i++;
        }
    }

    /*
    OPTION -corrections-gui
    */
    else if (command_line_option=="-corrections-gui")
    {
        ParameterFile paramFile;

        if( !paramFile.load(argv[argc-2]) )//parameter file defined in gui
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

        std::string path_image=get_path(argv[argc-6]);
        std::string path_watershed=argv[argc-5];
        std::string path_rf_predictions=argv[argc-4];

        std::string folder=argv[argc-1];
        if (folder=="no") folder="";
        path_watershed.append(folder.c_str());
        path_rf_predictions.append(folder.c_str());

        int i=2;
        while(i<argc-5)
        {
            std::string filepath_boundary_features=argv[argc-3];
            filepath_boundary_features.append(folder.c_str());
            filepath_boundary_features.append(get_filename(argv[i]));
            filepath_boundary_features.append(".bin");
            correct_prediction(filepath_boundary_features, path_image, path_watershed, path_rf_predictions, argv[argc-2], paramFile, 2);
            i++;
        }
    }

    /*
    OPTION -statistics
    */
    else if (command_line_option=="-statistics")
    {
        std::cout<<"Option -statistics"<<std::endl;

        ParameterFile paramFile;

        if( !paramFile.load("parameters.txt") )//default parameter file
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

       	Parameter<std::string> p_watershed;
       	p_watershed.assign("", "path_watershed", "pixel-classification/watershed-segmentation/");
        p_watershed.load(paramFile,"config");
        std::string path_watershed=p_watershed;

       	Parameter<std::string> p_rf_predictions;
       	p_rf_predictions.assign("", "path_rf_predictions", "rf-predictions/");
        p_rf_predictions.load(paramFile,"config");
        std::string path_rf_predictions=p_rf_predictions;

       	Parameter<std::string> p_plots;
       	p_plots.assign("", "path_plots", "plots/");
        p_plots.load(paramFile,"config");
        std::string path_plots=p_plots;

       	Parameter<std::string> p_results;
       	p_results.assign("", "path_results", "statistics-results/");
        p_results.load(paramFile,"config");
        std::string path_results=p_results;

        std::string folder=argv[argc-1];
        if (folder=="no") folder="";
        path_watershed.append(folder.c_str());
        path_rf_predictions.append(folder.c_str());
        path_plots.append(folder.c_str());
        path_results.append(folder.c_str());

        do_statistics(argv[argc-2],path_watershed,path_rf_predictions,path_plots,path_results,"parameters.txt", paramFile);
    }

    /*
    OPTION -statistics-gui
    */
    else if (command_line_option=="-statistics-gui")
    {
        ParameterFile paramFile;

        if( !paramFile.load(argv[argc-2]) )//parameter file defined in gui
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

        Parameter<int> low_grain_size;
        low_grain_size.assign("", "low_grain_size", 0);
        low_grain_size=atoi(argv[argc-12]);
        low_grain_size.save(paramFile,"config");

       	Parameter<int> high_grain_size;
       	high_grain_size.assign("", "high_grain_size", 20000);
        high_grain_size=atoi(argv[argc-11]);
        high_grain_size.save(paramFile,"config");

       	Parameter<int> grain_size_step;
       	grain_size_step.assign("", "grain_size_step", 5000);
        grain_size_step=atoi(argv[argc-10]);
        grain_size_step.save(paramFile,"config");

       	Parameter<int> min_bubble_distance;
       	min_bubble_distance.assign("", "min_bubble_distance", 0);
        min_bubble_distance=atoi(argv[argc-9]);
        min_bubble_distance.save(paramFile,"config");

       	Parameter<int> grain_step;
       	grain_step.assign("", "grain_step", 100);
        grain_step=atoi(argv[argc-8]);
        grain_step.save(paramFile,"config");

        Parameter<int> grain_size_min;
        grain_size_min.assign("", "grain_size_min", 4);
        grain_size_min=4;
        grain_size_min.save(paramFile,"config");

        std::string path_watershed=argv[argc-7];
        std::string path_rf_predictions=argv[argc-6];
        std::string path_plots=argv[argc-5];
        std::string path_results=argv[argc-4];

        std::string folder=argv[argc-1];
        if (folder=="no") folder="";
        path_watershed.append(folder.c_str());
        path_rf_predictions.append(folder.c_str());
        path_plots.append(folder.c_str());
        path_results.append(folder.c_str());

        int i=2;
        while(i<argc-12)
        {
            std::string filepath_boundary_features=argv[argc-3];
            filepath_boundary_features.append(folder.c_str());
            filepath_boundary_features.append(get_filename(argv[i]));
            filepath_boundary_features.append(".bin");

            do_statistics(filepath_boundary_features,path_watershed,path_rf_predictions,path_plots,path_results,argv[argc-2],paramFile);
            i++;
        }
    }

    /*
    OPTION -depth-profile
    */
    else if (command_line_option=="-depth-profile")
    {
        std::cout<<"Option -depth-profile"<<std::endl;

        ParameterFile paramFile;

        if( !paramFile.load("parameters.txt") )//default parameter file
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

       	Parameter<std::string> path_results;
       	path_results.assign("", "path_results", "statistics-results/");
        path_results.load(paramFile,"config");

        //single_depth_profile(path_results, paramFile);
        depth_profile(path_results, paramFile, 22);
    }

    /*
    OPTION -depth-profile-gui
    */
    else if (command_line_option=="-depth-profile-gui")
    {
        ParameterFile paramFile;

        if( !paramFile.load(argv[argc-1]) )//parameter file defined in gui
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

       	Parameter<int> nr_depths;
       	nr_depths.assign("", "nr_depths", 113);
        nr_depths=atoi(argv[argc-8]);
        nr_depths.save(paramFile,"config");

        Parameter<int> low_grain_size;
        low_grain_size.assign("", "low_grain_size", 0);
        low_grain_size=atoi(argv[argc-7]);
        low_grain_size.save(paramFile,"config");

       	Parameter<int> high_grain_size;
       	high_grain_size.assign("", "high_grain_size", 20000);
        high_grain_size=atoi(argv[argc-6]);
        high_grain_size.save(paramFile,"config");

       	Parameter<int> grain_size_step;
       	grain_size_step.assign("", "grain_size_step", 5000);
        grain_size_step=atoi(argv[argc-5]);
        grain_size_step.save(paramFile,"config");

       	Parameter<int> min_bubble_distance;
       	min_bubble_distance.assign("", "min_bubble_distance", 0);
        min_bubble_distance=atoi(argv[argc-4]);
        min_bubble_distance.save(paramFile,"config");

       	Parameter<int> grain_step;
       	grain_step.assign("", "grain_step", 100);
        grain_step=atoi(argv[argc-3]);
        grain_step.save(paramFile,"config");

        depth_profile(argv[argc-2], paramFile, atof(argv[argc-9]));
    }

    /*
    OPTION -cloudy-bands
    */
    else if (command_line_option=="-cloudy-bands")
    {
        if (argc<17)
        {
            std::cout<<"Error: option -cloudy-bands needs more parameters!"<<std::endl;
            return 0;
        }
        std::string path_rf_predictions=argv[14];
        std::string path_plots=argv[15];

        std::string folder=argv[argc-1];
        if (folder=="no") folder="";
        path_rf_predictions.append(folder.c_str());
        path_plots.append(folder.c_str());

        cloudy_bands(argv[2], argv[8], argv[3], argv[9], argv[4], argv[10], argv[5], argv[11], argv[6], argv[12], argv[7], argv[13],
                     path_rf_predictions, path_plots);
    }

    /*
    OPTION -double-grain-size
    */
    else if (command_line_option=="-double-grain-size")
    {
        if (argc<8)
        {
            std::cout<<"Error: option -cloudy-bands needs more parameters!"<<std::endl;
            return 0;
        }

        int minimal_grain_size=atoi(argv[2]);
        int minimal_bubble_distance=atoi(argv[3]);

        std::string suffix1=argv[4];
        suffix1.append("/");
        std::string suffix2=argv[5];
        suffix2.append("/");

        std::string path_results=argv[6];
        std::string path_plots=argv[7];

        double_grain_size(path_results, suffix1, suffix2, minimal_bubble_distance, minimal_grain_size, path_plots);
    }

    /*
    OPTION -view
    */
    else if (command_line_option=="-view")
    {
        std::cout<<"Option -view"<<std::endl;

        ParameterFile paramFile;

        if( !paramFile.load("parameters.txt") )//default parameter file
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

        Parameter<std::string> p_image;
        p_image.assign("","path_image", "images/");
        p_image.load(paramFile,"config");
        std::string path_image=p_image;

       	Parameter<std::string> p_watershed;
       	p_watershed.assign("", "path_watershed", "pixel-classification/watershed-segmentation/");
        p_watershed.load(paramFile,"config");
        std::string path_watershed=p_watershed;

       	Parameter<std::string> p_rf_predictions;
       	p_rf_predictions.assign("", "path_rf_predictions", "rf-predictions/");
        p_rf_predictions.load(paramFile,"config");
        std::string path_rf_predictions=p_rf_predictions;

       	Parameter<std::string> p_plots;
       	p_plots.assign("", "path_plots", "plots/");
        p_plots.load(paramFile,"config");
        std::string path_plots=p_plots;

       	Parameter<std::string> p_results;
       	p_results.assign("", "path_results", "statistics-results/");
        p_results.load(paramFile,"config");
        std::string path_results=p_results;

        std::string folder=argv[argc-1];
        if (folder=="no") folder="";
        path_image.append(folder.c_str());
        path_watershed.append(folder.c_str());
        path_rf_predictions.append(folder.c_str());
        path_plots.append(folder.c_str());
        path_results.append(folder.c_str());

        int minimal_bubble_distance=0;
        int minimal_grain_size=0;
        int mode=1;
        int start_value=0;

        view(argv[argc-2],path_image,path_watershed,path_rf_predictions,path_plots,path_results,"parameters.txt", paramFile, minimal_grain_size,
             minimal_bubble_distance, mode, start_value);
    }

    /*
    OPTION -view-gui
    */
    else if (command_line_option=="-view-gui")
    {
        ParameterFile paramFile;

        if( !paramFile.load(argv[argc-2]) )//parameter file defined in gui
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

        std::string path_image=get_path(argv[2]);
        std::string path_watershed=argv[argc-11];
        std::string path_rf_predictions=argv[argc-10];
        std::string path_plots=argv[argc-9];
        std::string path_results=argv[argc-8];
        std::string filepath_boundary_features=argv[argc-7];

        std::string folder=argv[argc-1];
        if (folder=="no") folder="";
        path_watershed.append(folder.c_str());
        path_rf_predictions.append(folder.c_str());
        path_plots.append(folder.c_str());
        path_results.append(folder.c_str());
        filepath_boundary_features.append(folder.c_str());

        int low_grain_size=atoi(argv[argc-6]);
        int minimal_bubble_distance=atoi(argv[argc-5]);
        int mode=atoi(argv[argc-4]);
        int start_value=atoi(argv[argc-3]);

        filepath_boundary_features.append(get_filename(argv[2]));
        filepath_boundary_features.append(".bin");

        view(filepath_boundary_features, path_image, path_watershed, path_rf_predictions, path_plots, path_results, argv[argc-2], paramFile, low_grain_size,
             minimal_bubble_distance, mode, start_value);

    }

    /*
    OPTION -analyze
    */
    else if (command_line_option=="-analyze")
    {
        ParameterFile paramFile;

        if( !paramFile.load(argv[7]) )//parameter file defined in gui
        {
            std::cout<<"Error: Parameter file could not be found!"<<std::endl;
            return 0;
        }

        std::string filepath_list=argv[2];
        std::string path_rf_predictions=argv[3];
        std::string path_results=argv[4];
        int low_grain_size=atoi(argv[5]);
        int minimal_bubble_distance=atoi(argv[6]);
        std::string path_plots=argv[8];

        int mode=atoi(argv[9]);
        int combo1=atoi(argv[10]);
        int spin=atoi(argv[11]);
        int combo2=atoi(argv[12]);
        std::string path_watershed=argv[13];
        bool correct_suffix=atoi(argv[14]);

        analyze(filepath_list, path_rf_predictions, path_results, argv[7], paramFile, low_grain_size, minimal_bubble_distance, path_plots, mode, combo1, spin,
            combo2, path_watershed, correct_suffix);
    }

    /*
    OPTION -new-depth-profile
    */
    else if (command_line_option=="-new-depth-profile")
    {
        std::string filepath_list=argv[2];
        std::string path_rf_predictions=argv[3];
        std::string path_results=argv[4];
        int low_grain_size=atoi(argv[5]);
        int minimal_bubble_distance=atoi(argv[6]);
        bool correct_suffix=atoi(argv[8]);

        int high_grain_size=atoi(argv[9]);
        int grain_size_step=atoi(argv[10]);
        float depth_bin_width=atof(argv[11]);

        new_depth_profiles(filepath_list, path_rf_predictions, argv[7], minimal_bubble_distance, path_results, correct_suffix, low_grain_size,
            high_grain_size, grain_size_step, depth_bin_width);
    }

    return 0;
}
