/*! \file seg.h 
 * \brief Class representing the Segmentation structure.
 */
/*! 
 * The macro <tt>SEG_STANDALONE</tt> enables the use of this class in standalone mode. However, the header file <tt>boundary_data_structures.h</tt> must be available for this to work.
 */
// IceGrain: Extraction and parameterization of grain boundary networks of ice
//
// Copyright (c) 2013 Tobias Binder.
// 
// This software was developed at the University of Heidelberg by
// Tobias Binder, Bjoern Andres, Thorsten Beier and Arthur Kuehlwein.
// Enquiries shall be directed to tobias.binder@iwr.uni-heidelberg.de.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// - Redistributions in binary form must reproduce the above copyright notice, 
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// - All advertising materials mentioning features or use of this software must 
//   display the following acknowledgement: ``This product includes the IceGrain
//   package developed by Tobias Binder and others''.
// - The name of the author must not be used to endorse or promote products 
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#pragma once

#include <stdio.h>
//#include <omp.h>

#include "CImg.h"

#include <iostream>
#include <vector>
#include <cgp/cgp_hdf5.hxx>
#include "marray.hxx"

#ifdef SEG_STANDALONE
//This include becomes necessary only if this class is included seperately in a project
#include "boundary_data_structure.h"
#endif

/*Class representing the Segmentation structure (cf. code manual for CIS)*/

class seg 
{
    public:
    //Attributes *************************************************************
    vigra::BasicImage<unsigned int>    ws_region_image;                 /*!< Watershed region image.*/
    marray::Marray<unsigned int>       one_boundings;                   /*!< One boundings.*/
    marray::Marray<unsigned int>       two_boundings;                   /*!< Two boundings.*/
    std::vector< std::vector<point> >  arcs;                            /*!< Segmentation arcs.*/
    std::vector<point>                 junctions;                       /*!< Segmentation junctions.*/
    int                                dim_x;                           /*!< Image width.*/
    int                                dim_y;                           /*!< Image height.*/
    bool                               two_pixel_boundary;              /*!< Flag for drawing two pixel wide boundaries.*/
    size_t                             numberOfBinningGroups;           /*!< Number of binning groups.*/
    
    //Methods ****************************************************************
    /*! Basic constructor. Creates attributes with empty data and sets the number of binning groups to 1000.
     * \param twoPixelBoundary Two pixel boundary flag
     */
    //Basic constructor with two pixel boundary flag, creates attributes with emtpy data
    seg(bool twoPixelBoundary);     
    /*! Load a watershed region image from a given HDF5 file.
     * \param filepath_to_ws_region_image Filepath to the watershed region image stored inside a HDF5 file
     */                         
    //Loads data from a *.h5 file
    void load_ws_image_from_file(std::string filepath_to_ws_region_image);
    /*! Load the remaining attributes from a HDF5 file.
     * \param filepath_to_ws_region_image Filepath to the HDF5 file
     */
    //Loads data from a *.objects.h5 file
    void load_obj_from_file(std::string filepath_to_ws_region_image);
    /*! Loads all data from a given HDF5 file into the class. Calls both <tt>load_ws_image_from_file()</tt> and <tt>load_obj_from_file()</tt>.
     * \param filepath_to_ws_region_image Filepath to the HDF5 file
     */
    //Creates a cgp data structure from existing *.h5 and *.objects.h5 files. Uses the functions "load_ws_image_from_file" and "load_obj_from_file".
    void load_cgp_data_structure(std::string filepath_to_ws_region_image);
    /*! Loads arcs from an given HDF5 file and stores them into a specified vector.
     * \param filepath_to_ws_region_image Filepath to the HDF5 file
     * \param arcs Vector the loaded arcs will be stored in
     */
    //Loads arcs from an existing *.objects.h5 file
    void load_arcs(std::string filepath_to_ws_region_image, std::vector< std::vector<point> > & arcs);
    //TODO: Destructor?
};

seg::seg(bool twoPixelBoundary)
{
    two_pixel_boundary = twoPixelBoundary;
    numberOfBinningGroups = 1000;
}



void seg::load_ws_image_from_file(std::string filepath_to_ws_region_image)
{    
    std::cout<<"Importing results from file: "<<std::endl;
    std::cout<<filepath_to_ws_region_image<<std::endl;
    
    //Open existing file
    hid_t file_id0 = H5Fopen(filepath_to_ws_region_image.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    
    //Open an existing dataset
    hid_t dataset_id = H5Dopen(file_id0, "/ws_image", H5P_DEFAULT);
    
    //Get filespace handle
    hid_t filespace = H5Dget_space(dataset_id);
    
    //Get dimensions
    hsize_t dims[3];    
    H5Sget_simple_extent_dims(filespace, dims, NULL);

    dim_y=dims[0];
    dim_x=dims[1];

    ws_region_image.resize(dims[1], dims[0]);

    //Loop over rows in ws_image
    for(int i = 0; i < dims[0]; i++)
    {
        //Datastructure to read in
        unsigned int * row_values = new unsigned int[dims[1]];
        
        //Dataspace for one row    
        hsize_t row[3];
        row[0] = 1;
        row[1] = dims[1];
        row[2] = 1;
        hid_t mdataspace_id = H5Screate_simple(3, row, NULL);

        //Select file hyperslap
        hsize_t start_h[3]; //Start of hyperslab
        hsize_t count[3];   //Block count

        count[0] = 1;
        count[1] = dims[1];
        count[2] = 1;

        start_h[0] = i;
        start_h[1] = 0;
        start_h[2] = 0;
        
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
        H5Dread(dataset_id, H5T_NATIVE_UINT, mdataspace_id, filespace, H5P_DEFAULT, row_values);

        for(int x = 0; x < dims[1]; x++)
        {
            unsigned int value = row_values[x];
            ws_region_image(x, i) = value;
        }
    
        //Close the memory space
        H5Sclose(mdataspace_id);

        delete row_values;
    }
    
    //Close the file space
    H5Sclose(filespace);

    //Close the data set
    H5Dclose(dataset_id);

    //Close the file
    H5Fclose(file_id0);

    std::cout<<"...done"<<std::endl;
}

void seg::load_obj_from_file(std::string filepath_to_ws_region_image)
{   
    std::cout<<"Importing results from file: "<<std::endl;
    std::cout<<filepath_to_ws_region_image<<std::endl;

    //Open an existing file
    hid_t file_id = H5Fopen(filepath_to_ws_region_image.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    {
        
        //Open an existing data set
            hid_t dataset_id = H5Dopen(file_id, "/neighborhood-2", H5P_DEFAULT);
    
            //Get file space handle
            hid_t filespace = H5Dget_space(dataset_id);
    
            //Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);
    
            size_t size[] = {dims[1],dims[0]};
            two_boundings.resize(size,size+2);
    
            //Loop over columns in two_boundings
            for(int i = 0; i < dims[1]; i++)
            {
                //Data structure to read in
                unsigned int * column_values = new unsigned int[dims[1]];
    
                //Data space for one column
                hsize_t column[2];
                column[0] = dims[0];
                column[1] = 1;  
                hid_t mdataspace_id = H5Screate_simple(2, column, NULL);
    
                //Select file hyperslab 
                hsize_t start_h[2];  //Start of hyperslab
                hsize_t count[2];    //Block count
                
                count[0]  = dims[0]; 
                count[1]  = 1;
    
                start_h[0]  = 0; 
                start_h[1]  = i;
    
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
                H5Dread(dataset_id, H5T_NATIVE_UINT, mdataspace_id, filespace, H5P_DEFAULT, column_values);
    
                for (int x = 0; x < dims[0]; x++)
                {            
                    unsigned int value=(unsigned int)column_values[x];
                    two_boundings(i,x)=value;
                }
                //Close the memory space
                H5Sclose(mdataspace_id);
    
                delete column_values;
            }

            //Close the file space
            H5Sclose(filespace);
    
            //Close the data set
            H5Dclose(dataset_id);
        }

        {
            //Open an existing dataset
            hid_t dataset_id = H5Dopen(file_id, "/neighborhood-1", H5P_DEFAULT);

            //Get file space handle
            hid_t filespace = H5Dget_space(dataset_id);
    
            //Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);
    
            size_t size[] = {dims[1],dims[0]};
            one_boundings.resize(size,size+2);
    
            //Loop over columns in one_boundings
            for(int i = 0; i < dims[1]; i++)
            {
                //Data structure to read in
                unsigned int * column_values = new unsigned int[dims[1]];

                //Data space for one column
                hsize_t column[2];
                column[0] = dims[0];
                column[1] = 1;  
                hid_t mdataspace_id = H5Screate_simple(2, column, NULL);
    
                //Select file hyperslab 
                hsize_t start_h[2];  //Start of hyperslab
                hsize_t count[2];    //Block count
                
                count[0]  = dims[0]; 
                count[1]  = 1;
    
                start_h[0]  = 0; 
                start_h[1]  = i;
    
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
                H5Dread(dataset_id, H5T_NATIVE_UINT, mdataspace_id, filespace, H5P_DEFAULT, column_values);

                for (int x = 0; x < dims[0]; x++)
                {            
                    unsigned int value=(unsigned int)column_values[x];
                    one_boundings(i,x)=value;
                }
                //Close the memory space
                H5Sclose(mdataspace_id);
    
                delete column_values;
            }
    
            //Close the file space
            H5Sclose(filespace);
    
            //Close the data set
            H5Dclose(dataset_id);
        }

        unsigned int max_labels[4];

        {
            //Open an existing data set
            hid_t dataset_id = H5Dopen(file_id, "/max-labels", H5P_DEFAULT);
    
            H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, max_labels);
    
            //Close the data set
            H5Dclose(dataset_id);
        }    

        //std::cout<<max_labels[0]<<" "<<max_labels[1]<<" "<<max_labels[2]<<" "<<max_labels[3]<<std::endl;        

        //Open an existing group
        hid_t one_sets_group_id = H5Gopen(file_id, "/1-sets", H5P_DEFAULT);
    
        for (size_t i = 1; i < max_labels[1]+1; i++)
        {
            std::stringstream datasetNameStream;
            datasetNameStream << i << "-1";
            std::stringstream binningGroupName;
            binningGroupName << "bin-" << i % numberOfBinningGroups;
            hid_t h = H5Gopen(one_sets_group_id, binningGroupName.str().c_str(), H5P_DEFAULT);

            //Open an existing data set
            hid_t dataset_id = H5Dopen(h, datasetNameStream.str().c_str(), H5P_DEFAULT);
    
            //Data structure to read in
            int junction_pos[3];
    
            H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, junction_pos);

            point junction_p;
            junction_p.x=junction_pos[1]/2;
            junction_p.y=junction_pos[0]/2;
    
            junctions.push_back(junction_p);
    
            //Close the data set
            H5Dclose(dataset_id);

            H5Gclose(h);
        }
    
        H5Gclose(one_sets_group_id);

        long unsigned int * arc_part_counters= new long unsigned int [max_labels[2]];

        {
            //Open an existing data set
            hid_t dataset_id = H5Dopen(file_id, "/parts-counters-2", H5P_DEFAULT);

            H5Dread(dataset_id, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, arc_part_counters);

            //Close the data set
            H5Dclose(dataset_id);
        }

        //Open an existing group
        hid_t two_sets_group_id = H5Gopen(file_id, "/2-sets", H5P_DEFAULT);
    
        arcs.resize(max_labels[2]);

        for (size_t i = 1; i < max_labels[2]+1; i++)
        {
            std::vector<point> this_arc;

            std::stringstream binningGroupName;
            binningGroupName << "bin-" << i % numberOfBinningGroups;
            hid_t h = H5Gopen(two_sets_group_id, binningGroupName.str().c_str(), H5P_DEFAULT);
    
            for (size_t part=1; part<arc_part_counters[i-1]+1; part++)
            {
                std::stringstream datasetNameStream;
            //Check for -2 and -3 ... and combine arc pixels    
                datasetNameStream << i <<"-"<<part;
            
                //Open an existing data set
                hid_t dataset_id = H5Dopen(h, datasetNameStream.str().c_str(), H5P_DEFAULT);

                //Get file space handle
                hid_t filespace = H5Dget_space(dataset_id);

                //Get dimensions
                hsize_t dims[2];
                H5Sget_simple_extent_dims(filespace, dims, NULL);

                //Datastructure to read in
                int * data_out_x = new int[dims[1]];
                int * data_out_y = new int[dims[1]];

                //Data space for one coordinate
                hsize_t row[2];
                row[0] = 1;
                row[1] = dims[1];  
                hid_t mdataspace_id = H5Screate_simple(2, row, NULL);

                // Select file hyperslab 
                hsize_t start_h[2];// Start of hyperslab
                hsize_t count[2];// Block count
            
                count[0]  = 1; 
                count[1]  = dims[1];

                start_h[0]  = 0; 
                start_h[1]  = 0;

                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
                H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_y);

                start_h[0]  = 1;

                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
                H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_x);

            for (int j = 0; j < dims[1]; j++)
            {
                point push;
                push.x = data_out_x[j];
                push.y = data_out_y[j];
                this_arc.push_back(push);
            }

            // Close the memory space
            H5Sclose(mdataspace_id);

            // Close the file space
            H5Sclose(filespace);

            // Close the data set
            H5Dclose(dataset_id);

            delete data_out_x;
            delete data_out_y;

        }

        for(int k=0;k<(int)this_arc.size();k++)
        {
            if(two_pixel_boundary==true)
            {
                point p1;
                point p2;

                //y is 0.5
                if(this_arc[k].x% 2 ==0 && this_arc[k].y% 2 !=0 )
                {
                    p1.x=(this_arc[k].x)/2;
                    p1.y=(this_arc[k].y)/2;

                    p2.x=(this_arc[k].x)/2;
                    p2.y=((this_arc[k].y)/2)+1;
                }

                //x is 0.5
                if(this_arc[k].x% 2 !=0 && this_arc[k].y% 2 ==0 )
                {
                    p1.x=(this_arc[k].x)/2;
                    p1.y=(this_arc[k].y)/2;

                    p2.x=(this_arc[k].x/2)+1;
                    p2.y=(this_arc[k].y)/2;
                }

                arcs[i-1].push_back(p1);
                arcs[i-1].push_back(p2);
            }

            if(two_pixel_boundary==false)
            {
                point p;

                p.x=(this_arc[k].x)/2;
                p.y=(this_arc[k].y)/2;

                arcs[i-1].push_back(p);
            }
        }

        H5Gclose(h);

        }

        H5Gclose(two_sets_group_id);

        // Close the file
        H5Fclose(file_id);

        std::cout<<"...done"<<std::endl;
}

void seg::load_cgp_data_structure(std::string filepath_to_ws_region_image)
{
    std::string arg_1 = filepath_to_ws_region_image;
    std::string arg_2 = filepath_to_ws_region_image;
    arg_2.resize(arg_2.size()-3);
    arg_2.append(".objects.h5");
    
    //Load from .h5 file
    load_ws_image_from_file(arg_1);
    //Load from .objects.h5 file
    load_obj_from_file(arg_2);
}

void seg::load_arcs(std::string filepath_to_ws_region_image, std::vector< std::vector<point> > & arcs)
{
    const size_t numberOfBinningGroups = 1000;

    filepath_to_ws_region_image.resize(filepath_to_ws_region_image.size()-4);
    filepath_to_ws_region_image.append(".objects.h5");

    std::cout<<"Importing arcs from file: "<<std::endl;
    std::cout<<filepath_to_ws_region_image<<std::endl;

    // Open an existing file
    hid_t file_id = H5Fopen(filepath_to_ws_region_image.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
 
    unsigned int max_labels[4];

    {
        // Open an existing dataset
        hid_t dataset_id = H5Dopen(file_id, "/max-labels", H5P_DEFAULT);

        H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, max_labels);

        // Close the dataset
        H5Dclose(dataset_id);
    }

    long unsigned int * arc_part_counters= new long unsigned int [max_labels[2]];

    {
        // Open an existing dataset
        hid_t dataset_id = H5Dopen(file_id, "/parts-counters-2", H5P_DEFAULT);

        H5Dread(dataset_id, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, arc_part_counters);

        // Close the dataset
        H5Dclose(dataset_id);
    }

    // Open an existing group
    hid_t two_sets_group_id = H5Gopen(file_id, "/2-sets", H5P_DEFAULT);

    arcs.resize(max_labels[2]);

    for (size_t i=1; i<max_labels[2]+1; i++)
    {
        std::vector<point> this_arc;

        std::stringstream binningGroupName;
        binningGroupName << "bin-" << i % numberOfBinningGroups;
        hid_t h = H5Gopen(two_sets_group_id, binningGroupName.str().c_str(), H5P_DEFAULT);

        for (size_t part=1; part<arc_part_counters[i-1]+1; part++)
        {
            std::stringstream datasetNameStream;
            datasetNameStream << i <<"-"<<part;//check for -2 and -3 ... and combine arc pixels
            
            // Open an existing dataset
            hid_t dataset_id = H5Dopen(h, datasetNameStream.str().c_str(), H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            // Datastructure to read in
            int * data_out_x = new int[dims[1]];
            int * data_out_y = new int[dims[1]];

            //dataspace for one coordinate
            hsize_t row[2];
            row[0] = 1;
            row[1] = dims[1];  
            hid_t mdataspace_id = H5Screate_simple(2, row, NULL);

            // select file hyperslab 
            hsize_t start_h[2];// start of hyperslab
            hsize_t count[2];// block count
            
            count[0]  = 1; 
            count[1]  = dims[1];

            start_h[0]  = 0; 
            start_h[1]  = 0;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_y);

            start_h[0]  = 1;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_x);

            for (int j=0;j<dims[1];j++)
            {
                point push;
                push.x=data_out_x[j];
                push.y=data_out_y[j];
                this_arc.push_back(push);
            }

            // Close the memoryspace
            H5Sclose(mdataspace_id);

            // Close the filespace
            H5Sclose(filespace);

            // Close the dataset
            H5Dclose(dataset_id);

            delete data_out_x;
            delete data_out_y;

        }

        for(int k=0;k<(int)this_arc.size();k++)
        {
            point p;

            p.x=(this_arc[k].x)/2;
            p.y=(this_arc[k].y)/2;

            arcs[i-1].push_back(p);
        }

        H5Gclose(h);

    }

    H5Gclose(two_sets_group_id);

    // Close the file
    H5Fclose(file_id);

    std::cout<<"...done"<<std::endl;
}

