/*! \file param.h
 * \brief Class representing the Network-Parameters structure.
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

//Those includes become necessary only if this class is included seperately
#include "boundary_data_structure.h"

/*Class representing the Network-Parameters structure (cf. code manual for IceGrain)*/

class param
{
    public:
    //Attributes **************************************************************
    std::vector<int>                            grain_areas;              /*!< Indices of grain areas.*/
    std::vector<float>                          grain_roundness;          /*!< Grain roundness of each grain area.*/
    std::vector<float>                          grain_box_flattening;     /*!< Grain box flattening of each grain area.*/
    std::vector<float>                          grain_box_width;          /*!< Grain box width of each grain area.*/
    std::vector<float>                          grain_box_height;         /*!< Grain box height of each grain area.*/
    std::vector< std::vector<float> >           ellipse_params;           /*!< Ellipse parameters for each area.*/
    std::vector<float>                          ellipse_long_axis;        /*!< Ellipse long axis for each area.*/
    std::vector<float>                          ellipse_flattening;       /*!< Ellipse flatteing for each area.*/
    std::vector<float>                          ellipse_long_axis_angle;  /*!< Ellipse long axis angle for each area.*/
    std::vector<float>                          grain_area_width;         /*!< Area width for each grain area.*/
    std::vector<float>                          grain_area_height;        /*!< Area height for each grain area.*/
    std::vector<float>                          grain_area_number;        /*!< Number for each grain area.*/
    std::vector<float>                          grain_arc_number;         /*!< Number for each grain arc.*/
    std::vector<float>                          grain_neighbors;          /*!< Neighbors for each grain.*/
    std::vector<float>                          grain_longest_arc_length; /*!< Longest arc length for each grain.*/
    std::vector< std::vector<int> >             grain_boundary_index;     /*!< Grain indices for each grain boundary.*/ 
    std::vector< std::vector<float> >           grain_boundary_phis;      /*!< Grain orientation for each grain boundary.*/
    std::vector< std::vector<float> >           grain_boundary_curvs;     /*!< Grain curvs for each grain boundary.*/
    std::vector<int>                            grain_junctions;          /*!< Grain junctions.*/
    std::vector< std::vector<int> >             grain_junction_index;     /*!< Grain indices for each grain junction.*/
    std::vector<float>                          turning_point;            /*!< Turning point.*/
    std::vector<bool>                           grain_junction;           /*!< Indices of grain junctions.*/
    std::vector<float>                          grain_perimeter_ratio;    /*!< Grain perimeter ratio for each grain.*/
    std::vector< std::vector<unsigned int> >    grain_area_boundaries;    /*!< Grain area boundaries for each grain.*/
    
    //Methods *****************************************************************
    /*! Load data structures from a given HDF5 file.
     * \param filepath_parameters Filepath to the HDF5 file containing the Network-Parameter structure
     */
    //Taken from 'extracted_parameters.h'
    void load_extracted_parameters(std::string filepath_parameters);
    /*! Load selected data structures from a given HDF5 file.
     * \param filepath_parameters Filepath to the HDF5 file containing the Network-Parameter structure
     */
    void load_selected_extracted_parameters(std::string filepath_parameters);
	/*! Save data structures into a HDF5 file.
     * \param filepath_parameters Filepath to the destination HDF5 file
     */
    void save_extracted_parameters(std::string filepath_parameters);
};

void param::load_extracted_parameters(std::string filepath_parameters)
{
    std::cout<<"Importing results from file: "<<std::endl;
    std::cout<<filepath_parameters<<std::endl;
    
    // Open an existing file
    hid_t file_save = H5Fopen(filepath_parameters.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    {
        //dataset for grain areas
        hid_t dataset_id = H5Dopen(file_save, "grain_areas", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        int * areas = new int[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, areas);

        for (int i=0; i<dims[0]; i++)
        {
            grain_areas.push_back(areas[i]);
        }

        delete areas;

        H5Dclose(dataset_id);
    }

    {
        //dataset for grain area width
        hid_t dataset_id = H5Dopen(file_save, "grain_area_width", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * width = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, width);

        grain_area_width.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_area_width[grain_areas[i]]=width[i];
        }

        delete width;

        H5Dclose(dataset_id);
    }

    {
        //dataset for grain area height
        hid_t dataset_id = H5Dopen(file_save, "grain_area_height", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * height = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, height);

        grain_area_height.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_area_height[grain_areas[i]]=height[i];
        }

        delete height;

        H5Dclose(dataset_id);
    }

    {
        //dataset for grain roundness
        hid_t dataset_id = H5Dopen(file_save, "grain_roundness", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * roundness = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, roundness);

        grain_roundness.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_roundness[grain_areas[i]]=roundness[i];
        }

        delete roundness;

        H5Dclose(dataset_id);
    }

    {
        //dataset for grain box flattening
        hid_t dataset_id = H5Dopen(file_save, "grain_box_flattening", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * flattening = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flattening);

        grain_box_flattening.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_box_flattening[grain_areas[i]]=flattening[i];
        }

        delete flattening;

        H5Dclose(dataset_id);
    }

    {
        htri_t check = H5Lexists(file_save, "grain_box_width", H5P_DEFAULT);
        if(!check)
        {
            std::cout<<"Repeat parameter extraction. Old version of parameter file found!"<<std::endl;
            exit(-1);
        }
        else
        {
            //dataset for grain box width
            hid_t dataset_id = H5Dopen(file_save, "grain_box_width", H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            float * width = new float[dims[0]];

            H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, width);

            grain_box_width.resize(grain_areas.back()+1);

            for (int i=0; i<dims[0]; i++)
            {
                grain_box_width[grain_areas[i]]=width[i];
            }

            delete width;

            H5Dclose(dataset_id);
        }
    }

    {
        //dataset for grain box height
        hid_t dataset_id = H5Dopen(file_save, "grain_box_height", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * height = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, height);

        grain_box_height.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_box_height[grain_areas[i]]=height[i];
        }

        delete height;

        H5Dclose(dataset_id);
    }

    {
        // Open an existing dataset
        hid_t dataset_id = H5Dopen(file_save, "/grain_ellipse_parameters", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        // Datastructure to read in
        float ** params = new float * [6];
        for(int i=0; i<6; i++)
        {
            params[i] = new float[dims[0]];
        }

        //dataspace for one param
        hsize_t column[2];
        column[0] = dims[0];
        column[1] = 1;  
        hid_t mdataspace_id = H5Screate_simple(2, column, NULL);

        // select file hyperslab 
        hsize_t start_h[2];// start of hyperslab
        hsize_t count[2];// block count
        
        count[0]  = dims[0]; 
        count[1]  = 1;

        start_h[0]  = 0; 

        for(int i=0; i<6; i++)
        {
            start_h[1]  = i;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_FLOAT, mdataspace_id, filespace, H5P_DEFAULT, params[i]);
        }

        ellipse_params.resize(grain_areas.back()+1);

        for(int i=0; i<6; i++)
        {
            for (int j=0;j<dims[0];j++)
            {
                ellipse_params[grain_areas[j]].push_back(params[i][j]);
            }  

            delete params[i];
        }

        // Close the memoryspace
        H5Sclose(mdataspace_id);

        // Close the filespace
        H5Sclose(filespace);

        delete params;

        // Close the dataset
        H5Dclose(dataset_id);
    }

    {
        //dataset for ellipse long axis
        hid_t dataset_id = H5Dopen(file_save, "grain_ellipse_long_axis", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * long_axis = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, long_axis);

        ellipse_long_axis.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            ellipse_long_axis[grain_areas[i]]=long_axis[i];
        }

        delete long_axis;

        H5Dclose(dataset_id);
    }

    {
        //dataset for ellipse flattening
        hid_t dataset_id = H5Dopen(file_save, "grain_ellipse_flattening", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * flattening = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flattening);

        ellipse_flattening.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            ellipse_flattening[grain_areas[i]]=flattening[i];
        }

        delete flattening;

        H5Dclose(dataset_id);
    }

    {
        //dataset for ellipse_long_axis_angle
        hid_t dataset_id = H5Dopen(file_save, "grain_ellipse_long_axis_angle", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * angle = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, angle);

        ellipse_long_axis_angle.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            ellipse_long_axis_angle[grain_areas[i]]=angle[i];
        }

        delete angle;

        H5Dclose(dataset_id);
    }

    {
        //dataset for grain arc number
        hid_t dataset_id = H5Dopen(file_save, "grain_arc_number", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        int * number = new int[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, number);

        grain_arc_number.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_arc_number[grain_areas[i]]=(float)number[i];
        }

        delete number;

        H5Dclose(dataset_id);
    }

    {
        htri_t check = H5Lexists(file_save, "grain_perimeter_ratio", H5P_DEFAULT);
        if(!check)
        {
            std::cout<<"Old version of parameter file: no grain perimeter ratio dataset found!"<<std::endl;
        }
        else
        {
            //dataset for grain perimeter ratio
            hid_t dataset_id = H5Dopen(file_save, "grain_perimeter_ratio", H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            float * number = new float[dims[0]];

            H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, number);

            grain_perimeter_ratio.resize(grain_areas.back()+1);

            for (int i=0; i<dims[0]; i++)
            {
                grain_perimeter_ratio[grain_areas[i]]=(float)number[i];
            }

            delete number;

            H5Dclose(dataset_id);
        }   
    }

    {
        //dataset for grain neighbors
        hid_t dataset_id = H5Dopen(file_save, "grain_neighbors", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        int * neighbors = new int[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, neighbors);

        grain_neighbors.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_neighbors[grain_areas[i]]=(float)neighbors[i];
        }

        delete neighbors;

        H5Dclose(dataset_id);
    }

    {
        //dataset for grain longest arc length
        hid_t dataset_id = H5Dopen(file_save, "grain_longest_arc_length", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * longest = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, longest);

        grain_longest_arc_length.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_longest_arc_length[grain_areas[i]]=longest[i];
        }

        delete longest;

        H5Dclose(dataset_id);
    }

    {
        //dataset for number turning points
        hid_t dataset_id = H5Dopen(file_save, "number_turning_points", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * number = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, number);

        for (int i=0; i<dims[0]; i++)
        {
            turning_point.push_back(number[i]);
        }

        delete number;

        H5Dclose(dataset_id);
    }

    grain_boundary_index.resize(turning_point.size());
    grain_boundary_phis.resize(turning_point.size());
    grain_boundary_curvs.resize(turning_point.size());

    {
        //Open an existing group
        hid_t group_id = H5Gopen(file_save, "/grain_boundary_index", H5P_DEFAULT);

        for(int j=1; j<grain_boundary_index.size()+1; j++)
        {
            char buffer [6];
            sprintf (buffer, "%d", j);

            //dataset for arcs
            hid_t dataset_id = H5Dopen(group_id, buffer, H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            //array to store arcs
            int * arcs=new int[dims[0]];
            
            H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arcs);

            for(int x=0;x<dims[0];x++)
            {
                grain_boundary_index[j-1].push_back(arcs[x]);
            }

            delete arcs;

            H5Dclose(dataset_id);
        }

        //Close the group
        H5Gclose(group_id);
    }

    {
        //Open an existing group
        hid_t group_id = H5Gopen(file_save, "/grain_boundary_phis", H5P_DEFAULT);

        for(int j=1; j<grain_boundary_phis.size()+1; j++)
        {
            char buffer [6];
            sprintf (buffer, "%d", j);

            //dataset for phis
            hid_t dataset_id = H5Dopen(group_id, buffer, H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            //array to store phis
            float * phis=new float[dims[0]];
            
            H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, phis);

            for(int x=0;x<dims[0];x++)
            {
                grain_boundary_phis[j-1].push_back(phis[x]);
            }

            delete phis;

            H5Dclose(dataset_id);
        }

        //Close the group
        H5Gclose(group_id);
    }

    {
        //Open an existing group
        hid_t group_id = H5Gopen(file_save, "/grain_boundary_curvs", H5P_DEFAULT);

        for(int j=1; j<grain_boundary_curvs.size()+1; j++)
        {
            char buffer [6];
            sprintf (buffer, "%d", j);

            //dataset for curvs
            hid_t dataset_id = H5Dopen(group_id, buffer, H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            //array to store curvs
            float * curvs=new float[dims[0]];
            
            H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, curvs);

            for(int x=0;x<dims[0];x++)
            {
                grain_boundary_curvs[j-1].push_back(curvs[x]);
            }

            delete curvs;

            H5Dclose(dataset_id);
        }

        //Close the group
        H5Gclose(group_id);
    }

    {
        htri_t check = H5Lexists(file_save, "grain_junction_index", H5P_DEFAULT);

        if(check)
        {
            // Open an existing dataset
            hid_t dataset_id = H5Dopen(file_save, "/grain_junction_index", H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            // Datastructure to read in
            int * data_out_1 = new int[dims[0]];
            int * data_out_2 = new int[dims[0]];
            int * data_out_3 = new int[dims[0]];
            int * data_out_4 = new int[dims[0]];

            //dataspace for one coordinate
            hsize_t column[2];
            column[0] = dims[0];
            column[1] = 1;  
            hid_t mdataspace_id = H5Screate_simple(2, column, NULL);

            // select file hyperslab 
            hsize_t start_h[2];// start of hyperslab
            hsize_t count[2];// block count
            
            count[0]  = dims[0]; 
            count[1]  = 1;

            start_h[0]  = 0; 
            start_h[1]  = 0;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_1);

            start_h[1]  = 1;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_2);

            start_h[1]  = 2;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_3);

            start_h[1]  = 3;

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_4);

            if(grain_junctions.size()!=dims[0])
            {
                std::cout<<"Error: Number of grain boundary junctions in network/parameter structure is different!"<<std::endl;
                exit(-1);
            }

            for (int j=0;j<dims[0];j++)
            {
                if (grain_junctions[j]+1>grain_junction_index.size()) grain_junction_index.resize(grain_junctions[j]+1);
                if (grain_junctions[j]+1>grain_junction.size()) grain_junction.resize(grain_junctions[j]+1,false);

                if (data_out_1[j]>0) grain_junction_index[grain_junctions[j]].push_back(data_out_1[j]-1);
                if (data_out_2[j]>0) grain_junction_index[grain_junctions[j]].push_back(data_out_2[j]-1);
                if (data_out_3[j]>0) grain_junction_index[grain_junctions[j]].push_back(data_out_3[j]-1);
                if (data_out_4[j]>0) grain_junction_index[grain_junctions[j]].push_back(data_out_4[j]-1);

                grain_junction[grain_junctions[j]]=true;
            }

            // Close the memoryspace
            H5Sclose(mdataspace_id);

            // Close the filespace
            H5Sclose(filespace);

            delete data_out_1;
            delete data_out_2;
            delete data_out_3;
            delete data_out_4;

            // Close the dataset
            H5Dclose(dataset_id);
        }
    }

    {
        if(!H5Lexists(file_save, "grain_area_boundaries", H5P_DEFAULT))
        {
            std::cout<<"Old version of parameter file: no grain area boundaries dataset found!"<<std::endl;
        }
        else
        {
            //Open the existing group for grain area boundaries
            hid_t group_id = H5Gopen(file_save, "/grain_area_boundaries", H5P_DEFAULT);

            for(int i = 0; i < grain_area_boundaries.size(); i++)
            {            
                char buffer[5];
                sprintf(buffer, "%d", i);
                htri_t check = H5Lexists(group_id, buffer, H5P_DEFAULT);

                if(check)
                {
                    //Data set for grain area boundaries
                    hid_t dataset_id = H5Dopen(group_id, buffer, H5P_DEFAULT);

                    //Get filespace handle
                    hid_t filespace = H5Dget_space(dataset_id);

                    //Get dimensions
                    hsize_t dims[2];
                    H5Sget_simple_extent_dims(filespace, dims, NULL);

                    //Array to store grain area boundaries
                    unsigned int * vals = new unsigned int[dims[0]];

                    H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);

                    for(int j = 0; j < dims[0]; j++)
                    {
                        grain_area_boundaries[i].push_back(vals[j]);
                    }

                    delete vals;

                    H5Dclose(dataset_id);
                }
            }

            //Close the group
            H5Gclose(group_id);
        }
    }

    //Close the file
    H5Fclose(file_save);

    std::cout<<"File closed"<<std::endl;
}

void param::load_selected_extracted_parameters(std::string filepath_parameters)
{
    std::cout<<"Importing results from file: "<<std::endl;
    std::cout<<filepath_parameters<<std::endl;
    
    // Open an existing file
    hid_t file_save = H5Fopen(filepath_parameters.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    {
        //dataset for grain areas
        hid_t dataset_id = H5Dopen(file_save, "grain_areas", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        int * areas = new int[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, areas);

        for (int i=0; i<dims[0]; i++)
        {
            grain_areas.push_back(areas[i]);
        }

        delete areas;

        H5Dclose(dataset_id);
    }

    {
        //dataset for grain roundness
        hid_t dataset_id = H5Dopen(file_save, "grain_roundness", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * roundness = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, roundness);

        grain_roundness.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_roundness[grain_areas[i]]=roundness[i];
        }

        delete roundness;

        H5Dclose(dataset_id);
    }

    {
        //dataset for grain box flattening
        hid_t dataset_id = H5Dopen(file_save, "grain_box_flattening", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * flattening = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flattening);

        grain_box_flattening.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_box_flattening[grain_areas[i]]=flattening[i];
        }

        delete flattening;

        H5Dclose(dataset_id);
    }

    {
        //dataset for ellipse_long_axis_angle
        hid_t dataset_id = H5Dopen(file_save, "grain_ellipse_long_axis_angle", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * angle = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, angle);

        ellipse_long_axis_angle.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            ellipse_long_axis_angle[grain_areas[i]]=angle[i];
        }

        delete angle;

        H5Dclose(dataset_id);
    }

    {
        htri_t check = H5Lexists(file_save, "grain_perimeter_ratio", H5P_DEFAULT);
        if(!check)
        {
            std::cout<<"Old version of parameter file: no grain perimeter ratio dataset found!"<<std::endl;
        }
        else
        {
            //dataset for grain perimeter ratio
            hid_t dataset_id = H5Dopen(file_save, "grain_perimeter_ratio", H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            float * number = new float[dims[0]];

            H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, number);

            grain_perimeter_ratio.resize(grain_areas.back()+1);

            for (int i=0; i<dims[0]; i++)
            {
                grain_perimeter_ratio[grain_areas[i]]=(float)number[i];
            }

            delete number;

            H5Dclose(dataset_id);
        }   
    }

    {
        //dataset for grain area width
        hid_t dataset_id = H5Dopen(file_save, "grain_area_width", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * width = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, width);

        grain_area_width.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_area_width[grain_areas[i]]=width[i];
        }

        delete width;

        H5Dclose(dataset_id);
    }

    {
        //dataset for grain area height
        hid_t dataset_id = H5Dopen(file_save, "grain_area_height", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        float * height = new float[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, height);

        grain_area_height.resize(grain_areas.back()+1);

        for (int i=0; i<dims[0]; i++)
        {
            grain_area_height[grain_areas[i]]=height[i];
        }

        delete height;

        H5Dclose(dataset_id);
    }

    //get dimensions for grain_boundary_curvs
    {
        //dataset for number turning points
        hid_t dataset_id = H5Dopen(file_save, "number_turning_points", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        grain_boundary_curvs.resize(dims[0]);

        H5Dclose(dataset_id);
    }

    {
        //Open an existing group
        hid_t group_id = H5Gopen(file_save, "/grain_boundary_curvs", H5P_DEFAULT);

        for(int j=1; j<grain_boundary_curvs.size()+1; j++)
        {
            char buffer [6];
            sprintf (buffer, "%d", j);

            //dataset for curvs
            hid_t dataset_id = H5Dopen(group_id, buffer, H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            //array to store curvs
            float * curvs=new float[dims[0]];
            
            H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, curvs);

            for(int x=0;x<dims[0];x++)
            {
                grain_boundary_curvs[j-1].push_back(curvs[x]);
            }

            delete curvs;

            H5Dclose(dataset_id);
        }

        //Close the group
        H5Gclose(group_id);
    }

    //Close the file
    H5Fclose(file_save);

    std::cout<<"File closed"<<std::endl;
}

void param::save_extracted_parameters(std::string filepath_parameters)
{
    //EXPORT RESULTS TO A HDF5 file
    filepath_parameters.append(".h5");

    std::cout<<"Exporting parameters to file: "<<std::endl;
    std::cout<<filepath_parameters<<std::endl;
    
    //creating the file
    hid_t file_save=H5Fcreate(filepath_parameters.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

    {
        //dataspace for grain areas
        hsize_t dims[2];
        dims[0] = grain_areas.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        int * areas = new int[grain_areas.size()];
        for (int i=0; i<grain_areas.size(); i++)
        {
            areas[i]=grain_areas[i];
        }

        //dataset for grain areas
        hid_t dataset_id = H5Dcreate1(file_save, "grain_areas", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, areas);

        delete areas;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for grain roundness
        hsize_t dims[2];
        dims[0] = grain_roundness.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * roundness = new float[grain_roundness.size()];
        for (int i=0; i<grain_roundness.size(); i++)
        {
            roundness[i]=grain_roundness[i];
        }

        //dataset for grain roundness
        hid_t dataset_id = H5Dcreate1(file_save, "grain_roundness", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, roundness);

        delete roundness;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for grain box flattening
        hsize_t dims[2];
        dims[0] = grain_box_flattening.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * flattening = new float[grain_box_flattening.size()];
        for (int i=0; i<grain_box_flattening.size(); i++)
        {
            flattening[i]=grain_box_flattening[i];
        }

        //dataset for grain box flattening
        hid_t dataset_id = H5Dcreate1(file_save, "grain_box_flattening", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flattening);

        delete flattening;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for grain box width
        hsize_t dims[2];
        dims[0] = grain_box_width.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * width = new float[grain_box_width.size()];
        for (int i=0; i<grain_box_width.size(); i++)
        {
            width[i]=grain_box_width[i];
        }

        //dataset for grain box width
        hid_t dataset_id = H5Dcreate1(file_save, "grain_box_width", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, width);

        delete width;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for grain box height
        hsize_t dims[2];
        dims[0] = grain_box_height.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * height = new float[grain_box_height.size()];
        for (int i=0; i<grain_box_height.size(); i++)
        {
            height[i]=grain_box_height[i];
        }

        //dataset for grain box height
        hid_t dataset_id = H5Dcreate1(file_save, "grain_box_height", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, height);

        delete height;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for ellipse params
        hsize_t dims[2];
        dims[0] = ellipse_params.size(); 
        dims[1] = 6; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        //dataset for ellipse params
        hid_t dataset_id = H5Dcreate1(file_save, "grain_ellipse_parameters", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        float ** params = new float * [6];
        for(int i=0; i<6; i++)
        {
            params[i] = new float[ellipse_params.size()];
            for (int j=0; j<ellipse_params.size(); j++) 
            {
                params[i][j]=ellipse_params[j][i];
            }
        }

        //dataspace for one param
        hsize_t column[2];
        column[0] = ellipse_params.size();
        column[1] = 1;  
        hid_t mdataspace_id = H5Screate_simple(2, column, NULL);

        // select file hyperslab 
        hsize_t start_h[2];// start of hyperslab
        hsize_t count[2];// block count
        
        count[0]  = ellipse_params.size(); 
        count[1]  = 1;

        start_h[0]  = 0; 

        for(int i=0; i<6; i++)
        {
            start_h[1]  = i;

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, mdataspace_id, dataspace_id, H5P_DEFAULT, params[i]);

            delete params[i];
        }

        // Close the memoryspace
        H5Sclose(mdataspace_id);

        // Close the filespace
        H5Sclose(dataspace_id);

        delete params;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for ellipse long axis
        hsize_t dims[2];
        dims[0] = ellipse_long_axis.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * long_axis = new float[ellipse_long_axis.size()];
        for (int i=0; i<ellipse_long_axis.size(); i++)
        {
            long_axis[i]=ellipse_long_axis[i];
        }

        //dataset for ellipse long axis
        hid_t dataset_id = H5Dcreate1(file_save, "grain_ellipse_long_axis", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, long_axis);

        delete long_axis;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for ellipse flattening
        hsize_t dims[2];
        dims[0] = ellipse_flattening.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * flattening = new float[ellipse_flattening.size()];
        for (int i=0; i<ellipse_flattening.size(); i++)
        {
            flattening[i]=ellipse_flattening[i];
        }

        //dataset for ellipse flattening
        hid_t dataset_id = H5Dcreate1(file_save, "grain_ellipse_flattening", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flattening);

        delete flattening;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for ellipse long axis angle
        hsize_t dims[2];
        dims[0] = ellipse_long_axis_angle.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * angle = new float[ellipse_long_axis_angle.size()];
        for (int i=0; i<ellipse_long_axis_angle.size(); i++)
        {
            angle[i]=ellipse_long_axis_angle[i];
        }

        //dataset for ellipse long axis angle
        hid_t dataset_id = H5Dcreate1(file_save, "grain_ellipse_long_axis_angle", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, angle);

        delete angle;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for grain area width
        hsize_t dims[2];
        dims[0] = grain_area_width.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * width = new float[grain_area_width.size()];
        for (int i=0; i<grain_area_width.size(); i++)
        {
            width[i]=grain_area_width[i];
        }

        //dataset for grain area width
        hid_t dataset_id = H5Dcreate1(file_save, "grain_area_width", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, width);

        delete width;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for grain area height
        hsize_t dims[2];
        dims[0] = grain_area_height.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * height = new float[grain_area_height.size()];
        for (int i=0; i<grain_area_height.size(); i++)
        {
            height[i]=grain_area_height[i];
        }

        //dataset for grain are height
        hid_t dataset_id = H5Dcreate1(file_save, "grain_area_height", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, height);

        delete height;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for grain arc number
        hsize_t dims[2];
        dims[0] = grain_arc_number.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        int * number = new int[grain_arc_number.size()];
        for (int i=0; i<grain_arc_number.size(); i++)
        {
            number[i]=(int)grain_arc_number[i];
        }

        //dataset for grain arc number
        hid_t dataset_id = H5Dcreate1(file_save, "grain_arc_number", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, number);

        delete number;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for grain perimeter ratio
        hsize_t dims[2];
        dims[0] = grain_perimeter_ratio.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * number = new float[grain_perimeter_ratio.size()];
        for (int i=0; i<grain_perimeter_ratio.size(); i++)
        {
            number[i]= grain_perimeter_ratio[i];
        }

        //dataset for grain perimeter ratio
        hid_t dataset_id = H5Dcreate1(file_save, "grain_perimeter_ratio", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, number);

        delete number;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for grain neighbors
        hsize_t dims[2];
        dims[0] = grain_neighbors.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        int * neighbors = new int[grain_neighbors.size()];
        for (int i=0; i<grain_neighbors.size(); i++)
        {
            neighbors[i]=(int)grain_neighbors[i];
        }

        //dataset for grain neighbors
        hid_t dataset_id = H5Dcreate1(file_save, "grain_neighbors", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, neighbors);

        delete neighbors;

        H5Dclose(dataset_id);
    }

    {
        //dataspace for grain longest arc length
        hsize_t dims[2];
        dims[0] = grain_longest_arc_length.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * longest = new float[grain_longest_arc_length.size()];
        for (int i=0; i<grain_longest_arc_length.size(); i++)
        {
            longest[i]=grain_longest_arc_length[i];
        }

        //dataset for grain longest arc length
        hid_t dataset_id = H5Dcreate1(file_save, "grain_longest_arc_length", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, longest);

        delete longest;

        H5Dclose(dataset_id);
    }

    {
        //Create group grain boundary index
        hid_t group_id = H5Gcreate1(file_save, "/grain_boundary_index", 0);

        for(int j=1; j<grain_boundary_index.size()+1; j++)
        {
            //array to store arcs
            int * arcs=new int[grain_boundary_index[j-1].size()];

            for(int x=0;x<grain_boundary_index[j-1].size();x++)
            {
                arcs[x]=grain_boundary_index[j-1][x];
            }

            //dataspace for arcs
            hsize_t dims[2];
            dims[0] = grain_boundary_index[j-1].size();
            dims[1] = 1; 
            hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

            char buffer [6];
            sprintf (buffer, "%d", j);

            //dataset for arcs
            hid_t dataset_id = H5Dcreate1(group_id, buffer, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
            
            H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arcs);

            delete arcs;

            H5Dclose(dataset_id);
        }

        //Close the group
        H5Gclose(group_id);
    }

    {
        //Create group grain boundary phis
        hid_t group_id = H5Gcreate1(file_save, "/grain_boundary_phis", 0);

        for(int j=1; j<grain_boundary_phis.size()+1; j++)
        {
            //array to store phis
            float * phis=new float[grain_boundary_phis[j-1].size()];

            for(int x=0;x<grain_boundary_phis[j-1].size();x++)
            {
                phis[x]=grain_boundary_phis[j-1][x];
            }

            //dataspace for phis
            hsize_t dims[2];
            dims[0] = grain_boundary_phis[j-1].size();
            dims[1] = 1; 
            hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

            char buffer [6];
            sprintf (buffer, "%d", j);

            //dataset for phis
            hid_t dataset_id = H5Dcreate1(group_id, buffer, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
            
            H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, phis);

            delete phis;

            H5Dclose(dataset_id);
        }

        //Close the group
        H5Gclose(group_id);
    }

    {
        //Create group grain boundary curvs
        hid_t group_id = H5Gcreate1(file_save, "/grain_boundary_curvs", 0);

        for(int j=1; j<grain_boundary_curvs.size()+1; j++)
        {
            //array to store curvs
            float * curvs=new float[grain_boundary_curvs[j-1].size()];

            for(int x=0;x<grain_boundary_curvs[j-1].size();x++)
            {
                curvs[x]=grain_boundary_curvs[j-1][x];
            }

            //dataspace for curvs
            hsize_t dims[2];
            dims[0] = grain_boundary_curvs[j-1].size();
            dims[1] = 1; 
            hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

            char buffer [6];
            sprintf (buffer, "%d", j);

            //dataset for curvs
            hid_t dataset_id = H5Dcreate1(group_id, buffer, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
            
            H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, curvs);

            delete curvs;

            H5Dclose(dataset_id);
        }

        //Close the group
        H5Gclose(group_id);
    }

    if (grain_junction.size()>0)
    {
        int nr_grain_junctions=0;

        for (int j=0; j<grain_junction.size(); j++)
            if(grain_junction[j]) nr_grain_junctions++;

        if (nr_grain_junctions>0)
        {
            //dataspace for grain_junction_index
            hsize_t dims[2];
            dims[0] = nr_grain_junctions; 
            dims[1] = 4; 
            hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

            //dataset for grain_junction_index    
            hid_t dataset_id = H5Dcreate1(file_save, "grain_junction_index", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);

            // Datastructure to write
            int * data_out_1 = new int[nr_grain_junctions];
            int * data_out_2 = new int[nr_grain_junctions];
            int * data_out_3 = new int[nr_grain_junctions];
            int * data_out_4 = new int[nr_grain_junctions];

            int index=0;

            for (int j=0; j<grain_junction_index.size(); j++)
            {
                if(grain_junction[j])
                {
                    if(grain_junction_index[j].size()>0) data_out_1[index]=grain_junction_index[j][0]+1;
                    else data_out_1[index]=0;
                    if(grain_junction_index[j].size()>1) data_out_2[index]=grain_junction_index[j][1]+1;
                    else data_out_2[index]=0;
                    if(grain_junction_index[j].size()>2) data_out_3[index]=grain_junction_index[j][2]+1;
                    else data_out_3[index]=0;
                    if(grain_junction_index[j].size()>3) data_out_4[index]=grain_junction_index[j][3]+1;
                    else data_out_4[index]=0;

                    index++;
                }
            }

            //dataspace for one coordinate
            hsize_t column[2];
            column[0] = nr_grain_junctions;
            column[1] = 1;  
            hid_t mdataspace_id = H5Screate_simple(2, column, NULL);

            // select file hyperslab 
            hsize_t start_h[2];// start of hyperslab
            hsize_t count[2];// block count
            
            count[0]  = nr_grain_junctions; 
            count[1]  = 1;

            start_h[0]  = 0; 
            start_h[1]  = 0;

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dwrite(dataset_id, H5T_NATIVE_INT, mdataspace_id, dataspace_id, H5P_DEFAULT, data_out_1);

            start_h[1]  = 1;

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dwrite(dataset_id, H5T_NATIVE_INT, mdataspace_id, dataspace_id, H5P_DEFAULT, data_out_2);

            start_h[1]  = 2;

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dwrite(dataset_id, H5T_NATIVE_INT, mdataspace_id, dataspace_id, H5P_DEFAULT, data_out_3);

            start_h[1]  = 3;

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_h, NULL, count, NULL);
            H5Dwrite(dataset_id, H5T_NATIVE_INT, mdataspace_id, dataspace_id, H5P_DEFAULT, data_out_4);

            // Close the memoryspace
            H5Sclose(mdataspace_id);

            // Close the filespace
            H5Sclose(dataspace_id);

            delete data_out_1;
            delete data_out_2;
            delete data_out_3;
            delete data_out_4;

            H5Dclose(dataset_id);
        }
    }

    {
        //dataspace for number turning points
        hsize_t dims[2];
        dims[0] = turning_point.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        float * number = new float[turning_point.size()];
        for (int i=0; i<turning_point.size(); i++)
        {
            number[i]=turning_point[i];
        }

        //dataset for number turning points
        hid_t dataset_id = H5Dcreate1(file_save, "number_turning_points", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, number);

        delete number;

        H5Dclose(dataset_id);
    }        

    {
        //Create a group for grain area boundaries
        hid_t group_id = H5Gcreate1(file_save, "/grain_area_boundaries", 0);       
        for(int i = 0; i < grain_area_boundaries.size(); i++)
        {            
            //write only non-empty entries to file
            if(grain_area_boundaries[i].size() > 0)
            {
                //Array to store the values
                unsigned int * vals = new unsigned int[grain_area_boundaries[i].size()];
                for(int j = 0; j < grain_area_boundaries[i].size(); j++)
                {
                    vals[j] = grain_area_boundaries[i][j]; 
                }

                //Data space for grain area boundaries
                hsize_t dims[2];
                dims[0] = grain_area_boundaries[i].size();
                dims[1] = 1;                
                hid_t dataspace_id  = H5Screate_simple(2, dims, NULL);

                char buffer[5];
                sprintf(buffer, "%d", i);

                //Data set for grain area boundaries
                hid_t dataset_id = H5Dcreate1(group_id, buffer, H5T_NATIVE_UINT, dataspace_id, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);

                delete vals;

                H5Dclose(dataset_id);
            }
        }

        //Close the group
        H5Gclose(group_id);            
    }
    
    //Close the file
    H5Fclose(file_save);

    std::cout<<"File closed"<<std::endl;
}

