/*! \file gbn.h 
 *	\brief Class representing the Grain-Boundary-Network structure.
 */
/*! 
 * The macro <tt>GBN_STANDALONE</tt> enables the use of this class in standalone mode. However, the header files 
 * <tt>boundary_probabilities.h</tt> and <tt>boundary_data_structures.h</tt> must be available for this to work.
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

#ifdef GBN_STANDALONE
//This include becomes necessary only if this class is included seperately in a project
#include "boundary_probabilities.h"
#include "boundary_data_structure.h"
#endif

/*Class representing the Grain-Boundary-Network structure (cf. code manual for IceGrain)*/

class gbn
{
    public:
    //Attributes *************************************************************
    std::vector<bool>               DEFAULT_VECTOR;                     /*!< Default empty vector.*/
    size_t                          nr_new_areas;						/*!< Number of areas.*/
    long *                          bubble_area_size;                   /*!< Size in pixels of each bubble area.*/
    long *                          grain_area_size;                    /*!< Size in pixels of each grain area.*/
    std::vector< std::vector<int> > grain_arc_index;                    /*!< Arc inidices of each grain arc.*/
    std::vector< std::vector<int> > bubble_arc_index;                   /*!< Arc indices of each bubble arc.*/
    std::vector<size_t>             region_labels;                      /*!< Region labels.*/
    std::vector<int>                grain_junctions;                    /*!< Grain junctions.*/
    std::vector<int>                grain_bubble_junctions;             /*!< Grain-bubble junctions.*/
    std::vector<bool>               grain_arc;                          /*!< Indices of grain arcs.*/
    std::vector<point>              grain_area_center_mass;             /*!< Center of mass positions for each grain area.*/
    std::vector<point>              bubble_area_center_mass;            /*!< Center of mass positions for each bubble area.*/
    std::vector<bool>               subgrain_arcs;                      /*!< Indices of bubble arcs.*/

    std::vector<bool>               found_bubble_arcs;                  /*!< Indices of found bubble arcs.*/
    std::vector<int>                found_bubble_areas;                 /*!< Indices of found grain areas.*/
    std::vector<float>              values;                             /*!< Unused as of now.*/

    ParameterFile *                 paramFileGBN;                       /*!< Pointer to parameter file used for gbn structure.*/

    //Methods ****************************************************************
    /*! Save the GBN structure into a HDF5 file.
     * \param dest_path Path to the target HDF5 file
     */
    void save_final_structure(std::string dest_path);
	/*! Load the GBN structure from a HDF5 file.
	 * \param filepath_new_classification Filepath to the HDF5 file
	 */
    void load_final_structure(std::string filepath_new_classification);
	/*! Load the grain sizes from a HDF5 file.
	 * \param filepath_new_classification Filepath to the HDF5 file
	 */
    void load_grain_sizes(std::string filepath_new_classification);
};

void gbn::save_final_structure(std::string dest_path)
{
    //EXPORT RESULTS TO A HDF5 file
    dest_path.append(".h5");

    std::cout<<"Exporting results to file: "<<std::endl;
    std::cout<<dest_path<<std::endl;
    
    //creating the file
    hid_t file_save=H5Fcreate(dest_path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

    {
        //dataspace for grain size
        hsize_t dims[2];
        dims[0] = nr_new_areas; 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        //dataset for grain size    
        hid_t dataset_id = H5Dcreate1(file_save, "grain_size", H5T_NATIVE_LONG, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, grain_area_size);

        H5Dclose(dataset_id);
    }

    {
        //dataspace for bubble size
        hsize_t dims[2];
        dims[0] = nr_new_areas; 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        //dataset for bubble size    
        hid_t dataset_id = H5Dcreate1(file_save, "bubble_size", H5T_NATIVE_LONG, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, bubble_area_size);

        H5Dclose(dataset_id);
    }

    {
        //dataspace for region labels
        hsize_t dims[2];
        dims[0] = region_labels.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        int * labels = new int[region_labels.size()];
        for (int i=0; i<region_labels.size(); i++)
        {
            labels[i]=(int)region_labels[i];
        }

        //dataset for region labels
        hid_t dataset_id = H5Dcreate1(file_save, "region_labels", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, labels);

        delete labels;

        H5Dclose(dataset_id);
    }

    if (grain_junctions.size()>0)
    {
        //dataspace for grain junctions
        hsize_t dims[2];
        dims[0] = grain_junctions.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        int * junctions = new int[grain_junctions.size()];
        for (int i=0; i<grain_junctions.size(); i++)
        {
            junctions[i]=grain_junctions[i];
        }

        //dataset for grain junctions    
        hid_t dataset_id = H5Dcreate1(file_save, "grain_junctions", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, junctions);

        delete junctions;

        H5Dclose(dataset_id);
    }

    if (grain_bubble_junctions.size()>0)
    {
        //dataspace for grain bubble junctions
        hsize_t dims[2];
        dims[0] = grain_bubble_junctions.size(); 
        dims[1] = 1; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        int * junctions = new int[grain_bubble_junctions.size()];
        for (int i=0; i<grain_bubble_junctions.size(); i++)
        {
            junctions[i]=grain_bubble_junctions[i];
        }

        //dataset for grain bubble junctions    
        hid_t dataset_id = H5Dcreate1(file_save, "grain_bubble_junctions", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, junctions);

        delete junctions;

        H5Dclose(dataset_id);
    }

    {
        //Create group grain_arcs_index
        hid_t group_id = H5Gcreate1(file_save, "/grain_arcs_index", 0);

        for(int j=1; j<grain_arc_index.size()+1; j++)
        {
            if (grain_arc_index[j-1].size()==0) continue;

            //array to store arcs
            int * arcs=new int[grain_arc_index[j-1].size()];

            for(int x=0;x<grain_arc_index[j-1].size();x++)
            {
                arcs[x]=grain_arc_index[j-1][x]+1;
            }

            //dataspace for arcs
            hsize_t dims[2];
            dims[0] = grain_arc_index[j-1].size();
            dims[1] = 1; 
            hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

            char buffer [5];
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
        //dataspace for grain_center_of_mass
        hsize_t dims[2];
        dims[0] = grain_arc_index.size(); 
        dims[1] = 2; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        //dataset for grain_center_of_mass    
        hid_t dataset_id = H5Dcreate1(file_save, "grain_center_of_mass", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);

        // Datastructure to write
        int * data_out_x = new int[grain_arc_index.size()];
        int * data_out_y = new int[grain_arc_index.size()];

        for (int j=0;j<grain_arc_index.size();j++)
        {
            if (grain_arc_index[j].size()==0)
            {
                data_out_x[j]=0;
                data_out_y[j]=0;
            }
            else
            {
                data_out_x[j]=grain_area_center_mass[j].x;
                data_out_y[j]=grain_area_center_mass[j].y;
            }
        }

        //dataspace for one coordinate
        hsize_t column[2];
        column[0] = grain_arc_index.size();
        column[1] = 1;  
        hid_t mdataspace_id = H5Screate_simple(2, column, NULL);

        // select file hyperslab 
        hsize_t start_h[2];// start of hyperslab
        hsize_t count[2];// block count
        
        count[0]  = grain_arc_index.size(); 
        count[1]  = 1;

        start_h[0]  = 0; 
        start_h[1]  = 0;

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_h, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_INT, mdataspace_id, dataspace_id, H5P_DEFAULT, data_out_x);

        start_h[1]  = 1;

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_h, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_INT, mdataspace_id, dataspace_id, H5P_DEFAULT, data_out_y);

        // Close the memoryspace
        H5Sclose(mdataspace_id);

        // Close the filespace
        H5Sclose(dataspace_id);

        delete data_out_x;
        delete data_out_y;

        H5Dclose(dataset_id);

    }

    {
        //dataspace for bubble_center_of_mass
        hsize_t dims[2];
        dims[0] = bubble_arc_index.size(); 
        dims[1] = 2; 
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

        //dataset for bubble_center_of_mass    
        hid_t dataset_id = H5Dcreate1(file_save, "bubble_center_of_mass", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);

        // Datastructure to write
        int * data_out_x = new int[bubble_arc_index.size()];
        int * data_out_y = new int[bubble_arc_index.size()];

        for (int j=0;j<bubble_arc_index.size();j++)
        {
            if (bubble_arc_index[j].size()==0)
            {
                data_out_x[j]=0;
                data_out_y[j]=0;
            }
            else
            {
                data_out_x[j]=bubble_area_center_mass[j].x;
                data_out_y[j]=bubble_area_center_mass[j].y;
            }
        }

        //dataspace for one coordinate
        hsize_t column[2];
        column[0] = bubble_arc_index.size();
        column[1] = 1;  
        hid_t mdataspace_id = H5Screate_simple(2, column, NULL);

        // select file hyperslab 
        hsize_t start_h[2];// start of hyperslab
        hsize_t count[2];// block count
        
        count[0]  = bubble_arc_index.size(); 
        count[1]  = 1;

        start_h[0]  = 0; 
        start_h[1]  = 0;

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_h, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_INT, mdataspace_id, dataspace_id, H5P_DEFAULT, data_out_x);

        start_h[1]  = 1;

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_h, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_INT, mdataspace_id, dataspace_id, H5P_DEFAULT, data_out_y);

        // Close the memoryspace
        H5Sclose(mdataspace_id);

        // Close the filespace
        H5Sclose(dataspace_id);

        delete data_out_x;
        delete data_out_y;

        H5Dclose(dataset_id);

    }

    {
        //Create group bubble_arcs
        hid_t group_id = H5Gcreate1(file_save, "/bubble_arcs", 0);

        for(int j=1; j<bubble_arc_index.size()+1; j++)
        {
            if (bubble_arc_index[j-1].size()==0) continue;

            //array to store arcs
            int * arcs=new int[bubble_arc_index[j-1].size()];

            for(int x=0;x<bubble_arc_index[j-1].size();x++)
            {
                arcs[x]=bubble_arc_index[j-1][x]+1;
            }

            //dataspace for arcs
            hsize_t dims[2];
            dims[0] = bubble_arc_index[j-1].size();
            dims[1] = 1; 
            hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

            char buffer [5];
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

    if (subgrain_arcs.size()!=DEFAULT_VECTOR.size());
    {
        int size=0;
        for (int i=0; i<subgrain_arcs.size(); i++)
        {
            if(subgrain_arcs[i]==true) size++;
        }

        if (size>0)
        {
            //dataspace for subgrain arcs
            hsize_t dims[2];
            dims[0] = size; 
            dims[1] = 1; 
            hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

            int * arcs = new int[size];

            int j=0;
            for (int i=0; i<subgrain_arcs.size(); i++)
            {
                if(subgrain_arcs[i]==true)
                {
                    arcs[j]=i+1;
                    j++;
                }
            }

            //dataset for subgrain arcs
            hid_t dataset_id = H5Dcreate1(file_save, "subgrain_arcs", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);

            H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arcs);

            delete arcs;

            H5Dclose(dataset_id);
        }
    }

    {
        //Create group Parameters
        hid_t group_id = H5Gcreate1(file_save, "/Parameters", 0);

        //Write the parameters into the group
        Parameter<int>::writeParamToHDF5("feature1", 0, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("feature2", 0, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("feature3", 0, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("feature4", 0, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("feature5", 0, group_id, paramFileGBN);

        Parameter<int>::writeParamToHDF5("curvature", 0, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("probmap", 0, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("arcsize", 0, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("region", 0, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("cross_section", 0, group_id, paramFileGBN);

        Parameter<std::string>::writeParamToHDF5("path_boundary_rf", "", group_id, paramFileGBN);

        Parameter<float>::writeParamToHDF5("bubble_boundary_threshold", (float)0.51, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("max_bubble_arc_length", 9001, group_id, paramFileGBN);
        Parameter<float>::writeParamToHDF5("subgrain_boundary_threshold", (float)0.91, group_id, paramFileGBN);
        Parameter<float>::writeParamToHDF5("no_boundary_threshold", (float)0.81, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("minimal_grain_size", 51, group_id, paramFileGBN);
        Parameter<int>::writeParamToHDF5("minimal_region_size", 50013, group_id, paramFileGBN);

        Parameter<int>::writeParamToHDF5("close_bubble_grain_size", 5001, group_id, paramFileGBN);
        //Close the group
        H5Gclose(group_id);

    }

    //Close the file
    H5Fclose(file_save);

    std::cout<<"File closed"<<std::endl;

}

void gbn::load_final_structure(std::string filepath_new_classification)
{
    std::cout<<"Importing results from file: "<<std::endl;
    std::cout<<filepath_new_classification<<std::endl;
    
    // Open an existing file
    hid_t file_save = H5Fopen(filepath_new_classification.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    {
        //dataset for grain size
        hid_t dataset_id = H5Dopen(file_save, "grain_size", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        nr_new_areas=dims[0];
        grain_area_size = new long[nr_new_areas];

        H5Dread(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, grain_area_size);

        H5Dclose(dataset_id);
    }

    {
        //dataset for bubble size    
        hid_t dataset_id = H5Dopen(file_save, "bubble_size", H5P_DEFAULT);

        bubble_area_size = new long[nr_new_areas];

        H5Dread(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, bubble_area_size);

        H5Dclose(dataset_id);
    }

    {
        //dataset for region labels
        hid_t dataset_id = H5Dopen(file_save, "region_labels", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        int * labels = new int[dims[0]];

        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, labels);

        for (int i=0; i<dims[0]; i++)
        {
            region_labels.push_back((size_t)labels[i]);
        }

        delete labels;

        H5Dclose(dataset_id);
    }

    {
        htri_t check = H5Lexists(file_save, "grain_junctions", H5P_DEFAULT);

        if(check)
        {
            //dataset for grain junctions
            hid_t dataset_id = H5Dopen(file_save, "grain_junctions", H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            int * junctions = new int[dims[0]];

            H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, junctions);

            for (int i=0; i<dims[0]; i++)
            {
                grain_junctions.push_back(junctions[i]);
            }

            delete junctions;

            H5Dclose(dataset_id);
        }
    }

    {
        htri_t check = H5Lexists(file_save, "grain_bubble_junctions", H5P_DEFAULT);

        if(check)
        {
            //dataset for grain_bubble junctions
            hid_t dataset_id = H5Dopen(file_save, "grain_bubble_junctions", H5P_DEFAULT);

            // Get filespace handle
            hid_t filespace = H5Dget_space(dataset_id);

            // Get dimensions
            hsize_t dims[2];
            H5Sget_simple_extent_dims(filespace, dims, NULL);

            int * junctions = new int[dims[0]];

            H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, junctions);

            for (int i=0; i<dims[0]; i++)
            {
                grain_bubble_junctions.push_back(junctions[i]);
            }

            delete junctions;

            H5Dclose(dataset_id);
        }
        else grain_bubble_junctions.resize(0);
    }

    {
        //Open an existing group
        hid_t group_id = H5Gopen(file_save, "/grain_arcs_index", H5P_DEFAULT);

        grain_arc_index.resize(nr_new_areas);
        grain_arc.resize(1,false);//an array that enables grain arc indeces

        for(int j=1; j<nr_new_areas+1; j++)
        {
            if (grain_area_size[j-1]==0) continue;

            char buffer [5];
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
                grain_arc_index[j-1].push_back(arcs[x]-1);
                if (arcs[x]>grain_arc.size()) grain_arc.resize(arcs[x],false);
                grain_arc[arcs[x]-1]=true;
            }

            delete arcs;

            H5Dclose(dataset_id);
        }

        //Close the group
        H5Gclose(group_id);

    }

    {
        //Open an existing group
        hid_t group_id = H5Gopen(file_save, "/bubble_arcs", H5P_DEFAULT);

        bubble_arc_index.resize(nr_new_areas);

        for(int j=1; j<nr_new_areas+1; j++)
        {
            if (bubble_area_size[j-1]==0) continue;
            else found_bubble_areas.push_back(j);

            char buffer [5];
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
                bubble_arc_index[j-1].push_back(arcs[x]-1);
                if (arcs[x]>grain_arc.size()) grain_arc.resize(arcs[x],false);
                grain_arc[arcs[x]-1]=false;
            }

            delete arcs;

            H5Dclose(dataset_id);
        }

        //Close the group
        H5Gclose(group_id);

    }

    {
        // Open an existing dataset
        hid_t dataset_id = H5Dopen(file_save, "/grain_center_of_mass", H5P_DEFAULT);

        grain_area_center_mass.resize(nr_new_areas);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        // Datastructure to read in
        int * data_out_x = new int[dims[0]];
        int * data_out_y = new int[dims[0]];

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
        H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_x);

        start_h[1]  = 1;

        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
        H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_y);

        for (int j=0;j<dims[0];j++)
        {
            grain_area_center_mass[j].x=data_out_x[j];
            grain_area_center_mass[j].y=data_out_y[j];
        }

        // Close the memoryspace
        H5Sclose(mdataspace_id);

        // Close the filespace
        H5Sclose(filespace);

        delete data_out_x;
        delete data_out_y;

        // Close the dataset
        H5Dclose(dataset_id);
    }

    {
        // Open an existing dataset
        hid_t dataset_id = H5Dopen(file_save, "/bubble_center_of_mass", H5P_DEFAULT);

        bubble_area_center_mass.resize(nr_new_areas);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        // Datastructure to read in
        int * data_out_x = new int[dims[0]];
        int * data_out_y = new int[dims[0]];

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
        H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_x);

        start_h[1]  = 1;

        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_h, NULL, count, NULL);
        H5Dread(dataset_id, H5T_NATIVE_INT, mdataspace_id, filespace, H5P_DEFAULT, data_out_y);

        for (int j=0;j<dims[0];j++)
        {
            bubble_area_center_mass[j].x=data_out_x[j];
            bubble_area_center_mass[j].y=data_out_y[j];
        }

        // Close the memoryspace
        H5Sclose(mdataspace_id);

        // Close the filespace
        H5Sclose(filespace);

        delete data_out_x;
        delete data_out_y;

        // Close the dataset
        H5Dclose(dataset_id);
    }

    {
        htri_t check = H5Lexists(file_save, "subgrain_arcs", H5P_DEFAULT);

        if(check)
        {
            //dataset for subgrain arcs
            hid_t dataset_id = H5Dopen(file_save, "subgrain_arcs", H5P_DEFAULT);

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
                if (arcs[x]>subgrain_arcs.size()) subgrain_arcs.resize(arcs[x],false);
                subgrain_arcs[arcs[x]-1]=true;
            }

            delete arcs;

            H5Dclose(dataset_id);
        }
        subgrain_arcs.resize(grain_arc.size(),false);
    }

    //Close the file
    H5Fclose(file_save);

    std::cout<<"File closed"<<std::endl;
}

void gbn::load_grain_sizes(std::string filepath_new_classification)
{
    std::cout<<"Importing results from file: "<<std::endl;
    std::cout<<filepath_new_classification<<std::endl;
    
    // Open an existing file
    hid_t file_save = H5Fopen(filepath_new_classification.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    {
        //dataset for grain size
        hid_t dataset_id = H5Dopen(file_save, "grain_size", H5P_DEFAULT);

        // Get filespace handle
        hid_t filespace = H5Dget_space(dataset_id);

        // Get dimensions
        hsize_t dims[2];
        H5Sget_simple_extent_dims(filespace, dims, NULL);

        nr_new_areas=dims[0];
        grain_area_size = new long[nr_new_areas];

        H5Dread(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, grain_area_size);

        H5Dclose(dataset_id);
    }

    //Close the file
    H5Fclose(file_save);

    std::cout<<"File closed"<<std::endl;
}


