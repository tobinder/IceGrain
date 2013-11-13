/*! \file remove_subgrain_arc.cpp
 * \brief Used for subgrain_extraction.
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
void remove_subgrain_arc(int arc_count,
                         marray::Marray<unsigned int> two_boundings,
                         size_t & nr_areas,
                         long * & bubble_area_size,
                         long * & grain_area_size,
                         std::vector< std::vector<int> > & old_grain_arc_index,
                         std::vector< std::vector<int> > & bubble_arc_index,
                         std::vector<size_t> & region_labels,
                         std::vector<point> & grain_area_center_mass,
                         std::vector<bool> & grain_arc,
                         std::vector<bool> & subgrain,
                         std::vector< std::vector<point> > & areas,
                         std::vector<int> & found_bubble_areas,
                         std::vector<int> & arc_state,
                         std::vector< std::vector<int> > & region_labels_reverse)
{
    if(arc_state[arc_count]==1)//grain arc
    {
        bool outside=false;
        int outside_label;

        //these two labels are selected for merging
        int label_merged = std::min(region_labels[two_boundings(arc_count,0)-1]-1,
                                    region_labels[two_boundings(arc_count,1)-1]-1);
        int label_removed = std::max(region_labels[two_boundings(arc_count,0)-1]-1,
                                     region_labels[two_boundings(arc_count,1)-1]-1);

        //look for arcs of these two labels -> no longer grain arcs
        for(int y=0;y<(int)two_boundings.shape(0);y++)
        {
            if ((region_labels[two_boundings(y,0)-1]-1==label_merged &&
                    region_labels[two_boundings(y,1)-1]-1==label_removed) ||
                (region_labels[two_boundings(y,1)-1]-1==label_merged &&
                    region_labels[two_boundings(y,0)-1]-1==label_removed))
            {
                arc_state[y]=0;
                grain_arc[y]=false;
                subgrain[y]=true;
                std::cout<<"remove subgrain arc "<<y<<" as grain arc"<<std::endl;
            }
        }

        //grain belongs now to outside -> remove as grain
        if(bubble_area_size[label_merged]==grain_area_size[label_merged])
        {
            outside=true;
            old_grain_arc_index[label_merged].clear();

            outside_label=label_merged+1;

            std::cout<<"label "<<label_removed+1<<" will be removed and merged with outside label "<<outside_label
                <<std::endl;
        }
        else if(bubble_area_size[label_removed]==grain_area_size[label_removed])
        {
            outside=true;
            old_grain_arc_index[label_removed].clear();

            outside_label=label_removed+1;
            label_removed=label_merged;

            std::cout<<"label "<<label_removed+1<<" will be removed and merged with outside label "<<outside_label
                <<std::endl;
            outside_label--;
        }
        else
        {
            //combine arc lists
            std::list<int> combined_grain_arc_list;

            for (int a=0; a<old_grain_arc_index[label_merged].size(); a++)
                combined_grain_arc_list.push_back(old_grain_arc_index[label_merged][a]);
            for (int a=0; a<old_grain_arc_index[label_removed].size(); a++)
                combined_grain_arc_list.push_back(old_grain_arc_index[label_removed][a]);

            //sort combined list und remove double entries
            combined_grain_arc_list.sort();
            combined_grain_arc_list.unique();

            //replace arc list of label_merged
            old_grain_arc_index[label_merged].clear();
            for (std::list<int>::iterator a=combined_grain_arc_list.begin(); a!=combined_grain_arc_list.end();++a) 
            {
                if (arc_state[*a]>0) old_grain_arc_index[label_merged].push_back(*a);//check if arc is bubble or grain arc
            }

            std::cout<<"label "<<label_removed+1<<" will be removed and merged with label "<<label_merged+1<<std::endl;
        }

        //erase label_removed -> decrease all labels higher than label_removed by one
        old_grain_arc_index.erase(old_grain_arc_index.begin()+label_removed);
        bubble_arc_index.erase(bubble_arc_index.begin()+label_removed);

        for(int l=0; l<region_labels.size(); l++)
            if (region_labels[l]>label_removed) region_labels[l]--;

        //update areas and area size
        if (!outside)
        {
            for (int p=0; p<areas[label_removed].size(); p++)
                areas[label_merged].push_back(areas[label_removed][p]);
            grain_area_size[label_merged]=areas[label_merged].size();
        }

        //areas.erase(areas.begin()+label_removed);
        //"erase" doesn't free memory, so swaping is necessary
        std::vector< std::vector<point> > areas2;

        for (int area=label_removed+1; area<areas.size(); area++)
            areas2.push_back(areas[area]);

        areas.resize(label_removed);

        for (int area=0; area<areas2.size(); area++)
            areas.push_back(areas2[area]);

        areas2.clear();

        grain_area_center_mass.erase(grain_area_center_mass.begin()+label_removed);
        found_bubble_areas.clear();
        nr_areas--;

        delete grain_area_size;
        delete bubble_area_size;
        grain_area_size = new long[nr_areas];
        bubble_area_size = new long[nr_areas];

        for (int area=0; area<nr_areas; area++)
        {
            if (old_grain_arc_index[area].size()>0) grain_area_size[area]=areas[area].size();
            else grain_area_size[area]=0;

            if (bubble_arc_index[area].size()>0)
            {
                bubble_area_size[area]=areas[area].size();
                found_bubble_areas.push_back(area);
            }
            else bubble_area_size[area]=0;
        }

        if (!outside)
        {
            //calculate new center of mass of merged grain
            long grain_area_x_sum=0;
            long grain_area_y_sum=0;

            //loop over all pixels with label
            for (int pixel=0;pixel<areas[label_merged].size();pixel++)
            {
                point p=areas[label_merged][pixel];
                grain_area_x_sum+=p.x;
                grain_area_y_sum+=p.y;
            }

            int xx=grain_area_x_sum/grain_area_size[label_merged];
            int yy=grain_area_y_sum/grain_area_size[label_merged];
            grain_area_center_mass[label_merged].x=xx;
            grain_area_center_mass[label_merged].y=yy;

            //update region_labels and region_labels_reverse
            for(int l=0; l<region_labels_reverse[label_removed].size(); l++)
            {
                region_labels_reverse[label_merged].push_back(region_labels_reverse[label_removed][l]);
                region_labels[region_labels_reverse[label_removed][l]-1]=label_merged+1;
            }

            region_labels_reverse.erase(region_labels_reverse.begin()+label_removed);
        }
        else
        {
            //update region_labels and region_labels_reverse
            for(int l=0; l<region_labels_reverse[label_removed].size(); l++)
            {
                region_labels[region_labels_reverse[label_removed][l]-1]=outside_label;
            }

            region_labels_reverse.erase(region_labels_reverse.begin()+label_removed);
        }
    }
}
