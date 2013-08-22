/*! \file ConvexPerimeter.h
 *  \brief Class representing a convex perimeter around a grain.
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
class ConvexPerimeter
{
    public:
    //Attributes
    std::vector<point> points;

    //Methods
    //Fit a convex perimeter to a set of pixels using Jarvis' March and return its length
    float fit(std::vector<point> pixels);

    private:
    //Methods
    //Sort by x-coordinates, then sort by y-coordinate 
    std::vector<point> sortByX(std::vector<point> points);    
};

float ConvexPerimeter::fit(std::vector<point> pixels)
{
    //Sort pixels by their x-coordinates.  
    std::vector<point> pixelsX;
    pixelsX = sortByX(pixels);

    //Find the points with the lowest and highest y-coordinates and their indices 
    point y_low;
    y_low.y = 100000;
    int y_low_index = 0;

    point y_high;
    y_high.y = 0;
    int y_high_index = 0;

    for(int a = 0; a < pixelsX.size(); a++)
    {
        if(pixelsX[a].y < y_low.y)
        {
            y_low = pixelsX[a];
            y_low_index = a;
        }
        if(pixelsX[a].y > y_high.y)
        {
            y_high = pixelsX[a];
            y_high_index = a;
        }
    }
    

    //Calculate the slopes of the edges of the resulting tetragon
    float m_leftHigh = (float)((float)y_high.y - (float)pixelsX[0].y)/((float)y_high.x - (float)pixelsX[0].x);
    float m_leftLow  = (float)((float)y_low.y - (float)pixelsX[0].y)/((float)y_low.x - (float)pixelsX[0].x);

    float m_rightHigh = (float)((float)pixelsX.back().y - (float)y_high.y)/((float)pixelsX.back().x - (float)y_high.x);
    float m_rightLow = (float)((float)pixelsX.back().y - (float)y_low.y)/((float)pixelsX.back().x - (float)y_low.x);
    
    //Disregard all points which are inside the area of this tetragon
    point empty;
    empty.x = -1;
    empty.y = -1;
    //Make sure the left- and rightmost points do not get deleted
    for(int b = 1; b < pixelsX.size()-1; b++)
    {
        //x- and y-coordinates of the current point
        int xb = pixelsX[b].x;
        int yb = pixelsX[b].y;

        //The y-coordinate of the point has to be inbetween these four values to be inside the tetragon
        int y_leftHigh = (float)pixelsX[0].y + m_leftHigh*((float)xb - (float)pixelsX[0].x);
        int y_leftLow = (float)pixelsX[0].y + m_leftLow*((float)xb - (float)pixelsX[0].x);
        int y_rightHigh = (float)y_high.y + m_rightHigh*((float)xb - (float)y_high.x);
        int y_rightLow = (float)y_low.y + m_rightLow*((float)xb - (float)y_low.x);

        if(yb <= y_leftHigh && yb >= y_leftLow && yb <= y_rightHigh && yb >= y_rightLow)
        {
            pixelsX[b] = empty;            
        }
    }

    //Use Jarvis' March algortihm to find the convex hull
    int i = 0;    

    //Length of the convex hull
    float length = 0.0f;
    //Old index, used for calculating the length of the convex hull
    int old_i = 0;

    int x0 = pixelsX[0].x;          //x-coordinate of the start point
    int xback = pixelsX.back().x;   //x-coordinate of the end point    
    int xi = pixelsX[i].x;          //x- and y-coordinates of the i-th point
    int yi = pixelsX[i].y;    
    std::vector<bool> already_used (pixelsX.size(),false);

    bool flag = false;              //Direction of the algorithm

    do
    {
        if(i == pixelsX.size()-1) flag = true;
        if(xi == -1) //Ignore points which are inside the tetragon
        {
            i++;
            xi = pixelsX[i].x;
            yi = pixelsX[i].y;
        }
        else
        {
            float low = 2*PI;
            point convex_point;
            int j_save = 0;      //j only exists in the for-loop
            float angle = 0.0f;
            
            for(int j = 0; j < pixelsX.size(); j++)
            {
                int xj = pixelsX[j].x;
                int yj = pixelsX[j].y;      
                if(xj != xi && xj != -1)
                {                    
                    point d;
                    d.x = xj - xi;
                    d.y = yj - yi;
                    if(flag == false)
                    {
                        if(d.x > 0)
                        {
                            angle = PI/2 - atan((double)(float)d.y/(float)d.x);
                        }
                        else if(d.x < 0)
                        {
                            angle = 3*PI/2 - atan((double)(float)d.y/(float)d.x);
                        }
                        else if(d.x == 0)
                        {
                            if(d.x == x0 || d.x == xback) angle = 0;
                            else angle = 2*PI;
                        }                        
                    }
                    if(flag == true)
                    {
                        if(d.x < 0)
                        {
                            angle = PI/2 - atan((double)(float)d.y/(float)d.x);
                        }
                        else if(d.x > 0)
                        {
                            angle = 3*PI/2 - atan((double)(float)d.y/(float)d.x);
                        }
                        else if(d.x == 0)
                        {
                            if(d.x == x0 || d.x == xback) angle = 0;
                            else angle = 2*PI;
                        }                        
                    }
                    if(angle < low && !already_used[j])
                    {
                       convex_point = pixelsX[j];
                       low = angle;
                       j_save = j;
                    }
                }
            }
            
            points.push_back(convex_point);
            old_i = i;
            i = j_save;           
            xi = pixelsX[i].x;
            yi = pixelsX[i].y;
            already_used[i]=true;
            
            //Calculate length of the convex hull's edges
            if(points.size() >= 1)
            {
                int x_old = pixelsX[old_i].x;
                int y_old = pixelsX[old_i].y;

                length = length + (float)sqrt((float)((float)xi-(float)x_old)*((float)xi-(float)x_old) + (float)((float)yi-(float)y_old)*((float)yi-(float)y_old));    
            }
        }      
    }
    while(xi != x0);
    
    
    /* //Debugging stuff
    std::cout << "Calculated convex perimeter length to " << length << std::endl; //Debug
    float ratio = length/(float)pixels.size();
    std::cout << "=> ratio =  " << ratio << std::endl; //Debug    
    */
    
    //points.clear();
    return length;
}

//2-dimensional sorting.
std::vector<point> ConvexPerimeter::sortByX(std::vector<point> points)
{
    std::vector<point> sorted = points;
    //Sort by x-coordinate
    for(int i = 0; i < sorted.size(); i++)
    {
        for(int j = i+1; j < sorted.size(); j++)
        {
            if(sorted[j].x < sorted[i].x)
            {
                point temp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = temp;
            }
        }
    }

    //Sort sorted x-coordinates by their y-coordinates
    for(int k = 0; k < sorted.size(); k++)
    {
        for(int l = k+1; l < sorted.size(); l++)
        {
            if(sorted[k].x == sorted[l].x)
            {
                if(sorted[l].y < sorted[k].y)
                {
                    point temp = sorted[k];
                    sorted[k] = sorted[l];
                    sorted[l] = temp;
                }
                else if(sorted[l].y == sorted[k].y)
                {
                    sorted.erase(sorted.begin()+l);
                    l--;
                }
            }
        }
    }    
    return sorted;
}
