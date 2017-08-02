/**
 * COBIMG, Constrained Optimization-based BIM Generator
 * Copyright (c) 2016-2017 The University of Hong Kong
 * Author: Fan Xue <xuef@hku.hk; fanxue@outlook.com>
 *
 * This file is part of cobimg.
 *
 * COBIMG is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * COBIMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with COBIMG.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#ifndef COBIMG_MEASUREMENT_H
#define COBIMG_MEASUREMENT_H

namespace cobimg
{
    class Measurement
    {
    public:
        enum Type
        {
            PHOTO_2D,
            POINT_CLOUD_3D, // TODO
            POLYGON_3D, // TODO
        };
    private:
        Type type;
        
    public:
        virtual bool readFile (char* filename, bool isColor = true) = 0;
        virtual double evaluate (char* filename) = 0;
        virtual void clear();
        void setType(Type t)
        {
            type = t;
        };
        Type getType()
        {
            return type;
        };
    };
}

#endif