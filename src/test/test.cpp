/**
 * OMG, Optimization-based Model Generator
 * Copyright (c) 2016-2017 The University of Hong Kong
 * Author: Fan Xue <xuef@hku.hk; fanxue@outlook.com>
 *
 * This file is part of omg.
 *
 * OMG is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with OMG.  If not, see <http://www.gnu.org/licenses/>.
 */
 
 #include "../core/photomeasure.h"
 
int main(int argc, char** argv)
{
    omg::PhotoMeasurement pm(4);
    char* fname = (char*)"1.bmp";
    char* fname1 = (char*)"2.bmp";
    bool read = pm.readFile(fname, false);
    double ret = pm.evaluate(fname1);
    std::cout << ret << std::endl;
}