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


#ifndef COBIMG_PHOTO_MEASUREMENT_H
#define COBIMG_PHOTO_MEASUREMENT_H

#include "measurement.h"

#define rb_w32_pow pow
#include <opencv/cv.h>
 
namespace cobimg
{
    class PhotoMeasurement : public Measurement
    {
    public:
        enum Metric
        {
            SSIM,
            MSSSIM,
            MSE,
            PSNR,
            VIFP,
        };
        
    private:
        Metric metric;
        IplImage* imgData;
        
        IplImage *img1 = NULL, *img2 = NULL, *img1_img2 = NULL, *img1_sq = NULL, *img2_sq = NULL;
		IplImage *mu1 = NULL, *mu2 = NULL, *mu1_mu2 = NULL, *mu1_sq = NULL, *mu2_sq = NULL;
		IplImage *sigma1_sq = NULL, *sigma2_sq = NULL, *sigma12 = NULL, *ssim_map = NULL;
		IplImage *temp1 = NULL, *temp2 = NULL, *temp3 = NULL;
        
        IplImage* ensureMem(IplImage* cvmat);
        void prepareMats();
        double evaluateSSIM (const IplImage* solution);
        double evaluateMSE (const IplImage* solution);
        double evaluatePSNR (const IplImage* solution);
        double evaluateMSSSIM (const IplImage* solution);
        double evaluateVIFP (const IplImage* solution);
        cv::Scalar computeSSIM(const cv::Mat& img1, const cv::Mat& img2);
        void applyGaussianBlur(const cv::Mat& src, cv::Mat& dst, int ksize, double sigma);
        void computeVIFP(const cv::Mat& ref, const cv::Mat& dist, int N, double& num, double& den);
    public:
        PhotoMeasurement(int m = 0);
        ~PhotoMeasurement();
        void setMetric(Metric m);
        Metric getMetric();
        IplImage* getImgData();
        bool readFile (char* filename, bool isColor = true);
        double evaluate (char* filename);
        void clear();
    };
}

#endif