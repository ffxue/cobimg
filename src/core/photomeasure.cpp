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

#include <limits>
 
// opencv 2.4.13 tested
#include <opencv/cv.h>	
#include <opencv/highgui.h>

#include "photomeasure.h"

namespace cobimg
{
    PhotoMeasurement::PhotoMeasurement(int metric_id)
    {
        setType(Type::PHOTO_2D);
        metric = Metric::SSIM;
        imgData = NULL;
        img1 = NULL;
        
        if (metric_id == 1)
            metric = Metric::MSSSIM;
        else if (metric_id == 2)
            metric = Metric::MSE;
        else if (metric_id == 3)
            metric = Metric::PSNR;
        else if (metric_id == 4)
            metric = Metric::VIFP;
    };
    PhotoMeasurement::~PhotoMeasurement()
    {
        clear();
    };
    
    void PhotoMeasurement::setMetric(Metric m)
    {
        metric = m;
    };
    
    PhotoMeasurement::Metric PhotoMeasurement::getMetric()
    {
        return metric;
    };
    bool PhotoMeasurement::readFile (char* filename, bool isColor)
    {
        int cvLoadOption = isColor ? CV_LOAD_IMAGE_COLOR : CV_LOAD_IMAGE_GRAYSCALE;
        imgData = cvLoadImage(filename, cvLoadOption);
    };
    IplImage* PhotoMeasurement::getImgData()
    {
        return imgData;
    };
    double PhotoMeasurement::evaluate (char* filename)
    {
        if (imgData == NULL)
            return std::numeric_limits<double>::max();
        PhotoMeasurement* solution = new PhotoMeasurement(metric);
        solution->readFile(filename, imgData->nChannels > 1);
        if (solution->getImgData() == NULL)
            return std::numeric_limits<double>::max();
        double value = std::numeric_limits<double>::max();
        // allocate memory and prelimianry computing for the first time use
        if (img1 == NULL)
            prepareMats();
        switch (metric)
        {
            case Metric::SSIM:
                value = evaluateSSIM(solution->getImgData());
                break;
            case Metric::MSSSIM:
                value = evaluateMSSSIM(solution->getImgData());
                break;
            case Metric::MSE:
                value = evaluateMSE(solution->getImgData());
                break;
            case Metric::PSNR:
                value = evaluatePSNR(solution->getImgData());
                break;
            case Metric::VIFP:
                value = evaluateVIFP(solution->getImgData());
                break;
                
        }
        solution->clear();
        return value;
    };
    IplImage* PhotoMeasurement::ensureMem(IplImage* cvmat)
    {
        if (cvmat != NULL)
            return cvmat;
        if (imgData == NULL)
            return NULL;
        CvSize size = cvSize(imgData->width, imgData->height);
        int depth = IPL_DEPTH_32F;
        return cvCreateImage(size, depth, imgData->nChannels);
    };
    void PhotoMeasurement::prepareMats()
    {
        if (imgData == NULL)
            return;
        img1    = ensureMem(img1);
        img1_sq = ensureMem(img1_sq);
        mu1     = ensureMem(mu1);
        mu1_sq  = ensureMem(mu1_sq);
        sigma1_sq = ensureMem(sigma1_sq);
        cvConvert(imgData, img1);
        cvPow( img1, img1_sq, 2 );
        cvSmooth( img1, mu1, CV_GAUSSIAN, 11, 11, 1.5 );
        cvPow( mu1, mu1_sq, 2 );
        cvSmooth( img1_sq, sigma1_sq, CV_GAUSSIAN, 11, 11, 1.5 );
        cvAddWeighted( sigma1_sq, 1, mu1_sq, -1, 0, sigma1_sq );
    }
    double PhotoMeasurement::evaluateSSIM (const IplImage* solution)
    {
        /**
         * The equivalent of Zhou Wang's SSIM matlab code using OpenCV.
         * from http://www.cns.nyu.edu/~zwang/files/research/ssim/index.html
         * The measure is described in :
         * "Image quality assessment: From error measurement to structural similarity"
         */
          
        // allocate memory for the first time use
        img2    = ensureMem(img2);
        img1_img2 = ensureMem(img1_img2);
        img2_sq = ensureMem(img2_sq);
        mu2     = ensureMem(mu2);
        mu1_mu2 = ensureMem(mu1_mu2);
        mu2_sq  = ensureMem(mu2_sq);
        sigma12 = ensureMem(sigma12);
        sigma2_sq = ensureMem(sigma2_sq);
        ssim_map = ensureMem(ssim_map);
        temp1   = ensureMem(temp1);
        temp2   = ensureMem(temp2);
        temp3   = ensureMem(temp3);
        
        // default parameters of SSIM
        const double C1 = 6.5025, C2 = 58.5225;
        
        cvConvert(solution, img2);
        cvPow( img2, img2_sq, 2 );
        cvMul( img1, img2, img1_img2, 1 );
        
        
        // PRELIMINARY COMPUTING
        cvSmooth( img2, mu2, CV_GAUSSIAN, 11, 11, 1.5 );        
        cvPow( mu2, mu2_sq, 2 );
        cvMul( mu1, mu2, mu1_mu2, 1 );
  
        cvSmooth( img2_sq, sigma2_sq, CV_GAUSSIAN, 11, 11, 1.5 );
        cvAddWeighted( sigma2_sq, 1, mu2_sq, -1, 0, sigma2_sq );
        cvSmooth( img1_img2, sigma12, CV_GAUSSIAN, 11, 11, 1.5 );
        cvAddWeighted( sigma12, 1, mu1_mu2, -1, 0, sigma12 );
        
        // FORMULA

        // (2*mu1_mu2 + C1)
        cvScale( mu1_mu2, temp1, 2 );
        cvAddS( temp1, cvScalarAll(C1), temp1 );
        // (2*sigma12 + C2)
        cvScale( sigma12, temp2, 2 );
        cvAddS( temp2, cvScalarAll(C2), temp2 );
        // ((2*mu1_mu2 + C1).*(2*sigma12 + C2))
        cvMul( temp1, temp2, temp3, 1 );
        // (mu1_sq + mu2_sq + C1)
        cvAdd( mu1_sq, mu2_sq, temp1 );
        cvAddS( temp1, cvScalarAll(C1), temp1 );
        
        // (sigma1_sq + sigma2_sq + C2)
        cvAdd( sigma1_sq, sigma2_sq, temp2 );
        cvAddS( temp2, cvScalarAll(C2), temp2 );
        // ((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2))
        cvMul( temp1, temp2, temp1, 1 );
        // ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2))
        cvDiv( temp3, temp1, ssim_map, 1 );

        CvScalar index_scalar = cvAvg( ssim_map );
        double ssim = 0.0;
        for (int i = 0; i < imgData->nChannels; i++)
            ssim += index_scalar.val[i];
        ssim /= imgData->nChannels;
        
        // through observation, there is approximately 
        // 1% error max with the original matlab program
        return 1.0 - ssim;
    };
    
    
    double PhotoMeasurement::evaluateMSSSIM (const IplImage* solution)
    {
        const double WEIGHT[] = {0.0448, 0.2856, 0.3001, 0.2363, 0.1333};
        const int NLEVELS = 5;
        
        double mssim[NLEVELS];
        double mcs[NLEVELS];

        cv::Mat im1[NLEVELS];
        cv::Mat im2[NLEVELS];
        
        int w = img1->width;
        int h = img1->height;
        
        img2    = ensureMem(img2);
        cvConvert(solution, img2);
        
        cv::Mat refMat(img1);
        cv::Mat solMat(img2);
        refMat.copyTo(im1[0]);
        solMat.copyTo(im2[0]);
        
        for (int l=0; l<NLEVELS; l++)
        {
            cv::Scalar res = computeSSIM(im1[l], im2[l]);
            mssim[l] = res.val[0];
            mcs[l] = res.val[1];

            if (l < NLEVELS-1)
            {
                w /= 2;
                h /= 2;
                im1[l+1] = cv::Mat(h, w, IPL_DEPTH_32F);
                im2[l+1] = cv::Mat(h, w, IPL_DEPTH_32F);
                
                cv::resize(im1[l], im1[l+1], cv::Size(w,h), 0, 0, cv::INTER_LINEAR);
                cv::resize(im2[l], im2[l+1], cv::Size(w,h), 0, 0, cv::INTER_LINEAR);
            }
        }

        // overall_mssim = prod(mcs_array(1:level-1).^weight(1:level-1))*mssim_array(level);
        double msssim = fabs(mssim[NLEVELS-1]);
        for (int l=0; l < NLEVELS-1; l++)
            msssim *= std::pow(fabs(mcs[l]), WEIGHT[l]);

        return 1.0 - msssim;
    }
    
    cv::Scalar PhotoMeasurement::computeSSIM(const cv::Mat& mimg1, const cv::Mat& mimg2)
    {
        const float C1 = 6.5025f;
        const float C2 = 58.5225f;
        int ht = mimg1.rows;
        int wt = mimg1.cols;
        int w = wt - 10;
        int h = ht - 10;
        
        cv::Mat mmu1(h,w,IPL_DEPTH_32F), mmu2(h,w,IPL_DEPTH_32F);
        cv::Mat mmu1_sq(h,w,IPL_DEPTH_32F), mmu2_sq(h,w,IPL_DEPTH_32F), mmu1_mmu2(h,w,IPL_DEPTH_32F);
        cv::Mat mimg1_sq(ht,wt,IPL_DEPTH_32F), mimg2_sq(ht,wt,IPL_DEPTH_32F), mimg1_mimg2(ht,wt,IPL_DEPTH_32F);
        cv::Mat msigma1_sq(h,w,IPL_DEPTH_32F), msigma2_sq(h,w,IPL_DEPTH_32F), msigma12(h,w,IPL_DEPTH_32F);
        cv::Mat tmp1(h,w,IPL_DEPTH_32F), tmp2(h,w,IPL_DEPTH_32F), tmp3(h,w,IPL_DEPTH_32F);
        cv::Mat ssim_map(h,w,IPL_DEPTH_32F), cs_map(h,w,IPL_DEPTH_32F);

        applyGaussianBlur(mimg1, mmu1, 11, 1.5);
        applyGaussianBlur(mimg2, mmu2, 11, 1.5);

        cv::pow(mmu1, 2, mmu1_sq);
        cv::pow(mmu2, 2, mmu2_sq);
        cv::multiply(mmu1, mmu2, mmu1_mmu2);
        cv::multiply(mimg1, mimg2, mimg1_mimg2);
        cv::pow(mimg1, 2, mimg1_sq);
        cv::pow(mimg2, 2, mimg2_sq);

        // msigma1_sq = filter2(window, mimg1.*mimg1, 'valid') - mmu1_sq;
        applyGaussianBlur(mimg1_sq, msigma1_sq, 11, 1.5);
        msigma1_sq -= mmu1_sq;
        // msigma2_sq = filter2(window, mimg2.*mimg2, 'valid') - mmu2_sq;
        applyGaussianBlur(mimg2_sq, msigma2_sq, 11, 1.5);
        msigma2_sq -= mmu2_sq;
        // msigma12 = filter2(window, mimg1.*mimg2, 'valid') - mmu1_mmu2;
        applyGaussianBlur(mimg1_mimg2, msigma12, 11, 1.5);
        msigma12 -= mmu1_mmu2;

        // cs_map = (2*msigma12 + C2)./(msigma1_sq + msigma2_sq + C2);
        tmp1 = 2*msigma12 + C2;
        tmp2 = msigma1_sq + msigma2_sq + C2;
        cv::divide(tmp1, tmp2, cs_map);
        // ssim_map = ((2*mmu1_mmu2 + C1).*(2*msigma12 + C2))./((mmu1_sq + mmu2_sq + C1).*(msigma1_sq + msigma2_sq + C2));
        tmp3 = 2*mmu1_mmu2 + C1;
        cv::multiply(tmp1, tmp3, tmp1);
        tmp3 = mmu1_sq + mmu2_sq + C1;
        cv::multiply(tmp2, tmp3, tmp2);
        cv::divide(tmp1, tmp2, ssim_map);

        // mssim = mean2(ssim_map);
        double mssim = cv::mean(ssim_map).val[0];
        // mcs = mean2(cs_map);
        double mcs = cv::mean(cs_map).val[0];

        cv::Scalar res(mssim, mcs);
        
        // free mem
        mmu1.release(); mmu2.release();
        mmu1_sq.release(); mmu2_sq.release(); mmu1_mmu2.release();
        mimg1_sq.release(); mimg2_sq.release(); mimg1_mimg2.release();
        msigma1_sq.release(); msigma2_sq.release(); msigma12.release();
        tmp1.release(); tmp2.release(); tmp3.release();
        ssim_map.release(); cs_map.release();
        
        return res;
    }
    
    void PhotoMeasurement::applyGaussianBlur(const cv::Mat& src, cv::Mat& dst, int ksize, double sigma)
    {
        int invalid = (ksize-1)/2;
        cv::Mat tmp(src.rows, src.cols, IPL_DEPTH_32F);
        cv::GaussianBlur(src, tmp, cv::Size(ksize,ksize), sigma);
        tmp(cv::Range(invalid, tmp.rows-invalid), cv::Range(invalid, tmp.cols-invalid)).copyTo(dst);
        tmp.release();
    }
    
    double PhotoMeasurement::evaluateMSE (const IplImage* solution)
    {
        img2    = ensureMem(img2);
        cvConvert(solution, img2);
        cv::Mat mat1(img1);
        cv::Mat mat2(img2);
        cv::Mat tmpMat(imgData->height, imgData->width, IPL_DEPTH_32F);
        cv::subtract(mat1, mat2, tmpMat);
        cv::multiply(tmpMat, tmpMat, tmpMat);
        return double(cv::mean(tmpMat).val[0]/(255*255));
    };
    double PhotoMeasurement::evaluatePSNR (const IplImage* solution)
    {
        img2    = ensureMem(img2);
        cvConvert(solution, img2);
        cv::Mat mat1(img1);
        cv::Mat mat2(img2);
        cv::Mat tmpMat(imgData->height, imgData->width, IPL_DEPTH_32F);
        cv::subtract(mat1, mat2, tmpMat);
        cv::multiply(tmpMat, tmpMat, tmpMat);
        return -double(10*log10(255*255/cv::mean(tmpMat).val[0]));
    };
    // H.R. Sheikh and A.C. Bovik, "Image information and visual quality," IEEE Transactions on
    // Image Processing, vol. 15, no. 2, pp. 430-444, February 2006.
    double PhotoMeasurement::evaluateVIFP (const IplImage* solution)
    {
        img2    = ensureMem(img2);
        cvConvert(solution, img2);
        cv::Mat mat1(img1);
        cv::Mat mat2(img2);
        //cv::Mat tmpMat(imgData->height, imgData->width, IPL_DEPTH_32F);
        
        const int NLEVS = 4;
        
        double num = 0.0;
        double den = 0.0;
        
        cv::Mat ref[NLEVS];
        cv::Mat dist[NLEVS];
        cv::Mat tmp1, tmp2;
        
        int w = imgData->width;
        int h = imgData->height;
        
        // for scale=1:4
        for (int scale=0; scale<NLEVS; scale++) 
        {
            // N=2^(4-scale+1)+1;
            int N = (2 << (NLEVS-scale-1)) + 1;
            
            if (scale == 0) 
            {
                mat1.copyTo(ref[scale]);
                mat2.copyTo(dist[scale]);
            }
            else 
            {
                // ref=filter2(win,ref,'valid');
                applyGaussianBlur(ref[scale-1], tmp1, N, N/5.0);
                // dist=filter2(win,dist,'valid');
                applyGaussianBlur(dist[scale-1], tmp2, N, N/5.0);
                
                w = (w-(N-1)) / 2;
                h = (h-(N-1)) / 2;
                
                ref[scale] = cv::Mat(h,w,IPL_DEPTH_32F);
                dist[scale] = cv::Mat(h,w,IPL_DEPTH_32F);
                
                // ref=ref(1:2:end,1:2:end);
                cv::resize(tmp1, ref[scale], cv::Size(w,h), 0, 0, cv::INTER_NEAREST);
                // dist=dist(1:2:end,1:2:end);
                cv::resize(tmp2, dist[scale], cv::Size(w,h), 0, 0, cv::INTER_NEAREST);
            }
            
            computeVIFP(ref[scale], dist[scale], N, num, den);
        }
        
        return 1.0 - double(num/den);
    };
    
    void PhotoMeasurement::computeVIFP(const cv::Mat& ref, const cv::Mat& dist, int N, double& num, double& den)
    {
        const float SIGMA_NSQ = 2.0f;
        
        int w = ref.cols - (N-1);
        int h = ref.rows - (N-1);
        
        cv::Mat tmp(h,w,IPL_DEPTH_32F);
        cv::Mat mu1(h,w,IPL_DEPTH_32F), mu2(h,w,IPL_DEPTH_32F), mu1_sq(h,w,IPL_DEPTH_32F), mu2_sq(h,w,IPL_DEPTH_32F), mu1_mu2(h,w,IPL_DEPTH_32F), sigma1_sq(h,w,IPL_DEPTH_32F), sigma2_sq(h,w,IPL_DEPTH_32F), sigma12(h,w,IPL_DEPTH_32F), g(h,w,IPL_DEPTH_32F), sv_sq(h,w,IPL_DEPTH_32F);
        cv::Mat sigma1_sq_th, sigma2_sq_th, g_th;
        
        // mu1 = filter2(win, ref, 'valid');
        applyGaussianBlur(ref, mu1, N, N/5.0);
        // mu2 = filter2(win, dist, 'valid');
        applyGaussianBlur(dist, mu2, N, N/5.0);
        
        const float EPSILON = 1e-10f;

        // mu1_sq = mu1.*mu1;
        cv::multiply(mu1, mu1, mu1_sq);
        // mu2_sq = mu2.*mu2;
        cv::multiply(mu2, mu2, mu2_sq);
        // mu1_mu2 = mu1.*mu2;
        cv::multiply(mu1, mu2, mu1_mu2);		
        
        // sigma1_sq = filter2(win, ref.*ref, 'valid') - mu1_sq;
        cv::multiply(ref, ref, tmp);
        applyGaussianBlur(tmp, sigma1_sq, N, N/5.0);
        sigma1_sq -= mu1_sq;
        // sigma2_sq = filter2(win, dist.*dist, 'valid') - mu2_sq;
        cv::multiply(dist, dist, tmp);
        applyGaussianBlur(tmp, sigma2_sq, N, N/5.0);
        sigma2_sq -= mu2_sq;
        // sigma12 = filter2(win, ref.*dist, 'valid') - mu1_mu2;
        cv::multiply(ref, dist, tmp);
        applyGaussianBlur(tmp, sigma12, N, N/5.0);
        sigma12 -= mu1_mu2;
        
        // sigma1_sq(sigma1_sq<0)=0;
        cv::max(sigma1_sq, 0.0f, sigma1_sq);
        // sigma2_sq(sigma2_sq<0)=0;
        cv::max(sigma2_sq, 0.0f, sigma2_sq);
        
        // g=sigma12./(sigma1_sq+1e-10);
        tmp = sigma1_sq + EPSILON;
        cv::divide(sigma12, tmp, g);
        
        // sv_sq=sigma2_sq-g.*sigma12;
        cv::multiply(g, sigma12, tmp);
        sv_sq = sigma2_sq - tmp;
        
        cv::threshold(sigma1_sq, sigma1_sq_th, EPSILON, 1.0f, cv::THRESH_BINARY);

        // g(sigma1_sq<1e-10)=0;
        cv::multiply(g, sigma1_sq_th, g);
        
        // sv_sq(sigma1_sq<1e-10)=sigma2_sq(sigma1_sq<1e-10);
        cv::multiply(sv_sq, sigma1_sq_th, sv_sq);
        cv::multiply(sigma2_sq, 1.0f - sigma1_sq_th, tmp);
        sv_sq += tmp;
        
        // sigma1_sq(sigma1_sq<1e-10)=0;
        cv::threshold(sigma1_sq, sigma1_sq, EPSILON, 1.0f, cv::THRESH_TOZERO);
        
        cv::threshold(sigma2_sq, sigma2_sq_th, EPSILON, 1.0f, cv::THRESH_BINARY);

        // g(sigma2_sq<1e-10)=0;
        cv::multiply(g, sigma2_sq_th, g);
        
        // sv_sq(sigma2_sq<1e-10)=0;
        cv::multiply(sv_sq, sigma2_sq_th, sv_sq);
        
        cv::threshold(g, g_th, 0.0f, 1.0f, cv::THRESH_BINARY);
        
        // sv_sq(g<0)=sigma2_sq(g<0);
        cv::multiply(sv_sq, g_th, sv_sq);
        cv::multiply(sigma2_sq, 1.0f - g_th, tmp);
        cv::add(sv_sq, tmp, sv_sq);
        
        // g(g<0)=0;
        cv::max(g, 0.0f, g);
        
        // sv_sq(sv_sq<=1e-10)=1e-10;
        cv::max(sv_sq, EPSILON, sv_sq);
        
        // num=num+sum(sum(log10(1+g.^2.*sigma1_sq./(sv_sq+sigma_nsq))));
        sv_sq += SIGMA_NSQ;
        cv::multiply(g, g, g);
        cv::multiply(g, sigma1_sq, g);
        cv::divide(g, sv_sq, tmp);
        tmp += 1.0f;
        cv::log(tmp, tmp);
        num += cv::sum(tmp)[0] / log(10.0f);
        
        // den=den+sum(sum(log10(1+sigma1_sq./sigma_nsq)));
        tmp = 1.0f + sigma1_sq / SIGMA_NSQ;
        cv::log(tmp, tmp);
        den += cv::sum(tmp)[0] / log(10.0f);
    };
    
    void PhotoMeasurement::clear()
    {
        if (imgData != NULL) cvReleaseImage(&imgData);
        if (img1 != NULL) cvReleaseImage(&img1);
        if (img2 != NULL) cvReleaseImage(&img2);
        if (img1_img2 != NULL) cvReleaseImage(&img1_img2);
        if (img1_sq != NULL) cvReleaseImage(&img1_sq);
        if (img2_sq != NULL) cvReleaseImage(&img2_sq);
        if (mu1 != NULL) cvReleaseImage(&mu1);
        if (mu2 != NULL) cvReleaseImage(&mu2);
        if (mu1_mu2 != NULL) cvReleaseImage(&mu1_mu2);
        if (mu1_sq != NULL) cvReleaseImage(&mu1_sq);
        if (mu2_sq != NULL) cvReleaseImage(&mu2_sq);
        if (sigma1_sq != NULL) cvReleaseImage(&sigma1_sq);
        if (sigma2_sq != NULL) cvReleaseImage(&sigma2_sq);
        if (sigma12 != NULL) cvReleaseImage(&sigma12);
        if (ssim_map != NULL) cvReleaseImage(&ssim_map);
        if (temp1 != NULL) cvReleaseImage(&temp1);
        if (temp2 != NULL) cvReleaseImage(&temp2);
        if (temp3 != NULL) cvReleaseImage(&temp3);
    };
}