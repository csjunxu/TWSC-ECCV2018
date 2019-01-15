The code in this package implements the Trilateral Weighted Sparse Coding Scheme for real color image denoising as described in the following paper:

> @article{TWSC_ECCV2018,        
>          author = {Jun Xu and Lei Zhang and David Zhang},        
>          title = {A Trilateral Weighted Sparse Coding Scheme for Real-World Image Denoising},        
>          journal = {ECCV},       
>          year = {2018}     
> }

Please cite the paper if you feel this code useful in your research.
Please see the file License.txt for the license governing this code.

  Version:       1.0 (13/07/2018), see ChangeLog.txt
  Contact:       Jun Xu <csjunxu@comp.polyu.edu.hk, nankaimathxujun@gmail.com>


Code
------------
Code will be updated on https://github.com/csjunxu/TWSC-ECCV2018

Data
------------
Please download the data from corresponding addresses.
1. cleanimages: 20 high quality commonly used natural gray scale images

2. nc: real noisy images with no ''ground truth''
                        This dataset can be found at http://demo.ipol.im/demo/125/
3. cc: 15 cropped real noisy images from CC [1]. 
                        This dataset can be found at  http://snam.ml/research/ccnoise
                        The smaller 15 cropped images can be found on in the directory 
                        ''Real_ccnoise_denoised_part'' of 
                        https://github.com/csjunxu/MCWNNM_ICCV2017
                                                The *real.png are noisy images;
                                                The *mean.png are "ground truth" images;
                                                The *ours.png are images denoised by CC.
4. dnd: The Darmstadt Noise Dataset [2] consists of 50 pairs of real noisy images, 
             each images provides 50 crops, resulting overall 1,000 crops provided on
             https://noise.visinf.tu-darmstadt.de/

[1] A Holistic Approach to Cross-Channel Image Noise Modeling and its Application to Image Denoising. 
     Seonghyeon Nam*, Youngbae Hwang*, Yasuyuki Matsushita, Seon Joo Kim. CVPR 2016.

[2] Benchmarking Denoising Algorithms with Real Photographs. Tobias Plötz and Stefan Roth. CVPR 2017.

Dependency
------------
This code is implemented purely in Matlab2014b and doesn't depends on any other toolbox.

Contact
------------
If you have any questions or suggestions with the code, or find a bug, please let us know. 
Contact Jun Xu at csjunxu@comp.polyu.edu.hk or nankaimathxujun@gmail.com.
