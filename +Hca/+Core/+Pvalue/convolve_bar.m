function [ barOut ] = convolve_bar(bar,  sigma, hsize )
    if sigma > 0
        % kernel
        ker = circshift(images.internal.createGaussianKernel(sigma, hsize),round(hsize/2));   

        % conjugate of kernel in phase space
        multF=conj(fft(ker'));

        % convolved with sequence ->
        barOut  = ifft(fft(bar).*multF); 
    else
        barOut = bar;
    end
    
end

