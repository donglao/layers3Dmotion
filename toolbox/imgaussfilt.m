function Output = imgaussfilt(A, sigma)
if nargin<2
    sigma=0.5;
end
I = double(A);

%Design the Gaussian Kernel
%Standard Deviation
%Window size
sz = 5;
[x,y]=meshgrid(-sz:sz,-sz:sz);

M = size(x,1)-1;
N = size(y,1)-1;
Exp_comp = -(x.^2+y.^2)/(2*sigma*sigma);
Kernel= exp(Exp_comp)/(2*pi*sigma*sigma);
 

%Initialize
Output=zeros(size(I));
%Pad the vector with zeros
I = padarray(I,[sz sz]);

%Convolution
for i = 1:size(I,1)-M
    for j =1:size(I,2)-N
        Temp = I(i:i+M,j:j+M).*Kernel;
        Output(i,j)=sum(Temp(:));
    end
end
%imshow(Output);
end