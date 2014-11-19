%Program to remove false contour By Dithering.
x=input('enter any image :');
s=imread(x);
s=rgb2gray(s);
s=imresize(s,[256 256],'bicubic');
s=double(s);
figure,imshow(uint8(s)),title('Original image');
noise = rand(256,256); %generates pseudorandom noise.
d=s+noise;
% quantization to 4 bits.
q4 = floor(d/16);
figure,imshow((q4-noise)),title('dithered 4 bit image');
%quantization to 2 bits.
q2 = floor(d/64);
figure,imshow((q2-noise)),title('dithered 2 bit image');
%quantization to 1bit.
q1 = floor(d/128);
figure, imshow((q1-noise)),title('dithered 1 bit image');
