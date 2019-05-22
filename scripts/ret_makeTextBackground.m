%% scipt that makes a page of text
% ret_makeTextBackground is the original version

clear all; clc; close all; 

%% modify here

% CHANGE THIS. path to the .txt file with words in a column
txtfile = 'words_list.txt'; 

% background color
bg_color        = [255];

% word color
word_color      = [0];                      

% CHANGE THIS. font name
word_font       = 'Helvetica';

% font size
word_fontSize   = 26;                       

% spacing between words of the same line, in units of pixels
word_spaceSize  = 8;                       

% assumes square. dim of object image
res             = 768;               

% whether or not we want the word to be bold
word_bold       = false;                    

% default is 2. With smaller font sizes we want larger (5 or larger)
word_sampsPerPt = 2;                        

% number of images
numImages       = 100;  

% CHANGE THIS. what to save the image matrix as. 
nameSave = 'Words_Spanish.mat';                

% CHANGE THIS. directory to save the matrix.
dirSave = './';

%% list of words to use
fid = fopen(txtfile);
L = textscan(fid,'%s');
fclose(fid);

% reads into struct, make into cell
L = L{1}; 


%% about the images
% Kendrick's version 768 x 768 x 3 x 100
% Specifically, 100 RGB images of resolution 768 x 768 


%% about some functions we'll use

% img = renderText(someText, fontName, fontSize, sampsPerPt, [antiAlias], ...
% [fractionalMetrics], [bold])


%% create the grascale images first

% master storage
II = uint8(zeros(res,res,3,numImages)); 

% blank (white) background
% we don't assign it to the color of the background because
% renderText makes an image of 1s and 0s, where 0 is the background and 1 is the word
% it is easier to go back and change this later. 
tem  = uint8(zeros(res, res)); 

Iheight = size(tem,1);
Ilength = size(tem,2);


for ii = 1: numImages

% initialize starting position
% origin (1,1) is upper left
wx = 1; 
wy = 1; 

while wy < Iheight
    
    while (wx < Ilength) && (wy < Iheight)

        % pick a random word
        wordIndR    = randi(length(L)); 
        wordImg     = renderText([L{wordIndR} '  '], word_font, word_fontSize, word_sampsPerPt, [], [], word_bold); 
        wordHeight  = size(wordImg,1);
        wordLength  = size(wordImg,2); 

        % add this word to the canvas
        tem(wy:wy+wordHeight-1,wx:wx+wordLength-1) = wordImg; 


        % update horizontal position
        wx = wx + wordLength; 

        % go to next line once we get to the end of this one
        if wx > Ilength; 
            wy = wy + wordHeight; 
            wx = 1; 
        end

    end
    
end

% crop the images so to fit
% (matrix may become enlarged because of last word in line)
tem = tem(1:res,1:res); 

% see what the image looks like
fig=figure('visible','off');
imagesc(tem);
axis off; axis square; 
title(['fontsize: ' num2str(word_fontSize) '. sampsPerPt: ' num2str(word_sampsPerPt)])
colormap gray

I = uint8(zeros(size(tem)));

% renderText makes an image of 1s and 0s, where 0 is the background and 1 is the word
% change the 1s to be word color and 0s to be background color
I(tem==0) = bg_color; 
I(tem==1) = word_color;

% make multi channel
Ithree = cat(3, I,I,I);

% store in master 
II(:,:,:,ii) = Ithree; 

% track progress
imshow(Ithree); title(ii)
ii

end


%% convert images to RGB and into class uint8
% maybe not necessary, but trying to prevent kendrick's code from crashing


%% convert to cell
images      = cell(1,1);
images{1}   = II; 


%% saving
% save('words_small4_2.mat','images'); 
save(fullfile(dirSave, [nameSave]), 'images')
