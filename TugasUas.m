function ImageProccesing ()
   MainFrm = figure ( ...
    'position', [250, 50, 900 600], ...
    'color', [0.5, 0.5, 0.5], ...
    'resize', 'off'); 
    
     TitleFrm1 = axes ( ... 
    'position', [0, 0.93, 0.5, 0.2], ... 
    'color',    [0.5, 0.5, 1], ...
    'xtick',    [], ... 
    'ytick',    [], ...  
    'xlim',     [0, 1], ... 
    'ylim',     [0, 1] );
    
    Title1 = text (0.25, 0.2, 'Gambar sebelumnya', 'fontsize', 20);
    
      ImgFrm1 = axes ( ...
    'position', [0, 0.38, 0.5, 0.55], ... 
    'xtick',    [], ... 
    'ytick',    [], ...
    'xlim',     [0, 1], ... 
    'ylim',     [0, 1]);
    
    TitleFrm2 = axes ( ... 
    'position', [0.5, 0.93, 0.5, 0.2], ... 
    'color',    [0.5, 0.5, 1], ...
    'xtick',    [], ... 
    'ytick',    [], ...  
    'xlim',     [0, 1], ... 
    'ylim',     [0, 1] );
    
    
    Title2 = text (0.25, 0.2, 'Gambar sesudah', 'fontsize', 20);
    
    TitleFrm3 = axes ( ... 
    'position', [0.18, 0.09, 0.16, 0.13], ... 
    'color',    [0.6, 0.6, 0.6], ...
    'xtick',    [], ... 
    'ytick',    [], ...  
    'xlim',     [0, 1], ... 
    'ylim',     [0, 1] );
    
    ImgFrm2 = axes ( ...
    'position', [0.5, 0.38, 0.5, 0.55], ... 
    'xtick',    [], ... 
    'ytick',    [], ...
    'xlim',     [0, 1], ... 
    'ylim',     [0, 1]);
    
    BukaGambar = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'BUKA GAMBAR', ...
    'units',    'normalized', ...
    'position', [0.01, 0.3, 0.16, 0.06], ...
    'callback', { @previewImage, ImgFrm1 });
    
    GrayAverage = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Grayscale(Average)', ...
    'units',    'normalized', ...
    'position', [0.18, 0.3, 0.16, 0.06], ...
    'callback', { @GrayscaleAverage, ImgFrm2});
    
    GrayLuminos = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Grayscale(Luminocity)', ...
    'units',    'normalized', ...
    'position', [0.35, 0.3, 0.16, 0.06], ...
    'callback', { @GreyLuminocity, ImgFrm2});
    
    GrayBiner = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Grayscale to Biner', ...
    'units',    'normalized', ...
    'position', [0.52, 0.3, 0.16, 0.06], ...
    'callback', { @GreytoBiner, ImgFrm2});
    
    PutarGambar = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'PUTAR GAMBAR', ...
    'units',    'normalized', ...
    'position', [0.01, 0.23, 0.16, 0.06], ...
    'callback', { @putargambar, ImgFrm2});
    
    Filterbatas = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Filter Batas', ...
    'units',    'normalized', ...
    'position', [0.18, 0.23, 0.16, 0.06], ...
    'callback', { @filterBatas, ImgFrm2});
    
    Filtermedian = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Filter Median', ...
    'units',    'normalized', ...
    'position', [0.35, 0.23, 0.16, 0.06], ...
    'callback', { @filterMedian, ImgFrm2});
    
    Filterpereratan = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Filter Pereratan', ...
    'units',    'normalized', ...
    'position', [0.52, 0.23, 0.16, 0.06], ...
    'callback', { @filterPererataan, ImgFrm2});
    
    Citranegatif = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'CITRA NEGATIF', ...
    'units',    'normalized', ...
    'position', [0.01, 0.16, 0.16, 0.06], ...
    'callback', { @citranegatif, ImgFrm2});
    
    Hapusgambar = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'HAPUS GAMBAR', ...
    'units',    'normalized', ...
    'position', [0.01, 0.09, 0.16, 0.06], ...
    'callback', { @HapusGambar, ImgFrm1});
    
    brightnessLabel = uicontrol (MainFrm, ...
    'style',    'text', ... 
    'string',   'Brightness', ...
    'units',    'normalized', ...
    'horizontalalignment', 'left', ...
    'position', [0.225, 0.16, 0.07, 0.04]);
    
    brightness = uicontrol (MainFrm, ...
    'style',    'edit', ... 
    'string',   '', ...
    'units',    'normalized', ...
    'position', [0.21, 0.11, 0.1, 0.04], ...
    'callback', {@tambahBrightness, ImgFrm2});
    
    cerminHorizontal = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Horizontal', ...
    'units',    'normalized', ...
    'position', [0.35, 0.16, 0.08, 0.06], ...
    'callback', {@CerminHorizontal, ImgFrm2});
    
    cerminVertikal = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Vertikal', ...
    'units',    'normalized', ...
    'position', [0.43, 0.16, 0.08, 0.06], ...
    'callback', {@CerminVertikal, ImgFrm2});
    
    translasi = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Translasi', ...
    'units',    'normalized', ...
    'position', [0.35, 0.09, 0.16, 0.06], ...
    'callback', {@Translasi, ImgFrm2});
    
    penyekalan = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Penyekalan', ...
    'units',    'normalized', ...
    'position', [0.35, 0.02, 0.16, 0.06], ...
    'callback', {@PenyekalanDimensi, ImgFrm2});
    
    affine = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Transformasi Affine', ...
    'units',    'normalized', ...
    'position', [0.52, 0.09, 0.16, 0.06], ...
    'callback', {@TransformasiAffine, ImgFrm2});
    
    rotasiLabel = uicontrol (MainFrm, ...
    'style',    'text', ... 
    'string',   'Rotasi', ...
    'units',    'normalized', ...
    'horizontalalignment', 'center', ...
    'position', [0.52, 0.16, 0.09, 0.06]);
    
    rotasi = uicontrol (MainFrm, ...
    'style',    'edit', ... 
    'string',   '', ...
    'units',    'normalized', ...
    'position', [0.615, 0.16, 0.065, 0.06], ...
    'callback', {@tambahRotasi, ImgFrm2});
    
    fourier1 = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'FFT', ...
    'units',    'normalized', ...
    'position', [0.69, 0.3, 0.16, 0.06], ...
    'callback', {@Fft1, ImgFrm2});
    
    fourier2 = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'log - FFT', ...
    'units',    'normalized', ...
    'position', [0.69, 0.23, 0.16, 0.06], ...
    'callback', {@Fft2, ImgFrm2});
    
    fourier3 = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'centered - log - FFT', ...
    'units',    'normalized', ...
    'position', [0.69, 0.16, 0.16, 0.06], ...
    'callback', {@Fft3, ImgFrm2});
    
    fourier4 = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'IFFT', ...
    'units',    'normalized', ...
    'position', [0.69, 0.09, 0.16, 0.06], ...
    'callback', {@Fft4, ImgFrm2});
       
    HighPass = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'High Pass', ...
    'units',    'normalized', ...
    'position', [0.52, 0.02, 0.16, 0.06], ...
    'callback', {@highPass, ImgFrm2});
    
    LowPass = uicontrol (MainFrm, ...
    'style',    'pushbutton', ... 
    'string',   'Low Pass', ...
    'units',    'normalized', ...
    'position', [0.69, 0.02, 0.16, 0.06], ...
    'callback', {@lowPass, ImgFrm2});
    
end
    
function previewImage (hObject, eventdata, ImageFrame)
  [fname, fpath] = uigetfile();
  a = imread (fullfile(fpath, fname));
  axes(ImageFrame);
  imshow(a, []);
  axis image off
  imwrite(a, 'test.png');
end

function Brightness = tambahBrightness (hObject, eventdata, ImageFrame)
  img = imread ('test.png');
  v = get (gcbo, "string");
  b = str2num(v);
  hasil = img + b;
axes(ImageFrame);
imshow(hasil);
axis image off
end

function GrayscaleAverage (hObject, eventdata, ImageFrame)
  img = imread ('test.png');
  abu = (img(:,:,1)+img(:,:,2)+img(:,:,3))/3;
axes(ImageFrame);
imshow(abu);
axis image off
end

function GreyLuminocity (hObject, eventdata, ImageFrame)
  img = imread('test.png');
  abulumi = (img(:,:,1) * 0.2989) + (img(:,:,2) * 0.5870) + (img(:,:,3) * 0.1141)
axes(ImageFrame);
imshow(abulumi, []);
axis image off
imwrite(abulumi, 'greyluminocity.png');
end

function GreytoBiner (hObject, eventdata, ImageFrame)
  img = imread('greyluminocity.png');
  [tinggi, lebar] = size(img);
  ambang = 100
  biner = zeros(tinggi, lebar);
  for baris = 1 : tinggi
    for kolom = 1 : lebar
      if img(baris, kolom) >= ambang
        biner(baris, kolom) = 0;
      else biner(baris, kolom) = 1;
      endif
    endfor
  endfor
axes(ImageFrame);
imshow(biner, []);
axis image off
end   

function putargambar (hObject, eventdata, ImageFrame)
  gambar = imread('test.png');
  [tinggi, lebar] = size(gambar);
  a = imcrop(gambar, [1 1 tinggi/2 lebar/2]);
  b = imcrop(gambar, [1 lebar/2 tinggi/2 lebar/2]);
  c = imcrop(gambar, [tinggi/2 1 tinggi/2 lebar/2]);
  d = imcrop(gambar, [tinggi/2 lebar/2 tinggi/2 lebar/2]);
  e = horzcat(b,a);
  f = horzcat(d,c);
  g = vertcat(e,f);
axes(ImageFrame);
imshow(g, []);
axis image off
endfunction

function citranegatif (hObject, eventdata, ImageFrame)
  gambar = imread('test.png');
  g = 255 - gambar;
axes(ImageFrame);
imshow(g, []);
axis image off
endfunction

function filterBatas (hObject, eventdata, ImageFrame)
  F = imread('test.png');
  ukuran = size(F);
  tinggi = ukuran(1);
  lebar = ukuran(2);
  G = F;
  for baris=2 : tinggi-1
    for kolom=2 : lebar-1
      minPiksel = min([F(baris-1, kolom-1) 
      F(baris, kolom+1)
      F(baris-1, kolom)
      F(baris, kolom-1)
      F(baris, kolom+1)
      F(baris+1, kolom-1)
      F(baris+1, kolom) 
      F(baris+1, kolom+1)]);
      
      maksPiksel = max([F(baris-1, kolom-1)
      F(baris-1, kolom) 
      F(baris, kolom+1)
      F(baris, kolom-1)
      F(baris, kolom+1) 
      F(baris+1, kolom-1)
      F(baris+1, kolom) 
      F(baris+1, kolom+1)]);
      
      if F(baris, kolom) < minPiksel
         G(baris, kolom) = minPiksel;
      else
         if F(baris, kolom) > maksPiksel
            G(baris, kolom) = maksPiksel;
         else 
            G(baris, kolom) = F(baris, kolom);
         endif
      endif
    endfor
  endfor

axes(ImageFrame);
imshow(G, []);
axis image off
clear;
endfunction

function filterMedian (hObject, eventdata, ImageFrame)
  F = imread('test.png');
  [tinggi,lebar] = size(F);
  for baris=2 : tinggi-1
    for kolom=2 : lebar-1
      data = [F(baris-1, kolom-1) ...
      F(baris-1, kolom) ...
      F(baris-1, kolom+1) ...
      F(baris, kolom-1) ...
      F(baris, kolom) ...
      F(baris, kolom+1) ...
      F(baris+1, kolom-1) ...
      F(baris+1, kolom) ...
      F(baris+1, kolom+1)];
      
      for i=1 : 8
        for j = i+1 : 9
          if data(i) > data(j)
            tmp = data (i);
            data(i) = data(j);
            data(j) = tmp;
          endif
        endfor
      endfor
      G(baris, kolom) = data(5);
    endfor
  endfor
 
axes(ImageFrame); 
imshow(G, []);
axis image off
clear;
endfunction

function filterPererataan (hObject, eventdata, ImageFrame)
    F = imread('test.png');
  [tinggi, lebar] = size(F);
  F2= double(F);
  for baris = 2 : tinggi-1
    for kolom = 2 : lebar-1
      jum = F2(baris-1, kolom-1)+ ...
            F2(baris-1, kolom) + ...
            F2(baris-1, kolom-1)+ ...
            F2(baris, kolom-1) + ...
            F2(baris, kolom) + ...
            F2(baris, kolom+1) + ...
            F2(baris+1, kolom-1)+ ...
            F2(baris+1, kolom)+ ...
            F2(baris+1, kolom+1);
      G(baris, kolom)= uint8(1/9 * jum);
    endfor
  endfor
axes(ImageFrame);
imshow(G, []);
axis image off
clear;
endfunction

function HapusGambar()
  cla reset;
endfunction  

function CerminHorizontal (hObject, eventdata, ImageFrame)
  img = imread ('test.png');
  [tinggi, lebar] = size(img);
  for y = 1 : tinggi
    for x = 1 : lebar
      x2 = lebar - x + 1;
      y2 = y;
      G(y, x) = img(y2, x2);
    endfor
  endfor
G = uint8(G);
axes(ImageFrame);
imshow(G);
axis image off
end

function CerminVertikal (hObject, eventdata, ImageFrame)
  img = imread ('test.png');
  [x, y, z] = size(img);
   
% Reverse elements of each column
% in each image plane (dimension)
  for plane = 1 : z
      len = x;
      for i = 1 : x 
          for j = 1 : y
   
            % To reverse the order of the element
            % of a column we can swap the 
            % topmost element of the row with
            % its bottom-most element  
              
              if i < x/2 
                  temp = img(i, j, plane);
                  img(i, j, plane) = img(len, j, plane);
                  img(len, j, plane) = temp; 
              end
          end
          len = len - 1;
      end
  end

axes(ImageFrame);
imshow(img);
axis image off
end

function Translasi (hObject, eventdata, ImageFrame)
  img = imread ('test.png');
  [tinggi, lebar] = size(img);
  G = zeros(size(img));
  G = uint8(G);
  for y = 1 : tinggi
    for x = 1 : lebar
      if (y + 20 >= 1) && (y + 20 <= tinggi) && ...
        (x + 20 >= 1) && (x + 20 <= lebar)
              G(y + 20, x + 20) = img(y, x); 
      endif
    endfor
  endfor
axes(ImageFrame);
imshow(G);
axis image off
end

function Rotasi = tambahRotasi (hObject, eventdata, ImageFrame)
  img = imread ('test.png');
  v = get (gcbo, "string");
  b = str2num(v);
  
  [tinggi, lebar] = size(img);
  G = zeros(size(img));
  G = uint8(G);
  rad = pi * b/180;
  cosa = cos(rad);
  sina = sin(rad);
  for y = 1 : tinggi
    for x = 1 : lebar
      x2 = round(x * cosa + y * sina);
      y2 = round(y * cosa - x * sina);
      if (y2 >= 1) && (y2 <= tinggi) && ...
        (x2 >= 1) && (x2 <= lebar)
        G(y, x) = img(y2, x2);
      else 
        G(y,x) = 0;
      endif
    endfor
  endfor
axes(ImageFrame);
imshow(G);
axis image off
end

function TransformasiAffine (hObject, eventdata, ImageFrame)
  img = imread ('test.png');
  [tinggi, lebar] = size(img);
  for y = 1 : tinggi
    for x = 1 : lebar
      x2 = 1 * x + 0.15 * y + 0;
      y2 = 0 * x + 1 * y + 0;
      if (x2 >= 1) && (x2 <= lebar) && ...
          (y2 >= 1) && (y2 <= tinggi)
          p = floor(y2);
          q = floor(x2);
          a = y2 - p;
          b = x2 - q;
          if (floor(x2) == lebar) || ...
            (floor(y2) == tinggi)
            G(y, x) = img(floor(y2), floor(x2));
          else 
            intensitas = (1 - a) * ((1 - b) * img(p, 1) + ...
            b * img(p, q + 1)) +       ...
            a * ((1 - b) * img(p + 1, q) + ...
            b * img(p + 1, q + 1));
            G(y, x) = intensitas;
          endif
      else 
          G (y, x) = 0;
      endif
    endfor
  endfor
  G = uint8(G);
axes(ImageFrame);
imshow(G);
axis image off
end

function PenyekalanDimensi (hObject, eventdata, ImageFrame)
img = imread ('test.png');

F = img;
sy = 2;
sx = 2;

[tinggi, lebar] = size(F);

tinggi_baru = tinggi * sy;
lebar_baru = lebar * sx;

F2 = double(F);
for y=1 : tinggi_baru
    y2 = ((y-1)/sy) + 1;
    for x=1 : lebar_baru
        x2 = ((x-1) / sx) + 1;
        G(y, x) = F(floor(y2), floor(x2));
    endfor
endfor
G = uint8(G);
axes(ImageFrame);
imshow(G);
axis image off
end

function Fft1 (hObject, eventdata, ImageFrame)
img = imread('test.png');
F = fft2(img);
axes(ImageFrame);
imshow(abs(F),[]);
axis image off
endfunction

function Fft2 (hObject, eventdata, ImageFrame)
img = imread('test.png');
axes(ImageFrame);
F = fft2(img);
imshow(log(abs(F)),[]);
axis image off
endfunction

function Fft3 (hObject, eventdata, ImageFrame)
img = imread('test.png');
axes(ImageFrame);
F = fft2(img);
F = fftshift(F);
imshow(log(abs(F)),[]);
axis image off
endfunction

function Fft4 (hObject, eventdata, ImageFrame)
img = imread('test.png');
axes(ImageFrame);
F = ifft2(img);
imshow(abs(F),[]);
axis image off
endfunction

function highPass (hObject, eventdata, ImageFrame)
a = imread('test.png');

input_image = a;

[M, N] = size(input_image);
  
FT_img = fft2(double(input_image));
n = 2; % one can change this value accordingly
  
% Assign Cut-off Frequency
D0 = 10; % one can change this value accordingly
  
% Designing filter
u = 0:(M-1);
v = 0:(N-1);
idx = find(u > M/2);
u(idx) = u(idx) - M;
idy = find(v > N/2);
v(idy) = v(idy) - N;
  
[V, U] = meshgrid(v, u);
D = sqrt(U.^2 + V.^2);
H = 1./(1 + (D0./D).^(2*n));
G = H.*FT_img;
output_image = real(ifft2(double(G)));

axes(ImageFrame);
imshow(output_image, [ ]);
axis image off
  
endfunction

function lowPass (hObject, eventdata, ImageFrame)
a = imread('test.png');

input_image = a;

[M, N] = size(input_image);
  
FT_img = fft2(double(input_image));
   
D0 = 30; 
  
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;
  
[V, U] = meshgrid(v, u);
  
D = sqrt(U.^2+V.^2);
  
H = double(D <= D0);
  
G = H.*FT_img;
   
output_image = real(ifft2(double(G)));

axes(ImageFrame);
imshow(output_image, [ ]);
axis image off
  
endfunction 
