function out = gauss_filter(img,sx,sy)
  if(nargin<3)
    sy = sx;
  end;
  hx = sum(fspecial('gaussian',ceil(8*sx),sx),1);
  hy = sum(fspecial('gaussian',ceil(8*sy),sy),1);
  nimg  = ones(size(img));
  img   = conv2(hy,hx,img,'same');
  nimg  = conv2(hy,hx,nimg,'same');
  out   = img./nimg;
%end function
