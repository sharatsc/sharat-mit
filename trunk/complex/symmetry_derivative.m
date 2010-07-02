%---------------------------------------------------------------
%
%
%sharat@mit.edu
%---------------------------------------------------------------
function g  = symmetry_derivative(n,sigma)
    [x,y]   = meshgrid(-3*sigma:0.5:3*sigma,-3*sigma:0.5:3*sigma);
    g       = exp(-(x.^2+y.^2)/(2*sigma.^2))/(2*pi*sigma.^2);
    for i = 1:abs(n)
      rd       =complex_image_gradient(real(g));
      id       =complex_image_gradient(imag(g));
      if(n>0)
	g        =complex(real(rd)-imag(id),imag(rd)+real(id));
      else
	g        =complex(real(rd)+imag(id),-imag(rd)+real(id));
      end;
    end;
    g        = imresize(g,ceil([6*sigma 6*sigma]),'bilinear');
%end function 
