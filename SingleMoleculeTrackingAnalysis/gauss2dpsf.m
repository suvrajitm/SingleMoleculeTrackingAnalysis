function psf = gauss2dpsf(param,data)
x = data(:,1); 
y = data(:,2);

Ibg = param(1);
I0 = param(2); 
x0 = param(3);
y0 = param(4);
sx = param(5);
sy = param(6);

psf = Ibg + I0*exp(-((x-x0).^2/(2*(sx)^2))-((y-y0).^2)/(2*(sy)^2));
