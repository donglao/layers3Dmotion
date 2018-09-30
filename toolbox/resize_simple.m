function output = resize_simple(input,h_mod,w_mod,fill)
[h,w,channels]=size(input);
a = round((h_mod - h)/2)+1;
b = round((w_mod - w)/2)+1;
output = ones(h_mod,w_mod,channels)*nan;
output(a:a+h-1,b:b+w-1,:) = input;
output(isnan(output)) = fill;
output = cast(output, class(input));

end