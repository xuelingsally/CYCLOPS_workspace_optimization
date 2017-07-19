function max_length = x_space_length(a, B, p)

max_length = 0;

for i=1:6
    a_b = a(:,i) + p;
    b = B(:,i);
    x_length = a_b(1) - b(1);
    
    if max_length < x_length
        max_length = x_length;
    end
    
end

end